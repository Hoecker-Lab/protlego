from typing import Tuple
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import parmed
import os
import shutil
import htmd.builder.amber as amber
import htmd.builder.charmm as charmm
from glob import glob

from moleculekit.tools.preparation import proteinPrepare
from moleculekit.tools.autosegment import autoSegment
from protlego.builder.chimera import Chimera
from protlego.builder.builder import Builder, NotCorrectPDBError
from moleculekit.molecule import Molecule

import logging

logger = logging.getLogger('protlego')
logging.getLogger('htmd').setLevel(logging.ERROR)

temperature = 300 * unit.kelvin
friction = 1.0 / unit.picoseconds
error_tolerance = 2.0 * unit.femtoseconds
distance_tolerance = 0.00001


def prepare_protein(chimera: Chimera) -> Molecule:
    """
    Builds a water box around the chimera
    :param chimera:
    :return:
    """
    non_standards = Builder.find_nonstandards(chimera)
    if non_standards:
        raise NotCorrectPDBError(f"PDB presents the non_standards residues {non_standards}."
                                 f" Call remove_residue() or Builder.mutate_nonstandards()"
                                 f" if you wish to minimize.")

    mol = proteinPrepare(chimera)
    if 'AR0' in mol.resname:
        mol.set("resname", "ARG", "resname AR0")
    mol = autoSegment(mol)
    mol.center()
    # D = maxDistance(mol,'all') + 6.0
    # smol=solvate(mol,minmax=[[-D, -D, -D], [D, D, D]])
    return mol


def _prepareMMsim(smol, ff: str, output: str = "/tmp/build"):
    if ff == "amber":
        mol = amber.build(smol, ionize=False, outdir=output)
        structure_prm = app.AmberPrmtopFile(f'{output}/structure.prmtop')
        inpcrd = app.AmberInpcrdFile(f'{output}/structure.crd')
        system = structure_prm.createSystem()
    elif ff == "charmm":
        mol = charmm.build(smol, ionize=False, outdir=output, psfgen='/sw/sci/app/namd/NAMD_2.9-x86_64-CUDA/psfgen')
        structure_prm = parmed.charmm.CharmmPsfFile(f'{output}/structure.psf')
        prmfiles = glob(os.path.join(output, 'topologies', '*.rtf'))
        os.rename(f"{output}/parameters", f"{output}/structure.prm")
        prmfiles.append(f"{output}/structure.prm")
        inpcrd = app.PDBFile(f'{output}/structure.pdb')
        params = parmed.charmm.CharmmParameterSet(*prmfiles)
        system = structure_prm.createSystem(params)
    else:
        raise ValueError("The specified forcefield is wrong. Please select amber or charmm")
    platform = mm.Platform.getPlatformByName('CPU')
    cpu_integrator = mm.VariableLangevinIntegrator(temperature, friction, error_tolerance)
    simulation = app.Simulation(structure_prm.topology, system, cpu_integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    state = simulation.context.getState(getEnergy=True, getForces=True)

    return state, system, cpu_integrator, simulation, structure_prm, inpcrd, mol


def score_pot_energy(chimera: Chimera, ff: str, output: str = "/tmp/build", keep_output_files: bool = False) \
        -> Tuple[Chimera, unit.quantity.Quantity]:
    """
    Score the potential energy of a chimera with the amber or charmm forcefields
    :param chimera: A chimera object to compute the pot. E on.
    :param ff: Forcefield to use. Amber or charmm. The latest version with be used.
    :param output: A folder where to keep the files. If not provided they will be stored in
    the /tmp folder and later removed.
    :param keep_output_files: Whether to keep the produced folders. Default is remove.

    :return: The chimera object that whose pot.E was computed and the potential energy value.
    """
    mol = prepare_protein(chimera)
    state, _, _, _, _, _, mol = _prepareMMsim(mol, ff, output)
    energy = state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    if keep_output_files is False:
        shutil.rmtree(output)
    print(f"Potential Energy: {energy}")
    return mol, energy


def _mol_chimera_wrapper(molecule: Molecule, chimera: Chimera) -> Chimera:
    molecule.write("/tmp/molecule.pdb")
    new_chimera = Chimera(filename="/tmp/molecule.pdb")
    os.remove("/tmp/molecule.pdb")

    return new_chimera


def minimize_potential_energy(chimera, ff: str,
                              output: str = "/tmp/build", keep_output_files=False, cuda=False,
                              restraint_backbone: bool = True) -> Tuple[unit.quantity.Quantity, Chimera]:
    """
    :param chimera: A chimera object where to perform the minimization
    :param forcefield: The forcefield to use for the minimization. Select between "amber" and "charmm"
    :param output: A folder where to keep the files. If not provided they will be stored in the /tmp folder and later removed.
    :param cuda: Whether to use GPU acceleration
    :param restraint_backbone: Keep the backbone atoms constraint in space

    :return: The chimera object that was minimized and the potential energy value.
    """
    smol = prepare_protein(chimera)
    state, system, integrator, simulation, structure_prm, inpcrd, mol = _prepareMMsim(smol, ff, output)
    pre_energy = state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    print("Energy before minimization", pre_energy)
    mol = _mol_chimera_wrapper(mol, chimera)
    pre_min_coords = simulation.context.getState(getPositions=True)
    # Standard values for the integrator and tolerance constraint
    if not cuda:
        # Default value
        tolerance = 10 * unit.kilojoule / unit.mole
    else:
        # High tolerance so the CPU only pre-minimizes
        tolerance = 1e6

    if restraint_backbone:
        # Applies an external force on backbone atoms
        # This allows the backbone to stay rigid, while severe clashes can still be resolved
        force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", 5.0 * unit.kilocalories_per_mole / unit.angstroms ** 2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for i, atom_crd in enumerate(inpcrd.positions):
            if mol.name[i] in ('CA', 'C', 'N'):
                force.addParticle(i, atom_crd.value_in_unit(unit.nanometers))
        system.addForce(force)

    # Setup CPU minimization
    integrator.setConstraintTolerance(distance_tolerance)
    simulation.minimizeEnergy(tolerance=tolerance)
    post_position = simulation.context.getState(getPositions=True).getPositions()
    post_state = simulation.context.getState(getEnergy=True, getForces=True)
    if cuda:
        min_coords = simulation.context.getState(getPositions=True)
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        gpu_integrator = mm.VariableLangevinIntegrator(temperature, friction, error_tolerance)
        gpu_integrator.setConstraintTolerance(distance_tolerance)
        gpu_min = app.Simulation(structure_prm.topology, system, gpu_integrator, platform, properties)
        gpu_min.context.setPositions(min_coords.getPositions())
        gpu_min.minimizeEnergy()
        post_position = gpu_min.context.getState(getPositions=True).getPositions()
        post_state = gpu_min.context.getState(getEnergy=True, getForces=True)

    post_energy = post_state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    print("Energy after minimization", post_energy)

    app.PDBFile.writeFile(structure_prm.topology, post_position, open(f"{output}/structure_minimized.pdb", 'w'),
                          keepIds=True)
    min_mol = Chimera(filename=f"{output}/structure_minimized.pdb")
    mol.coords = min_mol.coords
    mol.write(f"{output}/structure_minimized.pdb")
    if keep_output_files is False:
        shutil.rmtree(output)

    return post_energy, mol
