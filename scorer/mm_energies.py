from typing import Tuple
from simtk.openmm.app import *
from parmed import load_file
from simtk.openmm.app.modeller import Modeller
import simtk.openmm as mm
from simtk import unit
import os
import shutil

from moleculekit.tools.preparation import proteinPrepare
from moleculekit.tools.autosegment import autoSegment
from protlego.builder.chimera import Chimera
from protlego.builder.builder import Builder, NotCorrectPDBError
from moleculekit.molecule import Molecule

import logging

logger = logging.getLogger('protlego')

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

    if not os.path.exists(output):
        os.mkdir(output)

    smol = prepare_protein(chimera)
    smol.write(f"{output}/protein.pdb")
    pdb = PDBFile(f"{output}/protein.pdb")
    parm = load_file(f"{output}/protein.pdb")
    modeller = Modeller(pdb.topology, pdb.positions)

    if ff == 'amber':
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    if ff == 'charmm':
        forcefield = ForceField('charmm36.xml', 'charmm36/tip3p-pme-b.xml')

    modeller.addSolvent(forcefield, padding=1.0 * unit.nanometer)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
                                     nonbondedCutoff=1 * unit.nanometer, constraints=HBonds)
    if restraint_backbone:
        # Applies an external force on backbone atoms
        # This allows the backbone to stay rigid, while severe clashes can still be resolved
        force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", 5.0 * unit.kilocalories_per_mole / unit.angstroms ** 2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        atoms = [atom for atom in modeller.topology.atoms()]
        for idx, atom_crd in enumerate(modeller.positions):
            if idx >= len(parm.atoms): continue
            if parm.atoms[idx] in ('CA', 'C', 'N'):
                force.addParticle(idx, atom_crd.value_in_unit(unit.nanometers))
        system.addForce(force)

    integrator = mm.LangevinIntegrator(temperature, friction, error_tolerance)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Get pre-minimization energy (scoring)
    state = simulation.context.getState(getEnergy=True, getForces=True)
    pre_energy = state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    logger.info("Energy before minimization", pre_energy)

    # Standard values for the integrator and tolerance constraint
    if not cuda:
        # Default value
        tolerance = 10 * unit.kilojoule / unit.mole
    else:
        # High tolerance so the CPU only pre-minimizes
        tolerance = 1e6

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
        gpu_min = Simulation(modeller.topology, system, gpu_integrator, platform, properties)
        gpu_min.context.setPositions(min_coords.getPositions())
        gpu_min.minimizeEnergy()
        post_position = gpu_min.context.getState(getPositions=True).getPositions()
        post_state = gpu_min.context.getState(getEnergy=True, getForces=True)

    post_energy = post_state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    logger.info("Energy after minimization", post_energy)

    PDBFile.writeFile(modeller.topology, post_position, open(f"{output}/structure_minimized.pdb", 'w'), keepIds=True)
    min_mol = Chimera(filename=f"{output}/structure_minimized.pdb")

    if keep_output_files is False:
        shutil.rmtree(output)

    return post_energy, min_mol