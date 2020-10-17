import mdtraj as md
from protlego.builder.chimera import Chimera
from moleculekit.projections.metricdistance import MetricSelfDistance
from moleculekit.projections.metricdistance import contactVecToMatrix
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import re

import logging

logger = logging.getLogger('protlego')


def calc_sasa(chimera: Chimera = None, filename: str = None,
              probe_radius: float = 0.14, n_sphere_points: int = 960, sasa_type='total'):
    """

    Computes the Solvent Accessible Surface Area of the protein.
    This funcion uses the MDtraj shrake_rupley implementation as a basis.
    :param chimera: A Chimera object.
    :param filename: Path to a pdb file
    :param probe_radius: The radius of the probe, in nm.
    :param n_sphere_points: the number of points representing the sufrace of each atom. Higher values lead to more accuracy.
    :param sasa_type: Type of calculation to perform. To select from polar, apolar, or total.
    :return: areas: np.array containing the area of the chimera in Angstrom^2
    """
    sasa_types = ["polar", "apolar", "total"]
    if sasa_type not in sasa_types:
        raise ValueError(f"Invalid type. Expected one of {sasa_types}")
    if chimera and filename:
        raise ValueError("Only a Chimera object or the path to a pdb file must be specified")
    if not chimera and not filename:
        raise ValueError("At least a Chimera object or the path to a pdb file must be specified")
    if chimera:
        filename = "/tmp/structure.pdb"
        chimera.write(filename)

    polars = ['SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS']
    apolars = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'TRP', 'PHE', 'PRO']
    structure = md.load(filename)
    if sasa_type == 'polar':
        indices = [index for index, residue in enumerate(structure.topology.residues) if residue.name in polars]
    elif sasa_type == 'apolar':
        indices = [index for index, residue in enumerate(structure.topology.residues) if residue.name in apolars]
    else:
        indices = [index for index, residue in enumerate(structure.topology.residues)]

    sasa = md.shrake_rupley(structure, probe_radius=probe_radius, n_sphere_points=n_sphere_points, mode="residue")
    area = sasa[0][indices].sum()
    logger.info(f"Area is {area} (nm)^2")
    return area


def calc_dssp(chimera: Chimera = None, filename: str = None, simplified: bool = True):
    """
    Compute Dictionary of protein secondary structure (DSSP) secondary structure assignments.
    This funcion uses the MDtraj compute_dssp implementation as a basis.
    :param chimera: A Chimera object.
    :param filename: path to a pdb file
    :param simplified: Use the simplified 3-category assignment scheme. Otherwise the original 8-category scheme is used.
    :return: assignments np.ndarray. The secondary structure assignment for each residue
    """
    if chimera and filename:
        raise ValueError("Only a Chimera object or the path to a pdb file must be specified")
    if not chimera and not filename:
        raise ValueError("At least a Chimera object or the path to a pdb file must be specified")
    if chimera:
        filename = "/tmp/structure.pdb"
        chimera.write(filename)
    structure = md.load(filename)
    dssp = md.compute_dssp(structure, simplified=simplified)
    return dssp


def calc_dist_matrix(chimera: Chimera = None, filename: str = None, selection: str = 'residue', type='contacts',
                     plot=False):
    """
    Returns a matrix of C-alpha distances for a given pdb
    :param chimera: A Chimera object with n residues.
    :param filename: path to a pdb file
    :param selection: How to compute the distance. 'residue' (the closest two
    :param type: between contacts (contact map when distances are below 8 armstrongs) or distances atoms between two residues) or 'alpha' distance of the alpha carbons.
    :param plot: whether to plot the distance matrix. Default is False
    :return: matrix. np.array. An n by n distance matrix.
    """
    if chimera and filename:
        raise ValueError("Only a Chimera object or the path to a pdb file must be specified")
    if not chimera and not filename:
        raise ValueError("At least a Chimera object or the path to a pdb file must be specified")
    if filename:
        chimera = Chimera(filename=filename)

    if selection == 'residue':
        metr = MetricSelfDistance("protein", groupsel="residue", metric="distances", pbc=False)
        mapping = metr.getMapping(chimera)
        a = metr.project(chimera)
        matrix, _, _ = contactVecToMatrix(a[0], mapping.atomIndexes)
    elif selection == 'alpha':
        metr = MetricSelfDistance("protein and name CA", metric="distances", pbc=False)
        a = metr.project(chimera)
        mapping = metr.getMapping(chimera)
        matrix, _, _ = contactVecToMatrix(a, mapping.atomIndexes)
    else:
        raise ValueError("Specify a selection type: 'residue' or 'atom'")
    if type == "contacts":
        matrix = matrix < 8
    elif type != "contacts" and type != "distances":
        raise ValueError("Please select contact type between 'contacts' or distances")

    if plot:
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(111)
        cmap = 'binary'
        cax = ax.imshow(matrix, cmap=matplotlib.cm.get_cmap(cmap), interpolation='nearest', origin="lower")
        if type == 'distances':
            cmap = 'gist_rainbow'
            cax = ax.imshow(matrix, cmap=matplotlib.cm.get_cmap(cmap), interpolation='nearest', origin="lower")
            cbar = fig.colorbar(cax, cmap=matplotlib.cm.get_cmap(cmap))
        plt.xlabel('xlabel', fontsize=24)
        plt.ylabel('ylabel', fontsize=24)
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel("Residue index")
        plt.ylabel("Residue index")

    return matrix


# def calc_FCR(chimera:Chimera=None,pH:float=None) -> float:
#     """
#     Get the fraction of charged residues in the sequence. (pH keyword allows for a pH specific value)
#     This uses the localcider toolkit as a basis
#     Please refer to it for further documentation:
#     http://pappulab.github.io/localCIDER/
#     :param chimera: A chimera object
#     :param pH: The pH at which compute the charge seggregation
#     :return: Fraction of charged residues  (FCR)
#     """
#     sequence =chimera.sequence()['0']
#     SeqOb = SequenceParameters(sequence)
#     return SeqOb.get_FCR(pH=pH)
#
# def calc_kappa(chimera:Chimera) -> float:
#     """
#     Get the sequence’s kappa value.
#     κ is a parameter to describe the extent of charged amino acid mixing in a sequence.
#     (http://pappulab.wustl.edu/CIDER/help/)
#     This uses the localcider toolkit as a basis
#     Please refer to it for further documentation:
#     http://pappulab.github.io/localCIDER/
#     :param chimera: A chimera object
#     :return: kappa
#     """
#     sequence =chimera.sequence()['0']
#     SeqOb = SequenceParameters(sequence)
#     return SeqOb.get_kappa()


def calc_contact_order(chimera: Chimera = None, filename: str = None, diss_cutoff: int = 8):
    """
    The contact order of a protein is a measure of the locality of the inter-amino acid contacts in the
    native folded state. It is computed as the average seqeuence distance between residues that form contacts
    below a threshold in the folded protein divided by the total length of the protein"
    :param chimera: A Chimera object with n residues.
    :param filename: path to a pdb file
    :param diss_cutoff: The maximum distance in Armstrong between two residues to be in contact, default 8 Angstroms
    :return: the contact order (%)
    """
    if chimera and filename:
        raise ValueError("Only a Chimera object or the path to a pdb file must be specified")
    if not chimera and not filename:
        raise ValueError("At least a Chimera object or the path to a pdb file must be specified")
    if filename:
        chimera = Chimera(filename=filename)
    chimera.renumberResidues()
    metr = MetricSelfDistance("protein and noh", groupsel="residue", metric="contacts", threshold=diss_cutoff,
                              pbc=False)
    a = metr.project(chimera)
    mapping = metr.getMapping(chimera)
    matrix, _, _ = contactVecToMatrix(a[0], mapping.atomIndexes)
    triang = np.triu(matrix)
    idx1, idx2 = np.where(triang)
    total_contacts = len(idx1)
    total_residues = chimera.numResidues
    summation = np.sum(idx2 - idx1)
    co = 1 / (total_contacts * total_residues) * summation
    print(f"Contact order is {co*100} %")
    return co * 100


def hhbond_plot(chimera: Chimera = None, filename: str = None):
    """
    Computes a hhbond plot of a chimera object or a file. One of the two inputs
    must be provided.
    :param chimera: the chimera from where to compute the hydrogen bond plot
    :param filename: a path where to find the pdb file.
    :return: A contact map with hydrogen bond as metric
    """
    if chimera and filename:
        raise ValueError("Only a Chimera object or the path to a pdb file must be specified")
    if not chimera and not filename:
        raise ValueError("At least a Chimera object or the path to a pdb file must be specified")
    if chimera:
        filename = "/tmp/structure.pdb"
        chimera.write(filename)
    mol = md.load_pdb(filename)
    hbonds = md.baker_hubbard(mol, periodic=False)
    for hbond in hbonds:
        a = mol.topology.atom(hbond[0])
        residue1 = re.findall("\d+", str(a))[0]
        resid1 = int(residue1)
        b = mol.topology.atom(hbond[2])
        residue2 = re.findall("\d+", str(b))[0]
        resid2 = int(residue2)
        plt.plot(resid1, resid2, 'b.')
    plt.show()
