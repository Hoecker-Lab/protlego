from moleculekit.molecule import Molecule
from moleculekit.vmdviewer import viewer
import numpy as np
from scipy.spatial.distance import cdist
from collections import Counter
import sys, getopt
from graph_tool.all import *
import math
from typing import Tuple, Dict

import logging
import moleculekit

logging.getLogger(moleculekit.molecule.__name__).setLevel(logging.WARNING)
logger = logging.getLogger('hydrophobic_clusters')
from typing import NamedTuple


class Cluster(NamedTuple):
    area: np.array
    residues: list
    contacts: int
    ratio_contacts_residue: float
    ratio_area_residue: float


# Paths and variables
sel = "protein and chain A and not backbone and noh and resname ILE VAL LEU"
_ATOMIC_RADII = {'C': 1.88}
water_radius = 1.4
sphere_radius_carbon = _ATOMIC_RADII['C'] + water_radius
sphere_points = 610
sphere_area_const = 4.0 * math.pi * (sphere_radius_carbon ** 2) / sphere_points


class Atom:
    """
    Class to handle ILV atoms
    index: the vmd atom index
    mol: an htmd.htmd.molecule object
    """
    radius = _ATOMIC_RADII['C']  # All distances in Angstrom

    def __init__(self, index: int, mol: Molecule) -> None:
        self.index = index
        self.coords = mol.coords[index][:, 0]
        self.resid = mol.resid[index]
        self.point_coords = generate_sphere_points(self.coords, sphere_points, self.radius)
        self.neighbor_indices = self.get_neighbors(mol)

    def get_neighbors(self, mol: Molecule) -> np.array:
        """
        Provides all indices of atoms within 6.56 A of this atom.
        6.56 is the upper bound of a possible neighbor 1.88 (C) + 1.4 + 1.4 + 1.88 (C).
        """
        neighbor_indices = mol.get("index", sel=f"protein and chain A\
        and noh and not resid '{self.resid}' and within 6.56 of index '{self.index}'")
        return neighbor_indices


def generate_sphere_points(coords: np.array, n: int = 610, radius: float = 1.88) -> np.array:
    """

    :param coords: The coordinates of the center of the atom
    :param n: number of points to sample
    :param radius: the radius of the atom
    :return: a nx3 vector with the point coordinates
    """
    total_radius = radius + water_radius
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = k * inc
        points.append([math.cos(phi) * r, y, math.sin(phi) * r])
    vec = np.asarray(points)
    vec *= total_radius
    vec += coords
    return vec


def retrieve_neighbor_positions(atom: Atom, mol: Molecule) -> Tuple[np.array, Dict[int, int]]:
    """
    :param atom: an Atom object
    :param mol: a htmd.htmd.molecule object
    :return: A tuple object with the positions of the neighboring atoms. A dictionary indexing column positions to resid positions
    """
    positions = mol.coords[atom.neighbor_indices][:, :, 0]
    position_index_to_resid = {index: mol.resid[neighbor_indice] for index, neighbor_indice in
                               enumerate(atom.neighbor_indices)}
    return positions, position_index_to_resid


def retrieve_indices(matrix: np.array, coords: np.array, neighborpositions: np.array, radius: float = 1.88) -> np.array:
    """

    Computes if each of the n sphere points are penetrating neighboring spheres
    :param matrix: n x m Distance matrix where n is the number of sphere points and m the number of neighbors
    :param coords: the coordinates of the atom
    :param neighborpositions: Coordinates of the neighbors
    :param radius: radius of the atom
    :return: The atoms that are in closest with each n sphere points.
    """

    ## When a row contains several atoms in contact,
    #  select the one with the closest center to the center of atom A
    sphere_radius = water_radius + radius
    #When a point is within the sphere of other atoms, the atom that is closest is chosen.
    # We need to compute the distances between all atom-neighbor pairs.
    dist_center_atoms = cdist(np.reshape(coords, (1, 3)), neighborpositions)
    ranking = np.argsort(dist_center_atoms)
    valid = matrix <= sphere_radius
    idx2 = []
    for row in valid:
        if row.any():
            # 1. ordering the row with Falses and Trues according to the ranking order
            # 2. Getting the first element that is true, which will be the first in the
            # raking list, thus the index of the closest atom.
            idx2.append(ranking[0][np.where(np.isin(ranking, np.where(row)))[1][0]])
    return idx2


def fill_matrices(atom: Atom, mol: Molecule, atom_matrix: np.array,
                  resid_matrix: np.array, indices: np.array, resids: np.array) -> Tuple[np.array, np.array]:
    """
    :param atom: An Atom class
    :param mol: an htmd.htmd.molecule object
    :param atom_matrix: the index x index area matrix
    :param resid_matrix: the ILVresid x ILVresid area matrix
    :param indices: the indices that belong to ILV sidechain heavy atoms
    :param resids: the resids that belong to ILV sidechain heavy atoms
    :return: Updated atom_matrix and resid_matrix
    """
    neighbor_positions, position_index_to_resid = retrieve_neighbor_positions(atom, mol)
    # Compute distances between all sphere points and the neighbors.
    # THe shape will be 610 (rows) x nr. of neighbors (columns).
    distances = cdist(atom.point_coords, neighbor_positions)
    column_indices = retrieve_indices(distances, atom.coords, neighbor_positions)
    colpos_occurrences = Counter(column_indices)
    for colpos, occurrences in colpos_occurrences.items():
        if atom.neighbor_indices[colpos] in indices:
            area = sphere_area_const * occurrences
            atom_matrix[atom.index, atom.neighbor_indices[colpos]] = area
            index_i = np.where(resids == atom.resid)[0][0]
            index_j = np.where(resids == position_index_to_resid[colpos])[0][0]
            resid_matrix[index_i, index_j] += area
    return atom_matrix, resid_matrix


def create_graph(resid_matrix: np.array, resid_list: np.array, cutoff_area: float = 10.0) -> Graph:
    """

    :param resid_matrix: the ILVresid x ILVresid area matrix
    :param resid_list: the index x index area matrix
    :return: A Graph object where each component is a ILV cluster
    """
    g = Graph()
    g.vp.resid = g.new_vertex_property("int")
    g.ep.area = g.new_edge_property("float")
    # 1. Create all vertices
    for v in range(len(resid_matrix)):
        v1 = g.add_vertex()
        g.vp.resid[v1] = resid_list[v]

    # 2. Create edges and fill its values with areas
    for row_index, row in enumerate(resid_matrix):
        v1 = g.vertex(row_index)
        for column_index, area in enumerate(row):
            v2 = g.vertex(column_index)
            if area > cutoff_area:
                ae = g.add_edge(v1, v2)
                g.ep.area[ae] = area
    return g


def filter_mol(inputmolfile: str) -> Molecule:
    """
    Loads, filters the object to chain A and writes filtered pdb out.
    :param inputmolfile: path to the pdb file
    :return: the Molecule object
    """
    mol = Molecule(inputmolfile)
    mol.filter("protein and chain A")
    mol.write(f"{inputmolfile[:-4]}-chainA.pdb")
    return mol


def write_clusters(g: Graph, components: PropertyArray, inputmolfile: str, outputname: str) -> None:
    """
    This function prints the clusters to the terminal and outputs them
    into a VMD session
    :param g: A Graph object
    :param outputname: The pdb filename.
    :param outputname: The file name to output the VMD session to.
    :return:
    """
    f = open(outputname[:-4] + ".txt", "w")
    mol = Molecule(f"{inputmolfile[:-4]}-chainA.pdb")
    mol.reps.add(sel='protein', style='NewCartoon', color=8)
    for cluster_index in range(max(components) + 1):
        cluster = [i for i, x in enumerate(components) if x == cluster_index]
        if len(cluster) < 2: continue
        vfilt = g.new_vertex_property('bool')
        for i in cluster: vfilt[i] = True
        sub = GraphView(g, vfilt)
        area = np.sum([g.ep.area[edge] for edge in sub.edges()])
        f.write(f"Cluster index {cluster_index}:\t"
                f"Residues {len(cluster)}. Total area {area} A^2.\t "
                f"Number of contacts: {sub.num_edges()}.\t "
                f"Contacts / Residue: {sub.num_edges() / len(cluster) }.\t "
                f"Area / Residue {area / sub.num_edges()}\n")  # sum reverse and inverse
        resid_cluster = [g.vp.resid[i] for i in cluster]
        mol.reps.add('chain A and noh and not backbone and resid %s' % ' '.join(map(str, resid_cluster)),
                     style="VDW", color=cluster_index)

    vmdobject = viewer(dispdev='text')
    vmdobject.loadMol(mol, name=inputmolfile)
    vmdobject.send("save_state " + outputname)
    vmdobject.close()
    f.close()


def add_clusters(mol, g: Graph, components: PropertyArray):
    """

    :param mol:
    :param g:
    :param components:
    :return: Molecule representations
            A list with Cluster objects
    """
    clusters = []
    mol.reps.add(sel='protein', style='NewCartoon', color=8)
    for cluster_index in range(max(components) + 1):
        cluster = [i for i, x in enumerate(components) if x == cluster_index]
        if len(cluster) < 2: continue
        vfilt = g.new_vertex_property('bool')
        for i in cluster: vfilt[i] = True
        sub = GraphView(g, vfilt)
        area = np.sum([g.ep.area[edge] for edge in sub.edges()])
        logger.info(f"Cluster index {cluster_index}:"
                    f"Residues {len(cluster)}. Total area {area} A^2. "
                    f"Number of contacts: {sub.num_edges()}. "
                    f"Contacts / Residue: {sub.num_edges() / len(cluster) }. "
                    f"Area / Residue {area / sub.num_edges()}")  # sum reverse and inverse
        resid_cluster = [g.vp.resid[i] for i in cluster]
        clusters.append(Cluster(
            area=area,
            residues=resid_cluster,
            contacts=sub.num_edges(),
            ratio_contacts_residue=sub.num_edges() / len(cluster),
            ratio_area_residue=area / sub.num_edges()
        ))
        mol.reps.add('chain A and noh and not backbone and resid %s' % ' '.join(map(str, resid_cluster)),
                     style="VDW", color=cluster_index)
    return clusters


def write_largest_cluster(g: Graph, inputmolfile: str, outputname: str) -> None:
    """
    This function prints the clusters to the terminal and outputs them
    into a VMD session
    :param g: A Graph object
    :param outputname: The pdb filename.
    :param outputname: The file name to output the VMD session to.
    :return:
    """
    mol = Molecule(f"{inputmolfile[:-4]}-chainA.pdb")
    mol.reps.add(sel='protein', style='NewCartoon', color=8)
    l = label_largest_component(g)
    sub = GraphView(g, vfilt=l)
    resid_cluster = [g.vp.resid[i] for i in sub.vertices()]
    mol.reps.add('chain A and noh and not backbone and resid %s' % ' '.join(map(str, resid_cluster)),
                 style="VDW", color=1)  # red
    vmdobject = viewer(dispdev='text')
    vmdobject.loadMol(mol, name=inputmolfile)
    vmdobject.send("save_state " + outputname)
    vmdobject.close()


def postprocess_session(inputmolfile: str, outputname: str) -> None:
    """
    Modifies the VMD session to not include tmp files
    :param outputname: The vmd session (output file)
    :param inputmolfile: Path to the pdb file already processed (filtered and or protonated)
    :return: None. Modifies the file inline
    """
    f = open(outputname, "r")
    lines = f.readlines()
    f.close()
    f = open(f"{outputname}", "w")
    for line in lines:
        if "mol new /tmp/" in line:
            line = f"mol new {inputmolfile}" \
                   f" type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
            f.write(line)
        elif "mol addfile /tmp/" in line:
            continue
        else:
            f.write(line)
    f.close()


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError("usage:"):
        print('hydroClusters.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    # Initialize molecules
    logger.info("Filtering and writing PDB")
    mol = filter_mol(inputfile)
    resids = np.unique(mol.get("resid", sel=sel))
    dims = len((resids))
    indices = mol.get("index", sel=sel)
    dims_indices = len(mol.get("index", sel="protein and chain A"))

    # Initializing Final output
    logger.info("Initializing final output")
    contacts = np.zeros((dims, dims))
    atoms_to_atoms = np.zeros((dims_indices, dims_indices))

    # Filling matrices
    logger.info("Computing clusters")
    for index in indices:
        a = Atom(index, mol)
        if not a.neighbor_indices.any(): continue
        _, contacts = fill_matrices(a, mol, atoms_to_atoms, contacts, indices, resids)

    # Rendering clusters
    graph = create_graph(contacts, resids)
    comp, _ = label_components(graph)
    if comp.a.any():
        write_clusters(graph, comp.a, inputfile, outputfile)
        inputfile_processed = f"{inputfile[:-4]}-chainA.pdb"
        postprocess_session(inputfile_processed, outputfile)
        summary_output = outputfile[:-4] + "-major.vmd"
        write_largest_cluster(graph, inputfile, summary_output)
        postprocess_session(inputfile, summary_output)
        logger.info("Saving VMD sessions")
    else:
        logging.warning("The ILV residues are not in contact")


if __name__ == "__main__":
    main(sys.argv[1:])
