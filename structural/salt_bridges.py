from moleculekit.molecule import Molecule
from moleculekit.projections.metricdistance import MetricDistance
from pandas import DataFrame
from moleculekit.vmdviewer import viewer
import numpy as np
import sys, getopt
from protlego.structural.clusters import filter_mol, postprocess_session
from graph_tool.all import *

import logging
logging.getLogger('moleculekit').setLevel(logging.ERROR)
logger = logging.getLogger('salt_bridges')



def write_salt_bridges(data: np.ndarray, mapping:DataFrame, mol: Molecule, outputname: str) -> None:
    """
    This function outputs the HH networks into a VMD session
    :param g: A Graph object
    :param outputname: The pdb filename.
    :param outputname: The file name to output the VMD session to.
    :return:
    """
    f=open(outputname[:-4]+"-salt-bridges.txt","w")
    mol.reps.add(sel='protein' , style='NewCartoon', color=8)
    f.write(f"resid\tresid\n")
    for bond in mapping[data].atomIndexes.values:
        mol.reps.add(f"chain A and same residue as index {bond[0]}", style="Licorice", color="1")
        mol.reps.add(f"chain A and same residue as index {bond[1]}", style="Licorice", color="0")
        f.write(f"{bond[0]}\t"
                f"{bond[1]}\n")

    vmdobject = viewer(dispdev='text')
    vmdobject.loadMol(mol)
    vmdobject.send("save_state " + outputname)
    vmdobject.close()
    f.close()

def make_graph_salts(salts) -> Graph:
    """
    :param hbonds: a np.array indicating indexes of atoms that have an h-bond
    :param trajectory: an mdtraj trajectory. It can also be a pdb file processed with mdtraj as trajectory
    :return: A graph-tool graph object
    """
    g = Graph(directed=False)
    g.vp.resid = g.new_vertex_property("int")
    g.vp.chain = g.new_vertex_property("string")
    for salt in salts:
        donors = find_vertex(g, g.vp.resid, salt["residues"][0])
        acceptors = find_vertex(g, g.vp.resid, salt["residues"][1])

        if donors and g.vp.chain[donors[0]] == salt["chain"][0]:
            v1 = donors[0]
        else:
            v1 = g.add_vertex()
            g.vp.resid[v1] = salt["residues"][0]
            g.vp.chain[v1] = salt["chain"][0]

        if acceptors and g.vp.chain[acceptors[0]] == salt["chain"][1]:
            v2 = acceptors[0]
        else:
            v2 = g.add_vertex()
            g.vp.resid[v2] = salt["residues"][1]
            g.vp.chain[v2] = salt["chain"][1]

        g.add_edge(v1, v2)
    return g

def add_networks_salts(g, components):
    """
    :param g: a graph-tool graph object with the networks property
    :param components: the components from the graph
    :return: Includes the representations in the Chimera object.
    """
    salts = []

    for network_index in range(max(components) + 1):
        network = [i for i, x in enumerate(components) if x == network_index]  # these are the vertices
        resid_network = [g.vp.resid[i] for i in network]  # these are the resids
        chain_network = [g.vp.chain[i] for i in network]  # these are the chains

        salts.append({"residues": resid_network, "chain": chain_network})

    return salts


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError("usage:"):
        print ('salt_bridges.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('salt_bridges.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    #1. Load molecule
    logger.info("Filtering and writing PDB")
    mol = filter_mol(inputfile)

    #2. Compute distances
    logger.info("Computing distances among all polar residues")
    metr = MetricDistance('chain A and sidechain and acidic and element O',
                          'chain A and sidechain and basic and element N', metric="contacts", threshold=3.2, pbc=False)
    try:
        data = metr.project(mol)
    except:
        logger.error("Molecule has no basic or acidic residues")
        raise

    if len(np.shape(data)) > 1: data = data[0].copy()  # handling NMR structures
    mapping = metr.getMapping(mol)

    #3. Write txt and vmd session out
    write_salt_bridges(data, mapping, mol, outputfile)
    inputfile_processed = f"{inputfile[:-4]}-chainA.pdb"
    postprocess_session(inputfile_processed, outputfile)
    logger.info("Saving VMD session")

if __name__ == "__main__":
    main(sys.argv[1:])
