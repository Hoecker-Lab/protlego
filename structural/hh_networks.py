from moleculekit.molecule import Molecule
from moleculekit.tools.preparation import proteinPrepare
from moleculekit.vmdviewer import viewer
import mdtraj as md
from graph_tool.all import *
import numpy as np
import sys, getopt
from protlego.structural.clusters import postprocess_session
import logging

logging.getLogger('moleculekit').setLevel(logging.ERROR)
logger = logging.getLogger('protlego')


def protonate_mol(inputmolfile: str) -> Molecule:
    """
    Loads, filters the object to chain A and protonates the selection.
    It then writes the filtered protonated pdb out.
    :param inputmolfile: path to the pdb file
    :return: a Molecule object with the correct protonation state
    """
    mol1 = Molecule(inputmolfile)
    mol1.filter("protein and chain A")
    mol = proteinPrepare(mol1)
    mol.write(f"{inputmolfile[:-4]}-chainA_protonated.pdb")
    return mol


def make_graph_hh(hbonds: np.array, trajectory: md.Trajectory) -> Graph:
    """

    :param hbonds: a np.array indicating indexes of atoms that have an h-bond
    :param trajectory: an mdtraj trajectory. It can also be a pdb file processed with mdtraj as trajectory
    :return: A graph-tool graph object
    """
    g = Graph(directed=False)
    g.vp.resid = g.new_vertex_property("int")
    g.vp.atom = g.new_vertex_property("int")
    for i in hbonds:
        resid_d = trajectory.topology.atom(i[0]).residue.resSeq
        resid_a = trajectory.topology.atom(i[2]).residue.resSeq
        donors = find_vertex(g, g.vp.resid, resid_d)
        acceptors = find_vertex(g, g.vp.resid, resid_a)
        if donors:
            v1 = donors[0]
        else:
            v1 = g.add_vertex()
        if acceptors:
            v2 = acceptors[0]
        else:
            v2 = g.add_vertex()
        g.vp.resid[v1] = resid_d
        g.vp.resid[v2] = resid_a
        g.vp.atom[v1] = i[0]
        g.vp.atom[v2] = i[2]
        g.add_edge(v1, v2)
    return g


def write_networks(g: Graph, components: PropertyArray, inputmolfile: str, outputname: str) -> None:
    """
    :param g:  A graph-tool graph object
    :param components: the components of the graph
    :param inputmolfile: the pdb of the input protein
    :param outputname: a path wehre to write the output
    :return: Writes a file named outputname + "hh-networks.txt"
    """
    f = open(outputname[:-4] + "hh-networks.txt", "w")
    mol = Molecule(inputmolfile)
    mol.reps.add(sel='protein', style='NewCartoon', color=8)
    f.write(f"Atom index\tAtom index\tresid\tresid\n")
    for network_index in range(max(components) + 1):
        f.write(f"Network index {network_index}\n")

        network = [i for i, x in enumerate(components) if x == network_index]  # these are the vertices
        resid_network = [g.vp.resid[i] for i in network]  # these are the resids
        vfilt = g.new_vertex_property('bool')
        for i in network: vfilt[i] = True
        sub = GraphView(g, vfilt)
        # print for each edge the atoms and resid of the two pairs
        for edge in sub.edges():
            f.write(f"{sub.vp.atom[edge.source()]}\t"
                    f"{sub.vp.atom[edge.target()]}\t"
                    f"{sub.vp.resid[edge.source()]}\t"
                    f"{sub.vp.resid[edge.target()]}\n")

        mol.reps.add('chain A and noh and resid %s' % ' '.join(map(str, resid_network)),
                     style="Licorice", color=network_index)
        mol.reps.add('chain A and noh and resid %s' % ' '.join(map(str, resid_network)),
                     style="HBonds", color=network_index)

    vmdobject = viewer(dispdev='text')
    vmdobject.loadMol(mol, name=inputmolfile)
    vmdobject.send("save_state " + outputname)
    vmdobject.close()
    f.close()


def add_networks(mol, g: Graph, components: PropertyArray):
    """

    :param mol: the Chimera where to add the representations to.
    :param g: a graph-tool graph object with the networks property
    :param components: the components from the graph
    :return: INcludes the representations in the Chimera object.
    """
    bonds = []
    mol.reps.add(sel='protein', style='NewCartoon', color=8)
    for network_index in range(max(components) + 1):
        network = [i for i, x in enumerate(components) if x == network_index]  # these are the vertices
        resid_network = [g.vp.resid[i] for i in network]  # these are the resids
        bonds.append(resid_network)
        vfilt = g.new_vertex_property('bool')
        for i in network: vfilt[i] = True
        sub = GraphView(g, vfilt)

        mol.reps.add('chain A and noh and resid %s' % ' '.join(map(str, resid_network)),
                     style="Licorice", color=network_index)
        mol.reps.add('chain A and noh and resid %s' % ' '.join(map(str, resid_network)),
                     style="HBonds", color=network_index)
    return bonds


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError("usage:"):
        print('hh_newtorks.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('hh_newtorks.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    # 1. Protonate molecule
    _ = protonate_mol(inputfile)
    inputfile_protonated = f"{inputfile[:-4]}-chainA_protonated.pdb"

    t = md.load(inputfile_protonated)
    hbonds = md.baker_hubbard(t, sidechain_only=True)
    graph = make_graph_hh(hbonds, t)
    comp, _ = label_components(graph)
    if comp.a.any():
        write_networks(graph, comp, inputfile_protonated, outputfile)
        postprocess_session(inputfile_protonated, outputfile)
        logger.info("Saving VMD sessions")
    else:
        logger.warning("No Hydrogen Bonds found")


if __name__ == "__main__":
    main(sys.argv[1:])
