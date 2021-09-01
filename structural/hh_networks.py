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
    :return: the Molecule object
    """
    mol1 = Molecule(inputmolfile)
    mol1.filter("protein and chain A")
    mol = proteinPrepare(mol1)
    mol.write(f"{inputmolfile[:-4]}-chainA_protonated.pdb")
    return mol


def make_graph_hh(hbonds: np.array, trajectory: md.Trajectory, resi_to_chain: dict) -> Graph:
    """
    This function creates a graph with the hbond networks. Each component in the network will be a
    separate network in the final analysis.
    :param hbonds: a np.array indicating indexes of atoms that have an h-bond
    :param trajectory: an mdtraj trajectory. It can also be a pdb file processed with mdtraj as trajectory
    :return: A graph-tool graph object
    """
    
    g = Graph(directed=False)
    g.vp.resid = g.new_vertex_property("int")
    g.vp.chain = g.new_vertex_property("string")
    g.vp.atom = g.new_vertex_property("vector<int>")

    for i in hbonds:
    
        resid_d = trajectory.topology.atom(i[0]).residue.resSeq
        resid_a = trajectory.topology.atom(i[2]).residue.resSeq

        chain_d = resi_to_chain[resid_d]
        chain_a = resi_to_chain[resid_a]
        
        donors = find_vertex(g, g.vp.resid, resid_d)
        acceptors = find_vertex(g, g.vp.resid, resid_a)

        if donors and g.vp.chain[donors[0]] == chain_d:
            v1 = donors[0]
            g.vp.atom[v1].append(i[0])
        else:
            v1 = g.add_vertex()
            g.vp.resid[v1] = resid_d
            g.vp.chain[v1] = chain_d
            g.vp.atom[v1] = []
            g.vp.atom[v1].append(i[0])
        
        if acceptors and g.vp.chain[acceptors[0]] == str(chain_a):
            v2 = acceptors[0]
            g.vp.atom[v2].append(i[2])
        else:
            v2 = g.add_vertex()
            g.vp.resid[v2] = resid_a
            g.vp.chain[v2] = chain_a
            g.vp.atom[v2] = []
            g.vp.atom[v2].append(i[2])

        g.add_edge(v1, v2) 
    return g

def write_networks(g: Graph, components: PropertyArray, inputmolfile: str, outputname: str) -> None:
    """
    This function outputs the HH networks into a VMD session
    :param g: A Graph object
    :param outputname: The pdb filename.
    :param outputname: The file name to output the VMD session to.
    :return:
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
    bonds = []
    mol.reps.add(sel='protein', style='NewCartoon', color=8)
    
    for network_index in range(max(components) + 1):
        network = [i for i, x in enumerate(components) if x == network_index]  # these are the vertices
        resid_network = [g.vp.resid[i] for i in network]  # these are the resids
        chain_network = [g.vp.chain[i] for i in network]
        atoms_network = [[i for i in g.vp.atom[j]] for j in network]
        atoms_network = [atom for sublist in atoms_network for atom in sublist]
        bonds.append(resid_network)
        vfilt = g.new_vertex_property('bool')
        for i in network: vfilt[i] = True
        sub = GraphView(g, vfilt)

        mol.reps.add('noh and resid %s' % ' '.join(map(str, resid_network)),
                     style="Licorice", color=network_index)
        mol.reps.add('noh and resid %s' % ' '.join(map(str, resid_network)),
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
