from moleculekit.molecule import Molecule
from moleculekit.projections.metricdistance import MetricDistance
from pandas import DataFrame
from moleculekit.vmdviewer import viewer
import numpy as np
import sys, getopt
from protlego.structural.clusters import filter_mol, postprocess_session


import logging
logging.getLogger('moleculekit').setLevel(logging.ERROR)
logger = logging.getLogger('salt_bridges')



def write_salt_bridges(data: np.ndarray, mapping:DataFrame, mol: Molecule, outputname: str) -> None:
    """
    This function outputs the salt bridges into a VMD session
    :param data: A MetricDistance object
    :param mapping: A DataFrame object including the index - residue mapping
    :param mol: The pdb filename.
    :param outputname: The file prefix name to output the VMD session to. Example: "protein2"
    :return: A file with the VMD session named outputname+"-salt-bridges.txt"
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