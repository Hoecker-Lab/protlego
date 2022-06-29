import subprocess
from subprocess import CalledProcessError
from typing import Tuple
import numpy as np

from protlego.definitions import logger, TM_BIN

def get_tmalign_output(mobile: str, target:str, matrix_filename:str) -> np.ndarray:
    """Reads the output matrix of TM align

    :param mobile: path to the query pdb file
    :param target: path to the subject pdb file
    """
    try:
        subprocess.check_output([TM_BIN, mobile, target, '-m', matrix_filename, '-o', "TM.sup"])
    except Exception as e:
        logger.error(f"TMalign cannot align the molecules. Error follows {e}")
    matrot = [[None] * 4] * 3
    with open(matrix_filename, "r") as inputfile:
        for i, line in enumerate(inputfile.readlines()[2:5]):
            matrot[i] = [float(n) for n in line.strip().split()[1:]]

    # os.remove(matrix_filename)
    return np.array(matrot)


def tm2pymol(matrot: np.array) -> np.ndarray:
    """ This function takes a TM matrix and returns it in the correct format for Pymol

    :param matrot: A TM matrix (from get_tmalign_output)
    :return: pymolmat: Rotation matrix in the pymol format
    """
    pymolmat = [0] * 15 + [1]
    for index, i in enumerate(matrot):
        numbers = [float(n) for n in i]
        pymolmat[index * 4 + 3], pymolmat[index * 4 + 0], pymolmat[index * 4 + 1], pymolmat[index * 4 + 2] = numbers

    return np.array(pymolmat)


def tm2vmd(matrot: np.array) -> Tuple[np.array, np.array]:
    """ This function takes a TM matrix and returns it in the correct format for VMD
    :param matrot: A TM matrix (from get_tmalign_output)
    :return:
    vmdvec: a Rotation vector for VMD
    vmdmat: a translation matrix for VMD
    """
    vmdvec = []
    vmdmat = [[None] * 3] * 3
    for index, i in enumerate(matrot):
        vmdvec.append(i[0])
        vmdmat[index] = [float(n) for n in i[1:]]

    return np.array(vmdvec), np.array(vmdmat)

def get_new_resIDs(mol,index):
    """
    Produces an array which contains the new resid ID for each atom in the molecule
    :param mol: Molecule object to renumber
    :return: The residue index from which to start
    """
    length = len(mol.resid)
    c = int(index)
    seq = np.zeros(length, dtype=int)
    seq[0] = c
    for i in range(1, length):
        if mol.resid[i - 1] != mol.resid[i]:
            c += 1  # new sequence id
        seq[i] = c
    return seq
