from moleculekit.molecule import Molecule
from protlego.builder.chimera import Chimera
from subprocess import call
import numpy as np
import os

import logging
logging.getLogger('moleculekit').setLevel(logging.ERROR)
logger = logging.getLogger('protlego')

def rosettawrapper(runtype: str, relaxed_structures=None, output: str = "/tmp/build",chimera: Chimera = None,pdb: str = None):
    """

    :param runtype: takes score, minimize or relax as runtype
    :param relaxed_structures: set the number of relaxed structures for the runtype relax
    :param output: takes an output path (standard path is /tmp/build)
    :param chimera: takes a chimera object as inpuit structure
    :param pdb: takes the path of the pdb input structure
    :return:
    """

    if not chimera and not pdb:
        raise Exception('No chimera or pdb object as input given !')

    if chimera:
        chimera.write(output+"/structure.pdb")
        pdb = output+"/structure.pdb"
    if pdb:
        mol = Molecule(pdb)
        if mol.numFrames > 1:
            mol.dropFrames(keep=0)
            logger.info("Protein contains more than one model. Keeping only the first one")
            mol.write(output+'/structure.pdb')
            pdb = output+"/structure.pdb"

    if runtype == 'relax':
        call(['/agh/projects/software/rosetta/3.10/main/source/bin/relax.static.linuxgccrelease',
              f'-s {pdb}',
              '-out:suffix _relaxed',
              f'-nstruct {relaxed_structures}',
              '-ignore_zero_occupancy false',
              '-relax:default_repeats 5',
              f'-out:path:pdb {output}',
              f'-out:path:score {output}',
              '-in:ignore_waters'
              ])
    elif runtype == 'minimize':
        call(['/agh/projects/software/rosetta/3.10/main/source/bin/minimize.static.linuxgccrelease',
              f'-in:file:s {pdb}',
              '-run:min_type lbfgs_armijo_nonmonotone',
              '-run:min_tolerance 0.001',
              '-out:suffix _minimized',
              f'-out:path:score {output}',
              f'-out:path:pdb {output}'])
    elif runtype == 'score':
        call(['/agh/projects/software/rosetta/3.10/main/source/bin/score_jd2.static.linuxgccrelease',
              f'-in:file:s {pdb}',
              '-out:suffix _scored',
              f'-out:path:score {output}',
              f'-out:path:pdb {output}'])
    elif runtype == 'repack':
        resfile = output+"/resfile.res"
        counter = 1
        while os.path.isfile(resfile):
            resfile = output+f"/resfile{counter}.res"
            counter += 1
        with open(resfile, "w") as rf:
            rf.write("NATAA\nstart")
        call(['/agh/projects/software/rosetta/3.10/main/source/bin/fixbb.static.linuxgccrelease',
              f'-in:file:s {pdb}',
              '-out:suffix _repacked',
              f'-out:path:score {output}',
              f'-out:path:pdb {output}',
              '-minimize_sidechains',
              f'-resfile {resfile}'])
        if os.path.isfile(resfile):
            os.remove(resfile)

    aa = Molecule(pdb).numResidues

    if runtype == 'relax':
        file = open(f"{output}/score_relaxed.sc", "r")
    elif runtype == 'minimize':
        file = open(f"{output}/score_minimized.sc", "r")
    elif runtype == 'score':
        file = open(f"{output}/score_scored.sc", "r")
    elif runtype == 'repack':
        file = open(f"{output}/score_repacked.sc", "r")
    else:
        raise Exception('No corresponding score file could be found !')

    energies = []
    for i_line, line in enumerate(file):
        if i_line >=2:
            columns = line.split()
            energies.append(float(columns[1]))
    if runtype in ('score', 'minimize'):
        print(str(energies[0] / aa))

    elif runtype in ('repack', 'relax'):
        mean = np.mean(energies) / aa
        std = np.std(energies) / aa
        print ('Mean-value:' + str(mean) + ' Standard-deviation:' + str(std))
