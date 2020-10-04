"""
Pipeline to put the ligands on the SCOPe domains.
It maps every pdb to every domain, and then puts the ligands that present an atom at least 4 armstrongs
close to the domain.
"""
import os
from moleculekit.molecule import Molecule
import logging
import moleculekit

logging.getLogger(moleculekit.molecule.__name__).setLevel(logging.ERROR)

file = open("./dir.cla.scope.2.07-stable.txt", "r")


def handle_residues(resid):
    if len(resid.split('-')) == 3:
        resid_sel1 = '-' + resid.split('-')[1]
        resid_sel2 = resid.split('-')[2]
    elif len(resid.split('-')) == 4:
        resid_sel1 = '-' + resid.split('-')[1]
        resid_sel2 = '-' + resid.split('-')[3]
    else:
        raise ValueError(f"{resid}")
    return resid_sel1, resid_sel2


for i in range(4):
    file.readline()

for index, line in enumerate(file):

    columns = line.split()
    domain = columns[0]
    folder = domain[2:4]
    pdb = columns[1]

    output_file = f"/agh/db/scop/2.07/pdbstyle-2.07-ligands/{folder}/{domain}.pdb"
    if os.path.exists(output_file): continue
    pdb_file = f"/agh/db/pdb/data/{folder}/pdb{pdb}.ent.gz"
    if not os.path.exists(pdb_file):
        print(f"{pdb_file} not present")
        continue

    chains = columns[2].split(',')
    if chains == ['-']: continue
    prot_sel = ''
    for i, chain in enumerate(chains):
        if i > 0: prot_sel += " or "
        resids = chain.split(':')
        if not resids[1]:
            prot_sel += f"protein and chain {resids[0]}"
        else:
            resid_sel1, resid_sel2 = resids[1].split('-') if len(resids[1].split('-')) == 2 else handle_residues(
                resids[1])
            prot_sel += f"protein and chain {resids[0]} and resid '{resid_sel1}' to '{resid_sel2}'"
    selection = f"{prot_sel} or (not protein and same residue as within 4 of ({prot_sel}))"
    mol_pdb = Molecule(pdb_file, validateElements=False)
    mol_pdb.filter(selection)
    output_folder = f"/agh/db/scop/2.07/pdbstyle-2.07-ligands/{folder}"
    if not os.path.exists(output_folder): os.mkdir(output_folder)
    mol_pdb.write(output_file)
