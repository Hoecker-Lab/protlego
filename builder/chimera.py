from protlego.definitions import logger
from moleculekit.projections.metricdistance import MetricDistance
from protlego.structural.clusters import *
from protlego.structural.hh_networks import *
import os
from moleculekit.support import string_to_tempfile
import urllib.request

def get_SCOP_domain(domain):
    """
    :param domain: str. The SCOPe domain to download as pdb
    :return:
    """

    logger.info(f'Attempting to download domain {domain} from the SCOP server')
    url = f'https://scop.berkeley.edu/astral/pdbstyle/ver=2.07&id={domain}&output=text'
    connected = False
    while not connected:
        try:
            response = urllib.request.urlopen(url)
            text = response.read()
        except Exception as e:
            import time
            logger.warning(f'Failed to connect to SCOP with error {e}. Sleeping 5s and retrying.')
            time.sleep(5)
            continue
        connected = True
    filepath = string_to_tempfile(text.decode('ascii'), 'pdb')
    logger.info(f"File downloaded as {filepath}")

    return filepath

def get_FUZZLE_hhs(domain):
    """
    :param domain: str. The domain to download from Fuzzle as hhs

    :return: filepath: path where the file is located.
    """
    logger.info(f'Attempting to download hhs file for {domain} from the FUZZLE server')
    url = f'https://fuzzle.uni-bayreuth.de/hhs/scop95_2.07.psi.hhs/{domain}.hhs'
    connected = False
    while not connected:
        try:
            response = urllib.request.urlopen(url)
            text = response.read()
        except Exception as e:
            import time
            logger.warning(f'Failed to connect to FUZZLE with error {e}. Sleeping 5s and retrying.')
            time.sleep(5)
            continue
        connected = True
    filepath = string_to_tempfile(text.decode('ascii'), 'hhs')
    logger.info(f"File downloaded as {filepath}")	

    return filepath
class Chimera(Molecule):

    def add_crossover(self, crossover: int):
        start = self.resid[0]
        end = self.resid[-1]
        self.reps.add(sel=f"protein and resid '{start}' to '{crossover}'", style='NewCartoon', color=1)
        self.reps.add(sel=f"protein and resid '{crossover+1}' to '{end}'", style='NewCartoon', color=7)

    def __str__(self):  # TODO fix the naming of parents.
        return str(f"Protein with:"
                   f"\nNumber of residues: {self.numResidues}"
                   f"\nNumber of atoms: {self.numAtoms}")

    @property
    def sasa(self):
        return

    def calc_resid_dist(self, resid1: int, resid2: int) -> np.array:
        """
        Returns the distance between the alpha carbons of two residues"
        :param resid1:index of residue 1 in the PDB structure
        :param resid2:index of residue 2 in the PDB structure
        :return:
        """
        coord1 = self.get("coords", sel=f'resid {resid1} and name CA')
        coord2 = self.get("coords", sel=f'resid {resid2} and name CA')
        return cdist(coord1, coord2)

    def compute_hydrophobic_clusters(self, sel: str = "protein and not backbone and noh and resname ILE VAL LEU",
                                     cutoff_area: float = 10):
        """
        :param
            sel: VMD selection on which to compute the clusters.
            Default is every sidechain heavy atom ILE, VAL and LEU
            residues. "protein and not backbone and noh and resname ILE VAL LEU"
            cutoff_area: Minimum area between residues to be considered in contact.
	:return: A representation for each cluster
        """
        clusters = None

        # Removing previous visualizations
        [self.reps.remove(index) for index, rep in reversed(list(enumerate(self.reps.replist)))]

        resids = self.get("resid", sel=f"{sel}") # ILV cluster residues
        chains = self.get("chain", sel=f"{sel}")  # chains of ILV cluster residues 
        dims = len(resids) # length of iLV clusters
        indices = self.get("index", sel=f"{sel}") # the indices of the atoms in the ILV clusters

        # get a dictionary of atom index and resid position in resids list
        atom_to_residposition = {}
        for index in indices:
                resid = self.get("resid", sel=f"index {index}")[0]
                chain = self.get("chain", sel=f"index {index}")[0]
                index_residue = [j for j, residue in enumerate(resids) if (residue == resid and chains[j] == chain) ][0]
                atom_to_residposition[index] = index_residue


        logger.info("Initializing final output")
        contacts = np.zeros((dims, dims))
        
        logger.info("Computing clusters")
        for index in indices:
            a = Atom(index, self)
            if not a.neighbor_indices.any():
                continue
            contacts = fill_matrices(a, self, contacts, indices, atom_to_residposition)

        graph = create_graph(contacts, resids, chains, cutoff_area=cutoff_area)
        comp, _ = label_components(graph)
        if comp.a.any():
            clusters = add_clusters(self, graph, comp)
        else:
            logger.warning("There are not residues in contact for this selection")

        return clusters

    def compute_hydrogen_networks(self, sidechain_only: bool = True):
        """
        Analyzes the hydrogen bonds following the Baker Hubbard algorithm and
        clusters them. (Two residues are in the same cluster if there's a path
        of hydrogen networks between them). This function adds protons to the protein.
        between them
        :param sidechain_only: The whole residue or sidechain only. Note: Computing networks
        including backbone leads to very large clusters harder to visualize.
        :return: A representation for each network
        """
        # Remove previous structural representations
        bonds = None
        [self.reps.remove(index) for index, rep in reversed(list(enumerate(self.reps.replist)))]
        # [self.reps.remove(index) for index, rep in reversed(list(enumerate(self.reps.replist))) if
        #  rep.style == 'VDW']
        
        newmol = proteinPrepare(self)
        newmol.set("resname","HIS", "resname HID")
        newmol.set("resname", "HIS", "resname HIE")
        newmol.set("resname", "HIS", "resname HIP")

        newmol.write(f"/tmp/structure.pdb")
        
        resi_to_chain = dict(zip(self.resid,self.chain))

        bonds = ''
        t = md.load(f"/tmp/structure.pdb")
        os.remove(f"/tmp/structure.pdb")
        
        hbonds = md.baker_hubbard(t, sidechain_only=sidechain_only)
        graph = make_graph_hh(hbonds, t, resi_to_chain)
        comp, _ = label_components(graph)
        if comp.a.any():
            bonds = add_networks(self, graph, comp)
        else:
            logging.warning("No Hydrogen Bonds found")

        return bonds

    def compute_salt_bridges(self):

        salts = []
        [self.reps.remove(index) for index, rep in reversed(list(enumerate(self.reps.replist)))]

        metr = MetricDistance('sidechain and acidic and element O',
                              'sidechain and basic and element N', metric="contacts", threshold=3.2,
                              pbc=False)
        try:
            data = metr.project(self)
            mapping = metr.getMapping(self)
            
            if len(np.shape(data)) > 1: data = data[0].copy()  # handling NMR structures
            
            self.reps.add(sel='protein', style='NewCartoon', color=8)
            
            if mapping[data].atomIndexes.values.any():
                for salt in mapping[data].atomIndexes.values:
                    resid1 = self.get("resid", sel=f"same residue as index {salt[0]}")[0]
                    chain1 = self.get("chain", sel=f"same residue as index {salt[0]}")[0]
                    resid2 = self.get("resid", sel=f"same residue as index {salt[1]}")[0]
                    chain2 = self.get("chain", sel=f"same residue as index {salt[1]}")[0]
                
                    if [resid1, resid2] not in salts:
                        salts.append({"residues": [int(resid1), int(resid2)], "chain": [chain1, chain2]})
                        self.reps.add(f"protein and resid {resid1}", style="Licorice", color="1")
                        self.reps.add(f"protein and resid {resid2}", style="Licorice", color="0")
        except:
            logger.error("Molecule has no basic or acidic residues")
            raise
         
        graph = make_graph_salts(salts)
        comp, _ = label_components(graph)
        if comp.a.size != 0:
             salts = add_networks_salts(graph, comp)
        else:
            logger.warning('No salt bridges present in the structure')        
        return salts
