from difflib import SequenceMatcher
from typing import Dict
import os
from protlego.builder.alignments import *
from protlego.database.data import *
from protlego.builder.AAcodes import *
from protlego.builder.chimera import get_SCOP_domain, get_FUZZLE_hhs
from protlego.builder.chimera import Chimera

from csb.bio.io.hhpred import HHOutputParser
from csb.bio.hmm import HHpredHitAlignment
from moleculekit.molecule import Molecule
from moleculekit.projections.metricrmsd import MetricRmsd
from moleculekit.projections.metricsecondarystructure import MetricSecondaryStructure
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from protlego.definitions import TM_BIN, logger

aa_keys = special2standard.keys()
standard_aas = three2one.keys()

class NotDiverseChimeraError(Exception):
    pass


class BackboneClashError(Exception):
    pass


class NotCorrectPDBError(Exception):
    pass

class Builder():
    """
    Graph constructor and visualizer

    Examples
    --------
    myhit= fetch_id('10002347')
    a=Builder(myhit)
    aln=a.get_alignment(myhit.query,myhit.no)
    qPDB, sPDB = a.superimpose_structures(aln,partial_alignment=True)        
    chimeras=a.build_chimeras(partial_alignment=True)
    """
    def __init__(self, hit: Hit):
        qpdb_path = get_SCOP_domain(hit.query)
        spdb_path = get_SCOP_domain(hit.sbjct)
        logger.info(f'Loading {qpdb_path} as a chimera object') 
        self.qPDB = Chimera(qpdb_path, validateElements=False)
        os.remove(qpdb_path)
        if self.qPDB.numFrames > 1:
            self.qPDB.dropFrames(keep=0)
            logger.info("Query protein contains more than one model. Keeping only the first one")
        logger.info(f'Loading {spdb_path} as a chimera object')
        self.sPDB = Chimera(spdb_path, validateElements=False)
        os.remove(spdb_path)
        if self.sPDB.numFrames > 1:
            self.sPDB.dropFrames(keep=0)
            logger.info("Subject protein contains more than one model. Keeping only the first one")
        self.qaPDB, self.saPDB = {}, {}
        self.qpairs,self.spairs = [], []
        self.dst = []
        self.chim_positions = {}

    def get_alignment(self, query: str, no: str) -> HHpredHitAlignment:
        """ Obtain the HHS alignment 'no' for query 'query'.
        Only the alignment from the fragment region is retrieved.
        This implies that when the fragment is not located in the
        N-terminus hit.q_start and the position in the output won't
        be the same. For example, if q_start = 20, that aminoacid is
        in position 0 in aln.query.

        :param query: str. Domain query
        :param no: int. Specifies the position in the file (alignment with subject)

        :return: HHpredHitAlignment. Alignment between query and subject for the fragment region.
        """

        hhF = get_FUZZLE_hhs(query)
        try:
            hh = HHOutputParser().parse_file(hhF)
            pair = hh[int(no) - 1]
            aln = pair.alignment
            return aln
        except Exception as e:
            logger.error(f"Parsing of {hhF} failed. Error follows: {e}")

    @staticmethod
    def remove_residue(pdb: Molecule, resid):
        """
        Removes a specific residue by its name
        :param pdb: The pdb to mutate the residues to
        :param resid:
        :return: The pdb with the filtered residues
        """

        pdb.filter(f"not {resid}")
        return pdb

    @staticmethod
    def find_nonstandards(pdb: Molecule) -> list:
        """
        Finds non-standard aminoacids
        :param pdb: Molecule or Chimera object where to find non-standard aminoacids.
        :return: list of non-standard residues
        """
        non_standards = [aa for aa in np.unique(pdb.resname) if (aa in aa_keys or aa not in standard_aas)]
        if non_standards:
            for i in non_standards:
                if i != 'UNK':
                    logger.info(f"Found the following non-standar residue: {i}. "
                                f"Preserving in the original PDB")
                else:
                    logger.warning("Protein presents unknown residue UNK."
                                   " Call remove_residue() to remove it or provide parameters"
                                   " if you want to minimize it with AMBER or CHARMM.")
        return non_standards

    @staticmethod
    def mutate_nonstandards(pdb: Molecule) -> Molecule:
        """

        :param pdb: The pdb to mutate the residues to
        :return: The pdb with the mutated residues -if there are any-
        """
        non_standards = Builder.find_nonstandards(pdb)
        if non_standards:
            [pdb.mutateResidue(f"resname {i}", f"{special2standard[i]}") for i in non_standards]

        return pdb

    def seq2pdbaln(self, sequence: str, pdb: Molecule) -> List:
        """ Obtains the resid (positions in the PDB) of the aminoacids involved in the fragment.

        :param sequence: str. sequence of the fragment
        :param pdb: Molecule. PDB from where to obtain the residues
        :return: mapping. The the PDB Positions of the fragment. A List of tuples of the form: (sequence aminoacid, pdb resid)
        """

        # Mutate non standard residues
        copy_pdb = pdb.copy()  # making a copy to ensure the possible mutations don't affect the pdb
        copy_pdb = self.mutate_nonstandards(copy_pdb)

        # Mapping of residue to number
        try:
            pdb_sequence = copy_pdb.sequence()['0']
        except:
            raise NotCorrectPDBError(f"PDB structure from protein {pdb.viewname} not correct")

        copy_pdb = pdb.copy()  # second copy to preserve indexing
        pdb_indices = copy_pdb.get("index", sel="protein and name CA")

        # sequence from the HHS alignment (fragment)
        seq_positions = [(x,) for x in sequence]

        # obtain mapping
        matcher = SequenceMatcher(None, sequence, pdb_sequence, autojunk=False)
        for block in matcher.get_matching_blocks():
            i = 0
            for pos in range(block.a, block.a + block.size):
                seq_positions[pos] += (pdb_indices[block.b + i],)
                i += 1
        return seq_positions

    def _get_pairs(self, aln: HHpredHitAlignment):
        """ Obtain the positions of the pdb which are present in both alignments.
        Thus the number of atoms is the same

        :param aln: HHpredHitAlignment. The fragment sequence alignment
        """

        qpairs = []
        spairs = []
        preqpairs = []
        prespairs = []

        # Obtain the corresponding PDB residues to the fragment seq. alignment
        query_seq2pdb = self.seq2pdbaln(aln.query, self.qPDB)
        sbjct_seq2pdb = self.seq2pdbaln(aln.subject, self.sPDB)
        
        # Retrieving correct residue sequence numbers for each aligned position
        # i.e positions that both query and subject contain aminoacids
        for index, resi in enumerate(aln.columns):
            index += 1  # the first residue in HHpredHitAlignment is 1
            if aln.gap_at(index) is False:
                if len(query_seq2pdb[index - 1]) == 2 and len(sbjct_seq2pdb[index - 1]) == 2:
                    preqpairs.append(query_seq2pdb[index - 1][1])
                    prespairs.append(sbjct_seq2pdb[index - 1][1])
            else:
                # then this corresponds to a partial-alignment chunk and we can save it
                if preqpairs:
                    qpairs.append(preqpairs)
                    spairs.append(prespairs)
                    preqpairs = []
                    prespairs = []

        if preqpairs:
            qpairs.append(preqpairs)
            spairs.append(prespairs)
        self.qpairs = qpairs
        self.spairs = spairs

        self.global_qpairs = [[item for chunk in qpairs for item in chunk]]
        self.global_spairs = [[item for chunk in spairs for item in chunk]]



    def superimpose_structures(self, aln: HHpredHitAlignment, partial_alignment: bool = False):
        """ Moves the two molecules to the origin of coordinates.
            Aligns the full structures according the fragment and obtains RMSDs and distances.

        :param qpairs: List of CA to be aligned in the query structure
        :param spairs: List of CA to be aligned in the subject structure
        :return:
        """

        self._get_pairs(aln)

        # Re-align if command is called twice
        if self.qaPDB:
            self.qaPDB = {}
            self.saPDB = {}
            self.dst = []

        if partial_alignment is False:
            qpairs = self.global_qpairs
            spairs = self.global_spairs
        else:
            qpairs = self.qpairs
            spairs = self.spairs

        # Print info if the alignment was intended partial but there's only one chunk
        if len(qpairs) == 1 and partial_alignment is True:
            logger.info("The sequence alignment only contains one chunk. Performing global alignment")

        # We only need one query, center to the origin of coordinates.
        # It is the subject that aligns to this template.
        qmol = self.qPDB.copy()
        smol = self.sPDB.copy()
        qmol.center()
        smol.center()
        self.qaPDB[0] = qmol.copy()

        for index, qpair_chunk in enumerate(qpairs):
            logger.info(f"Performing alignment {index+1} with TMalign")
            saPDB, distance = self._superimpose_chunk(qpair_chunk, spairs[index], qmol, smol)
            self.dst.append(distance)
            self.saPDB[index] = saPDB.copy()

        self.global_dst = [[item for chunk in self.dst for item in chunk]]
        return self.qaPDB, self.saPDB

    def _superimpose_chunk(self, qpairs: list, spairs: list, qmol: Molecule, smol: Molecule) -> Tuple[
        Molecule, np.ndarray]:
        """
        Superimposes the part of the total sequence alignment defined
        by the indexes qpairs/spairs
        :param qpairs: list of indexes to align from the pdb query
        :param spairs: list of indexes to align from the pdb subject
        :param qmol: the query pdb
        :param smol: the subject pdb
        :return: smol: the subject molecule aligned.
                distances: a list of the distances between the alpha Carbons after alignment.
        """

        # Copying because we are going to cut the pdbs into the chunks
        copyq = qmol.copy()
        copys = smol.copy()

        copyq.filter('protein and backbone and same residue as index %s' % ' '.join(map(str, qpairs)))
        copys.filter('protein and backbone and same residue as index %s' % ' '.join(map(str, spairs)))
        copyq.write('/tmp/copyq.pdb')
        copys.write('/tmp/copys.pdb')

        # Matrix for VMD
        try:
            # We align subject onto the query
            tm_matrix = get_tmalign_output('/tmp/copys.pdb', '/tmp/copyq.pdb', "matrix")
        except Exception as e:
            raise ChildProcessError(f"TMalign cannot align the PDBs. Error follows: {e}")
        vectran, matrot = tm2vmd(tm_matrix)

        # remove copy files
        os.remove('/tmp/copyq.pdb')
        os.remove('/tmp/copys.pdb')

        # align the whole subject domain and fragment.
        # Copying so that the original smol does not lose the origin of coordinates
        s1mol = smol.copy()
        s1mol.rotateBy(matrot)
        s1mol.moveBy(vectran)
        copys.rotateBy(matrot)
        copys.moveBy(vectran)

        # Compute RMSD for the fragments
        rmsd = MetricRmsd(copyq, 'protein and name CA')
        data = rmsd.project(copys)
        logger.info(f"The RMSD between the fragments is {data} over {len(spairs)} alpha carbons")

        # Compute distances between the two selections
        bbq = copyq.get("coords", sel="protein and name CA")
        bbs = copys.get("coords", sel="protein and name CA")
        distances = np.diagonal(cdist(bbq, bbs))
        return s1mol, distances

    def _construct_chimera(self, qmol, smol, qstart, qend, sstart, send, combination):
        """
        :param qmol: Molecule. The query protein
        :param smol: Molecule. The subject protein in any of its positions
        :param qstart: int. Position to start the cut in the query
        :param qend: int. Position to end the cut in the query
        :param sstart: int. Position to start the cut in the sbjct
        :param send: int. Position to end the cut in the sbjct.
        :return: Molecule, DataFrame Objects.
        chim1: The resulting chimera
        mapping: The mapping from the old residue numbering to the new one
        """
        qmol_copy = qmol.copy()
        smol_copy = smol.copy()
        qmol_copy.filter(f"(protein and same residue as index '{qstart}' to '{qend}')\
         or (not protein and same residue as within 4 of protein and same residue as index '{qstart}' to '{qend}')")
        smol_copy.filter(f"(protein and same residue as index '{sstart}' to '{send}')\
         or (not protein and same residue as within 4 of protein and same residue as index '{qstart}' to '{qend}')")
        # Avoid chimeras that only have a few mutations from
        # one of the parents
        qmol_resid = qmol_copy.get("resid", sel="protein and name CA")
        smol_resid = smol_copy.get("resid", sel="protein and name CA")
        if qmol_resid.size < 10 or smol_resid.size < 10:
            raise NotDiverseChimeraError

        bbq = qmol_copy.get("coords", sel=f"protein and backbone")
        bbs = smol_copy.get("coords", sel=f"protein and backbone")

        distances = cdist(bbq, bbs)
        idx1, idx2 = np.where(distances < 1.3)
        if idx1.any() or idx2.any():
            raise BackboneClashError
        else:
            chim1 = Chimera()
            qmol_copy.renumberResidues()
            smol_copy.renumberResidues()

            if combination == 1:
                last_id = smol_resid[-1] + 1
                new_ids = get_new_resIDs(qmol_copy, last_id)
                qmol_copy.set("resid", new_ids)
                chim1.append(smol_copy)
                chim1.append(qmol_copy)
            else:
                last_id = qmol_resid[-1] + 1
                new_ids = get_new_resIDs(smol_copy, last_id)
                smol_copy.set("resid", new_ids)
                chim1.append(qmol_copy)
                chim1.append(smol_copy)
            chim1.set("chain", "A", "all")
        return chim1, last_id

    def build_chimeras(self, partial_alignment: bool = False, cutoff_distance: float = 1) -> Dict[str, Chimera]:
        """ Build all possible chimeras between the two proteins that fulfill
        these two criteria:
        1) That the distance between the fusion points is below the cutoff distance
        2) That the resulting chimera does not present any backbone clashes

        :return: A dictionary with all the possible chimeras
        """
        if self.dst is None:
            logger.error("You need to align the structures before building the chimeras")

        chimeras = {}
        outcomes = ['Query N-terminal', 'Subject N-terminal', 'Not enough mutations Query N-terminal',
                    'Not enough mutations Subject N-terminal', 'Backbone clash']
        self.chim_positions = dict(zip(outcomes, [[] for i in range(len(outcomes))]))

        q_indices = self.qPDB.get("index", sel="protein and name CA")
        qstart = min(q_indices)
        qend = max(q_indices)

        s_indices = self.sPDB.get("index", sel="protein and name CA")
        sstart = min(s_indices)
        send = max(s_indices)

        if partial_alignment is False:
            qpairs = self.global_qpairs
            spairs = self.global_spairs
            dst = self.global_dst
        else:
            qpairs = self.qpairs
            spairs = self.spairs
            dst = self.dst

        # Get the positions in the fragment closer than the cutoff
        for aln_index, chunk in enumerate(dst):
            if aln_index not in self.saPDB:
                logger.error(f"Alignment {aln_index+1} was not produced. Skipping to next alignment.")
                continue
            fusion_points = [index for index, distance in enumerate(chunk) if distance < cutoff_distance]
            # Build query-subject chimera
            for index in fusion_points:
                qMOL = self.qaPDB[0].copy()
                sMOL = self.saPDB[aln_index].copy()
                xo_query = qpairs[aln_index][index]
                xo_subject = spairs[aln_index][index]
                xo_index = [index for index, number in enumerate(self.global_qpairs[0]) if number == xo_query][0]
                residues = self.qPDB.get("resid", sel="index %s" % ' '.join(map(str, self.global_qpairs[0])))
                xo_resid = residues[xo_index]
                try:
                    xo_query_1 = qpairs[aln_index][index + 1]
                    xo_subject_1 = spairs[aln_index][index + 1]
                except:
                    # Position corresponds to C-terminus limit of the fragment
                    xo_query_1 = [i + 1 for i, qindex in enumerate(q_indices) if qindex == xo_query][0]
                    xo_subject_1 = [i + 1 for i, sindex in enumerate(s_indices) if sindex == xo_subject][0]

                # Combination query-subject
                try:
                    chimera1, xo = self._construct_chimera(qMOL, sMOL, qstart, xo_query,
                                                           xo_subject_1, send, 0)

                    self.chim_positions['Query N-terminal'].append(xo_query)
                    chimeras[f"comb1_{xo_resid}"] = chimera1
                    chimeras[f"comb1_{xo_resid}"].add_crossover(xo)
                except NotDiverseChimeraError:
                    self.chim_positions['Not enough mutations Query N-terminal'].append(xo_query)
                except BackboneClashError:
                    self.chim_positions['Backbone clash'].append(xo_query)

                # Combination subject-query
                try:
                    chimera2, xo = self._construct_chimera(qMOL, sMOL, xo_query_1, qend,
                                                           sstart, xo_subject, 1)

                    self.chim_positions['Subject N-terminal'].append(xo_query)
                    chimeras[f"comb2_{xo_resid}"] = chimera2
                    chimeras[f"comb2_{xo_resid}"].add_crossover(xo)
                except NotDiverseChimeraError:
                    self.chim_positions['Not enough mutations Subject N-terminal'].append(xo_query)
                except BackboneClashError:
                    self.chim_positions['Backbone clash'].append(xo_query)

        if not chimeras:
            logger.warning("No combination of query and subject produced a chimera that matched the criteria")

        return chimeras

    def plot_curves(self, query: str):
        """ Plots the distance between the alpha carbons in the structure for the fragment region along with the number of backbone clashes of the resulting chimera for each position

        :param dst: np.ndarray: an array containign the distance for each fragment position
        :param bbContacts1: an array containing the number of bb contacts for each fragment position for the resulting chimera of combination query-subject
        :param bbContacts2: an array containing the number of bb contacts for each fragment position for the resulting chimera of combination subject - query
        
        :return: a matplotlib.pyplot figure.
        """
        if self.chim_positions is None:
            logger.error("You need to build the chimeras before plotting. Call build_chimeras()")
            return

        dst = self.global_dst[0]

        residues = self.qPDB.get("resid", sel="index %s" % ' '.join(map(str, self.global_qpairs[0])))

        resids = {}
        distances = {}

        for key, value in self.chim_positions.items():
            if value:
                resids[key] = self.qPDB.get("resid", sel="index %s" % ' '.join(map(str, value)))
            else:
                resids[key] = np.zeros(0)
        for key, value in resids.items():
            if value.any():
                indices = [index for index, resi in enumerate(residues) if resi in value]
                distances[key] = [distance for index, distance in enumerate(dst) if index in indices]
            else:
                distances[key] = []

        color = [('#FCB711', "X"), ('#CC004C', "x"), ('gray', 'o'), ('gray', "o"), ('black', ".")]
        colors = dict(zip(resids.keys(), color))

        fig, ax = plt.subplots(figsize=(12, 9))
        ax.plot(residues, dst, '-', color='black', label='distance q-s')
        ax.set_xlabel(f"Residue in the fragment relative to domain {query}", fontsize=24)
        ax.set_ylabel(r'Distance ($\AA$)', fontsize=24)
        i = 0
        for key, value in sorted(resids.items()):
            if value.any():
                i += 1
                ax.plot(value, distances[key], colors[key][1], markersize=18, color=colors[key][0], label=key)

        ax.tick_params(labelsize=20)
        ax.tick_params(labelsize=20)
        ax.legend(loc=9, bbox_to_anchor=(0.5, 1.35), fontsize=18)
        plt.show()

    # def build_recombinations(self,partial_alignment:bool = False, cutoff_distance=1):
    #     """
    #     #TODO for now implemented only to be GLOBAL ALIGNMENT
    #     #TODO NEEDS TESTING - still working with resids
    #     #TODO perhaps is a good idea to follow this algorithm for the single crossover chimeras
    #     :return:
    #     """
    #     if self.dst is None:
    #         logger.error("You need to align the structures before building the chimeras")
    #         return
    #
    #     recombinations = {}
    #
    #     if partial_alignment is False:
    #         qpairs = self.global_qpairs
    #         spairs = self.global_spairs
    #         dst = self.global_dst
    #     else:
    #         qpairs = self.qpairs
    #         spairs = self.spairs
    #         dst = self.dst
    #
    #     combs = []
    #
    #     fusion_points = [index for index, distance in enumerate(dst[0]) if distance < cutoff_distance]
    #
    #     for pair in itertools.combinations(fusion_points, 2):
    #         if pair[0] < pair[1] and pair[1] - pair[0] > 10:
    #             combs.append(pair)
    #
    #     qMOL = self.qaPDB[0].copy()
    #     sMOL = self.saPDB[0].copy()
    #     for crossovers in combs:
    #         xo_query_1 = qpairs[0][crossovers[0]]
    #         xo_query_2 = qpairs[0][crossovers[1]]
    #         xo_subject_1 = spairs[0][crossovers[0]]
    #         xo_subject_2 = spairs[0][crossovers[1]]
    #         try:
    #             recomb1,recomb2= self._construct_recombination(qMOL,sMOL,crossovers[0],crossovers[1],
    #                                                             xo_query_1, xo_query_2,
    #                                                             xo_subject_1, xo_subject_2)
    #             recombinations[f"comb1_{crossovers[0]}_{crossovers[1]}"] = Chimera(recomb1)
    #             recombinations[f"comb2_{crossovers[0]}_{crossovers[1]}"] = Chimera(recomb2)
    #         except:
    #             continue
    #
    #     return recombinations
    #
    # def _construct_recombination(self, qmol,smol,index1,index2,
    #                              xo_query_1, xo_query_2,
    #                              xo_subject_1, xo_subject_2):
    #     qmol_1 = qmol.copy()
    #     qmol_2 = qmol.copy()
    #     qmol_3 = qmol.copy()
    #
    #     smol_1 = smol.copy()
    #     smol_2 = smol.copy()
    #     smol_3 = smol.copy()
    #
    #     qstart = min(qmol.resid)
    #     qend = max(qmol.resid)
    #
    #     sstart = min(qmol.resid)
    #     send = max(qmol.resid)
    #
    #     qmol_1.filter(f"protein and resid '{qstart}' to '{xo_query_1}'")
    #     qmol_2.filter(f"protein and resid '{xo_query_1+1}' to '{xo_query_2}'")
    #     qmol_3.filter(f"protein and resid '{xo_query_2+1}' to '{qend}'")
    #
    #     smol_1.filter(f"protein and resid '{sstart}' to '{xo_query_1}'")
    #     smol_2.filter(f"protein and resid '{xo_subject_1+1}' to '{xo_subject_2}'")
    #     smol_3.filter(f"protein and resid '{xo_subject_2+1}' to '{send}'")
    #
    #     # Combination query - subject
    #     chim1 = Molecule()
    #     chim1.append(qmol_1)
    #     chim1.append(smol_2)
    #     chim1.append(qmol_3)
    #     try:
    #         self.check_bb_clashes(chim1)
    #     except BackboneClashError:
    #         logger.warning(f"comb1_{index1}_{index2} present backbone clashes")
    #
    #
    #     # combination subject - query
    #     chim2 = Molecule()
    #     chim2.append(smol_1)
    #     chim2.append(qmol_2)
    #     chim2.append(smol_3)
    #     try:
    #         self.check_bb_clashes(chim1)
    #     except BackboneClashError:
    #         logger.warning(f"comb2_{index1}_{index2} present backbone clashes")
    #
    #     return chim1, chim2
    #
    # def check_bb_clashes(self,chim_mol):
    #     """
    #
    #     :param chim_mol:
    #     :return:
    #     """
    #     matrix = cdist(chim_mol.coords[:, :, 0], chim_mol.coords[:, :, 0])
    #     np.fill_diagonal(matrix, np.inf)
    #     idx1, idx2 = np.where(matrix < 1.3)
    #     if idx1.any() or idx2.any():
    #         raise BackboneClashError

