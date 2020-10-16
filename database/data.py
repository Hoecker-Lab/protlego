from protlego.definitions import ROOT_DIR, logger
from typing import List
import numpy as np
import sqlite3
from typing import NamedTuple

conn = sqlite3.connect(f'{ROOT_DIR}/fuzzle2.07.db')
cur = conn.cursor()

class Hit(NamedTuple):
    """
    Some of the documentation of this function was
    taken from the hhsuite python documentation:
    https://github.com/soedinglab/hh-suite/wiki
    as the sequence information from the Fuzzle hits
    come from HHsearch.
    The structural superimpositions were performed with
    TMalign:   https://zhanglab.ccmb.med.umich.edu/TM-align/

    =======================================================
    """
    id: int
    """ The database id for this hit """
    query: str
    """ The 7-letter SCOP95 code for the query domain"""
    q_scop_id: str
    """ The SCOP family the query belongs to"""
    no: int
    """ The HHsearch hit number for this query"""
    sbjct: str
    """ A 7-letter SCOP95 code for the subject domain """
    s_scop_id: str
    """ The SCOP family the subject belongs to"""
    s_desc: str
    """ Description of the subject domain."""
    prob: float
    """ HHsearch probability"""
    eval: float
    """ E-value"""
    pval: float
    """ p-value"""
    score: float
    """ The raw score is computed by the Viterbi HMM-HMM alignment excluding the secondary structure score. It is the sum of similarities of the aligned profile columns minus the position-specific gap penalties in bits."""
    ss: float
    """ The secondary structure score. This score tells you how well the PSIPRED-predicted (3-state) or actual DSSP-determined (8-state) secondary structure sequences agree with each other."""
    cols: int
    """ The number of aligned Match columns in the HMM-HMM alignment."""
    q_start: int
    """ Residue where the alignment starts for the query domain (sequence position)"""
    q_end: int
    """ Residue where the alignment ends for the query domain (sequence position) """
    s_start: int
    """ Residue where the alignment starts for the subject domain (sequence position)"""
    s_end: int
    """ Residue where the alignment ends for the subject domain (sequence position)"""
    hmm: int
    """ int """
    ident: float
    """ Identity %  for the sequence alignment"""
    q_sufam_id: str
    """ The superfamily the query belongs to"""
    s_sufam_id: str
    """ The superfamily the subject belongs to """
    q_fold_id: str
    """ The fold the query belongs to"""
    s_fold_id: str
    """ The fold the query belongs to"""
    rmsd_pair: float
    """ RMSD for the alignment between the two domains, strictly taking the alpha carbons from the structures that exactly appear in the HHsearch sequence alignment"""
    ca_pair: int
    """ The number of alpha carbon pairs that were used for the rmsd_pair calculation."""
    rmsd_tm_pair: float
    """RMSD for the TMalign alignment between the two domains, passing the sequence alignment as seed"""
    score_tm_pair: float
    """ TM-score for the rmsd_tm_pair superposition"""
    ca_tm_pair: int
    """ The number of alpha carbon pairs that were used for the rmsd_tm_pair calculation"""
    rmsd_tm: float
    """ RMSD for the TMalign alignment between the two domains without seed """
    score_tm: float
    """ TM-score for the rmsd_tm superposition"""
    ca_tm: int
    """ The number of alpha carbon pairs that were used for the rmsd_tm calculation"""
    q_tm_start: int
    """ Residue in the query structure where the rmsd_tm_pair alignment starts """
    q_tm_end: int
    """ Residue in the query structure where the rmsd_tm_pair alignment ends"""
    s_tm_start: int
    """ Residue in the subject structure where the rmsd_tm_pair alignment starts """
    s_tm_end: int
    """ Residue in the subject structure where the rmsd_tm_pair alignment ends"""
    q_cluster: str
    """ Query cluster where this fragment belongs to."""
    s_cluster: str
    """ Subject cluster where this fragment belongs to """

    def __repr__(self):
        return f"Hit between {self.query} and {self.sbjct} with probability {self.prob} %\n"


def parse_hit(line: List[str]) -> Hit:
    """

    :param line:
    :return:
    """
    return Hit(
        id=int(line[0]),
        query=line[1],
        q_scop_id=line[2],
        no=int(line[3]),
        sbjct=line[4],
        s_scop_id=line[5],
        s_desc=line[6],
        prob=float(line[7]),
        eval=float(line[8]),
        pval=float(line[9]),
        score=float(line[10]),
        ss=float(line[11]),
        cols=int(line[12]),
        q_start=int(line[13]),
        q_end=int(line[14]),
        s_start=int(line[15]),
        s_end=int(line[16]),
        hmm=int(line[17]),
        ident=int(line[18]),
        q_sufam_id=line[19],
        s_sufam_id=line[20],
        q_fold_id=line[21],
        s_fold_id=line[22],
        rmsd_pair=float(line[23]),
        ca_pair=int(line[24]),
        rmsd_tm_pair=float(line[25]) if line[25] else 'NULL',
        score_tm_pair=float(line[26]) if line[26] else 'NULL',
        ca_tm_pair=int(line[27]) if line[27] else 'NULL',
        rmsd_tm=float(line[28]),
        score_tm=float(line[29]),
        ca_tm=int(line[30]),
        q_tm_start=int(line[31]),
        q_tm_end=int(line[32]),
        s_tm_start=int(line[33]),
        s_tm_end=int(line[34]),
        q_cluster=None if (len(line) <= 35) else str(line[35]),
        s_cluster=None if (len(line) <= 35) else str(line[36])
    )


class Result:
    """ Class handling the data obtained from fuzzle"""

    def __init__(self, ahits: List[Hit]):
        self.hits = ahits

    def __repr__(self):
        return f"Query from Fuzzle with {len(self.hits)} hits" \
               f" belonging to {len(self.list_folds)} fold(s)"

    def __getitem__(self, key):
        return self.hits[key]

    def append(self, result):
        hits = self.hits + result.hits
        self.hits = np.unique(hits)

    @property
    def ids(self):
        """It returns all the hits IDs"""
        return [hit.id for hit in self.hits]

    @property
    def avg_len(self):
        """ It returns the Aminoacid average of the returned hits"""
        return np.mean([hit.cols for hit in self.hits])

    @property
    def std_len(self):
        """ It returns the Aminoacid standard deviation in the hits list"""
        return np.std([hit.cols for hit in self.hits])

    @property
    def unique_domains(self):
        """ It returns a list of unique domains"""
        unique_domains=np.unique([(hit.query, hit.sbjct) for hit in self.hits])
        logger.info(f"The query contains {len(unique_domains)} different domains")
        return unique_domains

    @property
    def unique_clusters(self):
        """ It returns a list of unique domains"""
        unique_clusters=np.unique([(hit.q_cluster, hit.s_cluster) for hit in self.hits])
        logger.info(f"The query contains {len(unique_clusters)} different clusters")
        return unique_clusters

    @property
    def list_folds(self):
        """ It returns the list of unique folds in the hits list"""
        return np.unique([(hit.q_fold_id, hit.s_fold_id) for hit in self.hits])

    @property
    def list_sufams(self):
        """ It returns the list of unique folds in the hits list"""
        return np.unique([(hit.q_sufam_id, hit.s_sufam_id) for hit in self.hits])

    @property
    def list_fams(self):
        """ It returns the list of unique folds in the hits list"""
        return np.unique([(hit.q_scop_id, hit.s_scop_id) for hit in self.hits])


def fetch_subspace(prob: int = 70, rmsd: float = 3.0,
                   ca_min: int = 10, ca_max: int = 200, score_tm_pair: float = 0.3, ratio: float = 1.25,
                   scop_q: str = None, diff_folds: bool = True) -> Result:
    """ Returns the entries in Fuzzle that satisfy the conditions:

    :param prob: Lower cutoff for the hit probability
    :param rmsd: Upper cutoff for rmsd_tm_pair
    :param ca_min: Lower cutoff for the number of AA (ca_tm_pair)
    :param ca_max: Upper cutoff for the number of AA (ca_tm_pair)
    :param score_tm_pair: Lower cutoff for the tm_score
    :param ratio: Proportion between cols/ca_tm_pair
    :param scop_q: A SCOP class. It will retrieve hits that contains domains from this class
    :return: A Result object
    """

    diff = 'q_fold_id != s_fold_id'
    if diff_folds is False:
        diff = None
    cond = None
    cond1 = None
    cond2 = None

    # handle SCOP class
    if scop_q:
        # SCOPe Level
        level = len(scop_q.split('.'))
        if level == 2:  # folds
            scop = "q_fold_id"
        elif level == 3:  # superfamilies
            scop = "q_sufam_id"
        elif level == 4:  # families
            scop = "q_scop_id"
        else:
            raise ValueError("The specified SCOP level does not exist")

        cond1 = scop + "='" + scop_q + "'"  # This prints: q_fold_id='c.23'. Not sure how to improve it.
        cond2 = scop.replace('q', 's') + "='" + scop_q + "'"

    ex = f"select * from hh207clusters where prob > {prob}" \
         f" and rmsd_tm_pair < {rmsd}" \
         f" and ca_tm_pair > {ca_min} and ca_tm_pair < {ca_max}" \
         f" and score_tm_pair > {score_tm_pair} " \
         f"and cast(cols as float)/ca_tm_pair < {ratio} "
    ex += f"and {diff} " if diff is not None else ""
    ex += f"and {cond} " if cond is not None else ""
    ex += f"and {cond1} " if cond1 is not None else ""
    ex += f"or {cond2} " if cond2 is not None else ""
    cur.execute(ex)

    ahits = []
    tmphits = cur.fetchall()
    for line in tmphits:
        ahits.append(parse_hit(line))
    return Result(ahits)


def _domain_from_PDB(pdb: str):
    """ Returns all domains included in a SCOP version that come from a PDB.

    :param pdb: A 4-letter pdb code
    :return: A list of domains that this PDB contains
    """
    pdb = pdb.lower()
    cur.execute(f"SELECT sdomain FROM scop_pdbref_scop WHERE pdbref like '{pdb}'")
    hits = cur.fetchall()
    domains = tuple(np.unique([x[0] for x in hits]))
    if domains:
        if len(domains) == 1:
            domains = str("('" + domains[0] + "')")
        cur.execute(f"SELECT ref FROM astral_astral95 WHERE id in {domains}")
        hits = cur.fetchall()
        domains = np.unique([x[0] for x in hits])
        logger.info(f"Domain(s) {domains} are present in the Fuzzle database"
                    f" as representatives for pdb {pdb.capitalize()}")
    else:
        from protlego.builder.builder import NotCorrectPDBError
        raise NotCorrectPDBError(f"PDB code {pdb.capitalize()} does not yet appear in the SCOP database")

    return domains


def fetch_byPDBs(pdb1: str, pdb2: str, prob: int = 70, rmsd: float = 3.0,
                 ca_min: int = 10, ca_max: int = 200, score_tm_pair: float = 0.3, ratio: float = 1.25,
                 diff_folds: bool = True):
    """ Includes all hits among the domains that belong to a pair of PDBs

    :param pdb1: The first PDB to check
    :param pdb2: The second PDB to check
    :param prob: the minimum allowed HHsearch probability
    :param rmsd: The maximum allowed RMSD (rmsd_tm_pair: "RMSD for the TMalign alignment between the two domains, passing the sequence alignment as seed)
    :param ca_min: The minimum allowed fragment length (for the TMalign alignment)
    :param ca_max: The maximun allowed fragment length (for the TMalign alignment)
    :param score_tm_pair: The minimum allowed TM-score (for the TMalign alignment)
    :param ratio: the maximum ratio for the sequence and structural alignment lengths (cols / ca_tm_pair)
    :param diff_folds: Whether to exclude hits from the same fold (True) or not (False)
    :return: A result class obtaining the hits that fulfill these criteria
    """
    domains1 = tuple(_domain_from_PDB(pdb1))
    domains2 = tuple(_domain_from_PDB(pdb2))
    if len(domains1) == 1: domains1 == domains1[0]
    if len(domains2) == 1: domains2 == domains2[0]
    domains = domains1 + domains2
    diff = 'q_fold_id != s_fold_id'
    if diff_folds is False:
        diff = None
    ex = f"select * from hh207clusters where prob > {prob}" \
         f" and rmsd_tm_pair < {rmsd}" \
         f" and ca_tm_pair > {ca_min} and ca_tm_pair < {ca_max}" \
         f" and score_tm_pair > {score_tm_pair}" \
         f" and cast(cols as float)/ca_tm_pair < {ratio}"
    ex += f" and query in {domains}"
    ex += f" and sbjct in {domains}"
    ex += f" and {diff} " if diff is not None else ""

    cur.execute(ex)
    ahits = []
    tmphits = cur.fetchall()
    for line in tmphits:
        ahits.append(parse_hit(line))
    return Result(ahits)


def fetch_byPDB(pdb: str, prob: int = 70, rmsd: float = 3.0,
                ca_min: int = 10, ca_max: int = 200, score_tm_pair: float = 0.3, ratio: float = 1.25,
                diff_folds: bool = True):
    """ Returns the entries in Fuzzle that contain the representative domains that correspond to that PDB

    :param prob: Lower cutoff for the hit probability
    :param rmsd: Upper cutoff for rmsd_tm_pair
    :param ca_min: Lower cutoff for the number of AA (ca_tm_pair)
    :param ca_max: Upper cutoff for the number of AA (ca_tm_pair)
    :param score_tm_pair: Lower cutoff for the tm_score
    :param ratio: Proportion between cols/ca_tm_pair
    :param scop_q: A SCOP class. It will retrieve hits that contains domains from this class
    :param query: A SCOP protein domain. It will retrieve hits that contain this query
    :return: A Result object
    """

    domains = _domain_from_PDB(pdb)
    hits = {}
    for query in domains:
        hits[query] = fetch_subspace(prob=prob, rmsd=rmsd,
                                     ca_min=ca_min, ca_max=ca_max, score_tm_pair=score_tm_pair, ratio=ratio,
                                     query=query, diff_folds=diff_folds)

    ahits = hits[domains[0]]
    for index, query in enumerate(domains):
        if index == 0: continue
        ahits.append(hits[query])

    return ahits


def fetch_by_domains(domain1: str, domain2: str, prob: int = 70, rmsd: float = 3.0,
                     ca_min: int = 10, ca_max: int = 200, score_tm_pair: float = 0.3, ratio: float = 1.25,
                     diff_folds: bool = True):
    """ Fetch all the hits between two parent domains

    :param domain1: The 7 letter code for one of the parents
    :param domain2: The 7 letter code for one of the parents
    :param prob: the minimum allowed HHsearch probability
    :param rmsd: The maximum allowed RMSD (rmsd_tm_pair: "RMSD for the TMalign alignment between the two domains, passing the sequence alignment as seed)
    :param ca_min: The minimum allowed fragment length (for the TMalign alignment)
    :param ca_max: The maximun allowed fragment length (for the TMalign alignment)
    :param score_tm_pair: The minimum allowed TM-score (for the TMalign alignment)
    :param ratio: the maximum ratio for the sequence and structural alignment lengths (cols / ca_tm_pair)
    :param diff_folds: Whether to exclude hits from the same fold (True) or not (False)
    :return: A result class obtaining the hits that fulfill these criteria
    """
    diff = 'q_fold_id != s_fold_id'
    if diff_folds is False:
        diff = None
    ex = f"select * from hh207clusters where prob > {prob}" \
         f" and (rmsd_tm_pair < {rmsd} or rmsd_tm_pair is NULL)" \
         f" and ca_tm_pair > {ca_min} and ca_tm_pair < {ca_max}" \
         f" and score_tm_pair > {score_tm_pair}" \
         f" and cast(cols as float)/ca_tm_pair < {ratio}" \
         f" and ((query='{domain1}' and sbjct='{domain2}')" \
         f" or (sbjct = '{domain1}' and query='{domain2}'))" \
         f" and query != sbjct"
    ex += f" and {diff} " if diff is not None else ""

    cur.execute(ex)
    ahits = []
    tmphits = cur.fetchall()
    if len(tmphits) == 1:
        hit = parse_hit(tmphits[0])
    else:
        for line in tmphits:
            ahits.append(parse_hit(line))
        hit = Result(ahits)
    return hit


def fetch_by_domain(domain: str, prob: int = 70, rmsd: float = 3.0,
                    ca_min: int = 10, ca_max: int = 200, score_tm_pair: float = 0.3, ratio: float = 1.25,
                    diff_folds: bool = True):
    """ Fetch all the hits that contain a specific domain

    :param domain: The 7 letter code for one of the parents
    :param prob: the minimum allowed HHsearch probability
    :param rmsd: The maximum allowed RMSD (rmsd_tm_pair: "RMSD for the TMalign alignment between the two domains, passing the sequence alignment as seed)
    :param ca_min: The minimum allowed fragment length (for the TMalign alignment)
    :param ca_max: The maximun allowed fragment length (for the TMalign alignment)
    :param score_tm_pair: The minimum allowed TM-score (for the TMalign alignment)
    :param ratio: the maximum ratio for the sequence and structural alignment lengths (cols / ca_tm_pair)
    :param diff_folds: Whether to exclude hits from the same fold (True) or not (False)
    :return: A result class obtaining the hits that fulfill these criteria
    """
    diff = 'q_fold_id != s_fold_id'
    if diff_folds is False:
        diff = None
    ex = f"select * from hh207clusters where prob > {prob}" \
         f" and (rmsd_tm_pair < {rmsd} or rmsd_tm_pair is NULL)" \
         f" and ca_tm_pair > {ca_min} and ca_tm_pair < {ca_max}" \
         f" and score_tm_pair > {score_tm_pair}" \
         f" and cast(cols as float)/ca_tm_pair < {ratio}" \
         f" and (sbjct = '{domain}' or query='{domain}')"
    ex += f" and {diff} " if diff is not None else ""

    cur.execute(ex)
    ahits = []
    tmphits = cur.fetchall()
    if len(tmphits) == 1:
        hit = parse_hit(tmphits[0])
    else:
        for line in tmphits:
            ahits.append(parse_hit(line))
        hit = Result(ahits)
    return hit


def validate_scopid(query: str) -> bool:
    """A SCOP domain is  A 7-character sid that consists of "d"
    followed by the 4-character PDB ID of the file of origin,
    the PDB chain ID ('_' if none, '.' if multiple as is the
    case in genetic domains), and a single character
    (usually an integer) if needed to specify the domain uniquely
    ('_' if not). Sids are currently all lower case,
    even when the chain letter is upper case.
    Examples:  include d4akea1, d1reqa2, and d1cph.1.
    :param query: The seven letter domain for the query
    """

    return len(query) == 7 and query[0] == 'd'


def fetch_id(fuzzle_id: int) -> Hit:
    """
    Returns the hit in fuzzle with that ID
    :param fuzzle_id: The Fuzzle HIT id to retrieve from hh207clusters
    :return: A Hit object
    """

    cur.execute("select * from hh207clusters where id = ?", (fuzzle_id,))

    tmphits = cur.fetchall()
    hit = parse_hit(list(tmphits[0]))
    return hit


def fetch_group(group1, group2=None, prob: int = 70, rmsd: float = 3.0,
                ca_min: int = 10, ca_max: int = 200, score_tm_pair: float = 0.3, ratio: float = 1.25,
                diff_folds: bool = True) -> Result:
    """ Fetching all hits between two specific groups (folds, superfamilies and families)
    or inside one specific group (group1)

    :param group1: The first group from where to search. E.g 'c.2'
    :param group2 (optional): The second group from where to search. E.g 'c.2'
    :param prob: the minimum allowed HHsearch probability
    :param rmsd: The maximum allowed RMSD (rmsd_tm_pair: "RMSD for the TMalign alignment between the two domains, passing the sequence alignment as seed)
    :param ca_min: The minimum allowed fragment length (for the TMalign alignment)
    :param ca_max: The maximun allowed fragment length (for the TMalign alignment)
    :param score_tm_pair: The minimum allowed TM-score (for the TMalign alignment)
    :param ratio: the maximum ratio for the sequence and structural alignment lengths (cols / ca_tm_pair)
    :return: A Result class with the hits that fulfill these criteria
    """

    diff = 'q_fold_id != s_fold_id'
    if diff_folds is False:
        diff = None

    level1 = len(group1.split('.'))
    if group2:
        level2 = len(group2.split('.'))
        if level1 != level2:
            raise ValueError("The two groups must belong to the same SCOP level")

    if level1 == 2:  # folds
        scop = "q_fold_id"
    elif level1 == 3:  # superfamilies
        scop = "q_sufam_id"
    elif level1 == 4:  # families
        scop = "q_scop_id"
    else:
        raise ValueError("The specified SCOP level does not exist")

    cond1 = scop + "='" + group1 + "'"
    cond_1 = scop.replace('q', 's') + "='" + group1 + "'"

    if group2:
        cond2 = scop.replace('q', 's') + "='" + group2 + "'"
        cond_2 = scop + "='" + group2 + "'"

    ex = f"select * from hh207clusters where prob > {prob}" \
         f" and rmsd_tm_pair < {rmsd}" \
         f" and ca_tm_pair > {ca_min} and ca_tm_pair < {ca_max}" \
         f" and score_tm_pair > {score_tm_pair}" \
         f" and cast(cols as float)/ca_tm_pair < {ratio}"
    ex += f" and (({cond1} and {cond2}) " if group2 else f" and ({cond1}"
    ex += f" or ({cond_1} and {cond_2})) " if group2 else f" or {cond_1})"
    ex += f" and {diff} " if diff is not None else ""

    cur.execute(ex)
    ahits = []
    tmphits = cur.fetchall()
    for line in tmphits:
        ahits.append(parse_hit(line))
    return Result(ahits)


def filter_hits_domain(ahits, domain):
    """Search all hits from a Result class where a certain domain appears

    :param ahits: An object Result
    :param domain: a SCOPe domain identifier
    :return: np.array. The starts and ends for the domains in all the hits it appears.
    """
    centersq = [(x.q_start, x.q_end, x.id) for x in ahits if x.query == domain]
    centerss = [(x.s_start, x.s_end, x.id) for x in ahits if x.sbjct == domain]
    return np.asarray(centersq + centerss)
