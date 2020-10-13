import csv, sqlite3
conn = sqlite3.connect('fuzzle2.07.db')

c = conn.cursor()
c.execute('''CREATE TABLE IF NOT EXISTS hh207clusters
            (id integer, query text, q_scop_id varchar, no integer, sbjct text,
            s_scop_id varchar, s_desc text, prob double precision,
            eval double precision, pval double precision, score double precision,
            ss double precision, cols integer, q_start integer, q_end integer,
            s_start integer, s_end integer, hmm integer, ident integer,
            q_sufam_id varchar, s_sufam_id varchar, q_fold_id varchar,
            s_fold_id varchar,rmsd_pair double precision, ca_pair integer,
            rmsd_tm_pair double precision, score_tm_pair double precision,
            ca_tm_pair integer, rmsd_tm double precision, score_tm double precision,
            ca_tm integer, q_tm_start integer, q_tm_end integer,s_tm_start integer,
            s_tm_end integer, q_cluster text, s_cluster text)''')

with open('./hh207clusters.csv','r') as fin:
    dr = csv.DictReader(fin)
    to_db = [(i['id'], i['query'],i['q_scop_id'],i['no'],
              i['sbjct'], i['s_scop_id'], i['s_desc'], i['prob'],
              i['eval'], i['pval'], i['score'], i['ss'],
              i['cols'], i['q_start'], i['q_end'], i['s_start'],
              i['s_end'], i['hmm'], i['ident'], i['q_sufam_id'],
              i['s_sufam_id'], i['q_fold_id'], i['s_fold_id'],
              i['rmsd_pair'], i['ca_pair'], i['rmsd_tm_pair'], i['score_tm_pair'],
              i['ca_tm_pair'], i['rmsd_tm'], i['score_tm'], i['ca_tm'],
              i['q_tm_start'], i['q_tm_end'], i['s_tm_start'], i['s_tm_end'],
              i['q_cluster'], i['s_cluster']) for i in dr]

c.executemany("INSERT INTO hh207clusters (id,query,q_scop_id,no,sbjct,s_scop_id,s_desc,prob,eval,pval,score,ss,cols,q_start,q_end,s_start,s_end,hmm,ident,q_sufam_id,s_sufam_id,q_fold_id,s_fold_id,rmsd_pair,ca_pair,rmsd_tm_pair,score_tm_pair,ca_tm_pair,rmsd_tm,score_tm,ca_tm,q_tm_start,q_tm_end,s_tm_start,s_tm_end,q_cluster,s_cluster) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);", to_db)

c.execute("CREATE INDEX query ON hh207clusters (query)")
c.execute("CREATE INDEX sbjct ON hh207clusters (sbjct)")
c.execute("CREATE INDEX query_subject ON hh207clusters (query,sbjct)")
c.execute("CREATE INDEX fuzzle_id ON hh207clusters (id)")
c.execute("CREATE INDEX folds ON hh207clusters (q_fold_id,s_fold_id)")
c.execute("CREATE INDEX sufams ON hh207clusters (q_sufam_id,s_sufam_id)")
c.execute("CREATE INDEX fams ON hh207clusters (q_scop_id,s_scop_id)")

conn.commit()
c.close()
