import os
import logging
import sqlite3


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
TM_BIN = f"{ROOT_DIR}/builder/TMalign"

logger = logging.getLogger('protlego')

conn = sqlite3.connect(f'{ROOT_DIR}/database/fuzzle2.07.db')
cur = conn.cursor()
