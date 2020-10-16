import os
import logging

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
TM_BIN = f"{ROOT_DIR}/TMalign"

print(TM_BIN)
logger = logging.getLogger('protlego')

