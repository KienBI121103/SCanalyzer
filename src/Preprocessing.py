import os 
import subprocess
import shutil
import pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from loguru import logger

from config import TOOLs

class Preprocessing:
    """Handel preprocessing of raw DNA sequence data"""
    def __init__(self, input_files: List[str], output_dir : str, thread_trim: float, )