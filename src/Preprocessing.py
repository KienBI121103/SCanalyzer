import os 
import subprocess
import shutil
import pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from loguru import logger

from config import TOOLs

class Preprocessing:
    """Handel preprocessing of raw RNA sequence data"""
    def __init__(self, input_dir: str, output_dir : str, threshold_trim: float, word_size: int, ref_path: str, num_threads: int = os.cpu_count() or 1):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.threshold_trim = threshold_trim
        self.word_size = word_size
        self.ref_path = Path(ref_path)
        self.num_threads = num_threads
        self.tools = TOOLS
        
        # Output file paths
        self.fastq_path: Optional[Path] = None
        self.fasta_trim_path: Optional[Path] = None
        
        
        # Validation tools and setup
        self._initialize()
        
    def _initialize(self) -> None:
        
        
        
        