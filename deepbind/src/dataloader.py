"""DeepBind truncating dataloader
"""
# python2, 3 compatibility
from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
from pyfaidx import Fasta
from concise.preprocessing import encodeDNA

from kipoi.data import Dataset
from kipoi.metadata import GenomicRanges
import numpy as np
import linecache
# --------------------------------------------
class BedToolLinecache(BedTool):
    """Faster BedTool accessor by Ziga Avsec
    Normal BedTools loops through the whole file to get the
    line of interest. Hence the access it o(n)
    Note: this might load the whole bedfile into memory
    """

    def __getitem__(self, idx):
        line = linecache.getline(self.fn, idx + 1)
        return pybedtools.create_interval_from_list(line.strip().split("\t"))


class TruncSeqDataset(Dataset):
    """
    Truncates the length of each input interval to 101bp
    
    Args:
        intervals_file: bed3 file containing intervals
        fasta_file: file path; Genome sequence
        target_file: file path; path to the targets in the csv format
    """

    def __init__(self, intervals_file, fasta_file, target_file=None, use_linecache=True):

        if use_linecache:
            linecache.clearcache()
            BT = BedToolLinecache
        else:
            BT = BedTool
        self.bt = BT(intervals_file)
        self.fasta_file = fasta_file
        self.fasta_extractor = Fasta(self.fasta_file)

        # Targets
        if target_file is not None:
            self.targets = pd.read_csv(target_file)
        else:
            self.targets = None

    def __len__(self):
        return len(self.bt)
    
    def _interval_to_fasta_id(self, interval):
        ''' convert a PyBedTools interval to Fasta string 
        id, for searching an associated .fasta file
        with PyFaidx '''
        return interval.chrom + ":" + interval.start + "-" + interval.stop
    
    def _compute_relative_coords(self, interval):
        ''' convert an absolute sequence interval into a relative one, and
        extract 101bp centered around the midpoint '''
        
        # compute new 101bp interval
        start = 0
        end = interval.end - interval.start
        
        midpoint = end - start // 2
        new_start = midpoint - 50
        new_stop = midpoint + 51
    
        # Intervals need to be 101bp wide
        assert new_stop - new_start == 101
        return new_start, new_stop

    def __getitem__(self, idx):
        if self.fasta_extractor is None:
            self.fasta_extractor = Fasta(self.fasta_file)
        interval = self.bt[idx]
        interval_fasta_id = self._interval_to_fasta_id(interval)

        if self.targets is not None:
            y = self.targets.iloc[idx].values
        else:
            y = {}

        # Run the fasta extractor
        start, end = self._compute_relative_coords(interval)
        record = self.fasta_extractor[interval_fasta_id]
        seq = record[start:end].seq
        
        return {
            "inputs": encodeDNA(seq),
            "targets": y,
            "metadata": {
                "ranges": GenomicRanges.from_interval(interval)
            }
        }