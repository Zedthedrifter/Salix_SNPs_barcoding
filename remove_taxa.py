#!/home/zchen/projects/rbge/zedchen/env/easy353/bin/python3

import subprocess as sbp
import sys
import argparse
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#list of functions
func_list=[
  'remove_taxa' #rename fqgz files to use array
]

#FUNCTIONS

#filter out bad taxa and bad seq
def remove_taxa(fasta,taxa):
  seqs=SeqIO.parse(fasta,'fasta')
  records=[SeqRecord(seq.seq,id=seq.description.replace(' ','').replace('.','_').split('salix')[0],name=seq.name,description=seq.description) for seq in seqs if taxa not in seq.description and len(set(str(seq.seq))) != 1]
  SeqIO.write(records,fasta,'fasta')

##################################################
def main(fasta,taxa):
  remove_taxa(fasta,taxa)
    

##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='remove the taxa that are not good for phylogeny construction', \
	  usage = 'remove_taxa.py <fasta> <taxa>')
  parser.add_argument('fasta', help='fasta file to be filtered', metavar='fasta')
  parser.add_argument('taxa', help='taxa to be removed', metavar='taxa')
  options = parser.parse_args()
  
  main(options.fasta,options.taxa)