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
def remove_taxa(fasta,bad_smp):
  for seq in SeqIO.parse(fasta,'fasta'): #parse can only be used as iterator once. you forgot that
    #print(seq.id)
    if seq.id.replace(' ','').replace('.','_') in bad_smp:
      print(seq.id)
  #the actual output
  seqs=SeqIO.parse(fasta,'fasta')
  records=[SeqRecord(seq.seq,id=seq.id.replace(' ','').replace('.','_'),name='',description='') 
  for seq in seqs 
  if seq.id.replace(' ','').replace('.','_') not in bad_smp and len(set(str(seq.seq))) != 1]
  SeqIO.write(records,fasta,'fasta')
  

##################################################
def main(fasta):  
  bad_smp=['S_phylicifolia_6'] #'P.gramineus_4','P.lucens_1','P.lucens_2'
  reformat=[i.replace(' ','').replace('.','_') for i in bad_smp] #just in case the formatting is wrong (i once put P.pusilus_6 and the '.' was not recognized later)
  remove_taxa(fasta,reformat)
    

##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='remove the taxa that are not good for phylogeny construction', \
	  usage = 'remove_taxa.py <fasta> <taxa>')
  parser.add_argument('fasta', help='fasta file to be filtered', metavar='fasta')
  parser.add_argument('taxa', help='taxa to be removed', metavar='taxa')
  options = parser.parse_args()
  
  main(options.fasta,options.taxa)