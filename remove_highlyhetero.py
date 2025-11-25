#!/home/zchen/projects/rbge/zedchen/env/easy353/bin/python3

import subprocess as sbp
import sys
from Bio import AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from collections import Counter
import numpy as np
import argparse

#remove sequences from an alignment that are simply too heterogeneous
#list of functions
func_list=[
  'get_consensus',
  'filter_divergent_sequences'
]

#FUNCTIONS

#make a concensus seq
def get_consensus(alignment):
  #"""Generate a simple consensus sequence"""
  consensus = []
  for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    most_common = Counter(column).most_common(1)[0][0]
    consensus.append(most_common)
  return ''.join(consensus)
  
  
def filter_divergent_sequences(input_file, output_file, threshold=0.5, method='majority'):
  # Read alignment
  alignment = AlignIO.read(input_file, 'fasta')
  # Calculate distance matrix
  calculator = DistanceCalculator('identity')
  dm = calculator.get_distance(alignment)
  if method == 'majority':
  # Method 1: Distance from consensus sequence
    consensus = get_consensus(alignment)
    distances = [calculator._pairwise(seq, consensus) for seq in alignment]
    ref_distance = 0  # Consensus is reference
  else:
  # Method 2: Distance from centroid (most central sequence)
    avg_distances = np.mean(dm.matrix, axis=1)
    centroid_idx = np.argmin(avg_distances)
    distances = dm.matrix[centroid_idx]
    ref_distance = avg_distances[centroid_idx]
  # Filter sequences
  kept_records = []
  for i, record in enumerate(alignment):
    if distances[i] - ref_distance <= threshold:
      kept_records.append(record)
    
    # Write output
  SeqIO.write(kept_records, output_file, 'fasta')
    
  print(f"Kept {len(kept_records)} of {len(alignment)} sequences "
        f"(removed {len(alignment)-len(kept_records)} outliers)")
##################################################
def main(input_file, output_file, threshold=0.5, method='majority'):
  filter_divergent_sequences(input_file, output_file, threshold=0.5, method='majority')
    

##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='Filter divergent sequences from FASTA alignment', \
	  usage = 'remove_highlyhetero.py <input_file> <output_file>')
  parser.add_argument('input', help='Input FASTA file')
  parser.add_argument('output', help='Output FASTA file')
  parser.add_argument('-t', '--threshold', type=float, default=0.4,
                       help='Maximum distance threshold (0-1)')
  parser.add_argument('-m', '--method', choices=['majority', 'cluster'], default='majority',
                       help='Filtering method (distance from consensus or cluster centroid)')
  args = parser.parse_args()
  
  main(args.input, args.output, args.threshold, args.method)