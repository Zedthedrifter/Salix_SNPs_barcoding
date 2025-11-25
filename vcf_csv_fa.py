#!/home/zchen/apps/env/easy353/bin/python3

import subprocess as sbp
import sys
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#import os
#os.environ['OMP_NUM_THREADS'] = '1'
#os.environ['MKL_NUM_THREADS'] = '1'
#os.environ['OPENBLAS_NUM_THREADS'] = '1'



#list of functions


#process in chunks instead for memory issue:
def parse_vcf(vcf, output, batch, step, acc):
	# just in case the vcf file is gz
	import gzip
	if vcf.endswith(".gz"):
		opener = gzip.open
	else:
		opener = open
	# process each line
	snps, snps_info, n = {}, {}, 0
	with opener(vcf, 'rt') as infile:
		for i in infile: 
			if i.startswith('#CHROM'):
				# VCF is tab-delimited, so split by tab
				heading = i.strip().split('\t')
				samples = heading[9:]
				samples = [s.split('/')[-1].replace('.bam', '') for s in samples]
				print(f"processing {len(samples)} samples: {samples}")
				sample_snp = {s: {} for s in samples}  
			else:
				if not i.startswith('#'): 
					n += 1
					if batch <= n and n < batch + step: #the batch selection method
						snps[f"SNP_{n}"] = {s: {} for s in samples}
						print(f"adding SNP_{n}") 
						# COLLECT SNP INFO - split by tab
						fields = i.strip().split('\t')
						snps_info[f"SNP_{n}"] = {heading_col: field_value for heading_col, field_value in zip(heading[:9], fields[:9])}
						# COLLECT SAMPLE INFO - samples start from 9th column
						sample_fields = fields[9:]
						genotypes = [gt.split(':')[0] for gt in sample_fields]
						prob = [gt.split(':')[1] for gt in sample_fields]
						for s, gt in zip(samples, genotypes):
							snps[f"SNP_{n}"][s]['GT'] = gt
							sample_snp[s][f"SNP_{n}"] = gt
						for s, pl in zip(samples, prob):
							snps[f"SNP_{n}"][s]['PL'] = pl.split(',')
	# OUTPUT SNP INFO
	#df = pd.DataFrame.from_dict(snps_info, orient='index')
	#df.to_csv(f'{output}/batch{acc}_all_SNPs_summary.csv', index=True, sep='\t')  # Make output tab-delimited too
	return (snps, snps_info, sample_snp)


#make a translation table for GT
def GT_translation_table(snps_info):
	vfa={k:{} for k in snps_info}
	for k,v in snps_info.items():
		vfa[k]['.']='-'
		if 'INDEL' in v['INFO']:
			total=max([len(v['REF'])]+[len(seq) for seq in v['ALT'].split(',')])
			vfa[k][0]=v['REF']+'-'*(total-len(v['REF']))
			for i,b in enumerate(v['ALT'].split(',')):
				vfa[k][i+1]=b+'-'*(total-len(b))
		else:
			vfa[k][0]=v['REF']
			for i,b in enumerate(v['ALT'].split(',')):
				vfa[k][i+1]=b
		#elif len(v['ALT'])>1:
		#    print(k,v)
		#    break
	return(vfa)

#keep homozygous for heterozygous with 0/* structure
def filter_snps(snps,samples,vfa):
	keep,fasta=[],{s:[] for s in samples}
	for snp,v in snps.items():
		GTs=[v2['GT'] for v2 in v.values()]
		#check if ambiguious
		check=[0 if '.' not in gt else 1 for gt in GTs]
		if sum(check)<=2: #at most two ambiguious samples
			#check for ambiguous bases (only one sample allowed to have ambiguous base)
			check=[0 if gt != './.' else 1 for gt in GTs]
			if sum(check)<=1:       
				#check if homozygous
				check=[0 if (gt.split('/')[1]==gt.split('/')[0] or '0' in gt) else 1 for gt in GTs]#might need to change this one. it requires no ambiguious SNPs
				if sum(check)==0:
					keep.append(snp)
	print(f"keep {len(keep)} SNPs that are homozygous or heterozygous with 0/* structure from {len(snps)} SNPs")
	return(keep)

#make bins to group snps
def group_into_bins(numbers, bin_size):
	min_val = min(numbers)
	max_val = max(numbers)
	
	# Adjust min_val to start from a multiple of bin_size
	start_bin = (min_val // bin_size) * bin_size
	
	# Create bins
	bins = {}
	current_bin_start = start_bin
	
	while current_bin_start <= max_val:
		current_bin_end = current_bin_start + bin_size - 1
		bin_key = f"{current_bin_start}-{current_bin_end}"
		bins[bin_key] = []
		current_bin_start += bin_size
	
	# Assign numbers to bins
	for num in numbers:
		# Calculate which bin this number belongs to
		bin_index = (num - start_bin) // bin_size
		bin_start = start_bin + (bin_index * bin_size)
		bin_end = bin_start + bin_size - 1
		
		# Handle edge case where num equals the upper bound
		if num == bin_end + 1:
			bin_start += bin_size
			bin_end += bin_size
		
		bin_key = f"{bin_start}-{bin_end}"
		
		# Add to bin if it exists (should always exist)
		if bin_key in bins:
			bins[bin_key].append(num)
	
	return bins

def make_bins(snps_info,keep,window=1000,minlen=0):
	chroms,chrom_bins={},{}
	for snp,v in snps_info.items():
		if snp in keep:
			if v['#CHROM'] not in chroms:
				chroms[v['#CHROM']]={}
				chroms[v['#CHROM']][int(v['POS'])]=snp
			else:
				chroms[v['#CHROM']][int(v['POS'])]=snp
	for chrom in chroms:
		numbers=list(chroms[chrom].keys())
		bins=group_into_bins(numbers, window)
		chrom_bins[chrom]={k:{chroms[chrom][p]:p for p in v} for k,v in bins.items() if len(v)>minlen} #remove bins that are too short
		print(f"making {len(chrom_bins[chrom])} bins for reference chromosome {chrom}")
	return(chrom_bins)

def summary_bins(chrom_bins):
	for chrom in chrom_bins:
		x=[int(i.split('-')[0]) for i in chrom_bins[chrom].keys()] #position
		y=[len(v) for v in chrom_bins[chrom].values()] #count of SNPs
		#print(x,y)
		fig, ax = plt.subplots(nrows=2, ncols=1,figsize=(60,40))
		ax[0].scatter(x,y)
		ax[0].set_title('SNP counts for 1000 bp windows along chromosome')
		for i in [50,100,150]:
		  ax[0].hlines(i,0,max(x),linestyles='dashed',color='k')
		ax[1].hist(y,bins=50) #distribution of SNP counts
		ax[1].set_title('number of windows that contain a certain number of SNPs')   
		plt.show()

def translate(snp,gts,vfa): #must design the function to be applied to series (snps is the same)
	output=[]
	for gt in gts:
		if '.' not in gt:
			allele=sum([int(i) for i in set(gt.split('/'))])
			#print(snp,gt,allele)
		else:
			allele='.'
		output.append(vfa[snp][allele])   
	#print(output,list(columns))
	row = pd.DataFrame([output])
	return(output)

#finally, make the fasta files from the vcfs
def make_fasta(df,chrom_bins,workdir,minlen,maxlen,acc):
	n,minlen,maxlen=1,int(minlen),int(maxlen)
	total=sum([sum([1 for snps in chrom_bins[chrom].values() if len(snps)>=minlen and len(snps)<=maxlen]) for chrom in chrom_bins ])
	print(f"making {total} fasta files as outputs")
	for chrom in chrom_bins:
		for Bin,snps in chrom_bins[chrom].items():
			if len(snps)>=minlen and len(snps)<=maxlen:
				tmp=df.loc[list(snps.keys())]
				records=[]
				for sample in tmp.columns:
					if sample != 'salix':
						fa=''.join(tmp[sample]).upper()
						records.append(SeqRecord(seq=Seq(fa),id=sample,name='',description=f'{chrom}_{Bin}'))
				SeqIO.write(records,f"{workdir}/batch{acc}_clip_{n}.fasta",'fasta')
				n+=1

#EXECUTE FUNCTIONS
def main(vcf,workdir,minlen=50,maxlen=100):
  print('start converting vcf to fasta')
  step=100000 #empirical decision: was 1000000 but way too slow
  vcfend=int(sbp.check_output(f"wc -l {vcf}|cut -f 1 -d ' '",shell=True).decode().strip('\n')) #the length of the file to parse all lines
  #step=100000
  #vcfend=100001
  for acc,batch in enumerate(range(0,vcfend,step)): #to avoid OOM: cannot read too much lines each time, split into batches
    print(f"processing batch{acc}: {batch}-{batch+step}")
    snps,snps_info,samples=parse_vcf(vcf,workdir,batch,step,acc)
    print('get all GTs')
    GTs=set([v2['GT'] for k,v in snps.items() for v2 in v.values() ])
    print(GTs)
    print('get translation table for GT')
    vfa=GT_translation_table(snps_info)
    #print(vfa)
    print('select snps to keep')
    keep=filter_snps(snps,samples,vfa)
    if len(keep)!=0: #if getting any output
      print('split the SNPs by position on the chromosome')
      chrom_bins=make_bins(snps_info,keep,window=1000)
      summary_bins(chrom_bins)
      #get the df summary of SNPs for all samples
      print('get all SNPs as df')
      df=pd.DataFrame.from_dict(samples,orient='columns')
      df=df.loc[keep] #take the good SNPs only
      #translate the GT to bases
      print('get all SNPs translated to nt')
      bps=df.apply(lambda x: translate(x.name,x,vfa), axis=1)
      bps = pd.DataFrame(bps.tolist(), index=bps.index)
      bps.columns = df.columns
      #output as fasta
      print(f'make fasta files for the snps, min seq len ={minlen}, max seq len={maxlen}')
      make_fasta(bps,chrom_bins,workdir,minlen,maxlen,acc)
##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='check references', \
	  usage = 'vcf_csv_fa.py --vcf <vcf> --workdir <outdir> --min_SNPs <100>')
  parser.add_argument('--vcf', help='vcf file to be processed', metavar='vcf')
  parser.add_argument('--workdir', help='output directory', metavar='workdir')
  parser.add_argument('--minlen', help='minimum SNP count per window for the clip to be recorded', metavar='minlen')
  parser.add_argument('--maxlen', help='maximum SNP count per window for the clip to be recorded', metavar='maxlen')
  options = parser.parse_args()
  
  main(options.vcf,options.workdir,options.minlen,options.maxlen)
