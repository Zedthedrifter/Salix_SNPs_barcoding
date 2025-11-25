#!/home/zchen/apps/env/easy353/bin/python3

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
  'rename_fqgz', #rename fqgz files to use array
  'rename_fasta_easy353', #rename fasta output from easy353 assembly to use array, each file is a gene
  'rename_fasta', #plastome assembly: extract contigs from a parent directory (each sample has its own assembly directory) and rename to the sample acc number
  'rename_contig', #rename the contigs from array name to actual species name based on a csv when generating the rename fastq.gz
  'gene_alignment_summary', #summary of the gene alignments: how many species, how many gaps in each alignment etc. a heatmap? 
  'easy353_salix_filter' #filtering out bad samples and select alignments with > 60% non-gaps for concatenation
]

#FUNCTIONS

#filtering for salix easy353
def easy353_salix_filter(indir,outdir):
  df=pd.read_csv(f"{indir}/alignment_summary.csv",index_col=0)
  summary_dict={}
  #filter out bad samples: Salix repens.
  df = df.drop(columns=[c for c in df.columns if 'repens' in c])
  #select for high quality alignments: at least % non gaps:
  for cut in [0,10,20,30,40,50,60,70,80,90]:
    f1=df[df.apply(lambda x:0 not in [1 if i>=cut else 0 for i in x], axis=1)] #selection for alignment longer than a certain cutoff
    summary_dict[cut]=list(f1.index)
    concat={k:'' for k in df.columns} #reset concatenation collector
    count=0
    #CONCATENATION
    for gene in summary_dict[cut]:
      seqs=SeqIO.parse(f"{indir}/{gene}fasta",'fasta') #somehow there is a . already?
      contigs={s.description.replace(' ','').replace('.','_').split('salix')[0]:s.seq for s in seqs if 'repens' not in s.description}
      for k,v in concat.items(): #add to each species
        concat[k]=v+contigs[k]
      count+=1
    print(f"added {count} orthogroups")
    lens=[len(v) for k,v in concat.items()]
    print('concatenated contig lenghts:\n',lens)
    records=[SeqRecord(v,id=k,name=k,description=f"easy353_{cut}%_alignment") for k,v in concat.items()]
    SeqIO.write(records,f"{outdir}/alignment_filter_{cut}.fasta",'fasta')
  with open(f'{outdir}/filtered_genes.csv','w') as f:
    w = csv.writer(f)
    w.writerows(summary_dict.items())
  
  #df.to_csv(f'{outdir}/alignment_summary.filtered.csv', index=True)  
#make summary of the alignments
def gene_alignment_summary(indir,outdir):
  #make a dataframe
  summary_dict={}
  genes=sbp.check_output(f"ls {indir}/*fasta",shell=True).decode().strip('\n').split('\n')
  for gene in genes:
    seqs=SeqIO.parse(gene,'fasta')
    gaps={s.description.replace(' ','').replace('.','_').split('salix')[0]:round(1-s.seq.count('-')/len(s.seq),2)*100 for s in seqs} #precentage of the alignment that's not gap
    summary_dict[gene.split('/')[-1].replace('fasta','')]=gaps
  df=pd.DataFrame.from_dict(summary_dict, orient='index')
  df=df.fillna(0)#showing %of none-gap bases in each gene alignment
  #make a heatmap showing the percentage of non-gap bases:
  fig, ax = plt.subplots(nrows=1, ncols=1)
  ax.pcolor(df)
  #ax.set_yticks(np.arange(0.5, len(df.index), 1), df.index)
  ax.set_yticks([],[])
  ax.set_xticks(np.arange(0.5, len(df.columns), 1), df.columns)
  ax.tick_params(axis='x', labelrotation=90)
  ax.set_xlabel('Salix samples')
  plt.show()
  df.to_csv(f'{outdir}/alignment_summary.csv', index=True)  #save CSV for future filtering
  

def rename_fqgz(indir,outdir,prefix):
  #make a dict for in and out files
  r1=sbp.check_output(f"ls {indir}/*_1.fq.gz",shell=True).decode().strip('\n').split('\n')
  r2=sbp.check_output(f"ls {indir}/*_2.fq.gz",shell=True).decode().strip('\n').split('\n')
  pairs=[(i,j) for i,j in zip(r1,r2)]
  test=[True if p[0].replace('_1.fq.gz','')==p[1].replace('_2.fq.gz','') else False for p in pairs ]
  print(test)
  if False not in test: #if all names are well matched
    print('all r1, r2 files are matched. Proceed to renaming')
    rename_dict={p:(f"{prefix}_{i+21}_1.fq.gz",f"{prefix}_{i+21}_2.fq.gz") for i,p in enumerate(pairs)} #starting from 21 instead of 1: 1-20 is done already
    for k,v in rename_dict.items(): #copy the files to new directory, renamed
      sbp.call(f"mv {k[0]} {outdir}/{v[0]}",shell=True)
      sbp.call(f"mv {k[1]} {outdir}/{v[1]}",shell=True)
    with open(f'{outdir}/renamed_fastq.csv','w') as f:
      w = csv.writer(f)
      w.writerows(rename_dict.items())
  else:
    print("detecting r1,r2 mismatches")

def rename_fasta_easy353(indir,outdir):
  #make a dict for in and out files
  prefix='gene'
  fas=sbp.check_output(f"ls {indir}/*fasta",shell=True).decode().strip('\n').split('\n')
  rename_dict={fa:f"{prefix}_{i+1}.fasta" for i,fa in enumerate(fas)}
  for k,v in rename_dict.items(): #copy the files to new directory, renamed
    sbp.call(f"mv {k} {outdir}/{v}",shell=True)
  with open(f'{outdir}/renamed.csv','w') as f:
    w = csv.writer(f)
    w.writerows(rename_dict.items())

#extract contigs from a parent directory (each sample has its own assembly directory) and rename to the sample acc number
def rename_fasta(indir,outdir,outfile):
  records=[]
  for dirc in sbp.check_output(f"ls {indir}",shell=True).decode().strip('\n').split('\n'):
    fa=sbp.check_output(f"ls {indir}/{dirc}/*graph1.1.path_sequence.fasta",shell=True).decode().strip('\n') #only get graph1.1.path_seq
    contig=SeqIO.read(fa,'fasta')
    records.append(SeqRecord(contig.seq,id=dirc,description='plastome',name=contig.name))
  SeqIO.write(records,f'{outdir}/{outfile}','fasta')

#rename contigs based on a csv index
def rename_contig(indir, outdir,infile,fcsv):
  df = pd.read_csv(fcsv, encoding="latin-1")
  d={k:v for k,v in zip(df['new index'],df['species'])}
  contigs=SeqIO.to_dict(SeqIO.parse(f"{indir}/{infile}",'fasta'))
  records=[SeqRecord(v.seq,id=d[k],name='',description='') for k,v in contigs.items()]
  SeqIO.write(records,f'{outdir}/{infile}','fasta')

def remove_taxa(infile,outfile,taxa):
  taxa=[sp.replace(' ','') for sp in taxa.split(',')]
  seqs=SeqIO.parse(infile,'fasta')
  records=[]
  for seq in seqs:
    hits=[1 if sp in seq.id else 0 for sp in taxa]
    if sum(hits)==0 and len(set(str(seq.seq))) != 1:
      records.append(SeqRecord(seq.seq,id=seq.id.replace('.','_'),name='',description=''))
  SeqIO.write(records,outfile,'fasta')
  print(taxa)

##################################################
def main():
  count=0
  #parser = argparse.ArgumentParser(description='rename files of a certain type within the directory to use slurm array later',  usage = 'rename_files.py -i <indir> -o <outdir> --prefix <pre>')
  parser = argparse.ArgumentParser(description="Python script with multiple executable functions")
  subparsers = parser.add_subparsers(dest='command', help='Available commands')
  
  # Common arguments that will be reused
  common_args = argparse.ArgumentParser(add_help=False)
  common_args.add_argument('-i','--indir', help='input directory', metavar='indir')
  common_args.add_argument('-o','--outdir', help='output directory', metavar='outdir')
  common_args.add_argument('--outfile', help='output fasta file name', metavar='outfile')
  common_args.add_argument('--fcsv', help='input csv file name', metavar='fcsv')
  common_args.add_argument('--infile', help='input fasta file name', metavar='infile')
  common_args.add_argument('--taxa', help='a list of taxa delimiter=,', metavar='taxa')
  
  #rename fqgz
  fqgz_parser=subparsers.add_parser('rename_fqgz', parents=[common_args],help='rename fq.gz files within the directory to use slurm array later', 
                                    usage = './rename_files.py rename_fqgz -i <indir> -o <outdir> --prefix <pre>')
  fqgz_parser.add_argument('--prefix', help='prefix', metavar='pre')
  fqgz_parser.set_defaults(func=rename_fqgz)
  
  #rename fasta from plastome assembly
  fasta_parser=subparsers.add_parser('rename_fasta', parents=[common_args],help='extract contigs from a parent directory and rename the contigs', 
                                      usage = './rename_files.py rename_fasta -i <indir> -o <outdir> --outfile <outfile>')
  fasta_parser.set_defaults(func=rename_fasta)
  
  #rename contigs in fasta file based on index wit csv
  rename_contig_parser=subparsers.add_parser('rename_contig', parents=[common_args], help='rename the contigs in a fasta file', 
                                      usage = './rename_files.py rename_contig -i <indir> --infile <indir> --csv <csv file> -o <outdir>')
  rename_contig_parser.set_defaults(func=rename_contig)
  
  #rename fasta of 353 gene assembly from easy353
  rename_fasta_parser=subparsers.add_parser('rename_fasta_easy353', parents=[common_args], help='rename the fasta gene files from easy353 assembly', 
                                      usage = './rename_files.py rename_fasta_easy353 -i <indir> -o <outdir>')
  rename_fasta_parser.set_defaults(func=rename_fasta_easy353)
  
  ##make summary of the alignments
  gene_alignment_summary_parser=subparsers.add_parser('gene_alignment_summary', parents=[common_args], help='make summary of the alignments', 
                                      usage = './rename_files.py gene_alignment_summary -i <indir> -o <outdir>')
  gene_alignment_summary_parser.set_defaults(func=gene_alignment_summary)
  
  ##make summary of the alignments
  easy353_salix_filter_parser=subparsers.add_parser('easy353_salix_filter', parents=[common_args], help='filter samples and gene alignments', 
                                      usage = './rename_files.py easy353_salix_filter -i <indir> -o <outdir>')
  easy353_salix_filter_parser.set_defaults(func=easy353_salix_filter)
  #remove unwanted taxa
  remove_taxa_parser=subparsers.add_parser('remove_taxa', parents=[common_args], help='remvoe records whose id contain selected taxa', 
                                      usage = './rename_files.py remove_taxa -infile <input fasta> -outfile <output fasta>')
  remove_taxa_parser.set_defaults(func=remove_taxa)
  
  
  #parse arguments
  args = parser.parse_args()

  if not args.command:
    parser.print_help()
    sys.exit(1)
  
  # Prepare common kwargs for the function call
  kwargs = {k: v for k, v in vars(args).items() 
              if k in ['infile', 'outfile', 'fcsv', 'genome','indir','outdir','surfix','species','taxa'] and v is not None}

  
  # Call the appropriate function
  args.func(**kwargs)
    

if __name__ == '__main__': 
  main()