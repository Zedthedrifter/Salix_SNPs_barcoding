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

#list of functions
func_list=[
  'process_tree',
  'process_treefiles',
  'is_monophyletic',
  'gene_summary'
]

#check each taxa
def is_monophyletic(tree, prefix):
  target_tips = [tip.name for tip in tree.get_terminals() if tip.name.startswith(prefix)]
  if not target_tips:
    return False
    # Find the MRCA of these tips
  try:
    mrca = tree.common_ancestor(target_tips)
  except:
    # This can happen if some tips are not found (shouldn't happen with our filtering)
    return False
  clade_tips = [tip.name for tip in mrca.get_terminals()]
  return all(tip.startswith(prefix) for tip in clade_tips) #check if all members start with the same prefix (if yes, monophyletic)

#check each tree
def process_tree(treefile,mono_dict):
  try:
    tree = Phylo.read(treefile, 'newick')
  except:
    print(f"Error reading tree file: {treefile}")
    return 0
  gene=treefile.split('/')[-1].replace('.fasta.treefile','')
  tip_names = [tip.name for tip in tree.get_terminals() if tip.name]
  #pres=set([''.join([i for i in name if i not in ['0','1','2','3','4','5','6','7','8','9']]) for name in tip_names])
  pres=[] #a better way to fix the name issues
  for name in tip_names: 
    for i in [0,1,2,3,4,5,6,7,8,9]:
      if str(i) in name:
        name=name.split(str(i))[0].strip(str(i)).strip('_')
    pres.append(name)
  pres=set(pres)
  for prefix in pres:
    if is_monophyletic(tree, prefix):
      mono_dict[gene][prefix]=1
    else:
      mono_dict[gene][prefix]=0

#iterate the whole directory
def process_treefiles(indir,outdir):
  treefiles=sbp.check_output(f"ls {indir}/*treefile",shell=True).decode().strip('\n').split('\n')
  mono_dict={treefile.split('/')[-1].replace('.fasta.treefile',''):{} for treefile in treefiles}
  for treefile in treefiles:
    process_tree(treefile,mono_dict) #add value to key
  df=pd.DataFrame.from_dict(mono_dict, orient='index')
  df=df.fillna(0)
  df.to_csv(f'{outdir}/gene_resolution.csv', index=True) 
  print(df)

#summary the power to resolve monophyletic group from each gene
def gene_summary(indir,outdir):
  df=pd.read_csv(f'{outdir}/gene_resolution.csv',index_col=0)
  names=pd.read_csv(f'{outdir}/sample_names.csv',index_col=0)
  #row sum: the number of taxa each gene can resolve
  df['total_taxa']=df.sum(axis=1, skipna=True)
  df = df.sort_values('total_taxa', ascending=False)
  print(df['total_taxa'].value_counts()) #a summary of genes with different resolving capacity
  #select for genes with high resolving capacity
  for cut in [3,2]:
    try:
      sbp.call(f"mkdir {outdir}/gene_{cut}taxa",shell=True)
    except:
      print(f"{outdir}/gene_{cut}taxa already exist")
    #concatenate the genes
    concat={k:'' for k in names.index } #reset concatenation collector
    count=0
    #CONCATENATION
    for gene in df[df['total_taxa']>=cut].index:#select for genes with at least xx resolving power
      print(gene)
      contigs=SeqIO.to_dict(SeqIO.parse(f"{indir}/{gene}/{gene}.fasta",'fasta'))
      genlen=[len(v.seq)for k,v in contigs.items()][0] #len of teh aligned gene
      for k,v in concat.items(): #add to each species
        try:
          concat[k]=v+contigs[k].seq
        except:
          concat[k]=v+Seq('-'*genlen)#if a certain sample doesn't have that gene, just fill with gap of the same length
      count+=1
    print(f"added {count} orthogroups")
    lens=[len(v) for k,v in concat.items()]
    print('concatenated contig lenghts:\n',lens)
    #print(concat)
    records=[SeqRecord(v,id=k,name=k,description=f"genes that can resolve {cut} taxa") for k,v in concat.items() if len(set(str(v))) != 1] #remove all gap cases
    SeqIO.write(records,f"{outdir}/gene_{cut}taxa/gene_{cut}taxa.fasta",'fasta')
    sbp.call(f'iqtree -s {outdir}/gene_{cut}taxa/gene_{cut}taxa.fasta -bb 1000 -redo -safe',shell=True)
  
#find a set of gene that resolves the most monophyletic groups
def gene_select(indir,outdir,cvg):
  f=open(f"{outdir}/astral_gene_selection.log",'a')
  f.write(f"the set of genes ensured {cvg} genes that resolve monophyletic group per taxa, when possible\n")
  df=pd.read_csv(f'{outdir}/gene_resolution.csv',index_col=0)
  df['total'] = df.sum(axis=1)
  #names=pd.read_csv(f'{outdir}/sample_names.csv',index_col=0)
  #step 1: for each column, find all rows with none-0 entry at this column, sort them by row sum, decending order
  gene_d={taxa:{k:v for k,v in zip(df[df[taxa]!=0].index,df[df[taxa]!=0]['total'])} for taxa in df.columns if taxa!='total'}
  gene_d={t: sorted(v.keys(), key=lambda k: v[k], reverse=True) for t,v in gene_d.items()}
  #step 2: select at least the designated number of genes
  genes,cvg=[],int(cvg)
  for k,v in gene_d.items():
    if len(v)<cvg:
      print(f"{k} can be resolved in a monophyletic group by only {len(v)} genes")
      f.write(f"{k} can be resolved in a monophyletic group by only {len(v)} genes\n")
    genes+=v[:cvg]
  #genes=sum([v[:cvg] for k,v in gene_d.items()],[])
  #step 3: take only unique genes
  genes=list(set(genes))
  print(genes)
  #step 4: collect treefiles for the genes
  try:
    sbp.call(f"rm {outdir}/rsl_{cvg}_genes.in.treefile",shell=True)
  except:
    print(f"no old compiled treefile to be removed")
  for gene in genes:  
    sbp.call(f"cat {indir}/{gene}.fasta.treefile >> {outdir}/rsl_{cvg}_genes.in.treefile", shell=True)
  f.write(f"output treefiles of {len(genes)} genes: {genes} \nto {outdir}/rsl_{cvg}_genes.treefile\n")
  
def gene_manual(indir,outdir,genels):
  #step 1: get the gene list
  for l in open(genels):
    if l.startswith('#'):
      continue #you can comment out a gene set
    else:
      l=l.split('#')[0].strip('#') #can add comments to explain the choice
      genes=l.split(',')[1:] #format: set1,gene_1, gene_10, etc. 
      genes=[g.replace(' ','') for g in genes]
      pre=l.split(',')[0]
      print(pre,genes)
    try:
      sbp.call(f"rm {outdir}/manual_{pre}.in.treefile",shell=True)
    except:
      print(f"no old compiled treefile to be removed")
    #step2: collect the trees of the selected genes
    for gene in genes:
      sbp.call(f"cat {indir}/{gene}.fasta.treefile >> {outdir}/manual_{pre}.in.treefile", shell=True)
  f=open(f"{outdir}/astral_gene_selection.log",'a')
  f.write(f"output treefiles of \n{genes}\n to {outdir}/rsl_{pre}_genes.treefile\n\n")
  

#FUNCTIONS
def main():
  #parser = argparse.ArgumentParser(description='rename files of a certain type within the directory to use slurm array later',  usage = 'rename_files.py -i <indir> -o <outdir> --prefix <pre>')
  parser = argparse.ArgumentParser(description="Python script with multiple executable functions")
  subparsers = parser.add_subparsers(dest='command', help='Available commands')
  
  # Common arguments that will be reused
  common_args = argparse.ArgumentParser(add_help=False)
  common_args.add_argument('-i','--indir', help='input directory', metavar='indir')
  common_args.add_argument('-o','--outdir', help='output directory', metavar='outdir')
  
  #process_treefiles 
  fqgz_parser=subparsers.add_parser('process_treefiles', parents=[common_args],help='go through all treefiles to generate a csv file to show how many monophyletic groups can each gene resolve', 
                                    usage = './count_mono.py process_treefiles -i <indir> -o <outdir>')
  fqgz_parser.set_defaults(func=process_treefiles)
  
  #gene_summary
  fqgz_parser=subparsers.add_parser('gene_summary', parents=[common_args],help='concatenate the sequences of selected genes', 
                                    usage = './count_mono.py gene_summary -i <indir> -o <outdir>')
  fqgz_parser.set_defaults(func=gene_summary)
  #gene_select
  fqgz_parser=subparsers.add_parser('gene_select', parents=[common_args],help='select a set of genes that results in at least n gene that resolves each taxa monophyletically, when possible, and compile the treefiles', 
                                    usage = './count_mono.py gene_summary -i <indir> -o <outdir> --cvg <threshold>')
  fqgz_parser.add_argument('--cvg', help='how many genes', metavar='cvg')
  fqgz_parser.set_defaults(func=gene_select)
  #gene_manual
  fqgz_parser=subparsers.add_parser('gene_manual', parents=[common_args],help='select a set of genes manually, and compile the treefiles', 
                                    usage = './count_mono.py gene_manual -i <indir> -o <outdir> --genels <genels>')
  fqgz_parser.add_argument('--genels', help='select the set of genes', metavar='genels')
  fqgz_parser.set_defaults(func=gene_manual)
  
  #######
  #parse arguments
  args = parser.parse_args()

  if not args.command:
    parser.print_help()
    sys.exit(1)
  
  # Prepare common kwargs for the function call
  kwargs = {'indir': args.indir,
            'outdir': args.outdir}
  
  # Add function-specific arguments
  if args.command == 'gene_select':
    kwargs['cvg'] = args.cvg
  if args.command == 'gene_manual':
    kwargs['genels'] = args.genels
  
  # Call the appropriate function
  args.func(**kwargs)
##################################################
if __name__ == '__main__': 
  main()