#!/home/zchen/apps/env/easy353/bin/python3
# 2025.08.12
# Zed ZChen@rbge.org.uk
# Input: vcf file, with GT:DP available (default)
# Output: csv file with genotype frequency for each species

from __future__ import division
import sys
import gzip
import copy
import argparse
import pandas as pd
import subprocess as sbp
from Bio import SeqIO
from Bio import Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#get samples vs species
def parse_names(names):
  df=pd.read_csv(names,index_col=0)
  id_sp=df.to_dict(orient='index')
  #id_sp={l.split(',')[0].split('/')[-1].replace('.bam',''):l.split(',')[1].strip('\n') for l in open(names)}
  return(id_sp)

#just one line of vcf file
def parse_single_record(l,samples,heading):
  snps={s:{} for s in samples}
  #COLLECT SNP INFO AND OUTPUT IN A CSV
  snps_info={i:j for i,j in zip(heading[:9],l.strip().split()[:9])}
  #COLLECT SAMPLE INFO
  l = l.strip().split()[9:] #the sample starts from teh 9th column
  genotypes=[gt.split(':')[0] for gt in l]
  for s,gt in zip (samples,genotypes):
    snps[s]['GT']=gt
  return(snps,snps_info)

#calcualte GT per species
# A MEMORY EFFICIENT METHOD, READ AND PROCESS ONE LINE EACH TIME
def species_specific_alleles(id_sp,sps,samples,snps,snps_info,min_freq,hqal,hqalgt,first):
  snp_al,snp_GT,sp_GT={},{},{}
  for sp in sps: #process each species
    gts=[snps[sample]['GT'] for sample,spp in id_sp.items() if sp == spp['species'] and sample in samples]
    sp_GT[sp]=gts#store the genotype 
    alleles=[]
    for gt in gts:
      alleles+=[int(i) for i in gt.split('/') if i != '.'] #remove the uncertain ones and collect all alleles
    #calculate allele frequency of the sample
    snp_al[sp]={al:round(alleles.count(al)/len(alleles)*100) for al in set(alleles)} #store the value
    #calculate GT frequency
    gts=[snps[sample]['GT'] for sample,spp in id_sp.items() if sp == spp['species'] and sample in samples]
    snp_GT[sp]={gt:round(100*gts.count(gt)/len(gts),2) for gt in set(gts)} #the dictionary with all GT frequencies for all SNPs
  tmp1={f"{species}_{al}": v.get(al,0) for species,v in snp_al.items() for al in [0,1]}
  tmp2={f"{species}_{gt}": v.get(gt,0) for species,v in snp_GT.items() for gt in ['0/0','0/1','1/1']}
  #print(snp_GT[k])
  #filter for SNPs with HQ species specific allele
  for al,typ in {0:'REF',1:'ALT'}.items(): #assume two alleles
    freqs={sp:snp_al[sp].get(al,0) for sp in snp_al}
    present=[1 for f in freqs.values() if f > 0 ]
    if sum(present)==1 : #the allele is present in only one species
      ssallele,target={},[s for s,f in freqs.items() if f >0][0]
      if target not in ['S.myrsinites','S.purpurea']: #remove polyploid S.myrsinites and single sample purpurea (decoy)
        ssallele['target sp']=target
        ssallele['target allele']=snps_info[typ]
        ssallele['target freq%']=freqs[target]
        
        #second filter: if the allele is present in all samples of the species:
        present_in_all=[1 if str(al) in gt else 0 for gt in sp_GT[target]]
        presence_freq=round(sum(present_in_all)/len(present_in_all)*100,2)
        ssallele['sample freq%']=presence_freq
        
        ssallele={**ssallele,**snps_info} #add snps info
        ssallele={**ssallele,**tmp1} #add allele freq for each species
        #print(k,al, sum(freqs),high) #make the first dictionary for collecting values
        if presence_freq>=int(min_freq):
          first+=1
          print(first,ssallele['target sp'],ssallele['target freq%'],ssallele['POS'],ssallele['#CHROM'])
          #print(tmp2)
          hqallele=ssallele #record if hq
          hqalGT={**ssallele,**tmp2} #add GT freq info
          
          #output results
          if first==1: #add heading once
            hqal.write(','.join(hqallele.keys()).replace(' ','_'))
            hqal.write('\n')
            hqalgt.write(','.join(hqalGT.keys()).replace(' ','_'))
            hqalgt.write('\n')
          #record the results
          hqal.write(','.join([str(v).replace(',',';') for v in hqallele.values()])) #make sure there's not comma to disrupt the csv format
          hqal.write('\n')
          hqalgt.write(','.join([str(v).replace(',',';') for v in hqalGT.values()]))
          hqalgt.write('\n')
          continue #no need to check the other allele
      else:
        continue
  return first
      
#CHECK WITH SLIDE WINDOW: WHICH WINDOW CONTAINS SNPS OF ALL/MANY SPECIES
def slide_window_check(chrm,pos_dict,ws=1000,step=100):
  chrm_dict=pos_dict[chrm]
  all_pos=list(chrm_dict.keys())
  bins=[(w_start,w_start+ws) for i,w_start in enumerate(range(int(min(all_pos)-step/2),int(max(all_pos)+step/2),step))]
  bin_dict={}
  for pos,sp in chrm_dict.items():
    for bin in bins:
      if pos >=bin[0]:
        if  pos < bin[1]:
          if bin not in bin_dict:
            bin_dict[bin]={pos:sp}
          else:
            bin_dict[bin][pos]=sp
      else:
        break
  stats=[len(set(v.values())) for k,v in bin_dict.items()]
  sps=set(chrm_dict.values())
  print(sps)
  print(f'{chrm} \nstep size={step} \nfound {len(bin_dict)} 1kb windows with ss SNPs \nmax={max(stats)} species per window in {stats.count(max(stats))} windows')
  #select windows with maximum species counts
  #or use species number -2 or -1? should include more windows if necessary
  selected={k:{sp:list(v.values()).count(sp) for sp in sps} for k,v in bin_dict.items() if len(set(v.values())) == max(stats)}
  for k,v in selected.items():
    print(k,v)
  print('######################################')
  #unoverlapping: #NEED TO REWRITE THE LOGIC
  keep,tested={},[]
  for k1,v1 in selected.items():
    if k1 not in tested:
      tested.append(k1)
      tmp={'k':k1,'v':v1}
      for k2,v2 in {k:v for k,v in selected.items() if k not in tested}.items():
        if tmp['k'][0]<k2[0] and tmp['k'][1]>k2[0]: #if the two windows overlaps
          tested.append(k2)
          #CHECK THE PRESENCE OF EACH SPECIES IN THE WINDOW
          print(tmp['k'],list(tmp['v'].values()).count(1),k2,list(v2.values()).count(1))
          if list(tmp['v'].values()).count(1) <= list(v2.values()).count(1): #IF LESS SPECIES ARE PRESENT WITH ONLY ONE SNP
            #CHECK THE NUMBER OF EACH SPECIES
            score=[1 if tmp['v'][k] > v2[k] else -1 if tmp['v'][k] < v2[k] else 0 for k in tmp['v']]
            if sum(score) >=0: #if the original one is better, keep it unchanged
              tmp['k']=tmp['k']
              tmp['v']=tmp['v']
            else:
              tmp['k']=k2
              tmp['v']=v2
          else:
            tmp['k']=k2
            tmp['v']=v2
      #collect the final winner
      keep[tmp['k']]=tmp['v']
  print('WINDOWS TO KEEP:')
  for k in keep:
    print(k)
  keep=[{**{'chrm':chrm,'start':k[0], 'end':k[1]},**v} for k,v in keep.items()] #unfold the nested dictionary and add chromosome and position information
  print('######################################')
  return(keep)

#============================================nested function for PARENT FUNCTION 3: EXTRACT SNPs from species specific SNPs enriched regions
#make a translation table for GT
def GT_translation_table(v):
  vfa={}
  if 'INDEL' in v['INFO']:
    total=max([len(v['REF'])]+[len(seq) for seq in v['ALT'].split(',')])
    vfa[0]=v['REF']+'-'*(total-len(v['REF']))
    for i,b in enumerate(v['ALT'].split(',')):
      vfa[i+1]=b+'-'*(total-len(b))
  else:
    vfa[0]=v['REF']
    for i,b in enumerate(v['ALT'].split(',')):
      vfa[i+1]=b
	#elif len(v['ALT'])>1:
	#    print(v)
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


def extract_snp_fa(gene,record,vcf,samples,heading,id_sp):
  print(f"PROCESSING {gene}")
  SNPs=sbp.check_output(f"bcftools view {vcf} -r 'ENA|{record['chrm']}|{record['chrm']}.1':{record['start']}-{record['end']}|grep -v '#'",shell=True).decode().strip('\n').split('\n')
  contigs={sample:'' for sample in id_sp}
  count=0
  #PROCESS SINGLE SNPS: memory efficient
  for i,l in enumerate(SNPs):
    snps,snps_info=parse_single_record(l,samples,heading)
    vfa=GT_translation_table(snps_info)
    #for each sample
    for s, v in snps.items():
      gt=v['GT'].replace('.','0') #TREAT ALL AMBIGUIOUS BASES AS IDENTICAL TO REF
      allele=sum([int(i) for i in set(gt.split('/'))])
      base=vfa[allele]
      contigs[s]+=base
    count+=len(base)
  #remove empty ones
  contigs={k:v for k,v in contigs.items() if len(v) ==count}
  print(f"PROCESSING {len(SNPs)} SNPs, contig length={count} in {len(contigs)} samples")
  return(contigs)

#DONE WITH THE TRAINING DATASET######################################################################################################################################################################################
#ONTO THE QUERY SET
#EXTRACT GOOD SNPS
def extract_SNPs(snps,snps_info,snp_out):
  #reset the keys for extraction
  print('number of SNPs before renaming keys',len(snps),len(snps_info))
  snps={f"{v['#CHROM']}_{v['POS']}":snps[k] for k,v in snps_info.items()}
  snps_info={f"{v['#CHROM']}_{v['POS']}":v for k,v in snps_info.items()}
  print('number of SNPs after renaming keys',len(snps),len(snps_info))
  #extract
  print('number of species specific SNPs',len(snp_out)) 
  snps={k:snps.get(f"{v['#CHROM']}_{v['POS']}",'NA') for k,v in snp_out.items()}
  snps_info={k:snps_info.get(f"{v['#CHROM']}_{v['POS']}",'NA') for k,v in snp_out.items()}
  print('number of SNPs after extraction',len(snps),len(snps_info))
  return(snps,snps_info)

#SHOW SAMPLE GT
def specific_SNPs_GT(samples,snps,snps_info,snp_out,output):
  GTs={}
  for k,v in snp_out.items():
    if snps_info[k]=='NA': #if the snp is not called in those samples
      print(f'{k} not called in all samples: identical to reference')
    else:
      #print(snps_info[k],v)
      if snps_info[k]['ALT']==v['ALT']:
        GTs[k]={s:snps[k][s]['GT'] for s in samples}
      else:
        for s in snps[k]:
          if snps[k][s]['GT']=='0/0':
            GTs[k][s]='0/0'
          else:
            GTs[k]='novo'
  df=pd.DataFrame.from_dict(GTs,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_GT.csv', index=True)  
  #convert GTs to list of alleles
  allele_table={0:'REF',1:'ALT'}
  alleles={snp:{} for snp in GTs}
  for snp,v in GTs.items():
    alleles[snp]={sample:[snps_info[snp][allele_table[int(i)]] for i in gt.split('/')] for sample,gt in v.items()}
    alleles[snp]['target sp']=snp_out[snp]['target sp']
    alleles[snp]['target allele']=snp_out[snp]['target allele']
  df=pd.DataFrame.from_dict(alleles,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_allele.csv', index=True)  
  return(GTs,alleles)

#TRANSLATE SAMPLE GT TO SPECIES
def GT_to_species(alleles,snp_out,output):
  species={}
  for k,v in alleles.items():
    species[k]={sample:snp_out[k]['target sp'] if snp_out[k]['target allele'] in al else 'NA' for sample,al in v.items()}
    species[k]['target sp']=snp_out[k]['target sp']
    species[k]['target allele']=snp_out[k]['target allele']
  df=pd.DataFrame.from_dict(species,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_allele_sample_species.csv', index=True)  
  return(species)

def output_all_snp_fasta(snps,snps_info,snpsfa):
  longer=max([len(snps_info['ALT']),len(snps_info['REF'])])
  REF=snps_info['REF']+'-'*(longer-len(snps_info['REF']))
  ALT=snps_info['ALT']+'-'*(longer-len(snps_info['ALT']))
  print(REF,ALT)
  for s in snpsfa:
    if '1' in snps[s]['GT']:
      snpsfa[s]+=ALT
    else:
      snpsfa[s]+=REF
#
################################################################################################################################################################################
#
#PARENT FUNCTION 1 : mem efficient loop version
#using samples with known species, identify species specific SNPs
def specifi_SNPs(vcf,names,min_freq,output):
  id_sp=parse_names(names) #sample names
  ssallele,hqallele,hqalGT={},{},{} #data collector
  hqal=open(f'{output}hq_specific_allele_freq.csv','w')
  hqalgt=open(f'{output}hq_specific_allele_GT_freq.csv','w')
  if vcf.endswith(".gz"):
  	opener = gzip.open
  else:
  	opener = open
  #process each line
  first=0
  for i in opener(vcf,'rt'):
    if i.startswith('#CHROM'):
      heading = i.strip().split()
      samples=heading[9:]
      samples=[s.split('/')[-1].replace('.bam','').replace('.sorted','') for s in samples]
      print(f"processing {len(samples)} samples: {samples}")
      sps=list(set([sp['species'] for sample,sp in id_sp.items() if sample in samples])) #only consider species present in the sample. the sample-species list can be longer/contain more species than the one analysed
      snpsfa={sp:'' for sp in samples}
      #print(snpsfa)
      sp_nb=len(sps)
    else:
      if not i.startswith('#'): 
        snps,snps_info=parse_single_record(i,samples,heading) #parse one line
        #put all snps into a fasta file
        first=species_specific_alleles(id_sp,sps,samples,snps,snps_info,int(min_freq),hqal,hqalgt,first) #output all species specific snps
  print('Complete calling species specific SNPs from training dataset')
  
  #output all SNPs as fasta
  #records=[SeqRecord(Seq(seq),id=s,name='',description='') for s,seq in snpsfa.items()]
  #SeqIO.write(records,f'{output}/SNPs_aligned.fasta','fasta')
  #output all HQ records
  #df=pd.DataFrame.from_dict(hqallele,orient='index') #turn to dataframe 
  #print(df)
  #df.to_csv(f'{output}/hq_specific_allele_freq.csv', index=True)  
  #df=pd.DataFrame.from_dict(hqalGT,orient='index') #output allele as well as GT freq
  #print(df)
  #df.to_csv(f'{output}/hq_specific_allele_GT_freq.csv', index=True)  

#MAKE DICT OF SPECIES SPECIFIC SNP POSITIONS
def parse_ssSNP_csv(csv):
  pos=sbp.check_output(f"cat {csv}|cut -d ',' -f 1,5,6",shell=True).decode().strip('\n').split('\n')
  pos_dict,keys={},['species','pos']
  for l in pos[1:]: #23844
    l=l.split(',')
    chrm=l[1].split('|')[1]
    if chrm not in pos_dict:
      pos_dict[chrm]={int(l[2]):l[0]}
    else:
      pos_dict[chrm][int(l[2])]=l[0]
  return(pos_dict)
#==========================================================
#PARENT FUNCTION 2: EXTRACT REGION WITH SPECIES SPECIFIC SNPS
def extract_ssSNP_regions(names,output):
  #get chromosome positions
  pos_dict=parse_ssSNP_csv(names)
  high_rsl_wds=[]
  n=0
  #check number of species per window
  for chrm in pos_dict: #ITERATE OVER EACH CHROMOSOME
    if n<100:
      keep=slide_window_check(chrm,pos_dict,ws=1000,step=100)
      n+=1
      high_rsl_wds+=keep
  #find regions enriched with ss SNPs
  high_rsl_wds={f"gene_{i}":v for i,v in enumerate(high_rsl_wds)} #add index for gene array later
  df=pd.DataFrame.from_dict(high_rsl_wds,orient='index')
  print(df)
  df.to_csv(f'{output}high_resolution_windows.csv', index=True) 
  print('Complete EXTRACTING REGION WITH SPECIES SPECIFIC SNPS')

#==========================================================
#PARENT FUNCTION 3: convert vcf of selected region to fasta
def vcf_to_fasta(vcf,csv,names,output):
  df=pd.read_csv(csv,index_col=0)
  pos_dict=df.to_dict(orient='index')
  id_sp=parse_names(names) #sample names
  print(df)
  #get the heading and samples for the vcf
  if vcf.endswith(".gz"):
  	opener = gzip.open
  else:
  	opener = open
  for i in opener(vcf,'rt'):
    if i.startswith('#CHROM'):
      heading = i.strip().split()
      samples=heading[9:]
      samples=[s.split('/')[-1].replace('.bam','').replace('.sorted','') for s in samples]
      print(f"processing {len(samples)} samples: {samples}")
      sps=list(set([sp['species'] for sample,sp in id_sp.items() if sample in samples])) #only consider species present in the sample. the sample-species list can be longer/contain more species than the one analysed
      snpsfa={sp:'' for sp in samples}
      sp_nb=len(sps)
      break
  #parse one WINDOW
  for k,v in pos_dict.items():
    contigs=extract_snp_fa(k,v,vcf,samples,heading,id_sp)
    records=[SeqRecord(seq=Seq(v),id=id_sp[k]['names'],name='',description='') for k,v in contigs.items()] #RENAME THE CONTIGS to species_idx
    SeqIO.write(records,f"{output}/{k}.fasta",'fasta')    
  print('Complete CONVERTING VCF TO FASTA (SNPS ONLY)')

#==========================================================
#PARENT FUNCTION 4 #haven't converted to loop version
#for samples without known species, find species/potential parents of hybrids based on species specific SNPs
def find_species(vcf,names,output):
  #read species specific SNPs from the previous step
  df=pd.read_csv(names,index_col=0)
  print(df)
  snp_out=df.to_dict(orient='index')
  #COLLECT SNPS
  snps,snps_info,samples=parse_vcf(vcf,output)
  #EXTRACT SPECIES SPECIFIC SNPS (IF EXIST)
  snps,snps_info=extract_SNPs(snps,snps_info,snp_out)
  GTs,alleles=specific_SNPs_GT(samples,snps,snps_info,snp_out,output)
  species=GT_to_species(alleles,snp_out,output)



#==========================================================
#now making the execution function  
def main():
  count=0
  #parser = argparse.ArgumentParser(description='rename files of a certain type within the directory to use slurm array later',  usage = 'rename_files.py -i <indir> -o <outdir> --prefix <pre>')
  parser = argparse.ArgumentParser(description="Python script with multiple executable functions")
  subparsers = parser.add_subparsers(dest='command', help='Available commands')
  
  # Common arguments that will be reused
  common_args = argparse.ArgumentParser(add_help=False)
  common_args.add_argument('-v','--vcf', help='vcf file with GT and DP info',metavar='vcf')
  common_args.add_argument('--csv', help='csv file, high resolution region, species specific SNPs, etc',metavar='csv')
  common_args.add_argument('-n','--names', help='sample ID to species names corresponding file, csv: sample ID/file path, species',metavar='names')
  common_args.add_argument('-o','--output',help='path to output directory',metavar='output')
  common_args.add_argument('--min_freq', help='minimum allele presence frequency in samples',metavar='min_freq')
  common_args.add_argument('--rename', help='rename the sample to species_index',metavar='rename')
  
  #find species specific SNPs
  specifi_SNPs_parser=subparsers.add_parser('specifi_SNPs', parents=[common_args],help='find and output species specific SNPs', 
                                    usage = './species_specific_allele.py main_parent -v <vcf file> -n <ID_species.txt> --min_freq <minimum allele presence frequency in samples%> -o <output_dir/prefix> ')
  specifi_SNPs_parser.set_defaults(func=specifi_SNPs)
  
  #assign species to unclassified samples based on SNPs
  find_species_parser=subparsers.add_parser('find_species', parents=[common_args],help='find and output hits on identified species specific SNPs', 
                                    usage = './species_specific_allele.py find_species -v <vcf file> -n <species specific SNPs.csv> -o <output_dir/prefix> ')
  find_species_parser.set_defaults(func=find_species)
  
  #extract_ssSNP_regions
  extract_ssSNP_regions_parser=subparsers.add_parser('extract_ssSNP_regions', parents=[common_args],help='find and output hits on identified species specific SNPs', 
                                    usage = './species_specific_allele.py extract_ssSNP_regions -n <species specific SNPs.csv> -o <output_dir/prefix> ')
  extract_ssSNP_regions_parser.set_defaults(func=extract_ssSNP_regions)
  
  #convert vcf to fasta
  vcf_to_fasta_parser=subparsers.add_parser('vcf_to_fasta', parents=[common_args],help='convert vcf of selected regions to fasta', 
                                    usage = './species_specific_allele.py vcf_to_fasta --csv <high_resolution_windows.csv> -n <ID_species.txt> -v <vcf files> -o <output_dir/prefix> ')
  vcf_to_fasta_parser.set_defaults(func=vcf_to_fasta)
  
  
  #parse arguments
  args = parser.parse_args()

  if not args.command:
    parser.print_help()
    sys.exit(1)
  
  # Prepare common kwargs for the function call
  kwargs = {k: v for k, v in vars(args).items() 
              if k in ['vcf', 'names', 'output','min_freq', 'csv', 'rename'] and v is not None}
   
  # Call the appropriate function
  args.func(**kwargs)
    

if __name__ == '__main__': 
  main()
