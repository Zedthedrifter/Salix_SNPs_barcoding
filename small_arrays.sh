#!/bin/bash
#SBATCH --job-name="sarray" 
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=20G
#SBATCH --array=13
#SBATCH --cpus-per-task=1

function fastp_fix {

INDIR=$1
OUTDIR=$2

fastp \
  -i $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz \
  -I $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz \
  -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_R1_filtered.fq.gz \
  -O $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_R2_filtered.fq.gz \
  --correction \  # enable base correction for overlapped regions
  --detect_adapter_for_pe \  # auto-detect adapters
  -h $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_fastp_report.html
}

#===================================================================
function main {

#UNIVERSAL VARIABLES=======================================================================================================================
env_name=captus #dnadiff is in captus.bowtie2: easy353; snps: after mapping
USER=zedchen 
WORKDIR=$SCRATCH/Barcoding_km/SNP_salix  #
prefix=salix
#SLURM_ARRAY_TASK_ID=21

####constants############
DATA=$HOME/projects/rbge/pholling/barcoding/salix_combined/results/00_reads/
CSV=$DATA/renamed.csv
#bcftools=$HOME/apps/manual/bcftools-1.22/bcftools #somehow it's broken and i had to install bcftools from source
astral=$HOME/apps/manual/ASTRAL/astral.5.7.8.jar 
REF=$WORKDIR/refs
Ptri=$REF/P_trichoides_chromosome.fasta #~50% alignment rate
Ppus=$REF/Ppusillus.chromosome.fasta #

#==========================================================================================================================================


#subdirectories

RESULT1=$WORKDIR/results/01_sorted_bam #Skip the sam stage--> pipe to bam directly
RESULT2=$WORKDIR/results/02_VCFs
RESULT3=$WORKDIR/results/03_SNP_fasta
RESULT4=$WORKDIR/results/04_IQTREE
RESULT5=$WORKDIR/results/05_ASTRAL
RESULT=$WORKDIR/results/0

#USAGES
#=============================================
#read count for each sample, much faster to run in array
zcat $DATA/${prefix}_${SLURM_ARRAY_TASK_ID}_R1_filtered.fq.gz|egrep '^\+$' -c

#how many reads are mapped
#samtools view -F 4 -c $RESULT1/REP/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

#FIX UNPAIRED READS IN SALIX_13
#fastp_fix $DATA $DATA
}

main
