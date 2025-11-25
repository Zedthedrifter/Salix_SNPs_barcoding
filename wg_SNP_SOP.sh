#!/bin/bash
#SBATCH --job-name=VCF
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --array=0
#SBATCH --mem=8G

#40-47,54-58,70-73,75,86
#1-31,39-47,54-58,70-73,75,86 #no myrsinites
#other salix diploids:
#caprea(54-58),repens(70-73,75)

function display_run_time {

SECONDS=$1
FUNCTION=$2

hours=$((SECONDS / 3600))
minutes=$(( (SECONDS % 3600) / 60 ))
seconds=$((SECONDS % 60))
printf "%s -- Duration: %02d:%02d:%02d\n" $FUNCTION $hours $minutes $seconds
}

###################################################################################
#STEP 0: SUBSAMPLING
function subsampling {

INDIR=$1
OUTDIR=$2
SUBSAMPLE=$3

#PARSE
sub=$(grep "${prefix}_${SLURM_ARRAY_TASK_ID}," $SUBSAMPLE|cut -d ',' -f 3)
frac=$(grep "${prefix}_${SLURM_ARRAY_TASK_ID}," $SUBSAMPLE|cut -d ',' -f 4)
echo $sub $frac

#SUBSAMPLING
seqtk sample -s 10 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz $frac > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq
seqtk sample -s 10 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz $frac > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq

#COMPRESS
gzip ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq
gzip ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq
#rm -f ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq #no need
echo "DONE SUBSAMPLING ${prefix}_${SLURM_ARRAY_TASK_ID}"
}


#STEP1: Indexing
#ENV: snps
function chrom_index {

REFDIR=$1
FASTA=$2
IDX=$3

cd $REFDIR
bowtie2-build -f $FASTA $IDX
}

#STEP2: MAPPING
function chrom_map {

INDIR=$1
IDX=$2
OUTDIR=$3

echo 'RUNNING BOWTIE2 MAPPING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
bowtie2-align-s --wrapper basic-0 \
                -x $IDX \
                -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz \
                -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz \
                -p 20 \
                -N 1 \
                -L 20 \
                --threads 8 \
                --phred33 \
                --very-sensitive-local \
                --no-discordant \
                --no-mixed \
                --no-unal \
                --time \
                --rg-id ${prefix}_${SLURM_ARRAY_TASK_ID} \
                --rg SM:${prefix}_${SLURM_ARRAY_TASK_ID} \
                --rg PL:'ILLUMINA' |\
                samtools view -Sbh -F 4 -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam #all mapped reads
#number of all mapped reads
#samtools view -F 4 -c $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam > $log

#only paired mapped reads with q>30
samtools view $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam -Sbh -F 4 -f 3 -q 30 -@ 8 |samtools sort  -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

#echo 'RUNNING SAMTOOLS COVERAGE ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
#samtools coverage $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $SORTED/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt
}
function bwa_map {

INDIR=$1
REF=$2
OUTDIR=$3

echo 'RUNNING BWA MAPPING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
bwa mem $REF \
        ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz \
        ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz \
        -t 8 \
        -k 20 |\
        samtools view -Sbh -F 4 -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam

#only paired mapped reads with q>30
samtools view $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam -Sbh -F 4 -f 3 -q 30 -@ 8 |samtools sort  -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
#INDEX
echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

#echo 'RUNNING SAMTOOLS COVERAGE ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
#samtools coverage $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $SORTED/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt
}
#a patch
function index_bam {

OUTDIR=$1

echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

echo 'RUNNING SAMTOOLS COVERAGE ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools coverage $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $SORTED/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt

}

#need to do this for sample1-16 once. i can do this in a simple for loop though; array is so much faster
function fix_RG_bam {

OUTDIR=$1

RG='@RG\tID:'${prefix}_${SLURM_ARRAY_TASK_ID}'\tSM:'${prefix}_${SLURM_ARRAY_TASK_ID}'\tPL:ILLUMINA'

samtools addreplacerg -r $RG \
                      -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fixed.bam \
                      $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
}

#STEP3: SNP calling : 
function snp_calling {

bamlist=$1
OUTDIR=$2
REF=$3
NAME=$4

echo 'calling SNPs'
#ls $INDIR/${prefix}*sorted.bam > bamlist.txt #need to change back to sorted.bam

bcftools mpileup -Ou -f $REF --bam-list $bamlist --threads 80 \
                 --annotate INFO/AD,FORMAT/DP,FORMAT/AD | \
bcftools call -Ou -mv | \
bcftools filter -s LowQual -e 'QUAL<20 || INFO/DP<100' > $OUTDIR/${NAME}.flt1.vcf #combined depth across samples>100 or quality >20
#this is a soft filter, which does not remove any snps. you can remove them later in the df easily

#rm bamlist.txt

}

#STEP4: SNP filtering
function filter_vcf {

INDIR=$1
OUTDIR=$2
NAME=$3

# Prep: Filter SNPs using vcftools

vcftools --vcf $INDIR/${NAME}.flt1.vcf \
         --out $OUTDIR/${NAME}.flt2 \
         --recode --recode-INFO-all \
         --minQ 30 \
         --max-missing 0.95 \
         --min-meanDP 5 \
         --max-meanDP 30 \
         --min-alleles 2 \
         --max-alleles 2 \
         --hwe 0.05

#make a compressed and indexed version for bcftools view -r (region) later
bgzip $OUTDIR/${NAME}.flt2.recode.vcf $OUTDIR/${NAME}.flt2.vcf.gz
tabix -p vcf $OUTDIR/${NAME}.flt2.vcf.gz
#-min-meanDP 5 : at least 5 reads per sample, otherwise not well supported
#--max-meanDP: if coverage > 30 (normalized to 15), the site is likely to be paralog (or i should drop to 25? filtering is rather quick)
#I removed --mac 5, which requires allele count to be >=5. Allele count is simply the number of times that allele appears over all individuals at that site. this varies for sample size etc. 
#--max-alleles 2: no max allele
#--remove-indels: keep indels? you can remove it from the csv later
}

#
function extract_loci {

VCF=$1
OUTDIR=$2
NAME=$3
min_freq=$4


#./calculate_snp.freq.py specifi_SNPs \
./species_specific_allele.py specifi_SNPs \
        -v $VCF \
        -n $NAME \
        --min_freq $min_freq \
        -o $OUTDIR #output directory

}
#

function rename_contigs {

INDIR=$1
OUTDIR=$2
CSV=$3

./rename_files.py rename_contig -i $INDIR --infile gene_${SLURM_ARRAY_TASK_ID}.fasta --fcsv $CSV -o $OUTDIR
}

function iqtree_per_gene {

INDIR=$1
OUTDIR=$2

mkdir $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
cp $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
#remove bad taxa (sample is bad) 'repens' and other all gap contigs
#./remove_taxa.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta #in this case 'repens' does not exist, just remove all gap is good enough 
#remove the outliers from an alignment
#./remove_highlyhetero.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -t 0.4 -m majority 
iqtree -s $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -bb 1000 -redo -safe
mv $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}*treefile $OUTDIR
rm -rf $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
}

function run_astral {

INDIR=$1
OUTDIR=$2
cvg=$3 #can resolve >= xx phylogeny

echo 'running astral'
./count_mono.py gene_select -i $INDIR -o $OUTDIR --cvg $cvg
#create astral phylogeny with selected genes
java -jar $astral -i $OUTDIR/rsl_${cvg}genes.in.treefile -o $OUTDIR/rsl_${cvg}genes.out.treefile 2>out.log

}

function run_astral_all_genes {

INDIR=$1
OUTDIR=$2
NAME=$3

echo 'running astral with all genes'
rm $OUTDIR/$NAME.in.treefile -f
cat $INDIR/*treefile > $OUTDIR/$NAME.in.treefile
#create astral phylogeny with selected genes
java -jar $astral -i $OUTDIR/$NAME.in.treefile -o $OUTDIR/$NAME.out.treefile 2>out.log
echo -e '\n\n\nDONE running astral with all genes\n\n\n'

}



#EXECUTION=================================================================================================================================

function main {

#UNIVERSAL VARIABLES=======================================================================================================================
env_name=snps_test #dnadiff is in captus.bowtie2: easy353; snps: after mapping
USER=zedchen 
WORKDIR=$SCRATCH/Barcoding_km/SNP_salix 
prefix=salix
#SLURM_ARRAY_TASK_ID=1

####constants############
DATA=$HOME/projects/rbge/pholling/barcoding/salix_combined/results/00_reads
CSV=/mnt/shared/projects/rbge/pholling/barcoding/salix_combined/results/00_reads/renamed.csv
bcftools=$HOME/apps/manual/bcftools-1.22/bcftools #somehow it's broken and i had to install bcftools from source
astral=$HOME/apps/manual/ASTRAL/astral.5.7.8.jar 
REF=$WORKDIR/refs
Scin=$REF/S_cinerea.chromosome.fasta
Sher=$REF/S_herbacea.chromosome.fasta
Scap=$REF/S_caprea.chromosome.fasta
Srep=$REF/S_repens.chromosome.fasta
Sret=$REF/S_reticulata.chromosome.fasta
bamlist2=/mnt/shared/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/bamlist.all_dip.CAP.txt
bamlist3=/mnt/shared/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/bamlist.all_dip.REP.txt
SUBSAMPLE=/home/zchen/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/subsampling.csv
min_freq=90
NAME=/mnt/shared/projects/rbge/pholling/barcoding/salix_combined/results/00_reads/sample_sp.csv
#==========================================================================================================================================


#subdirectories
REFR1=$WORKDIR/results/ref_01_mummer #REF RESULT
RESULT0=$WORKDIR/results/00_subsampling
RESULT1=$WORKDIR/results/01_sorted_bam #Skip the sam stage--> pipe to bam directly
RESULT2=$WORKDIR/results/02_VCFs
RESULT3=$WORKDIR/results/03_SNP_fasta
RESULT4=$WORKDIR/results/04_IQTREE
RESULT5=$WORKDIR/results/05_ASTRAL
RESULT=$WORKDIR/results/0

#USAGES
#=============================================
function setup_dir {

mkdir $WORKDIR
mkdir $WORKDIR/results
mkdir $WORKDIR/refs
mkdir $REFR1
mkdir $RESULT1
mkdir $RESULT0 -p
mkdir $RESULT1/CIN_bt2
mkdir $RESULT1/CAP_bt2
mkdir $RESULT1/REP_bt2
mkdir $RESULT1/CIN_bwa
mkdir $RESULT1/CAP_bwa
mkdir $RESULT1/REP_bwa
mkdir $RESULT2
mkdir $RESULT3/CAP -p
mkdir $RESULT3/REP -p
mkdir $RESULT4/CAP -p
mkdir $RESULT4/REP -p
mkdir $RESULT5/CAP -p
mkdir $RESULT5/REP -p
}

#SLURM_ARRAY_TASK_ID=1
setup_dir #WILL NOT OVERWRITE DIR, JUST LEAVE IT

#REFERENCE CURATION
#./process_ref.py $REF $REF
#chromosome_homolog: 10G mem
#chromosome_homolog $WORKDIR/refs  $WORKDIR/refs $REFR1 #

#STEP 0: SUBSAMPLING ARRAY, MEM=3G
#subsampling $DATA $RESULT0 $SUBSAMPLE

#===========================================================
#STEP_1: INDEXING: SINGLE
function index_ref {

chrom_index $REF $Scin CIN #polyploid
chrom_index $REF $Sher HER #not good, too different and hybridize a lot with others
chrom_index $REF $Scap CAP
chrom_index $REF $Sret RET
chrom_index $REF $Srep REP
#
bwa index $Scap
bwa index $Srep
}
#index_ref

#===========================================================
#STEP_2: env: snps_test; MAPPING: ARRAY 1-47, also add caprea (54-58) and repens (70-73,75)
#mem=9G is enough BEFORE SUBSAMPLING. none exceeded 8G when processing the big salix 47 files (must be less for other salix/potamogeton)
function master_map {

#BWA-MEM
bwa_map $RESULT0 $Scap $RESULT1/CAP_bwa
bwa_map $RESULT0 $Srep $RESULT1/REP_bwa
#BT2
#chrom_map $RESULT0 $REF/CAP $RESULT1/CAP_bt2 
#chrom_map $RESULT0 $REF/REP $RESULT1/REP_bt2 
}

#master_map

#index_bam #just to fix 1-13. there was a bug and they could not be indexed
#fix_RG_bam $RESULT1/CIN 

#===========================================================
#STEP_3: SNP CALLING: SINGLE, 36G
#CAPREA
SECONDS=0
#snp_calling $bamlist2 $RESULT2 $Scap CAP
display_run_time $SECONDS SNP_CALLING
#
SECONDS=0
#filter_vcf $RESULT2 $RESULT2 CAP
display_run_time $SECONDS SNP_FILTERING

#REPENS
SECONDS=0
#snp_calling $bamlist3 $RESULT2 $Srep REP
display_run_time $SECONDS SNP_CALLING
#
SECONDS=0
#filter_vcf $RESULT2 $RESULT2 REP
display_run_time $SECONDS SNP_FILTERING
#
#===========================================================
#STEP 4: EXTRACTING SPECIES SPECIFIC SNPS
SECONDS=0
#extract_loci $RESULT2/CAP.flt2.recode.vcf $RESULT2/CAP_ $NAME $min_freq
#extract_loci $RESULT2/REP.flt2.recode.vcf $RESULT2/REP_ $NAME $min_freq 
display_run_time $SECONDS species_specific_SNPs
#===========================================================

#STEP 5: FIND SNP ENRICHED REGIONS
#./species_specific_allele.py extract_ssSNP_regions -n $RESULT2/CAP_hq_specific_allele_GT_freq.csv -o $RESULT2/CAP_
#./species_specific_allele.py extract_ssSNP_regions -n $RESULT2/REP_hq_specific_allele_GT_freq.csv -o $RESULT2/REP_

#STEP_4: convert VCF to fasta: SINGLE
./species_specific_allele.py vcf_to_fasta --csv $RESULT2/CAP_high_resolution_windows.csv -n $NAME -v $RESULT2/CAP.flt2.vcf.gz -o $RESULT3/CAP
#./species_specific_allele.py vcf_to_fasta --csv $RESULT2/REP_high_resolution_windows.csv -n $NAME -v $RESULT2/REP.flt2.vcf.gz -o $RESULT3/REP

#STEP_5: IQTREE PER GENE:gene array: env=easy353 MEM=1G
#iqtree_per_gene $RESULT3/CAP $RESULT4/CAP #0-54
#iqtree_per_gene $RESULT3/REP $RESULT4/REP #0-71

#STEP_6: ASTRAL TREE: Run as single
function RUN_ASTRAL {

INDIR=$1
OUTDIR=$2

./count_mono.py process_treefiles -i $INDIR -o $OUTDIR
run_astral_all_genes $INDIR $OUTDIR SNPs_all_genes
run_astral $INDIR $OUTDIR 1 #get at least 2 genes for each mono taxa
run_astral $INDIR $OUTDIR 2 #get at least 2 genes for each mono taxa
}

#RUN_ASTRAL $RESULT4/CAP $RESULT5/CAP
#RUN_ASTRAL $RESULT4/REP $RESULT5/REP



#NO LONGER USEFUL STUFF
#collect the chromosome info into  csv
function extract_info {

INDIR=$1

rm $INDIR/clip_positions.csv -f
for fa in $(ls $INDIR/*fasta); do a=$(head $fa -n1|cut -d ' ' -f 2) ; b=$(echo $fa|rev|cut -d '/' -f 1|rev); echo $b , $a >> $INDIR/clip_positions.csv; done

}

#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN160 --minlen 160 --maxlen 200 #
#extract_info $RESULT3/CIN160
#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN150 --minlen 150 --maxlen 200 #634, probably the best
#extract_info $RESULT3/CIN150
#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN140_160 --minlen 140 --maxlen 160 #way too many clips selected if you set a wide range
#extract_info $RESULT3/CIN140_160
#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN40_50 --minlen 40 --maxlen 50 #try SNPs on a more conserved region
#extract_info $RESULT3/CIN40_50

#STEP_5: RENAME
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN160 -o $RESULT3/CIN160_renamed
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN150 -o $RESULT3/CIN150_renamed
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN140_160 -o $RESULT3/CIN140_160_renamed
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN40_50 -o $RESULT3/CIN40_50_renamed

#STEP_6:IQTREE PER GENE:gene array: env=easy353 MEM=1G
#rename_contigs $RESULT3/CIN160_renamed $RESULT3/CIN160_renamed $CSV #1-308
#rename_contigs $RESULT3/CIN150_renamed $RESULT3/CIN150_renamed $CSV #1-825
#iqtree_per_gene $RESULT3/CIN150_renamed $RESULT4/CIN150 #1-825 
#rename_contigs $RESULT3/CIN140_160_renamed $RESULT3/CIN140_160_renamed $CSV #1-978
#iqtree_per_gene $RESULT3/CIN140_160_renamed $RESULT4/CIN140_160 #1-978

#STEP_7: ASTRAL TREE: Run as single
#./count_mono.py process_treefiles -i $RESULT4/CIN150 -o $RESULT5/CIN150
#./count_mono.py process_treefiles -i $RESULT4/CIN160 -o $RESULT5/CIN160
#./count_mono.py process_treefiles -i $RESULT4/CIN140_160 -o $RESULT5/CIN140_160
#use astral method to make phylogeny: how many high-resolution genes per taxa do we need?

#RUN_ASTRAL $RESULT4/CIN160 $RESULT5/CIN160
#RUN_ASTRAL $RESULT4/CIN140_160 $RESULT5/CIN140_160
#run_astral_all_genes $RESULT4/CIN140_160 $RESULT5/CIN140_160 SNPs_all_genes

#should i look for species specific genotypes?
#or just turn vcf to fasta --> chop --> filter --> phylogeny --> astral

}

main