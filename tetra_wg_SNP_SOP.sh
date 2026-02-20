#!/bin/bash
#SBATCH --job-name="IQTREE"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=6
#SBATCH --array=0
#SBATCH --mem=14G

#list of tetraploid samples:
#35,52,63,67,78,109,129,28,33,79,80,118,31,68,74,85,86,40,64,70,71,72,48,73,91,92,93,94,116,23,59,62,95,24,34,41,47,117,133,100,101,102,103
#accidentally put 115
#need to remove 77
#just did same sample renaming

function display_run_time {

SECONDS=$1
FUNCTION=$2

hours=$((SECONDS / 3600))
minutes=$(( (SECONDS % 3600) / 60 ))
seconds=$((SECONDS % 60))
printf "%s -- Duration: %02d:%02d:%02d\n" $FUNCTION $hours $minutes $seconds
}

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
                -N 0 \
                -L 20 \
                --threads 8 \
                --phred33 \
                --sensitive-local \
                --no-discordant \
                --no-unal \
                --time \
                --rg-id ${prefix}_${SLURM_ARRAY_TASK_ID} \
                --rg SM:${prefix}_${SLURM_ARRAY_TASK_ID} \
                --rg PL:'ILLUMINA' |\
                samtools view -Sbh -F 4 -f 3 -@ 8|\
                samtools sort -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
#-q 30 is way too stringent
# --no-mixed #not allowing one in pair to be unmapped is too stringent
echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

echo 'RUNNING SAMTOOLS COVERAGE ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools coverage $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $SORTED/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt
}

function bwa_map {

INDIR=$1
REF=$2
OUTDIR=$3

echo 'RUNNING BWA MAPPING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
echo "${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz"
echo "${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz"

bwa mem $REF \
        ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz \
        ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz \
        -t 4 \
        -k 20 |\
        samtools view -Sbh -F 4 -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam - #Note the - at the end tells samtools view to read from stdin.

#only paired mapped reads with q>30/TOO HIGH, JUST USE Q10
samtools view $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam -Sbh -F 4 -f 3 -q 10 -@ 8 |samtools sort  -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
#samtools view $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam -Sbh -F 4 -f 3 -q 10 -@ 8 |\
#samtools addreplacerg -r "ID:${prefix}_${SLURM_ARRAY_TASK_ID} SM:${prefix}_${SLURM_ARRAY_TASK_ID} LB:${prefix}_${SLURM_ARRAY_TASK_ID} PL:ILLUMINA" |\
#samtools sort  -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
#INDEX
#echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
}

#baw doesnot have heading by default
function samtools_add_heading {

INDIR=$1
OUTDIR=$2

samtools addreplacerg \
         -r "ID:${prefix}_${SLURM_ARRAY_TASK_ID}" \
         -r "SM:${prefix}_${SLURM_ARRAY_TASK_ID}" \
         -r "LB:${prefix}_${SLURM_ARRAY_TASK_ID}" \
         -r "PL:ILLUMINA" \
         -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.bam \
         $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam 
#INDEX
#echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.bam
}

function samtools_sort {

INDIR=$1

samtools sort -@ 8 -o ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.bam
samtools index ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
}

#CHECK THE IDENTITY OF UNMAPPED READS
function unmapped_check {

INDIR=$1
REF=$2
OUTDIR=$3

echo 'RUNNING BWA MAPPING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
echo "${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz"
echo "${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz"

#bwa mem $REF ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz  ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz  -t 8  -k 20 |samtools view -Sbh -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam

#COLLECT UNMAPPED READS
echo 'EXTRACTING UNMAPPED READS FROM' ${prefix} ${SLURM_ARRAY_TASK_ID}
#samtools view $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam -Sbh -f 4 -q 30 -@ 8 |
#samtools view $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.unmapped.bam -Sbh -f 4 -f 3 -@ 8| samtools fastq -1 $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_unmapped_R1.fq.gz -2 $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_unmapped_R2.fq.gz -c 9
#ASSEMBLY
spades.py --only-assembler -1 $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_unmapped_cleaned_R1.fq.gz -2 $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_unmapped_cleaned_R2.fq.gz -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID} --careful -t 8 -m32
#echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
}

#a patch
function index_bam {

OUTDIR=$1

echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

echo 'RUNNING SAMTOOLS COVERAGE ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools coverage $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $SORTED/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt

}



#STEP3: SNP calling : 

#ALTERNATIVE CALLING FUNCTION FOR TETRAPLOID SAMPLES
#too memory heavy. run in parallel by chromosomes
function freebayes_calling_2_pipe {

bamlist=$1
OUTDIR=$2
REF=$3
NAME=$4

echo "doing freebayes"

while read bam; do
  samplename=$(basename $bam .bam)
  echo -e "$samplename"
done < $bamlist_tri_tetra > samples.tsv

#set up max DP
SAMPLE_COUNT=$(wc -l samples.tsv|cut -f 1 -d ' ')
MAX_DP=$((SAMPLE_COUNT * 40))
#--samples samples.tsv\
#--bam-list $bamlist \

#get chromosome: ARRAY 1-14
CHROM=$(cat $REF|grep '>' | sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 1 -d ' ') #use index to select chrome
CHROM_NAME=$(echo $CHROM|cut -f 2 -d '|')
echo "calling SNPs on $CHROM_NAME"


freebayes -f $REF $(cat $bamlist) \
          --region ${CHROM/>/} \
          --ploidy 4 \
          --min-alternate-count 2 \
          --min-alternate-fraction 0.01 \
          --use-best-n-alleles 4 \
          --haplotype-length 0 \
          --min-mapping-quality 20 \
          --min-base-quality 20 \
          --min-supporting-allele-qsum 20 \
          --min-coverage 100 \
          --read-mismatch-limit 4 \
          --read-snp-limit 4 \
  |  bcftools view  -O u \
  | bcftools filter -e "QUAL<30 || INFO/DP > $MAX_DP || F_MISSING > 0.2 || MQM < 30 || MQMR < 30 || SAP > 30 || INFO/RPP > 30" \
  | bcftools filter -i 'INFO/AC[1] >= 4' \
  | bcftools filter -i 'COUNT(FMT/DP<5) <= 2 && COUNT(FMT/GQ<20) <= 2 && COUNT(FMT/DP>40) <= 2 ' -Oz -o $OUTDIR/${NAME}_${CHROM_NAME}.flt3.vcf.gz

bcftools index $OUTDIR/${NAME}_${CHROM_NAME}.flt3.vcf.gz #need index for consensus extraction
 
#Allele frequency bounds - remove very rare/very common SNPs

#Quality/Depth ratio - normalized quality score

#Genotype quality filter - GQ < 20

#Minor allele count - require at least 2 samples with alt allele
#SAP=Strand bias p-value (Phred-scaled)
#RPP = Read Placement Probability for alternate allele (Phred-scaled)
#INFO/AC[1] >= 4: include allele >=4 copies: assume 1 copy per sample, at least one copy to be present in all samples of one species)
}

#another solution: gatk4
#conda install gatk4 --name snps_test
function gatk_SNP_calling {

bamlist=$1
OUTDIR=$2
REF=$3
NAME=$4

    echo "Running GATK tetraploid whole-genome calling for $NAME"
    
    # Calculate DP threshold
    SAMPLE_COUNT=$(wc -l < "$bamlist")
    MAX_DP=$((SAMPLE_COUNT * 40))  # 43×60 = ~2580
    
    # Direct joint calling (no per-sample GVCFs)
    INPUT_BAMS=$(while read bam; do echo "-I $bam"; done < "$bamlist")
    
    gatk HaplotypeCaller \
        -R "$REF" \
        $INPUT_BAMS \
        -ploidy 4 \
        --min-base-quality-score 20 \
        --minimum-mapping-quality 20 \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20 \
        -O "${OUTDIR}/${NAME}.raw.vcf.gz"
    

}


#STEP4: SNP filtering
function filter_vcf {

INDIR=$1
OUTDIR=$2
NAME=$3

# Prep: Filter SNPs using vcftools
#expected depth: 30
vcftools --vcf $INDIR/${NAME}.flt1.vcf \
         --out $OUTDIR/${NAME}.flt2 \
         --recode --recode-INFO-all \
         --min-alleles 2 \
         --max-alleles 2 \
         --hwe 0.05

# Then use bcftools for advanced INDEL filtering
bcftools view $OUTDIR/${NAME}.tmp.recode.vcf | \
bcftools filter \
    -e 'TYPE="indel" && (IDV<5 || IMF<0.25 || VDB<0.1)' \
    -Oz -o $OUTDIR/${NAME}.flt2.vcf.gz
    
#index
bcftools index $OUTDIR/${NAME}.flt2.vcf.gz
 
#PER CHROMOSOME
#vcftools --vcf /mnt/shared/scratch/zchen/Barcoding_km/SNP_potamogeton/results/02_VCFs/PUS_q10.flt1.vcf --chr 'ENA|OZ286224|OZ286224.1' --out /mnt/shared/scratch/zchen/Barcoding_km/SNP_potamogeton/results/02_VCFs/PUS_q10.flt2.chr1 --recode --recode-INFO-all          --minQ 20          --max-missing 0.9          --min-meanDP 5          --max-meanDP 30          --min-alleles 2          --max-alleles 2          --hwe 0.05

#-minDP 5 : at least 5 reads per sample, otherwise not well supported
#--maxDP: if coverage > 20 (normalized to 15), the site is likely to be paralog (or i should drop to 25? filtering is rather quick)
#I removed --mac 5, which requires allele count to be >=5. Allele count is simply the number of times that allele appears over all individuals at that site. this varies for sample size etc. 
#--max-alleles 2: no max allele
#--remove-indels: keep indels? you can remove it from the csv later
}

function extract_loci {

VCF=$1
OUTDIR=$2
NAME=$3
min_freq=$4


#./calculate_snp.freq.py specifi_SNPs \
./species_specific_allele.tetra.py specifi_SNPs \
        -v $VCF \
        -n $NAME \
        --min_freq $min_freq \
        -o $OUTDIR #output directory

}
#
function extract_loci_by_chromosome {

INDIR=$1
OUTDIR=$2
NAME=$3
min_freq=$4
REF=$5
PREFIX=$6

#GET CHROMOSOME NAME
CHROM=$(cat $REF|grep '>' | sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 1 -d ' ') #use index to select chrome
CHROM_NAME=$(echo $CHROM|cut -f 2 -d '|')
echo "calling SNPs on $CHROM_NAME"

#./calculate_snp.freq.py specifi_SNPs \
./species_specific_allele.tetra.py specifi_SNPs \
        -v $INDIR/${PREFIX}_${CHROM_NAME}.flt3.vcf.gz \
        -n $NAME \
        --min_freq $min_freq \
        -o $OUTDIR/${PREFIX}_${CHROM_NAME}_ #output directory

}

#
function extract_window_by_chromosome {

DI_VCF=$1
INDIR=$2
OUTDIR=$3
PREFIX=$4
REF=$5

#GET CHROMOSOME NAME
CHROM=$(cat $REF|grep '>' | sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 1 -d ' ') #use index to select chrome
CHROM_NAME=$(echo $CHROM|cut -f 2 -d '|')
echo "calling SNPs on $CHROM_NAME"

#GET VCF FOR TETRA
VCF2=$INDIR/${PREFIX}_${CHROM_NAME}_hq_specific_allele_freq.csv

#COMBINE HQ SPECIFIC SNPS
cat $DI_VCF $VCF2|grep "${CHROM_NAME}"|cut -d ',' -f 1-6 > ${OUTDIR}/${PREFIX}_${CHROM_NAME}_combined_SNPs.csv

#./calculate_snp.freq.py specifi_SNPs \
./species_specific_allele.tetra.py extract_ssSNP_regions \
        -n ${OUTDIR}/${PREFIX}_${CHROM_NAME}_combined_SNPs.csv \
        -o ${OUTDIR}/${PREFIX}_${CHROM_NAME}

}
#
#MAKING CONSENSUS SEQUENCES FROM VCF

function extract_consensus {

VCF=$1
REF=$2
WINDOWS=$3
OUTDIR=$4
TETRA_INDIR=$5
PREFIX=$6
REF_SP=$7

# 1. SET UP PARAMETERS
GENE=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 1 -d ',')
NAME=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 2 -d ',') #THE CHROMOSOME
START=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 3 -d ',')
END=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 4 -d ',')
REGION="ENA|${NAME}|${NAME}.1:${START}-${END}"

# 2. EXTRACT REGION from reference
samtools faidx $REF "$REGION" > ${OUTDIR}/${GENE}.fasta
# Add sample name to header
sed -i "1s/^>.*/>${REF_SP}_${NAME}_${START}_${END}/" ${OUTDIR}/${GENE}.fasta

# 3. List all DIPLOID samples in the VCF#need to split this into two steps: diploids and tetraploids, then combine results
mkdir ${OUTDIR}/${GENE} #a temporary directory so the files won't mix and interfere with each other
SAMPLES=$(bcftools query -l $VCF)

# 4. Generate consensus for each DIPLOID sample
for SAMPLE in $SAMPLES; do
    TMP=$(echo $SAMPLE|rev|cut -f 1 -d '/'|rev);SAMPLE_NAME=${TMP/.sorted.bam/} #get the last field and remove .sort.bam
    echo "Processing $SAMPLE_NAME..."
    
    # Generate consensus for this sample
    samtools faidx $REF $REGION | bcftools consensus $VCF -s $SAMPLE -I -o ${OUTDIR}/${GENE}/tmp.${GENE}.${SAMPLE_NAME}.fasta
          
    # Add sample name to header
    sed -i "1s/^>.*/>${SAMPLE_NAME}/" ${OUTDIR}/${GENE}/tmp.${GENE}.${SAMPLE_NAME}.fasta
done

# 3.1: LIST ALL TETRA PLOID SAMPLES

TETRA_VCF=$TETRA_INDIR/${PREFIX}_${NAME}.flt3.vcf.gz
SAMPLES=$(bcftools query -l $TETRA_VCF) #TATRAPLOID SAMPLE LIST

# 4.1: Generate consensus for each TETRAPLOID sample
for SAMPLE in $SAMPLES; do
    TMP=$(echo $SAMPLE|rev|cut -f 1 -d '/'|rev);SAMPLE_NAME=${TMP/.sorted.bam/} #get the last field and remove .sort.bam
    echo "Processing $SAMPLE_NAME..."
    
    # Generate consensus for this sample
    samtools faidx $REF $REGION | bcftools consensus $TETRA_VCF -s $SAMPLE -I -o ${OUTDIR}/${GENE}/tmp.${GENE}.${SAMPLE_NAME}.fasta
          
    # Add sample name to header
    sed -i "1s/^>.*/>${SAMPLE_NAME}/" ${OUTDIR}/${GENE}/tmp.${GENE}.${SAMPLE_NAME}.fasta
done

# 5. Combine all into one multi-FASTA
cat ${OUTDIR}/${GENE}/tmp.${GENE}.*.fasta >> ${OUTDIR}/${GENE}.fasta

rm -rf ${OUTDIR}/${GENE}
}

#

function easy353_mafft {

INDIR=$1
OUTDIR=$2
CSV=$3

FASTA=$(ls $INDIR/*fasta|sed -n "${SLURM_ARRAY_TASK_ID}p"|rev|cut -f 1 -d '/'|rev)
GENE=${FASTA/.fasta/}
#RENAME AND GET A COPY
./rename_files.py rename_contig -i $INDIR --infile $FASTA --fcsv $CSV -o $OUTDIR
mv $OUTDIR/$FASTA $OUTDIR/tmp.$FASTA
mafft --maxiterate 10000 $OUTDIR/tmp.$FASTA > $OUTDIR/$FASTA
rm -f $OUTDIR/tmp.$FASTA
#remove any space that might cause trouble later
sed -i 's/ //g' $OUTDIR/$FASTA
}

function iqtree_per_gene {

INDIR=$1
OUTDIR=$2

FASTA=$(ls $INDIR/*fasta|sed -n "${SLURM_ARRAY_TASK_ID}p"|rev|cut -f 1 -d '/'|rev)
GENE=${FASTA/.fasta/}

mkdir $OUTDIR/$GENE
cp $INDIR/$FASTA $OUTDIR/$GENE
sed -i 's/>P\.trichoides_OZ[0-9]*_[0-9]*_[0-9]*/>P_trichoides_ref/' $OUTDIR/$GENE/$FASTA #rename the referenece in the treefile
#remove bad taxa (sample is bad) and all gap contigs
./remove_taxa.py $OUTDIR/$GENE/$FASTA #'P_pusillus_6','P_berchtoldii_5','P_berchtoldii_6','Undetermined','S_filiformis_4','S_filiformis_5','P_perfoliatus_6','S_filiformis_45','P_pusillus_2'
iqtree -s $OUTDIR/$GENE/$FASTA -bb 1000 -redo -safe
mv $OUTDIR/$GENE/*treefile $OUTDIR
rm -rf $OUTDIR/$GENE


}

function run_astral {

INDIR=$1
OUTDIR=$2
cvg=$3 #can resolve >= xx phylogeny

echo 'running astral'
rm -f $OUTDIR/rsl_${cvg}genes.in.treefile
./count_mono.py gene_select -i $INDIR -o $OUTDIR --cvg $cvg
#create astral phylogeny with selected genes
java -jar $astral -i $OUTDIR/rsl_${cvg}genes.in.treefile -o $OUTDIR/rsl_${cvg}genes.out.treefile 2>out.log

}

function run_astral_all_genes {

INDIR=$1
OUTDIR=$2
NAME=$3

echo 'running astral with all genes'
cat $INDIR/*treefile > $OUTDIR/$NAME.in.treefile
#create astral phylogeny with selected genes
java -jar $astral -i $OUTDIR/$NAME.in.treefile -o $OUTDIR/$NAME.out.treefile 2>out.log
echo -e '\n\n\nDONE running astral with all genes\n\n\n'

}

#EXECUTION=================================================================================================================================

function main {

#UNIVERSAL VARIABLES=======================================================================================================================
prefix=potamogeton
#SLURM_ARRAY_TASK_ID=21

####constants############
DATA=$HOME/projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/
CSV=$DATA/renamed.csv
SUBSAMPLE=/mnt/shared/projects/rbge/zedchen/barcoding/potamogeton_20250714/SNP_scripts/subsampling.csv
TETRASUBSAMPLE=/mnt/shared/projects/rbge/zedchen/barcoding/potamogeton_20250714/SNP_scripts/subsampling_tetraploids.csv
#bcftools=$HOME/apps/manual/bcftools-1.22/bcftools #somehow it's broken and i had to install bcftools from source
astral=$HOME/apps/manual/ASTRAL/astral.5.7.8.jar 
REF=$HOME/projects/rbge/zedchen/References/Pota_ref
Ptri=$REF/P_trichoides_chromosome.fasta #~50% alignment rate
Ppus=$REF/Ppusillus.chromosome.fasta #
bamlist1=/mnt/shared/projects/rbge/zedchen/barcoding/potamogeton_20250714/SNP_scripts/bamlist_TRI.txt
bamlist2=/mnt/shared/projects/rbge/zedchen/barcoding/potamogeton_20250714/SNP_scripts/bamlist_PUS.txt
bamlist_tri_tetra=/mnt/shared/projects/rbge/zedchen/barcoding/potamogeton_20250714/SNP_scripts/TRI_bamlist_tetra.fofn
exclude=/mnt/shared/projects/rbge/zedchen/barcoding/potamogeton_20250714/SNP_scripts/exclude.txt
diploid_specificSNPs=/mnt/shared/scratch/zchen/Barcoding_km/SNP_potamogeton/results//02_VCFs/TRI_hq_specific_allele_freq.csv #combined the diploid results for hq window finding
DI_vcf=/mnt/shared/scratch/zchen/Barcoding_km/SNP_potamogeton/results/02_VCFs/TRI_q10.flt3.vcf.gz
min_freq=90
NAME=/mnt/shared/projects/rbge/zedchen/barcoding/potamogeton_20250714/SNP_scripts/sample_sp.csv

#==========================================================================================================================================


#subdirectories
WORKDIR=$SCRATCH/Barcoding_km/SNP_potamogeton/TETRA #
DI_DIR=$SCRATCH/Barcoding_km/SNP_potamogeton/results
RESULT0=$WORKDIR/00_subsampling
RESULT1=$WORKDIR/01_sorted_bam #Skip the sam stage--> pipe to bam directly
RESULT2=$WORKDIR/02_VCFs
RESULT2_1=$WORKDIR/02_GATK_calls
RESULT3=$WORKDIR/03_HQ_windows
RESULT4=$WORKDIR/04_SNP_fasta
RESULT5=$WORKDIR/05_aligned_fasta
RESULT6=$WORKDIR/06_IQTREE
RESULT7=$WORKDIR/07_ASTRAL
RESULT=$WORKDIR/0

#USAGES
#=============================================
function setup_dir {

mkdir $WORKDIR -p
mkdir $RESULT0
mkdir $RESULT1
mkdir $RESULT1/TRI_bwa
mkdir $RESULT1/PUS_bwa
mkdir $RESULT2
mkdir $RESULT2_1
mkdir $RESULT3
mkdir $RESULT4
mkdir $RESULT5
mkdir $RESULT6
mkdir $RESULT7
}

setup_dir #WILL NOT OVERWRITE DIR, JUST LEAVE IT
#

#REFERENCE CURATION
#./process_ref.py $REF $REF

#===========================================================
#STEP 0: SUBSAMPLING ARRAY, MEM=3G, env=snps_test
SECONDS=0
#subsampling $DATA $RESULT0 $TETRASUBSAMPLE
display_run_time $SECONDS SUBSAMPLE

#===========================================================
#STEP_1: INDEXING: SINGLE
function master_index {

#samtools dict "$Ptri" -o "${Ptri%.*}.dict"
#chrom_index $REF $Ptri TRI
#chrom_index $REF $Ppus PUS
#INDEXING: BWA
bwa index $Ptri
bwa index $Ppus
}

SECONDS=0
#master_index
display_run_time $SECONDS SUBSAMPLE

#===========================================================
#STEP_2: MAPPING: ARRAY-SAMPLE 21-126
#mem=10G (according to mem report: <=4.5G when using half the data)
function master_map {

#BWA-MEM
bwa_map $RESULT0 $Ptri $RESULT1/TRI_bwa
#bwa_map $RESULT0 $Ppus $RESULT1/PUS_bwa
#BT2
#chrom_map $DATA $REF/TRI $RESULT1/TRI #~50% alignment rate
#chrom_map $DATA $REF/PUS $RESULT1/PUS #~50% alignment rate
}

SECONDS=0
#master_map
#samtools_add_heading $RESULT1/TRI_bwa/no_heading $RESULT1/TRI_bwa
#samtools_sort $RESULT1/TRI_bwa
#unmapped_check $RESULT0 $Ppus $RESULT1/PUS_bwa #ENV: salix
display_run_time $SECONDS MAPPING

#===========================================================
#STEP_3: SNP CALLING: SINGLE, 128G
#PROCESS BAM FILES TO REMOVE BAD SAMPLES:
function master_call {

echo "CALLING SNPs"
#ls $RESULT1/TRI_bwa/*sorted.bam -1| grep -v -f exclude.txt > $bamlist1 #run this when the excluding file is changed
#ls $RESULT1/PUS_bwa/*sorted.bam -1| grep -v -f exclude.txt > $bamlist2

##PUS
#snp_calling $bamlist2 $RESULT2 $Ppus PUS_q10
#filter_vcf $RESULT2 $RESULT2 pus
#
}

#master_call

#TRI
SECONDS=0

#let's use this one. the filters are really good
#freebayes_calling_2_pipe $bamlist_tri_tetra $RESULT2 $Ptri TRI_q10
display_run_time $SECONDS SNP_CALLING

#===========================================================
#STEP 4: EXTRACTING SPECIES SPECIFIC SNPS
SECONDS=0
#extract_loci_by_chromosome $RESULT2 $RESULT2 $NAME $min_freq $Ptri TRI_q10
#extract_loci $RESULT2/PUS_q10.flt2.recode.vcf $RESULT2/PUS_ $NAME $min_freq
#display_run_time $SECONDS species_specific_SNPs

#STEP 5: FIND SNP ENRICHED REGIONS: combine diploid and tetraploid
#extract_window_by_chromosome $diploid_specificSNPs $RESULT2 $RESULT3 TRI_q10 $Ptri

#COMBINE ALL CSV
function combine_wd {

rm TRI_q10_combined_high_resolution_windows.csv
for f in $(ls $RESULT3/*high_resolution_windows.csv); do 
tail -n +2 $f >> $RESULT3/TRI_q10_combined_high_resolution_windows.csv; done

}

#combine_wd

#STEP_6: convert VCF to fasta: array 1-1216
#extract_consensus $DI_vcf $Ptri $RESULT3/TRI_q10_combined_high_resolution_windows.csv $RESULT4 $RESULT2 TRI_q10 P.trichoides

#STEP 7: COPY, RENAME, ALIGNMENT, IQTREE, env=easy353 , array: 1-1216
#easy353_mafft $RESULT4 $RESULT5 $CSV #MEM >=2G
#SLURM_ARRAY_TASK_ID=2
#iqtree_per_gene $RESULT5 $RESULT6 #MEM=4G


#STEP_7: ASTRAL TREE: Run as single

echo "PICK REGIONS WITH THE HIGHEST RESOLVING POWER"
cp $SAMPLE_SP_CSV $OUTDIR/sample_names.csv
#./count_mono.py process_treefiles -i $INDIR -o $OUTDIR

function RUN_ASTRAL {

INDIR=$1
OUTDIR=$2
SAMPLE_SP_CSV=$3

echo "RUNNING ASTRAL ANALYSIS"

run_astral $INDIR $OUTDIR 1 #get at least 1 genes for each mono taxa
run_astral $INDIR $OUTDIR 2 #get at least 2 genes for each mono taxa
run_astral_all_genes $INDIR $OUTDIR SNPs_all_genes
}

RUN_ASTRAL $RESULT6 $RESULT7 $NAME
 
}

main