#!/bin/bash
#SBATCH --job-name=cap_SNP_calling
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --array=1
#SBATCH --mem=8G

#40-47,54-58,70-73,75,86
#salix diploid: 1-31,39-47,54-58,70-73,75,86 #no myrsinites
#other salix diploids:
#caprea(54-58),repens(70-73,75)

#salix_tetra:48-52,65-69,76-85
#alba: 76-78
#aurita: 48-52 (53 is bad)
#cinerea: 79-82
#pentandra: 83-85
#phylicifolia: 65-69


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
seqtk sample -s 10 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz $frac | gzip > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz #> ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq
seqtk sample -s 10 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz $frac | gzip > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz #> ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq

#COMPRESS


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
OUTFILE=$2

RG='@RG\tID:'${prefix}_${SLURM_ARRAY_TASK_ID}'\tSM:'${prefix}_${SLURM_ARRAY_TASK_ID}'\tPL:ILLUMINA'

samtools addreplacerg -r $RG \
                      -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fixed.bam \
                      $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

ls $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fixed.bam >> $OUTFILE
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fixed.bam
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

#STEP 4:
#THIS IS A NEW FILTER. SAME PARAMETER AS POTAMOGETON TETRA AND DIPLOID
function bcftools_filter {

INVCF=$1
OUTDIR=$2
NAME=$3

#set up max DP
meanDP=20
SAMPLE_COUNT=$(bcftools query -l $INVCF|wc -l|cut -f 1 -d ' ')
MAX_DP=$((SAMPLE_COUNT * $meanDP))
echo "$SAMPLE_COUNT samples, max DP=${meanDP}x${SAMPLE_COUNT}=${MAX_DP}"

echo "Filtering diploid SNPs to tetraploid-like standards: $NAME"
    
     bcftools view $INVCF -O u 2>/dev/null \
 |   bcftools view -i '1==1' 2>/dev/null \
 |   bcftools +fill-tags -- -t AC,AN,AF,MAF,HWE  \
 |   bcftools filter   -e "QUAL<30 || INFO/DP > $MAX_DP" \
 |   bcftools filter   -i 'INFO/AC[1] >= 3' \
 |   bcftools filter   -i '(COUNT(FMT/DP<5) + COUNT(FMT/DP="."))<= 2 && COUNT(FMT/DP>30) <= 1' \
 |   bcftools filter   -e 'INFO/MAF[0] < 0.01'  \
 |   bcftools filter   -e 'INFO/MQ < 30'  \
 |   bcftools filter   -e 'INFO/VDB < 0.1'  -Oz -o $OUTDIR/${NAME}.flt3.vcf.gz
    
    echo "Indexing filtered VCF..."
    bcftools index $OUTDIR/${NAME}.flt3.vcf.gz -f

}

#tetraploid SNP calling
function freebayes_calling_2_pipe {

bamlist=$1
OUTDIR=$2
REF=$3
NAME=$4

echo "doing freebayes"

while read bam; do
  samplename=$(basename $bam .bam)
  echo -e "$samplename"
done < $bamlist > samples.tsv

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
  | bcftools filter -i 'INFO/AC[1] >= 3' \
  | bcftools filter -i '(COUNT(FMT/DP<5) + COUNT(FMT/DP="."))<= 2 && COUNT(FMT/GQ<20) <= 2 && COUNT(FMT/DP>40) <= 2 ' -Oz -o $OUTDIR/${NAME}_${CHROM_NAME}.flt3.vcf.gz

bcftools index $OUTDIR/${NAME}_${CHROM_NAME}.flt3.vcf.gz #need index for consensus extraction
 
#Allele frequency bounds - remove very rare/very common SNPs

#Quality/Depth ratio - normalized quality score

#Genotype quality filter - GQ < 20

#Minor allele count - require at least 2 samples with alt allele
#SAP=Strand bias p-value (Phred-scaled)
#RPP = Read Placement Probability for alternate allele (Phred-scaled)
#INFO/AC[1] >= 4: include allele >=4 copies: assume 1 copy per sample, at least one copy to be present in all samples of one species)
}
#
#

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
echo "extract species specific SNPs on $CHROM_NAME"

#./calculate_snp.freq.py specifi_SNPs \
./species_specific_allele.tetra.py specifi_SNPs \
        -v $INDIR/${PREFIX}_${CHROM_NAME}.flt3.vcf.gz \
        -n $NAME \
        --min_freq $min_freq \
        -o $OUTDIR/${PREFIX}_${CHROM_NAME}_ #output directory

}
##
#
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
echo "finding species specific SNP enriched windows on $CHROM_NAME"

#GET VCF FOR TETRA
VCF2=$INDIR/${PREFIX}_${CHROM_NAME}_hq_specific_allele_freq.csv

#COMBINE HQ SPECIFIC SNPS
cat $DI_VCF $VCF2|grep "${CHROM_NAME}"|cut -d ',' -f 1-6 > ${OUTDIR}/${PREFIX}_${CHROM_NAME}_combined_SNPs.csv

#./calculate_snp.freq.py specifi_SNPs \
./species_specific_allele.tetra.py extract_ssSNP_regions \
        -n ${OUTDIR}/${PREFIX}_${CHROM_NAME}_combined_SNPs.csv \
        -o ${OUTDIR}/${PREFIX}_${CHROM_NAME}

}
##

#CONVERT VCF TO FASTA BY CALLING CONSENSUS
function extract_consensus {

VCF=$1
REF=$2
OUTDIR=$3
PREFIX=$4
REF_SP=$5

echo "EXTRACT CONSENSUS SEQUENCE FROM $REF REF AND $VCF"

# 1. SET UP PARAMETERS
WINDOWS=$RESULT3/${PREFIX}_high_resolution_windows.csv

GENE=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 1 -d ',')
NAME=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 2 -d ',') #THE CHROMOSOME
START=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 3 -d ',')
END=$(cat $WINDOWS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -f 4 -d ',')
REGION="ENA|${NAME}|${NAME}.1:${START}-${END}"


mkdir ${OUTDIR}/${GENE}_tmp #a temporary directory so the files won't mix and interfere with each other

# 2. EXTRACT REGION from reference
samtools faidx $REF "$REGION" > ${OUTDIR}/${GENE}.fasta
# Add sample name to header
sed -i "1s/^>.*/>${REF_SP}_ref/" ${OUTDIR}/${GENE}.fasta

# 3. List all DIPLOID samples in the VCF#need to split this into two steps: diploids and tetraploids, then combine results
SAMPLES=$(bcftools query -l $VCF)

# 4. Generate consensus for each DIPLOID sample
for SAMPLE in $SAMPLES; do
    TMP=$(echo $SAMPLE|rev|cut -f 1 -d '/'|rev);SAMPLE_NAME=${TMP/.sorted.bam/} #get the last field and remove .sort.bam
    echo "Processing $SAMPLE_NAME..."
    
    # Generate consensus for this sample
    samtools faidx $REF $REGION | bcftools consensus $VCF -s $SAMPLE -I -o ${OUTDIR}/${GENE}_tmp/tmp.${GENE}.${SAMPLE_NAME}.fasta
          
    # Add sample name to header
    sed -i "1s/^>.*/>${SAMPLE_NAME}/" ${OUTDIR}/${GENE}_tmp/tmp.${GENE}.${SAMPLE_NAME}.fasta
done

# 5. Combine all into one multi-FASTA
cat ${OUTDIR}/${GENE}_tmp/tmp.${GENE}.*.fasta >> ${OUTDIR}/${GENE}.fasta

rm -rf ${OUTDIR}/${GENE}_tmp
}

#combine di and tetra vcf for extraction 

function extract_consensus_di_tetra {

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

mkdir $OUTDIR/tmp_$GENE
cp $INDIR/$FASTA $OUTDIR/tmp_$GENE
sed -i 's/>S\.caprea_OZ[0-9]*_[0-9]*_[0-9]*/>S.caprea_ref/' $OUTDIR/tmp_$GENE/$FASTA #rename the referenece in the treefile
sed -i 's/>S\.repens_OZ[0-9]*_[0-9]*_[0-9]*/>S.repens_ref/' $OUTDIR/tmp_$GENE/$FASTA #rename the referenece in the treefile
#sed -i 's/>${REF_SP}_OZ[0-9]*_[0-9]*_[0-9]*/>${REF_SP}_ref/' $OUTDIR/tmp_$GENE/$FASTA #rename the referenece in the treefile: this doesn't work
#remove bad taxa (sample is bad) and all gap contigs
#./remove_taxa.py $OUTDIR/$GENE/$FASTA 
iqtree -s $OUTDIR/tmp_$GENE/$FASTA -bb 1000 -redo -safe
mv $OUTDIR/tmp_$GENE/*treefile $OUTDIR
rm -rf $OUTDIR/tmp_$GENE


}
function run_astral {

INDIR=$1
OUTDIR=$2
cvg=$3 #can resolve >= xx phylogeny
RESOLUTION=$4

echo 'running astral'
./count_mono.py gene_select -i $INDIR -o $OUTDIR --cvg $cvg -r $RESOLUTION
#create astral phylogeny with selected genes
sed -i 's/S\.caprea_OZ[0-9]*_[0-9]*_[0-9]*/S.caprea_ref/' $OUTDIR/rsl_${cvg}_genes.in.treefile #rename reference seq
sed -i 's/S\.repens_OZ[0-9]*_[0-9]*_[0-9]*/S.repens_ref/' $OUTDIR/rsl_${cvg}_genes.in.treefile #rename reference in tree
java -jar $astral -i $OUTDIR/rsl_${cvg}_genes.in.treefile -o $OUTDIR/rsl_${cvg}_genes.out.treefile 2>out.log
echo -e "\n\n\nDONE running astral with ${cvg}-gene per species coverage\n\n\n"
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
echo -e "\n\n\nDONE running astral with all genes\n\n\n"

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
Scap=$REF/S_caprea.chromosome.fasta #size=365398294 bp
Srep=$REF/S_repens.chromosome.fasta #size=322687971 bp
Sret=$REF/S_reticulata.chromosome.fasta
bamlist2=/mnt/shared/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/bamlist.all_dip.CAP.txt
bamlist3=/mnt/shared/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/bamlist.all_dip.REP.txt
bamlist_tetra_CAP=/mnt/shared/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/bamlist.all_TETRA.CAP.txt
bamlist_tetra_REP=/mnt/shared/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/bamlist.all_TETRA.REP.txt
SUBSAMPLE=/home/zchen/projects/rbge/zedchen/barcoding/salix_combined/SNP_scripts/subsampling.csv
min_freq=90
NAME=/mnt/shared/projects/rbge/pholling/barcoding/salix_combined/results/00_reads/sample_sp.csv
#==========================================================================================================================================


#subdirectories
RESULT0=$WORKDIR/00_subsampling
RESULT1=$WORKDIR/01_sorted_bam #Skip the sam stage--> pipe to bam directly
RESULT2=$WORKDIR/02_VCFs
RESULT3=$WORKDIR/03_HQ_windows
RESULT4=$WORKDIR/04_SNP_fasta
RESULT5=$WORKDIR/05_aligned_fasta
RESULT6=$WORKDIR/06_IQTREE
RESULT7=$WORKDIR/07_ASTRAL
RESULT=$WORKDIR/results/0

#USAGES
#=============================================
function setup_dir {

mkdir $WORKDIR
mkdir $WORKDIR/refs
mkdir $RESULT1
mkdir $RESULT0 -p
mkdir $RESULT1/CIN_bwa
mkdir $RESULT1/CAP_bwa
mkdir $RESULT1/REP_bwa
mkdir $RESULT1/CAP_bwa_all
mkdir $RESULT1/REP_bwa_all
mkdir $RESULT2
mkdir $RESULT3 -p
mkdir $RESULT4/CAP -p
mkdir $RESULT4/REP -p
mkdir $RESULT4/CAP_tetra -p
mkdir $RESULT4/REP_tetra -p
mkdir $RESULT4/CAP_tetra2 -p
mkdir $RESULT4/REP_tetra2 -p
mkdir $RESULT5/CAP -p
mkdir $RESULT5/REP -p
mkdir $RESULT5/CAP_tetra -p
mkdir $RESULT5/REP_tetra -p
mkdir $RESULT5/CAP_tetra2 -p
mkdir $RESULT5/REP_tetra2 -p
mkdir $RESULT6/CAP -p
mkdir $RESULT6/REP -p
mkdir $RESULT6/CAP_tetra -p
mkdir $RESULT6/REP_tetra -p
mkdir $RESULT6/CAP_tetra2 -p
mkdir $RESULT6/REP_tetra2 -p
mkdir $RESULT7/CAP -p
mkdir $RESULT7/REP -p
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
fix_RG_bam $RESULT1/CAP_bwa bamlist.all_TETRA.CAP.txt
fix_RG_bam $RESULT1/REP_bwa bamlist.all_TETRA.REP.txt

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
#bcftools_filter $RESULT2/CAP.flt1.vcf $RESULT2 CAP
display_run_time $SECONDS SNP_FILTERING

#REPENS
SECONDS=0
#snp_calling $bamlist3 $RESULT2 $Srep REP
display_run_time $SECONDS SNP_CALLING
#
SECONDS=0
#bcftools_filter $RESULT2/REP.flt1.vcf $RESULT2 REP
display_run_time $SECONDS SNP_FILTERING
#
#filter_individual $RESULT2 $RESULT2 CAP
#===========================================================
#TETRA CALLING/EXTRACT LOCI/WINDOW
#Freebayes: SNP calling for tetraploid samples
#
#freebayes_calling_2_pipe $bamlist_tetra_CAP $RESULT2 $Scap CAP #1-19???
#freebayes_calling_2_pipe $bamlist_tetra_REP $RESULT2 $Srep REP #1-18
#memory: 
#64G:48-53,65-69,76-85

#EXTRACT: REP: 1-18, CAP: 1-19
#SLURM_ARRAY_TASK_ID=2
#extract_loci_by_chromosome $RESULT2 $RESULT2 $NAME $min_freq $Srep REP
#extract_window_by_chromosome /mnt/shared/scratch/zchen/Barcoding_km/SNP_salix/02_VCFs/REP_hq_specific_allele_freq.csv $RESULT2 $RESULT3 REP $Srep
#
#extract_loci_by_chromosome $RESULT2 $RESULT2 $NAME $min_freq $Scap CAP
#extract_window_by_chromosome /mnt/shared/scratch/zchen/Barcoding_km/SNP_salix/02_VCFs/CAP_hq_specific_allele_freq.csv $RESULT2 $RESULT3 CAP $Scap

function combine_wd {

PREFIX=$1

rm $RESULT3/${PREFIX}_combined_high_resolution_windows.csv
#
for f in $(ls $RESULT3/${PREFIX}_OZ*high_resolution_windows.csv); do 
head -1 $f; done |sort|uniq >> $RESULT3/${PREFIX}_combined_high_resolution_windows.csv
#
for f in $(ls $RESULT3/${PREFIX}_OZ*high_resolution_windows.csv); do 
tail -n +2 $f >> $RESULT3/${PREFIX}_combined_high_resolution_windows.csv; done

}

#combine_wd REP
#combine_wd CAP
#===========================================================


#STEP 4: EXTRACTING SPECIES SPECIFIC SNPS
SECONDS=0
#extract_loci $RESULT2/CAP.flt3.vcf.gz $RESULT3/CAP_ $NAME $min_freq
#extract_loci $RESULT2/REP.flt3.vcf.gz $RESULT3/REP_ $NAME $min_freq
display_run_time $SECONDS species_specific_SNPs
#===========================================================

#EXTRACT TETRA ALONG WITH DIPLOID ONES (using old diploid windows)
#extract_consensus_di_tetra  $RESULT2/CAP.flt3.vcf.gz $Scap $RESULT3/CAP_high_resolution_windows.csv $RESULT4/CAP_tetra $RESULT2 CAP S.caprea #array: 2-126
#extract_consensus_di_tetra  $RESULT2/REP.flt3.vcf.gz $Srep $RESULT3/REP_high_resolution_windows.csv $RESULT4/REP_tetra $RESULT2 REP S.repens #array: 2-405

#EXTRACT TETRA ALONG WITH DIPLOID ONES: using windows called with tetra and diploid SNPs together
#extract_consensus_di_tetra  $RESULT2/CAP.flt3.vcf.gz $Scap $RESULT3/CAP_high_resolution_windows.csv $RESULT4/CAP_tetra $RESULT2 CAP S.caprea #array: 2-126
#extract_consensus_di_tetra  $RESULT2/REP.flt3.vcf.gz $Srep $RESULT3/REP_combined_high_resolution_windows.csv $RESULT4/REP_tetra2 $RESULT2 REP S.repens #array: 2-322

#STEP 7: COPY, RENAME, ALIGNMENT, IQTREE, env=easy353 , array
#easy353_mafft $RESULT4/CAP $RESULT5/CAP $CSV #MEM >=2G
#easy353_mafft $RESULT4/REP $RESULT5/REP $CSV #MEM >=2G
#
#easy353_mafft $RESULT4/CAP_tetra $RESULT5/CAP_tetra $CSV #MEM >=2G #array: 2-126
#easy353_mafft $RESULT4/REP_tetra $RESULT5/REP_tetra $CSV #MEM >=2G : 2-405
#
#easy353_mafft $RESULT4/REP_tetra2 $RESULT5/REP_tetra2 $CSV #MEM >=2G : 1-321
#easy353_mafft $RESULT4/CAP_tetra2 $RESULT5/CAP_tetra2 $CSV #MEM >=2G : 1-?

#STEP_6: IQTREE PER GENE:gene array: env=easy353 MEM=4G
#iqtree_per_gene $RESULT5/CAP $RESULT6/CAP S.caprea #1-125
#iqtree_per_gene $RESULT5/REP $RESULT6/REP S.repens #1-404
#
#iqtree_per_gene $RESULT5/CAP_tetra $RESULT6/CAP_tetra #1-125
#./count_mono.py process_treefiles -i $RESULT6/CAP_tetra -o $RESULT7/CAP/di_tetra_
#iqtree_per_gene $RESULT5/REP_tetra $RESULT6/REP_tetra #1-404
#./count_mono.py process_treefiles -i $RESULT6/REP_tetra -o $RESULT7/REP/di_tetra_
#
#iqtree_per_gene $RESULT5/REP_tetra2 $RESULT6/REP_tetra2 S.repens #1-321
#./count_mono.py process_treefiles -i $RESULT6/REP_tetra2 -o $RESULT7/REP/di_tetra2_

#STEP_7: ASTRAL TREE: Run as single
#./count_mono.py process_treefiles -i $RESULT6/CAP -o $RESULT7/CAP
function RUN_ASTRAL {

INDIR=$1
OUTDIR=$2

#run_astral_all_genes $INDIR $OUTDIR SNPs_all_genes
run_astral $INDIR $OUTDIR 1 di_tetra_gene_resolution.csv #get at least 1 genes for each mono taxa
run_astral $INDIR $OUTDIR 2 di_tetra_gene_resolution.csv #get at least 2 genes for each mono taxa
run_astral $INDIR $OUTDIR 3 di_tetra_gene_resolution.csv #get at least 2 genes for each mono taxa
}

RUN_ASTRAL $RESULT6/CAP_tetra $RESULT7/CAP
RUN_ASTRAL $RESULT6/REP_tetra $RESULT7/REP


}

main