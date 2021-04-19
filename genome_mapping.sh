# 08.02.18, Denise
# adapting pipelines to realign whole genomes.


module load BWA
module load SAMtools
module load picard


# variables ##need to list here all of the files and directories that I am going to use
# example to fill:
ref=/nesi/nobackup/uoo02327/denise/Kea-Kaka_take2/realignment/N_meridionalis_pseudochr.fa
datadir=/nesi/nobackup/uoo02327/denise/Kea-Kaka_take2/realignment/sequences/
#next two lines are optional, they depend on the sample name, it can include only the extension of the sample files
#but they are used below to get the correct sample name in the readgroup
fq1=_R1_001.fastq
fq2=_R2_001.fastq

#the whole list needs to be in "" and the samples need to be separated by a white space
platform="Illumina"
picard=/opt/nesi/mahuika/picard/2.1.0/picard.jar

#index the reference fasta file
if [ ! -e $ref.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
bwa index $ref
else
echo "BWA Index file found"
fi

#create dictionary file of the reference if it doesn't exist
dictfile=${ref%.*}.dict
if [ ! -e "${ref%.*}.dict" ]; then
echo "Dictionary file of reference does not exist: creating dictionary"
java -jar $picard CreateSequenceDictionary R=$ref O=${ref%.*}.dict
else
echo "Dictionary file found"
fi

#index the reference if it is not indexed already
idxfile=${ref%.*}.fai
if [ ! -e "$ref.fai" ]; then
echo "Reference is not indexed: indexing with SAMTOOLS"
samtools faidx $ref
else
echo "Index for reference file found"
fi
#####################################################

samplist=$(awk "NR==$SLURM_ARRAY_TASK_ID" filelist.txt)

for samp in $samplist

do

echo "processing data for $samp"

#unzip data if necessary

if [ -e $datadir${samp}${fq1}.gz ]; then
echo "Unzipping fastq files"
gunzip -f $datadir${samp}${fq1}.gz
fi
if [ -e $datadir${samp}${fq2}.gz ]; then
echo "Unzipping fastq files"
gunzip -f $datadir${samp}${fq2}.gz
fi


#get readgroup info
##split the machine/lane info from the fq file
infoline=$(head -n 1 ${datadir}$samp$fq1)
instrument=`echo $infoline | cut -d ':' -f1`
instrument=$(echo $instrument | sed -e 's/ /_/g') # this is to fix issues with spaces in the instrument name
instrumentrun=`echo $infoline | cut -d ':' -f2`
flowcell=`echo $infoline | cut -d ':' -f3`
lane=`echo $infoline | cut -d ':' -f4`
index=`echo $infoline | cut -d ':' -f10`
##the next two lines come from the sample name itself (sample and library in my case)
sampname=`echo $samp | cut -d '_' -f1`
library=`echo $samp | cut -d '_' -f2`
#work that info into some
rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
rgpl="PL:${platform}"
rgpu="PU:${flowcell}.${lane}"
rglb="LB:${sampname}_${library}"
rgsm="SM:${sampname}"

#specify the paired-end reads
fastq1=${datadir}$samp$fq1
fastq2=${datadir}$samp$fq2

#transforming fastq to unsortedBAM
echo "producing uBAM for $samp"
java -jar $picard FastqToSam \
FASTQ=$fastq1 \
FASTQ2=$fastq2 \
OUTPUT=${samp}_fastqtosam.bam \
READ_GROUP_NAME=$rgid \
SAMPLE_NAME=$rgsm \
LIBRARY_NAME=$rglb \
PLATFORM_UNIT=$rgpu \
PLATFORM=$rgpl

echo "Finished processing $samp"

done
##############################

#merge bam files for same sample

echo "processing $sampname"

mkdir ${sampname}_process/

mv ${sampname}_*_fastqtosam.bam ${sampname}_process/

cd ${sampname}_process/

echo "merging $sampname"
ls *.bam > bamlist
samtools merge -f ${sampname}_fastqtosam.bam -b bamlist


echo "sorting $sampname"
java -jar $picard SortSam \
I=${sampname}_fastqtosam.bam \
O=${sampname}_sorted_fqtosam.bam \
SORT_ORDER=queryname \
TMP_DIR=./tmp

#marking IlluminaAdapters
echo "marking IlluminaAdapters for $sampname"
java -jar $picard MarkIlluminaAdapters \
I=${sampname}_sorted_fqtosam.bam \
O=${sampname}_markilluminaadapters.bam \
M=${sampname}_markilluminaadapters_metrics.txt \
TMP_DIR=./tmp

#cleaning up
rm ${sampname}_fastqtosam.bam

#Aligning and merging with uBAM all piped
echo "Aligning and merging with uBAM for $sampname"
java -jar $picard SamToFastq \
I=${sampname}_markilluminaadapters.bam \
FASTQ=/dev/stdout QUIET=true \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true \
TMP_DIR=./tmp | bwa mem -M -t 8 -p $ref /dev/stdin | java -jar $picard MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=${sampname}_markilluminaadapters.bam \
OUTPUT=${sampname}_piped.bam \
R=$ref CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS  \
TMP_DIR=./tmp

#cleaning up
rm ${sampname}_markilluminaadapters.bam

#Getting stats for the alignment
echo "Getting stats for $sampname"
samtools stats ${sampname}_piped.bam > ${sampname}_samtools_stats.txt

echo "Finished alignment $sampname"

# sorting .bam file and moving it out
echo "Sorting BAM file for $sampname"
java -jar $picard SortSam \
I=${sampname}_piped.bam \
O=${sampname}_sorted.bam \
SORT_ORDER=coordinate \
TMP_DIR=./tmp

# removing duplicates from .bam file and moving it out
echo "Remove duplicates from BAM file for $sampname"
java -jar $picard MarkDuplicates \
I=${sampname}_sorted.bam \
O=${sampname}_sorted_rmdup.bam \
M=${sampname}_marked_dup_metrics.txt \
REMOVE_DUPLICATES=TRUE \
ASSUME_SORTED=TRUE \
TMP_DIR=./tmp

mv ${sampname}_sorted_rmdup.bam ..

cd ..
