# Denise Martini
# Pipeline to realign whole genomes using BWA/SAMtools/Picard tools
# This is meant to be run on samples/individuals that are associated with
# multiple sets of paired-end sequence files. The script takes an argument that is
# a line number from a 'filelist.txt' (this file needs to be in the directory where the script is run,
# if in a different location this should be specified below, line 73).
# An example of how 'filelist.txt' should look like is at lines 60-66.
# To run the script on the files in the 2nd line of filelist.txt, run:
# bash genome_mapping.sh 2

# In case you are running this in a module-loading server type (e.g. NeSI):
# module load BWA
# module load SAMtools
# module load picard
# If using the lines above you can remove the variable sign ($) from $samtools and $bwa
# in the rest of the script (not from $picard), and also comment out lines 30 and 31.

# variables ##need to list here all of the files and directories that I am going to use
# example to fill:
ref=/path/to/reference/file.fa
datadir=/path/to/directory/containing/sequences/

#next two lines are optional, they depend on the sample name, they can include only the extension of the sample files
#but they are used below to get the correct sample name in the readgroup, so they should be adjusted
#ex. this is set for files named as: 'indv2_lib1_lane4_R1_001.fastq'
fq1=_R1_001.fastq
fq2=_R2_001.fastq

#the whole list needs to be in "" and the samples need to be separated by a white space
platform="Illumina"
picard=/path/to/picard.jar
bwa=/path/to/bwa
samtools=/path/to/samtools

#index the reference fasta file
if [ ! -e $ref.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
$bwa index $ref
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
$samtools faidx $ref
else
echo "Index for reference file found"
fi
#####################################################

# here you need a 'filelist.txt' which has every file associated with the same sample
# (e.g. from different lanes, or libraries) in a single line, white-space separated,
# one sample per line, example:
###
# indv1_file1 indv1_file2 indv1_file3
# indv2_lib1_lane4 indv2_lib1_lane5 indv2_lib2_lane4 indv2_lib2_lane5
# indv3_lib1_lane4 indv3_lib1_lane5 indv3_lib2_lane4 indv3_lib2_lane5
###
# when running this as an array job on a server (e.g. with a SLURM scheduler):
# samplist=$(awk "NR==$SLURM_ARRAY_TASK_ID" filelist.txt)
# otherwise the line number is provided as an argument to the script:
samplist=$(awk "NR==$1" filelist.txt)

# this section works through each of the white-space separated files (samp) in the
# chosen line of 'filelist.txt'
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

##the next two lines come from the names of files in 'filelist.txt' (sample_library here)
##they should be fixed to reflect whatever the user's 'filelist.txt' looks like
sampname=`echo $samp | cut -d '_' -f1`
library=`echo $samp | cut -d '_' -f2`

#work that info into some read group identifiers:
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

#merge bam files for same sample (merge files from the same line of 'filelist.txt')

echo "processing $sampname"

mkdir ${sampname}_process/

mv ${sampname}_*_fastqtosam.bam ${sampname}_process/

cd ${sampname}_process/

echo "merging $sampname"
ls *.bam > bamlist
$samtools merge -f ${sampname}_fastqtosam.bam -b bamlist


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
$samtools stats ${sampname}_piped.bam > ${sampname}_samtools_stats.txt

echo "Finished alignment $sampname"

# sorting .bam file
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
