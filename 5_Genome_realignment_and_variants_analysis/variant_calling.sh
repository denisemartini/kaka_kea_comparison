# 09.02.18, Denise
# outputting variant calling for both PSMC and SnpEff

module load GCC
module load SAMtools
module load BCFtools
module load BEDTools
module load picard

picard=/opt/nesi/mahuika/picard/2.1.0/picard.jar


# variables ##need to list here all of the files and directories that I am going to use
# example to fill:
ref=/nesi/nobackup/uoo02327/denise/Kea-Kaka_take2/realignment/N_meridionalis_pseudochr.fa

# the sample
samp=$(awk "NR==$SLURM_ARRAY_TASK_ID" samplist.txt)

echo 'Processing '${samp}
# index bam file and collect a couple of metrics
echo 'Indexing '${samp}
samtools index ${samp}_sorted_rmdup.bam
java -jar $picard CollectWgsMetrics \
I=${samp}_sorted_rmdup.bam \
O=${samp}_wgs_metrics.txt \
R=${ref}

# extract min and max depth values from the stats
meancov=$(grep -A1 'MEAN_COVERAGE' ${samp}_wgs_metrics.txt | tail -1 | awk '{print $2}')
min=$(echo $meancov/3 | bc -l | cut -d '.' -f1)
max=$(echo $meancov*2 | bc -l | cut -d '.' -f1)

# prepare a mask file to remove from the pileup the parts with low coverage, so that pileup goes faster
echo 'Preparing mask file for '${samp}
bedtools genomecov -ibam ${samp}_sorted_rmdup.bam -bga | awk '$4 < 5' > ${samp}_lowcov.bed
bedtools subtract -a ${ref}_genome.bed -b ${samp}_lowcov.bed > ${samp}_mask.bed

# call consensus for PSMC
echo 'Preparing PSMC input for '${samp}
samtools mpileup -uf $ref ${samp}_sorted_rmdup.bam -l ${samp}_mask.bed \
--min-MQ 20 --min-BQ 20 | bcftools call -c - | vcfutils.pl vcf2fq -d $min -D $max | gzip > ${samp}_diploid.fq.gz

# call variants for SnpEff
echo 'Calling variant for '${samp}
bcftools mpileup --min-MQ 20 --min-BQ 20 \
--output-type u \
--fasta-ref $ref ${samp}_sorted_rmdup.bam | bcftools call \
--multiallelic-caller --variants-only \
--output-type v - --output ${samp}.vcf

bgzip ${samp}.vcf
tabix -p vcf ${samp}.vcf.gz
