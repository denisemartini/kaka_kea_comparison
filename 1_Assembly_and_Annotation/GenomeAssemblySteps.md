# Whole Genome Assembly

### Preprocessing steps
Before proceeding with the genome assembly, we used kmergenie and sga-qc to get an assessment of the genome size and k-mer size to use for assembly. We then used Trimmomatic to cut adapter sequences and low-quality bases.

#### KMERGENIE (v1.6424)
To run Kmergenie on a set of sequences, you need to list them in a file, one line per sequence file:
```bash
ls /path/to/sequences/*.fastq.gz > list_reads.txt
```
Then you can run kmergenie on the list like:
```bash
/path/to/kmergenie list_reads.txt
# if kmergenie has been added to your PATH, you don't need to specify the path in the command
```
Kmergenie by default tests all kmers between 15 and 121, so to change those lower and upper limits you need to set the values of the -k and -l options. The -t options sets the number of threads to use (starts up several instances of kmergenie at the same time). An example of how to run it this way, also redirecting the output to a log file:
```bash
/path/to/kmergenie list_reads.txt -k 131 -l 41 -t 8 list_reads.txt &> kmergenie.log
```
A sample output report and the downloadable software are available at: http://kmergenie.bx.psu.edu/.

#### SGA-PREQC (v0.9.4)
To run SGA_PREQC on all sequences, they need to be unzipped and concatenated in one file:
```bash
cd /path/to/sequences/
gunzip *.fastq.gz
cat *R1*.fastq >> ./all_forward.fastq
cat *R2*.fastq >> ./all_reverse.fastq
```
Then, running all steps:
```bash
# setting this up as a variable and adding it to my PATH
$SGA_HOME=/path/to/SGA_directory
PATH=$SGA_HOME/bin:$PATH

# 1) preprocessing
sga preprocess --pe-mode 1 all_forward.fastq all_reverse.fastq > my_genome.fastq
# 2) index the data
sga index -a ropebwt -t 8 --no-reverse my_genome.fastq
# 3) run preqc
sga preqc -t 8 my_genome.fastq > my_genome.preqc							
# 4) make the report
$SGA_HOME/src/bin/sga-preqc-report.py my_genome.preqc $SGA_HOME/src/examples/preqc/*
```
The effective amount of time these steps take depends on the system, but the index and preqc take very long (24-32hrs in our case), while the report only takes a few seconds.
Further instructions and the software are available at: https://github.com/jts/sga.

#### TRIMMOMATIC (v0.38)
We set up Trimmomatic to trim sequences containing adapters and low average quality bases, finally removing sequences shorter than 36bp after trimming.
```bash

```

### Assembly
The main genome draft was produced with SOAPdenovo2, while the Satsuma2-Chromosemble pipeline was used to scaffold the assembly with synteny information from a reference species.

#### SOAPdenovo2 (v2.04-r240)
We ran [SOAPdenovo2](https://github.com/aquaskyline/SOAPdenovo2) with several different kmer sizes, starting from the values identified by kmergenie.
An example of how to run this (and output in a log file):
```bash
/path/to/SOAPdenovo-127mer all -s path/to/soap.config -K 53 -o my_genome_53 -p 8 -N 1157400000 &> SOAPdenovo53.log
/path/to/SOAPdenovo-127mer all -s path/to/soap.config -K 43 -o my_genome_43 -p 8 -N 1157400000 &> SOAPdenovo43.log
```
The parameter -s requires a path to a configuration file that includes information on the input files to use for the assembly (and the order in which they are to be used, more info and an example are available at: https://github.com/aquaskyline/SOAPdenovo2), -K is the kmer size to test, -o is the output file prefix, -p is the number of threads to use and -N is the estimated genome size (from SGA-PREQC in our case).

#### Satsuma (v3.1.0)
We used the Chromosemble pipeline in [Satsuma](http://satsuma.sourceforge.net/) to run two steps of super-scaffolding by synteny on our assembly, using the zebra finch (_Taeniopygia guttata_) genome as reference.
The command only requires a target file (-t, the reference genome), a query file (-q, the assembly to scaffold), an output folder (-o) and the number of threads to use (-n). As optional parameters, setting -thorough as true (1) runs a (slow) higher sensitivity alignment and setting -pseudochr as true (1) maps the scaffolds into chromosomes.
Example:
```bash
Chromosemble -t /path/to/reference.fa -q /path/to/assembly.fa \
-o /path/to/output/directory -n 8 -thorough 1 -pseudochr 1 |& tee -a log_thorough-satsuma.txt
# if needed to set SATSUMA2_PATH, first run:
export SATSUMA2_PATH=/path/to/satsuma/bin
```
More options and information available at: http://satsuma.sourceforge.net/.

### Final evaluation
We used [BUSCO v3](https://busco-archive.ezlab.org/v3/) to verify the completeness of the genome in the final assembly.
We ran a genome assessment (setting the genome mode, -m) using the avian database (aves_odb9 in lineage, -l), with 8 cores (-c) and setting "chicken" as the Augustus species gene finding parameter (-sp):
```bash
python /path/to/run_BUSCO.py -i /path/to/final_assembly.fa -o /path/to/output/directory -l /path/to/busco-v3-datasets/aves_odb9 -m genome -c 8 -sp chicken > BUSCO_final_assembly.log
# in case of trouble, it is probably because the Augustus config directory is in a non-writable path; the solution is to copy the contents of the directory somewhere you have write permission and run:
export AUGUSTUS_CONFIG_PATH=/path/to/writable/augustus-3.3/config
```
We then ran the assemblathon_stats.pl script from https://github.com/KorfLab/Assemblathon to collect a few quality statistics on the final assembly.
```bash
/path/to/assemblathon_stats.pl /path/to/final_assembly.fa > final_assembly.stats
```
