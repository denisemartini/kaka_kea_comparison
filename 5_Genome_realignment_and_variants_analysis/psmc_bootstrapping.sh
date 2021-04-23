#!/bin/bash -e
# Denise, 09.10.18
# Script to run PSMC in a 100 bootstrap fashion for my samples
# modified on 14.02.19 to take the sample name as an argument, not in a loop

samp=$1
psmcdir=/usr/local/psmc-r49

echo "pre-processing "$samp", starting at "$(date)
${psmcdir}/fq2psmcfa -q20 ${samp}_diploid.fq.gz > ${samp}.psmcfa
${psmcdir}/splitfa ${samp}.psmcfa > ${samp}_split.psmcfa

echo "running psmc on "$samp", starting at "$(date)
${psmcdir}/psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o ${samp}.psmc ${samp}.psmcfa

echo "bootstrapping "$samp", starting at "$(date)
parallel -j8 "/usr/local/psmc-r49/psmc -N30 -t5 -r5 -b -p "4+30*2+4+6+10" ${samp}_split.psmcfa -o ${samp}_round-{}.psmc" :::: <(seq 100)

cat ${samp}.psmc ${samp}_round-*.psmc > ${samp}_combined.psmc

echo $samp "done at "$(date)
