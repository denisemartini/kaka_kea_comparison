list="a b c d e f"
rep="1 2 3"
# the species list is an argument taken from the command line, it needs to be a list where elements are whitespace separated and all enveloped in ""
# eg: "Nnot Nmer Tgut"
# sbatch nesi_BEB_postprocess.sh "Nnot Nmer Tgut"
species=$1

for sp in $species
do
	for n in $rep
	do
		for f in $list
		do
			cat alt_${sp}_${n}/alt_${sp}_${n}_${f}.mlc | sed -n '/Data set/p;/Bayes Empirical Bayes/,/^The grid/p' | grep -v "^The grid" | grep -v "BEB" | grep -v "Positive" > BEB_sites_${sp}_${n}_${f}.txt
			cat aln_list_${f} | while read line
			do
				gene=$(echo ${line} | cut -d '_' -f1-2)
				nr=$(awk -v var="${line}" '$0~var{print NR}' aln_list_${f})
				sed -i "s|Data set ${nr}\$|${gene}|" BEB_sites_${sp}_${n}_${f}.txt
			done
		done
		cat BEB_sites_${sp}_${n}_*.txt > BEB_sites_${sp}_${n}_all.txt
		rm BEB_sites_${sp}_${n}_?.txt
	done
	for gene in $(awk '{print $1}' sig_PS_genes_${sp}.txt)
	do
		echo ${gene} >> sig_PSG_sites_${sp}.txt
		for n in $rep
		do
			sed -n "/${gene}\$/,/OG_/p" BEB_sites_${sp}_${n}_all.txt | grep -v "OG" | tr "\n" "\t" >> sig_PSG_sites_${sp}.txt
			echo "" >> sig_PSG_sites_${sp}.txt
		done
	done
done
