## script to postprocess outputs from HMMERscan and merge them into counts per family ID per species
## Denise, some fixes on 04.02.19
samplist="N_meridionalis F_peregrinus T_guttata N_notabilis M_undulatus G_gallus"

for samp in $samplist
do
  # extract only first result for each gene
  sort -k3,3 -u ${samp}_output | sort -k1,1 > ${samp}_unique_output.txt

  for ID in `awk '{print $1}' ${samp}_unique_output.txt | sort -u`
  do
    # write ID in output file
    echo -n `echo $ID`$'\t' >> ${samp}_family_count.txt
    # write count in output file
    echo `awk -v var="$ID" '$0~var{print $0}' ${samp}_unique_output.txt | awk 'END{print NR}'` >> ${samp}_family_count.txt
  done
  # because I will be curious as to how many families I have and how many genes have been assigned in total:
  echo $samp >> stats.txt
  echo "number of identified gene families: "`wc -l ${samp}_family_count.txt | cut -d ' ' -f1` >> stats.txt
  echo "number of genes assigned: "`awk '{sum += $2} END{print sum}' ${samp}_family_count.txt` >> stats.txt

done

# merging species counts in one file, first two by two
join -a1 -a2 -o 0,1.2,2.2 -e "0" N_meridionalis_family_count.txt N_notabilis_family_count.txt > join1.txt
join -a1 -a2 -o 0,1.2,2.2 -e "0" G_gallus_family_count.txt T_guttata_family_count.txt > join2.txt
join -a1 -a2 -o 0,1.2,2.2 -e "0" F_peregrinus_family_count.txt M_undulatus_family_count.txt > join3.txt

# then first four together
join -a1 -a2 -o 0,1.2,1.3,2.2,2.3 -e "0" join1.txt join2.txt > join4.txt
# finally, the last merge, including a header first
echo "ID Nmer Nnot Ggal Tgut Fper Mund" > final_treefam_output.txt
join -a1 -a2 -o 0,1.2,1.3,1.4,1.5,2.2,2.3 -e "0" join4.txt join3.txt >> treefam_output.txt

# final fixes to the output
grep -v "#" treefam_output.txt > final_treefam_output.txt
sed -i -e 's/TF/\tTF/' final_treefam_output.txt
sed -i -e 's/ /\t/g' final_treefam_output.txt
sed -i -e 's/ID/Description\tID/' final_treefam_output.txt
