list="a b c d e f"
# the species list is an argument taken from the command line, it needs to be a list where elements are whitespace separated and all enveloped in ""
# eg: "Nnot Nmer Tgut"
# the second argument is the file containing the list of genes tested
# sbatch nesi_lnL_postprocess.sh "Nnot Nmer Tgut" OGs_custom_tested_list.txt
species=$1

for sp in $species
do
  echo "lnL_1" >> alt_${sp}_1.lnl
  echo "lnL_2" >> alt_${sp}_2.lnl
  echo "lnL_3" >> alt_${sp}_3.lnl
  for f in $list
  do
    grep 'lnL' alt_${sp}_1/alt_${sp}_1_${f}.mlc | awk '{print ($5)}' >> alt_${sp}_1.lnl
    grep 'lnL' alt_${sp}_2/alt_${sp}_2_${f}.mlc | awk '{print ($5)}' >> alt_${sp}_2.lnl
    grep 'lnL' alt_${sp}_3/alt_${sp}_3_${f}.mlc | awk '{print ($5)}' >> alt_${sp}_3.lnl
  done

  paste $2 alt_${sp}_1.lnl alt_${sp}_2.lnl alt_${sp}_3.lnl > alt_${sp}_all.lnl
  rm alt_${sp}_1.lnl alt_${sp}_2.lnl alt_${sp}_3.lnl
done

for sp in $species
do
  echo "lnL_1" >> null_${sp}_1.lnl
  echo "lnL_2" >> null_${sp}_2.lnl
  echo "lnL_3" >> null_${sp}_3.lnl
  for f in $list
  do
    grep 'lnL' null_${sp}_1/null_${sp}_1_${f}.mlc | awk '{print ($5)}' >> null_${sp}_1.lnl
    grep 'lnL' null_${sp}_2/null_${sp}_2_${f}.mlc | awk '{print ($5)}' >> null_${sp}_2.lnl
    grep 'lnL' null_${sp}_3/null_${sp}_3_${f}.mlc | awk '{print ($5)}' >> null_${sp}_3.lnl
  done

  paste $2 null_${sp}_1.lnl null_${sp}_2.lnl null_${sp}_3.lnl > null_${sp}_all.lnl
  rm null_${sp}_1.lnl null_${sp}_2.lnl null_${sp}_3.lnl
done
