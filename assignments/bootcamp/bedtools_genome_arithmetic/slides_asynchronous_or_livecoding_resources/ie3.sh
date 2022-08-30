
#USAGE: ie3.sh input_vcf_file nucleotide_of_interest

nucoi=$2
grep -v "#" $1 | awk '{if ($4 == $nucoi) {print $5}}' | sort | uniq -c