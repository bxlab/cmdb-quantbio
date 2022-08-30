#USAGE: ie3.sh input_vcf_file
#broad goal: What is the most common alternate allele for a reference allele of C

#error free code

grep -v "#" $1 | awk '{if ($4 == "C") {print $5}}' | sort | uniq -c

#USAGE: ie3.sh input_vcf_file nucleotide_of_interest
#refine the goal to work for any nucleotide


#code with error
nucoi=$2
grep -v "#" $1 | awk '{if ($4 == $nucoi) {print $5}}' | sort | uniq -c

