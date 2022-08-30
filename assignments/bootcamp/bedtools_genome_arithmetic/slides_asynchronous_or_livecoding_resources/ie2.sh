
# goal: number of SNPs in the 100 genomes file that intersect each gene nad how many unique genes are represented

genefile=/Users/cmdb/data/bed_files/genes.bed
vcffile=/Users/cmdb/data/vcf_files/random_snippet.vcf

bedtools intersect -a $genefile -b $vcffile > intersect_out_ie2.bed

#Number of SNPs
wc -l intersect_out_ie2.bed

#number of unique genes
cut -f 4 intersect_out_ie2.bed | sort | uniq | wc -l



















