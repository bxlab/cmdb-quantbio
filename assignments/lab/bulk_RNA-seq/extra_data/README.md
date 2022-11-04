Note any pre-loaded data used in this README

```
wget "https://www.dropbox.com/s/hxjvua05abz8zqb/all_annotated.csv?dl=1" -O all_annotated.csv
cut -d, -f1,3-7,11-15 all_annotated.csv > dros_gene_expression.csv
sed -i '' 's/14A/14/g' dros_gene_expression.csv
```
