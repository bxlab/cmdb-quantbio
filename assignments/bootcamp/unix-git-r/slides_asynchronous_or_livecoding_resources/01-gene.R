library( "tidyverse" )

df <- read_tsv( "~/Data/hg38/hg38-gene-metadata-feature.tsv" )
df
View( df )

distinct( df, chromosome_name )
df %>%
  group_by( chromosome_name ) %>%
  summarize( count=n() ) %>%
  arrange( -count ) %>%
  View()

df %>%
  group_by( gene_biotype ) %>%
  summarize( count=n() ) %>%
  arrange( -count )

df_coding <- df %>%
  filter( gene_biotype=="protein_coding" )
df_coding

df_count <- df_coding %>%
  group_by( chromosome_name ) %>%
  summarize( count=n() ) %>%
  arrange( -count ) # %>%
  # View( "df_coding" )

setwd( "~/qbb2024-answers/day1-morning" )
write_tsv( df_count, "counts.tsv" )
