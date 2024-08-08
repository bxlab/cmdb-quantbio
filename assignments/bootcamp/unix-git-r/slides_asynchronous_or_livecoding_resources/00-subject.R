library( "tidyverse" )

df <- read_tsv( "~/Data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt" )
df

distinct( df, SEX )
distinct( df, AGE )

df %>%
  distinct( AGE )

df %>%
  group_by( AGE ) %>%
  summarize( n() )

df %>%
  filter( AGE=="70-79" )

df_oldest <-df %>%
  filter( AGE=="70-79" )
df_oldest

setwd( "~/qbb2024-answers/day1-morning" )
write_tsv( df_oldest, "oldest.tsv" )

df2 <- read_tsv( "oldest.tsv" )
df2
