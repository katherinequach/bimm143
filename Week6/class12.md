# Class 12 Pt. 1 (RNASeq Galaxy)
Katherine Quach (A18541014)

- [Section 1: Proportion of G/G in a
  population](#section-1-proportion-of-gg-in-a-population)
- [Section 4: Population Scale Analysis
  \[OPTIONAL\]](#section-4-population-scale-analysis-optional)

## Section 1: Proportion of G/G in a population

Downloaded a CSV file from Ensemble
<https://useast.ensembl.org/Homo_sapiens/Variation/Sample?db=core;r=17:39895045-39895145;v=rs8067378;vdb=variation;vf=959672880#373531_tablePanel>

> Q5: What proportion of the Mexican Ancestry in Los Angeles sample
> population (MXL) are homozygous for the asthma associated SNP (G\|G)?

\[HINT: You can filter the displayed genotypes by entering the
population code MXL. Then either

count those of interest or download a CVS file for this population and
use excel or the R functions read.csv(), and table() to answer this
question\]

Here we read this CSV file

``` r
mxl <-read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head (mxl)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
table(mxl$Genotype..forward.strand.)
```


    A|A A|G G|A G|G 
     22  21  12   9 

``` r
table(mxl$Genotype..forward.strand.) / nrow(mxl) * 100
```


        A|A     A|G     G|A     G|G 
    34.3750 32.8125 18.7500 14.0625 

> Q6. Back on the ENSEMBLE page, use the “search for a sample” field
> above to find the particular sample HG00109. This is a male from the
> GBR population group. What is the genotype for this sample?

Now let’s look at a different population. I picked the GBR.

``` r
gbr <- read.csv("373522-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

Find proportion of GIG

``` r
round(table(gbr$Genotype..forward.strand.) / nrow(gbr) * 100, 2)
```


      A|A   A|G   G|A   G|G 
    25.27 18.68 26.37 29.67 

This variant that is associated with childhood asthma is more frequent
in the GBR populatio than the BKL population.

Let’s now dig into this further.

## Section 4: Population Scale Analysis \[OPTIONAL\]

One sample is obviously not enough to know what is happening in a
population. You are interested in assessing genetic differences on a
population scale. So, you processed about ~230 samples and did the
normalization on a genome level. Now, you want to find whether there is
any association of the 4 asthma-associated SNPs (rs8067378…) on ORMDL3
expression.

How many samples do we have?

``` r
expr <- read.table("genotype.txt")
head(expr)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
nrow(expr)
```

    [1] 462

``` r
table(expr$geno)
```


    A/A A/G G/G 
    108 233 121 

``` r
library(ggplot2)
```

Let’s make a box pot

``` r
ggplot(expr) + 
  aes(geno, exp, fill=geno) +
  geom_boxplot(notch=TRUE)
```

![](class12_files/figure-commonmark/unnamed-chunk-10-1.png)
