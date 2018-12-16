Picking CCLE cell lines for pathway analysis
================

This is for an R01 application. We want to study cellular signaling using proteomics. The idea is that current knowledge about signaling is very “fuzzy” due to overlapping of signaling components between pathways, e.g. both EGFR and FGFR signal through the MAPK signaling cascade.

    We were criticized about how model cell lines are picked for the study, which was somewhat arbitrary. So I want to see if I can suggest a way to choose the cells based on the expression level of the pathway members. 

RNA-seq data were downloaded from the CCLE. The data are read counts, so I will have to normalize them. I will use the DESeq2 package to do so.

``` r
library(tidyverse)
library(ggthemes)
library(DESeq2)
library(robustbase)
```

``` r
# import the raw count data, there are two rows of unwanted materials, so I will have to skip it
raw.counts <- read.delim("../../CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct",
                         header = TRUE, stringsAsFactors = FALSE,
                         skip = 2)
as.tibble(raw.counts)
```

    ## # A tibble: 56,318 x 1,078
    ##    Name  Description X22RV1_PROSTATE X2313287_STOMACH X253JBV_URINARY~
    ##    <chr> <chr>                 <dbl>            <dbl>            <dbl>
    ##  1 ENSG~ DDX11L1                   0                2                0
    ##  2 ENSG~ WASH7P                 2316             1538             1094
    ##  3 ENSG~ MIR1302-11                5                1               19
    ##  4 ENSG~ FAM138A                   0                1               30
    ##  5 ENSG~ OR4G4P                    0                0                3
    ##  6 ENSG~ OR4G11P                   0                0                2
    ##  7 ENSG~ OR4F5                     0                0                0
    ##  8 ENSG~ RP11-34P13~              63               40               23
    ##  9 ENSG~ CICP27                  165               33               46
    ## 10 ENSG~ AL627309.1             1877             1650              302
    ## # ... with 56,308 more rows, and 1,073 more variables:
    ## #   X253J_URINARY_TRACT <dbl>, X42MGBA_CENTRAL_NERVOUS_SYSTEM <dbl>,
    ## #   X5637_URINARY_TRACT <dbl>, X59M_OVARY <dbl>,
    ## #   X639V_URINARY_TRACT <dbl>, X647V_URINARY_TRACT <dbl>,
    ## #   X697_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>, X769P_KIDNEY <dbl>,
    ## #   X786O_KIDNEY <dbl>, X8305C_THYROID <dbl>, X8505C_THYROID <dbl>,
    ## #   X8MGBA_CENTRAL_NERVOUS_SYSTEM <dbl>, A101D_SKIN <dbl>,
    ## #   A1207_CENTRAL_NERVOUS_SYSTEM <dbl>, A172_CENTRAL_NERVOUS_SYSTEM <dbl>,
    ## #   A204_SOFT_TISSUE <dbl>, A2058_SKIN <dbl>, A253_SALIVARY_GLAND <dbl>,
    ## #   A2780_OVARY <dbl>, A375_SKIN <dbl>,
    ## #   A3KAW_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>, A427_LUNG <dbl>,
    ## #   A498_KIDNEY <dbl>, A4FUK_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>,
    ## #   A549_LUNG <dbl>, A673_BONE <dbl>, A704_KIDNEY <dbl>, ABC1_LUNG <dbl>,
    ## #   ACCMESO1_PLEURA <dbl>, ACHN_KIDNEY <dbl>, AGS_STOMACH <dbl>,
    ## #   ALLSIL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>,
    ## #   AM38_CENTRAL_NERVOUS_SYSTEM <dbl>,
    ## #   AML193_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>,
    ## #   AMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>,
    ## #   AN3CA_ENDOMETRIUM <dbl>, ASPC1_PANCREAS <dbl>, AU565_BREAST <dbl>,
    ## #   BC3C_URINARY_TRACT <dbl>,
    ## #   BCP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>, BCPAP_THYROID <dbl>,
    ## #   BDCM_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>, BEN_LUNG <dbl>,
    ## #   BFTC905_URINARY_TRACT <dbl>, BFTC909_KIDNEY <dbl>,
    ## #   BHT101_THYROID <dbl>, BHY_UPPER_AERODIGESTIVE_TRACT <dbl>,
    ## #   BICR16_UPPER_AERODIGESTIVE_TRACT <dbl>,
    ## #   BICR18_UPPER_AERODIGESTIVE_TRACT <dbl>,
    ## #   BICR22_UPPER_AERODIGESTIVE_TRACT <dbl>,
    ## #   BICR31_UPPER_AERODIGESTIVE_TRACT <dbl>,
    ## #   BICR56_UPPER_AERODIGESTIVE_TRACT <dbl>,
    ## #   BICR6_UPPER_AERODIGESTIVE_TRACT <dbl>,
    ## #   BL41_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>,
    ## #   BL70_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>, BT12_SOFT_TISSUE <dbl>,
    ## #   BT20_BREAST <dbl>, BT474_BREAST <dbl>, BT483_BREAST <dbl>,
    ## #   BT549_BREAST <dbl>, BV173_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>,
    ## #   BXPC3_PANCREAS <dbl>, C2BBE1_LARGE_INTESTINE <dbl>, C32_SKIN <dbl>,
    ## #   C8166_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>,
    ## #   CA46_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>, CADOES1_BONE <dbl>,
    ## #   CAKI1_KIDNEY <dbl>, CAKI2_KIDNEY <dbl>, CAL120_BREAST <dbl>,
    ## #   CAL12T_LUNG <dbl>, CAL148_BREAST <dbl>,
    ## #   CAL27_UPPER_AERODIGESTIVE_TRACT <dbl>, CAL29_URINARY_TRACT <dbl>,
    ## #   CAL33_UPPER_AERODIGESTIVE_TRACT <dbl>, CAL51_BREAST <dbl>,
    ## #   CAL54_KIDNEY <dbl>, CAL62_THYROID <dbl>, CAL78_BONE <dbl>,
    ## #   CAL851_BREAST <dbl>, CALU1_LUNG <dbl>, CALU3_LUNG <dbl>,
    ## #   CALU6_LUNG <dbl>, CAMA1_BREAST <dbl>, CAOV3_OVARY <dbl>,
    ## #   CAOV4_OVARY <dbl>, CAPAN1_PANCREAS <dbl>, CAPAN2_PANCREAS <dbl>,
    ## #   CAS1_CENTRAL_NERVOUS_SYSTEM <dbl>,
    ## #   CCFSTTG1_CENTRAL_NERVOUS_SYSTEM <dbl>, CCK81_LARGE_INTESTINE <dbl>,
    ## #   CFPAC1_PANCREAS <dbl>, CH157MN_CENTRAL_NERVOUS_SYSTEM <dbl>,
    ## #   CHAGOK1_LUNG <dbl>, CHP126_AUTONOMIC_GANGLIA <dbl>,
    ## #   CHP212_AUTONOMIC_GANGLIA <dbl>,
    ## #   CI1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE <dbl>, CJM_SKIN <dbl>,
    ## #   CL11_LARGE_INTESTINE <dbl>, CL14_LARGE_INTESTINE <dbl>, ...

DESeq2 requires a matrix with gene name as row names, and a separated matrix for annotation containing the column names (which is the cell lines in this case).

``` r
# make the coutn data matrix
count.data <- as.matrix(raw.counts[, 3:1078])
row.names(count.data) <- raw.counts$Name         # the row names are the Ensembl IDs
count.data[1:3, 1:3]
```

    ##                   X22RV1_PROSTATE X2313287_STOMACH X253JBV_URINARY_TRACT
    ## ENSG00000223972.4               0                2                     0
    ## ENSG00000227232.4            2316             1538                  1094
    ## ENSG00000243485.2               5                1                    19

``` r
# annotation matrix, this contain a single column named condition, which is the name of the cell lines
coldata <- as.matrix(data.frame(condition = names(raw.counts[, 3:1078]))) 
row.names(coldata) <- names(raw.counts[, 3:1078])   # also make the cell lines as the col.names
head(coldata, 3)
```

    ##                       condition              
    ## X22RV1_PROSTATE       "X22RV1_PROSTATE"      
    ## X2313287_STOMACH      "X2313287_STOMACH"     
    ## X253JBV_URINARY_TRACT "X253JBV_URINARY_TRACT"

``` r
# create the DESeq object
dds <- DESeqDataSetFromMatrix(countData = count.data,
                              colData = coldata,
                              design = ~condition) # design is just the condition in this case
```

    ## converting counts to integer mode

``` r
# then create the normalization factor using estimateSizeFactors()
dds <- estimateSizeFactors(dds)

# retrieve the normalized count using counts(, normalized = TRUE)
nor.CCLE <- raw.counts[, 1:2]
colnames(nor.CCLE) <- c("Ensembl", "Gene.names")
nor.CCLE <- cbind(nor.CCLE, data.frame(counts(dds, normalized = TRUE), stringsAsFactors = FALSE), row.names = NULL)

nor.CCLE[1:4,1:4]
```

    ##             Ensembl Gene.names X22RV1_PROSTATE X2313287_STOMACH
    ## 1 ENSG00000223972.4    DDX11L1         0.00000         2.085369
    ## 2 ENSG00000227232.4     WASH7P      1753.92512      1603.648846
    ## 3 ENSG00000243485.2 MIR1302-11         3.78654         1.042685
    ## 4 ENSG00000237613.2    FAM138A         0.00000         1.042685

``` r
### signal transduction -> by RTK -> EGFR ###
egfr <- read.delim("../../signaling_genes/by_EGFR.tsv", header = TRUE, stringsAsFactors = FALSE)
egfr <- egfr %>% separate(MoleculeName, into = c("Molecule" , "Gene.name"), sep = " ") %>%
  filter(MoleculeType == "Proteins")
```

    ## Warning: Expected 2 pieces. Additional pieces discarded in 2 rows [4, 6].

``` r
egfr <- egfr$Gene.name[!duplicated(egfr$Gene.name)]

### response to stimulations -> stress -> hypoxia
hypoxia <- read.delim("../../signaling_genes/by_hypoxia.tsv", header = TRUE, stringsAsFactors = FALSE)
hypoxia <- hypoxia %>% separate(MoleculeName, into = c("Molecule" , "Gene.name"), sep = " ") %>%
  filter(MoleculeType == "Proteins")
```

    ## Warning: Expected 2 pieces. Additional pieces discarded in 3 rows [2, 3,
    ## 4].

``` r
hypoxia <- hypoxia$Gene.name[!duplicated(hypoxia$Gene.name)]
```

``` r
nor.CCLE <- read.delim("Normalized_CCLE_RNAseq.txt", header = TRUE, stringsAsFactors = FALSE)
```

To see what is the best way of picking the cell lines to study, I want to try Rank sum of the expression, mean expresssion, median and the sum of expression of all genes involved in a particular pathway / stimulation.

``` r
##################################################################
##################### the following parts work with EGFR 
##################################################################

### let's try using rank to see whether it is a good way to differentiate the cell lines ###

# copy the df 1st
rank.expression <- nor.CCLE

# then rank order all of the gene expressions, producing a matrix
rank.expression <- apply(rank.expression[, 3:1078], 2, 
                         FUN = rank, ties.method = "min") # rank will order ascendingly, 
                                                          # so lower expression get a smaller rank number
                                                          # ties.method = "min" put the ties to a smaller number

# add the Ensembl IDs and gene names back to the matrix, produce a df
rank.expression <- cbind(nor.CCLE[, 1:2], rank.expression)

# then select genes that are only in the egfr pathway
egfr.rank <- rank.expression[rank.expression$Gene.names %in% egfr, ]

# add up the rank of the egfr pathway components 
egfr.matrix <- data.frame(Ensembl = rank.expression$Ensembl, stringsAsFactors = FALSE)
egfr.matrix <- data.frame(Rank.sum = colSums(egfr.rank[, 3:1078]), stringsAsFactors = FALSE)
```

``` r
# make a df containig the normalized expression of egfr components
egfr.ccle <- nor.CCLE[nor.CCLE$Gene.names %in% egfr, ]

#####################
### median ##########
#####################

# calculate the medians and add to the egfr.matrix
egfr.matrix$Median <- robustbase::colMedians(as.matrix(egfr.ccle[, 3:1078]))

#####################
### mean ##########
#####################

egfr.matrix$Mean <- colMeans(as.matrix(egfr.ccle[, 3:1078]))

#####################
### median ##########
#####################

egfr.matrix$Sum <- colSums(as.matrix(egfr.ccle[, 3:1078]))
```

``` r
egfr.matrix$Cell.lines <- row.names(egfr.matrix)
egfr.matrix.long <- egfr.matrix %>% gather(key = "Methods", value = "Expressions", 1:4)

ggplot(egfr.matrix.long) +
  geom_point(mapping = aes(x = Cell.lines, y = Expressions)) +
  facet_wrap(~ Methods, scales = "free_y")
```

![](pickCells_files/figure-markdown_github/clean%20up%20and%20visualize-1.png)

``` r
# boxplot
ggplot(egfr.matrix.long) +
  geom_boxplot(mapping = aes(x = Methods, y = Expressions))
```

![](pickCells_files/figure-markdown_github/clean%20up%20and%20visualize-2.png)

``` r
# get the outliner value, this is the max value of the whisker
max.whisk <- quantile(egfr.matrix$Sum, 0.75) + 1.5 * IQR(egfr.matrix$Sum)

# get the cell lines in the outlying region
Cells <- data.frame(Cell.lines = egfr.matrix$Cell.lines[egfr.matrix$Sum >= max.whisk])
Cells$Signal <- "EGFR"
```
