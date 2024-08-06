# 3 Counting Reads
For counting the number of reads mapped to different exons, featureCounts from the Rsubread package was used. 

### 3.1 Load required libraries-
The following packages were used for this part of the analysis. Make sure to install the packages if they are not already installed. 

```{r}
packages <- c("DEXSeq","Rsubread","tibble", "dplyr", "EnhancedVolcano", "edgeR")

#Load libraries
invisible(lapply(packages, library, character.only = TRUE))
```
### 3.2 Load the processed annotation file- 

```{r}
# Load processed annotation file
load("/sc/arion/scratch/naskat01/AC_RNAseq/00_fastq/script/mm10_exon_anno.RData")

### 3.3 Read the bam files obtained in alignment step (after mapping) as input for featureCounts.  
Before running this step, make sure you have the folder containing bam files in the current working directory.
# Read all BAM files as input for featureCounts
filesToCount <- dir("/sc/arion/scratch/naskat01/AC_RNAseq/bam/sortedbam", pattern=".bam$", full.names=T)
countData <- featureCounts(files=filesToCount, annot.ext=anno, isGTFAnnotationFile=FALSE, minMQS=0, useMetaFeatures=FALSE, allowMultiOverlap=TRUE, largestOverlap = TRUE, countMultiMappingReads=FALSE, primaryOnly=TRUE, isPairedEnd=TRUE, nthreads = 12)

# Non-specific filtering: Remove the exons with low counts
isexpr <- rownames(countData$counts)[rowSums(cpm(countData$counts) > 1) >= 3]
countData$counts <- countData$counts[rownames(countData$counts) %in% isexpr, ]
anno <- anno %>% filter(GeneID %in% rownames(countData$counts))

# Remove genes with only 1 site and NA in geneIDs
dn <- anno %>% group_by(GeneID) %>% summarise(nsites=n()) %>% filter(nsites > 1 & !is.na(GeneID))
anno <- anno %>% filter(GeneID %in% dn$GeneID)
countData$counts <- countData$counts[rownames(countData$counts) %in% anno$GeneID, ]

# 4. Differential Splicing and Exon usage analysis

## 4.1 Using DEXSeq pipeline for differential exon analysis

### 4.1.1 Load library and create a sample table to define the experimental design

library(DEXSeq) %>% invisible
sampleTable = data.frame(
   row.names = c( "3454", "187D", "3A5C", 
                 "190E", "6F1A", "620D", "0C7A", "28", "4C24", "4D01", "6D6F", "6552", "1B11", "3425", "4960", "2665", "342E", "6B18", "366C", "6C61" ),
   Sex = rep(c("Male","Female"),c(10,10)),Group = rep(c("Wildtype", "FDD-Knockin", "Wildtype", "FDD-Knockin"),c(5,5,5,5)),condition = rep(c("Male_Wildtype", "Male_FDD-Knockin","Female_Wildtype", "Female_FDD-Knockin"),c(5,5,5,5)), libType = rep(c("paired-end")))


### 4.1.2 Prepare the exon information file

exoninfo = anno[anno$GeneID %in% rownames(countData$counts),]
exoninfo <- GRanges (seqnames=  anno$Chr, 
ranges= IRanges (start=anno$Start, end=anno$End, width = anno$Width), strand = Rle(anno$Strand))
mcols(exoninfo)$TranscriptIDs <- anno$TranscriptIDs
mcols(exoninfo)$Ticker <- anno$Ticker
mcols(exoninfo)$ExonID <- anno$ExonID
mcols(exoninfo)$n <- anno$n
mcols(exoninfo)$GeneID <- anno$GeneID

transcripts_l = strsplit(exoninfo$TranscriptIDs,  "\\,")

#Save the countData, sampleTable and exoninfo and transcripts in one object for further use 
save(countData, sampleTable, exoninfo, transcripts_l, file="AS_countdata.RData")
