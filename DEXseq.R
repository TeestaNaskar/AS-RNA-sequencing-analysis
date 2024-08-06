#this script is for differential expression for splice site variants
setwd("/sc/arion/scratch/naskat01/AC_RNAseq/00_fastq/script/")
packages <- c("DEXSeq","Rsubread","tibble", "dplyr", "EnhancedVolcano", "edgeR")

#Load libraries
invisible(lapply(packages, library, character.only = TRUE))

load("AS_countdata.RData")


sampledata = sampleTable[sampleTable$Sex=="Male",]
male_bam_files <- paste0(rownames(sampledata), "_sorted.bam")

library(DEXSeq) %>% invisible
dxd = DEXSeqDataSet(
   countData$counts[,colnames(countData$counts) %in% male_bam_files],
   sampleData=sampledata,
   design= ~ sample + exon + Group:exon,
   featureID = exoninfo$ExonID,
   groupID = exoninfo$GeneID,
   featureRanges = exoninfo,
   transcripts = transcripts_l)

# Explore/Inspect the DEXSeq object
head(counts(dxd), 5)
colData(dxd)
split(seq_len(ncol(dxd)), colData(dxd)$exon)
head(rowRanges(dxd),3)
sampleAnnotation(dxd)

### 4.1.4 Normalisation and Dispersion Estimation

dxd = estimateSizeFactors(dxd) 
dxd = estimateDispersions(dxd)

### 4.1.5 Testing for Differential exon Usage

dxd = testForDEU(dxd)
#Estimate fold changes
dxd = estimateExonFoldChanges(dxd, fitExpToVar="Group")
dxr = DEXSeqResults(dxd)
dxr
mcols(dxr)$description
table(tapply(dxr$padj, dxr$groupID, any))
dgene = data.frame(perGeneQValue=perGeneQValue(dxr)) %>% rownames_to_column("groupID")
dexon = dxr %>% data.frame() %>%  dplyr::select(-matches("dispersion|stat|countData|genomicData")) %>% inner_join(dgene) %>% arrange(perGeneQValue) %>% distinct()
writexl::write_xlsx(dexon,"DEXseq_significant_genes.xlsx", col_names = TRUE)
