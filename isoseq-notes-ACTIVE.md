## _NOTE: This IsoSeq analysis uses our newly assembeled and annotated genome_ ##

## Introduction
This workflow describes modified methods used for the [IsoSeq bulk workflow](https://isoseq.how/clustering/cli-workflow.html) on five multiplexed _Actias luna_ samples. The goal of this workflow is to align PacBio IsoSeq long reads back to my newly assembled and annotated genome for the Luna moth, to look at differential IsoForm expression across life stages, specific to silk production.

## **04/7/2022; Action Items for IsoSeq Analysis**

https://isoseq.how/clustering/high-level-workflow.html

1. Run lima seperately on all samples from Run02 (1_H06 on the cluster)
2. Run refine seperately on all samples from Run02 (1_H06 on the cluster)
3. Check lima and refine step for any glaring errors. (ie. in lima summary, these stats should not be over 40%)

```
Undesired 5p--5p pairs        : xxx (xx%) <- Only with --isoseq
Undesired 3p--3p pairs        : xxx (xx%) <- Only with --isoseq
Undesired single side         : xxx (xx%) <- Only with --isoseq
Undesired no hit              : xxx (xx%) <- Only with --isoseq
```

4. Repeat steps 1-3 for Run01 (1_E06 on the cluster)
5. Combine each sample from each run into a .fofn file (each fofn has 2 file names within it, from Run01-Run02)
6. Run clustering step, and check back in with CGS here to ask if each .fofn is run seperately on clustering step.
7. Run mapping with pbmm2 to assembled and annotated Luna genome 
8. Run collapse step with isoseq3 collapse function.
9. Retrieve collapse statistics with [this python script](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step#summarizing-and-plotting-transcriptexonintron-stats-after-collapse), and visualize statistics using this R script.

```
R SCRIPT FOR VISUALIZING ISOSEQ COLLAPSE STATS

# --------------------------------------
library(ggplot2)
library(dplyr)
library(ggthemes)
# --------------------------------------
# https://github.com/Magdoll/cDNA_Cupcake/blob/master/beta/plot_simple_stats_post_collapse.R
#args <- commandArgs(trailingOnly = TRUE)
#input.file1 <- args[1]
#input.file2 <- args[2]
#PREFIX <- args[3]

setwd("~/Dropbox (UFL)/GigaByte_SilkGenomes_Data/IsoSeq")
x.1 <- read.table("PiV3_iso.collapsed.simple_stats.txt", sep='\t',header=T)
x.2 <- read.table("PiV3_iso.collapsed.exon_stats.txt",sep='\t',header=T)

par(mfrow=c(3,3))
x <- x.1
# plot length histogram
ggplot(x, aes(length)) + geom_histogram(binwidth=200, fill ='lightblue', color='black') + xlab("Transcript Length (bp)") + ylab("Number of Unique Transcripts") + labs(title="Mapped Unique Transcript Lengths")  + theme_tufte() + xlim(c(0,10000))
ggsave("PiV3_Iso.Rplot.histogram.png", width=6, height=4.5)

# plot genomic length histogram
ggplot(x, aes(genomic_length/1000)) + geom_histogram(binwidth=100, fill='lightgreen', color='black') + xlab("Transcript Genomic Length (kb)") + ylab("Number of Unique Transcripts") + labs(title="Mapped Unique Transcript Genomic Lengths")  + theme_few()
ggsave("PiV3_Iso.Rplot.genomic_len_histogram.png", width=6, height=4.5)

# plot number of exons
x$exon_cat <- "1"
x[x$num_exon>=2,"exon_cat"] <- "2-5"
x[x$num_exon>=10,"exon_cat"] <- "10-20"
x[x$num_exon>20,"exon_cat"] <- ">20"
x$exon_cat <- factor(x$exon_cat, levels=c("1","2-5","10-20",">20"))
ggplot(x, aes(exon_cat)) + geom_bar(fill='darkgreen') + xlab("Number of Exons") + ylab("Number of Transcripts") + theme_tufte() + labs(title="Number of Exons Per Transcript")
ggsave("PiV3_Iso.Rplot.num_exons.png", width=6, height=4.5)

# number oof isoforms per locus
t <- x %>% group_by(locus) %>% summarise(count=n())
t$isocat <- "1"
t[t$count>=2, "isocat"] <- "2-5"
t[t$count>=6, "isocat"] <- "6-10"
t[t$count>10, "isocat"] <- "11-20"
t[t$count>20, "isocat"] <- ">20"
t$isocat <- factor(t$isocat, levels=c("1","2-5","6-10","11-20",">20"))
ggplot(t, aes(isocat)) + geom_bar(fill='pink') + xlab("Number of Isoforms") + ylab("Number of Locus") + theme_tufte() + labs(title="Number of Isoforms Per Locus")
ggsave("PiV3_Iso.Rplot.num_isoforms.png", width=6, height=4.5)

# ---------------
# read exon/intron size
# ---------------
x <- x.2
x$exon_size_cat <- "<100bp"
x[x$exon_size>=100,"exon_size_cat"] <- "100-200"
x[x$exon_size>200,"exon_size_cat"] <- "201-500"
x[x$exon_size>500,"exon_size_cat"] <- "501-1000"
x[x$exon_size>1000,"exon_size_cat"] <- "1001-10000"
x[x$exon_size>10000,"exon_size_cat"] <- ">10000"
x$exon_size_cat <- factor(x$exon_size_cat, levels=c("<100bp","100-200","201-500","501-1000","1001-10000",">10000"))
ggplot(x, aes(exon_size_cat)) + geom_bar() + theme_tufte() + xlab("Exon Lengths (bp)") + ylab("Count") + labs(title="Distribution of Exon Lengths")
ggsave("PiV3_Iso.Rplot.num_exon_sizes.png", width=6, height=4.5)

x1 <- subset(x, !is.na(intron_size))
x1$intron_size_cat <- "<100bp"
x1[x1$intron_size>100, "intron_size_cat"] <- "100-1000"
x1[x1$intron_size>1000, "intron_size_cat"] <- "1000-10000"
x1[x1$intron_size>10000, "intron_size_cat"] <- "10kb-100kb"
x1[x1$intron_size>100000, "intron_size_cat"] <- "100kb-1Mb"
x1[x1$intron_size>1000000, "intron_size_cat"] <- ">1Mb"
x1$intron_size_cat <- factor(x1$intron_size_cat, levels=c("<100bp","100-1000","1000-10000", "10kb-100kb", "100kb-1Mb", ">1Mb")
)
ggplot(x1, aes(intron_size_cat)) + geom_bar() + theme_tufte() + xlab("Intron Lengths (bp)") + ylab("Count") + labs(title="Distribution of Intron Lengths")
ggsave("PiV3_Iso.Rplot.num_intron_sizes.png", width=6, height=4.5)
```

## 04/10/2022; Run IsoSeq Lima Step 

First, we must make an array script which uses the following line of code in the format of lima input file, primer name file, and output name. Note, here we just use --isoseq flag. When making the array script, we should use 8 cores and 4 hours per sample, and utilize environmental variables for input, primer name file and output. See example code from workflow below.

```
lima movieX.hifi_reads.bam barcoded_primers.fasta movieX.fl.bam --isoseq
```



