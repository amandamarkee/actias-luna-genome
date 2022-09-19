
## **09/19/2022**

- Today I archived my assembly notes from HiPerGator computing cluster, into this repository as an archived text file.
- Moving forward, I will transcribe those methods here, and this document will serve ad the active workflow notes.
- I will link all scripts used in the README file once I've finished the full workflow.

- The summarized methods I use for my genome assembly and annotation follow [Kawahara et al., 2022](https://doi.org/10.46471/gigabyte.64). 
 The steps used in the paper are:
   1. Genome size estimations and genome profiling ([K-Mer Counter (KMC) v.3.1.1](https://github.com/refresh-bio/KMC)) with k-mer size=21.
   2. Sequence assembly and analysis
      - Contig assembling using Hifiasm v0.13-r307
      - Assembly statistics were calculated using [assembly_stats.py](https://github.com/MikeTrizna/assembly_stats)
      - Genome completeness was determined using BUSCO v.5.2.2 (bd10 reference Endopterygota)
      - Contamination was detected using BlobTools v1.0.
   3. Genome annotation
      - Structure Annotation: 
         - RepeatMasker (followed by RepeatModeler2) 
         - species-specific gene model training: BUSCO v.4.1.4
         - predicted genes with the homology-based gene prediction: GeMoMa v1.6.4
         - generate additional ab initio gene predictions: MAKER v3.01.03
      - Functional Annotation:
         - Add functional annotations to the predicted proteins: BlastP
         - Blast2GO for go term annotation

## **09/19/2022; Data Accumulation and File Explanations **

- Samples were sequenced at the University of Florida's Interdisciplinary Center for Biotechnology Research (ICBR)
- Once ICBR has completed HiFi sequencing, the core emailed me regarding my run reports (see run-reports directory)
- To initiate data delivery, I coordinated with ICBR's bioinformatics team following [these instructions](https://biotech.ufl.edu/wp-content/uploads/2022/04/Configuring-Globus-with-ICBR-as-a-HiPerGator-Account-Holder.pdf)
- All raw data for the genomic sequence data (NS2497) can be found in the following HiPerGator directory
  - /blue/kawahara/amanda.markee/ICBR/NS-2497
      - 0000000200 = folder that contains full completed hifi reads w/metadata
      - 0000000205 = folder that contains css / .json  
      - 1_A05 = folder that contains raw data, and ccs reports


