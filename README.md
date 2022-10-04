# _Actias luna_ Genome Project
Genome assembly and annotation workflow for the luna moth, _Actias luna_. This repository format was modified from YiMing Weng's [genome assembly and annotation repository](https://github.com/yimingweng/Kely_genome_project) for the tomato pinworm _Keiferia lycopersicella_.

## Background
- Moths and butterflies (Lepidoptera) are one of the dominant groups of silk producers in the insect world, but research is severely limited to model organisms such as the domestic silkworm moth (_Bombyx mori_) due to it's economic importance.
- The luna moth (_Actias luna_) is a charismatic native species found in Florida, known to produce silk during throuhgout it's lifespan.
- A comprehensive document about this species can be found on the [Featured Creatures website](https://entnemdept.ufl.edu/creatures/misc/moths/luna_moth.htm) published by University of Florida.
- For my Master's thesis project, my research aims to characterize silk production in the luna moth on a genomic and phenotypic scale. 
- The focus of this repository will be to provide the workflow for assembling and annotating a high-quality reference genome for this species, in order to resolve ecological/evolutionary questions pertaining to it's silk production.
<img width="2560" alt="methods" src="https://user-images.githubusercontent.com/56971761/190717558-06e504b8-9414-42d0-a1ba-eb2ffab90cc9.png">

<br />

## General Workflow
This section describes the methods I used for both assembling and annotating the luna moth genome from PacBio High Fidelity (HiFi) reads. Methods in this repository follow those described in [Kawahara et al., 2022](https://gigabytejournal.com/articles/64)

### Read Quality Assessment [(FastQC)](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md#09192022-raw-read-quality-assessment-with-fastqc)

### Contig Assembly [(hifiasm)](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md#09192022-genome-assembly-with-hifiasm)

### Assembly Quality Assessment [(BUSCO)](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md#10032022-genome-completeness-with-busco)

### Polishing and Purging Haplotigs [(assemblystats.py)](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md#10032022-genome-assembly-quality-assessment-with-assemblystatspy)

### Genome Size Estimation [(KMC and GenomeScope)](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md#10032022-genome-size-estimation-kmer-with-kmc-and-genomescope)

### Contamination Filtering [(BlobTools)](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md#10042022-contaminaiton-filtering-with-blobtools)

### Feature Annotation (TBD: Maker/Baker?)

### Functional Annotation (TBD: Maker/Baker?)



<br />

## References

Kawahara, A.Y., Storer, C.G. Markee, A., Heckenhauer, J., Powell, A., Plotkin, P., Hotaling, S., Cleland, T.P., Dikow, R.B., Dikow, T., Kuranishi, R.B., Messcher, R., Pauls, S.U., Stewart, R.J., Tojo, K., Frandsen, P.B. 2022. Long-read HiFi sequencing correctly assembles repetitive heavy fibroin silk genes in new moth and caddisfly genomes. Gigabyte, 2022 doi: 10.46471/gigabyte.64.
