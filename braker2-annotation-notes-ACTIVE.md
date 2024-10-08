## _STOP: Please see the [assembly](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md) note for the first half of the instructions before proceeding_ ##

## **10/9/2022; Annotation – RepeatModeler2**
Now that my genome assembly is complete, we have a high-quality and accurate assembly to begin building our annotation from. The ultimate goal is to perform both feature and functional annotation, but we will begin with feature annotation. For my annotation, I will be following [Dr. YiMing Weng's](https://github.com/yimingweng/Kely_genome_project/blob/main/note.md#09072022) following the BRAKER2 protocol.

First, we need to run [RepeatModeler2](https://www.pnas.org/doi/10.1073/pnas.1921046117) to identify and model for repeat elements in the genome.

Navigate to your working directory for annotation, and execute the following script on your genome assembly using the following code:
```
[amanda.markee@login1 aluna_annotation]$ pwd
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation
sbatch repeatmodeler.sh /blue/kawahara/amanda.markee/Aluna_genome/blobtools/aluna_final_assembly.fasta aluna_repeatmodeler
```

```
#!/bin/bash

#SBATCH --job-name=aluna_repeatmodeler2.sh
#SBATCH -o aluna_repeatmodeler2.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 30

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmodeler/2.0
module load seqkit/2.0.0

genome=${1} # full path to the genome assembly
outprefix=${2}

# build the RM2 database for the genome
BuildDatabase -name aluna_genome ${genome}

# run RepeatModeler with the database
RepeatModeler -database aluna_genome -pa 10 -LTRStruct >& ${outprefix}.out

cat aluna_genome-families.fa | seqkit fx2tab | awk '{ print "Aluna_1.0_"$0 }' | seqkit tab2fx > aluna_genome-families.prefix.fa
cat aluna_genome-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > aluna-genome-families.prefix.fa.known
cat aluna_genome-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > aluna-genome-families.prefix.fa.unknown
```
</br>

## **11/7/2022; Feature Annotation (Part 1) – RepeatMasker**

Once RepeatModeler2 is complete, we can move on to masking the genome with the RepeatModeler2 output information. Here, I will be following reccomendations from [Dr.Card](https://darencard.net/blog/2022-07-09-genome-repeat-annotation/) about how to comprehensively mask a genome with repeat sequences. There are four steps, and I am following [Dr. YiMing Weng's](https://github.com/yimingweng/Kely_genome_project/blob/main/note.md#09072022-1) four seperate scripts below.

We will use the output of RepeatModeler2 to run RepeatMasker. Including the repeats from RepeatModeler2, I will use 4 different peices of evidence to mask the repeat regions for the genome: (1) mask the simple and short repeats, detected by RepeatMasker; (2) mask repeats based on existing databases (Repbase, Lepidoptera database); (3) mask genome based on the output of RepeatModeler2; (4) calculate the percentage of hardmasking vs softmasking. Note that I used soft-masking for all the repeat regions.

</br>

## Step 1: Feature Annotation – Mask Simple Repeats

```
cd /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_masker
sbatch aluna_repeatmask_step1.slurm
```
```
###########################  script content  ###########################

#!/bin/bash
#SBATCH --job-name=aluna_repeatmask_step1
#SBATCH -o aluna_repeatmask_step1.out
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir aluna_repeatmasker_step1

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-xsmall \
-e RMBlast \
-gff \
-noint \
-no_is \
-dir aluna_repeatmasker_step1 \
/blue/kawahara/amanda.markee/Aluna_genome/blobtools/aluna_final_assembly.fasta &> aluna_repeatmasker_step1.out
```

## Step 2: Feature Annotation – Mask Repeats Based on Existing Databases

Once we mask the simple repeat elements, we will use the existing database [Repbase](https://www.girinst.org/repbase/) to continue building our gene model. First we navigate too our RepeatMasker directory, then run the following script in that directory.

```
cd /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_masker
sbatch aluna_repeatmask_step2.slurm
```

```
###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=aluna_repeatmask_step2
#SBATCH -o aluna_repeatmask_step2.out
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir aluna_repeatmasker_step2

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-species Lepidoptera \
-dir aluna_repeatmasker_step2 \
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_masker/aluna_repeatmasker_step1/aluna_final_assembly.fasta.masked &> aluna_repeatmasker_step2.out

mv aluna_repeatmask_step2.out out_files
```

## Step 3: Feature Annotation – Mask Genome Based on RepeatModeler2 Output

Next, we use the output from RepeatModeler2 (from the beginning of the annotation workflow) to continue masking repeat regions. 

```
cd /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_masker
sbatch aluna_repeatmask_step3.slurm
```

```
###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=aluna_repeatmask_step3
#SBATCH -o aluna_repeatmask_step3.out
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir aluna_repeatmasker_step3

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-lib /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_modeler/aluna-genome-families.prefix.fa.known \
-dir aluna_repeatmasker_step3 \
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_masker/aluna_repeatmasker_step2/aluna_final_assembly.fasta.masked.masked &> aluna_repeatmasker_step3.out

mv aluna_repeatmask_step3.out aluna_repeatmasker_step3
```

## Step 4: Feature Annotation – Calculate Masking Percentage
The final step is to calculate the masking percentage to determine the soft masking rate vs. the hard masking rate. I used the following script.

```
cd /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_masker/aluna_repeatmasker_step3
bash /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/repeat_masker/aluna_repeatmasker_step3/maskrate.sh aluna_final_assembly.fasta.masked.masked.masked
```

```
###########################  script content  ###########################
#!/bin/bash

input=${1}

total=$(cat ${input} | grep -v ">" | wc -c)
softmask=$(cat ${input} | grep -v ">"| tr -dc a-z | wc -c)
hardmask=$(cat ${input} | grep -v ">"| tr -dc "N" | wc -c)
softrate=$(echo 100*$softmask/$total | bc -l | grep -o  "^.....")
hardrate=$(echo 100*$hardmask/$total | bc -l | grep -o  "^.....")

if [[ $softmask == 0 ]]
then
echo -e "softmasking rate is 0%"
else
echo -e "softmasking rate is ${softrate}%"
fi

if [[ $hardmask == 0 ]]
then
echo -e  "hardmasking rate is 0%"
else 
echo -e  "hardmasking rate is ${hardrate}%"
fi
```

softmasking rate is **40.58%** | hardmasking rate is 0%


## **12/5/2022; Feature Annotation – Building the Gene Model in FunAnnotate**

During a feature annotation, we can use multiple pieces of evidence to build and support the strongest and most accurate gene model. It's best to to have RNAseq data from the same species as your genome when building your model, so that the splicing sites and the gene features can be better caught by the predictor. If you do not have RNAseq data available for your species, a suboptimal alternative is to predict the gene model based on the ab initio (signal sensors and content sensors) and the protein sequence data from other insect species. 

All FunAnnotate work will be done in this directory
```
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/funannotate
```

First, we must set up Funannotate using the script below to build the appropriate environment in HiPerGator.
```
#!/bin/bash
#SBATCH --job-name=funannotate_setup
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10gb
#SBATCH --time=12:00:00
#SBATCH --output=funannotate_%j.log
#SBATCH --account=plantpath
#SBATCH --qos=plantpath

pwd; hostname; date

module load funannotate
module load funannotate/1.8.13

if [[ ! -d config ]]; then
module purge; module load augustus
rsync -a $HPC_AUGUSTUS_CONF .
fi

export AUGUSTUS_CONFIG_PATH=$(pwd)/config
export FUNANNOTATE_DB=$(pwd)/funannotate_db

funannotate setup -i all --update -d ./funannotate_db

date
```

Next, I want to start prepping my assembly by training using existing A.luna RNAseq reads, and my IsoSeq transcriptome. HiPerGator support assisted in writing the script below to run the training step:

```
#!/bin/bash
#SBATCH --job-name=fannotate_train
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=120:00:00
#SBATCH --output=fannotate_%j.log
date;hostname;pwd

module purge; module load funannotate/1.8.13

funannotate train -i masked_genome.fasta -s 2017RNALibPool01-7_S7_L001_R2_001.fastq.gz --pacbio_isoseq instar4_isoseq.fa --species "Actias luna" -o out

date

```

## **12/7/2022; Feature Annotation – Building the Gene Model in FunAnnotate (cont.)**

FunAnnotate uses Evidence Modeler to combine ab initio gene model predictions with evidence (transcripts or proteins) aligned to the genome. Therefore, the evidence that you supply at runtime for --transcript_evidence and --protein_evidence are important. By default, funannotate will use the UniProtKb/SwissProt curated protein database for protein evidence. However, you can specify other forms of protein evidence, perhaps from a well-annotated closely related species, using the --protein_evidence option. Multiple files can be passed to both --transcript_evidence or --protein_evidence by separating the files by spaces, for example:

```
funannotate predict --input genome.fa --species "Awesome species" --transcript_evidence --pacbio_isoseq isoseq.fasta myESTs.fa \
    -o output --protein_evidence closely_related.fasta $FUNANNOTATE_DB/uniprot_sprot.fasta
```

For the pipeline I am using, FunAnnotate uses the following steps for annotation:
1. Align Transcript Evidence to genome using minimap2
2. Align Protein Evidence to genome using Diamond/Exonerate.
3. Parse BAM alignments generating hints file
4. Parse PASA gene models and use to train/run Augustus, snap, GlimmerHMM
5. Extract high-quality Augustus predictions (HiQ)
6. Run Stringtie on BAM alignments, use results to run CodingQuarry
7. Pass all data to Evidence Modeler and run
8. Filter gene models (length filtering, spanning gaps, and transposable elements)
9. Predict tRNA genes using tRNAscan-SE
10. Generate an NCBI annotation table (.tbl format)
11. Convert to GenBank format using tbl2asn
12. Parse NCBI error reports and alert user to invalid gene models


## **1/11/2023; Feature Annotation – Building the Gene Model in FunAnnotate (cont.)**

The first step to the FunAnnotate pipeline after setup is to train the assembly. I used the following script for training, with RNA short read sequences, and PacBio Isoseq sequences for evidence.

```
#!/bin/bash
#SBATCH --job-name=fannotate_train
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=120:00:00
#SBATCH --output=fannotate_%j.log
date;hostname;pwd

module purge; module load funannotate/1.8.13

funannotate train -i masked_genome.fasta -s 2017RNALibPool01-7_S7_L001_R2_001.fastq.gz --pacbio_isoseq instar4_isoseq.fa --species "Actias luna" -o out

date
```

The next step is the gene prediction step (1/11/2023), using the following script:

```
#!/bin/bash
#SBATCH --job-name=funannotate_predict
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128gb
#SBATCH --time=120:00:00
#SBATCH --output=funannotate_predict_%j.log

date;hostname;pwd

module purge; module load funannotate/1.8.13

funannotate predict -i masked_genome.fasta -o predict_out --species "Actias luna" \
    --transcript_evidence trinity.fasta.clean --rna_bam transcript.alignments.bam --pasa_gff Actias_luna_pasa.valid_blat_alignments.gff3

date
```

However this is where I am experiencing issues with the protein database that FunAnnotate uses, Diamond/Exonerate. For some reason, during the gene prediction portion, the pipeline is able to complete steps 1-2, but provides the following error during what I believe to be the protein evidence training step:

```
[amanda.markee@login1 archive]$ cat funannotate_predict_55173243.log 
Mon Jan  9 17:25:03 EST 2023
c0709a-s3.ufhpc
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/funannotate/predict
-------------------------------------------------------
[Jan 09 05:25 PM]: OS: Red Hat Enterprise Linux Server 7.9, 64 cores, ~ 527 GB RAM. Python: 3.8.13
[Jan 09 05:25 PM]: Running funannotate v1.8.13
[Jan 09 05:25 PM]: Parsed training data, run ab-initio gene predictors as follows:
  Program      Training-Method
  augustus     pasa           
  genemark     selftraining   
  glimmerhmm   pasa           
  snap         pasa           
[Jan 09 05:26 PM]: Loading genome assembly and parsing soft-masked repetitive sequences
[Jan 09 05:26 PM]: Genome loaded: 147 scaffolds; 533,189,247 bp; 41.40% repeats masked
[Jan 09 05:26 PM]: Aligning transcript evidence to genome with minimap2
[Jan 09 05:27 PM]: Found 38,028 alignments, wrote GFF3 and Augustus hints to file
[Jan 09 05:27 PM]: Extracting hints from RNA-seq BAM file using bam2hints
[Jan 09 05:27 PM]: Mapping 554,221 proteins to genome using diamond and exonerate
[Jan 09 05:50 PM]: Found 207,366 preliminary alignments with diamond in 0:21:24 --> generated FASTA files for exonerate in 0:00:57
[Jan 09 07:07 PM]: Exonerate finished in 1:17:13: found 4,201 alignments
[Jan 09 07:07 PM]: Running GeneMark-ES on assembly
[Jan 10 12:23 AM]: 60,834 predictions from GeneMark
[Jan 10 12:23 AM]: Filtering PASA data for suitable training set
[Jan 10 12:23 AM]: CMD ERROR: diamond blastp --query augustus.training.proteins.fa --db aug_training.dmnd --more-sensitive -o aug.blast.txt -f 6 qseqid sseqid pident --query-cover 80 --subject-cover 80 --id 80 --no-self-hits
[Jan 10 12:23 AM]: diamond v2.0.15.153 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org
Please cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)

#CPU threads: 64
Scoring parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)
Temporary directory: 
#Target sequences to report alignments for: 25
Opening the database...  [0s]
Error: Incomplete database file. Database building did not complete successfully.
```

Because I am having too many issues getting FunAnnotate to work correctly, I am switching to the BRAKER2 pipeline. This will include all annotation steps thus far until building the gene model (ie. RepeatMasker and RepeatModeler steps are still relevant). 

## 1/16/2023; Feature Annotation – Building the Gene Model in BRAKER2

To build the feature (also known as structural) annotation in BRAKER2, I will be following the pipeline provided by [Dr. Keating Godfrey's](https://github.com/rkeatinggodfrey/Hyles_lineata_genome/tree/main/Annotation) annotation of the Hyles liniata genome. I created my softmasked genome using RepeatMasker and RepeatModeler, and now I am moving forward with using hits from a protien database, RNA seq data for _A.luna_, as well as PacBio IsoSeq transcriptome data as evidence for building my gene model. Because Dr. Godfrey's annotation notes do not use long read IsoSeq data, I will first conduct her pipeline entirely without my IsoSeq data, and then use use a [slightly modified annotation protocol](https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md) which incorporates long-read data, at a later date.

![Screen Shot 2023-01-16 at 1 37 42 PM](https://user-images.githubusercontent.com/56971761/212746627-d47f8945-a360-4bc6-91a5-6293528f4751.png)

## 1/16/2023; Feature Annotation – BRAKER2 setup

Resources:

Running BRAKER with proteins of any evolutionary distance:
Running BRAKER with RNA-seq data: https://github.com/Gaius-Augustus/BRAKER#braker-with-rna-seq-data
Output files: https://github.com/Gaius-Augustus/BRAKER#output-of-braker

All BRAKER2 annotation will be done in the following directory:

```
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2
```

Below is a complete list of allinput files and their names used in this annotation:

    1. masked_genome.fasta - genome sequence in FASTA format that has been softmasked for repeats

    2. transcript.alignments.bam - spliced alignments of short-read RNA-Seq in BAM format (RNA-Seq hints can be used instead of the BAM file, see the  documentation of BRAKER for usage information)

    3. B_mori_protein.fasta - a large database of protein sequences in FASTA format (e.g. a suitable OrthoDB partition)

    4.  instar1.fastq,instar2.fastq,instar3.fastq,... - list of assembled subread libraries from long-read RNA-Seq as FASTQ files (protocol has only been tested with PacBio ccs reads)

Note: I will be using instar subreads (input file #4) that came off of PacBio Sequel IIe, which are in bam format. I will be using [this protocol](https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html) to convert bam to fastq. 


## Step 1 of BRAKER2 Feature Annotation – Running with protein sequences
## (a) Retrieve protein sequences

### Retrieve protein sequences from a well-annotated, closely related species

I downloaded B. mori protein sequences into the folder where I am running BRAKER2

Script to download from NCBI:

```
#!/bin/bash
#SBATCH --job-name=ncbi_download
#SBATCH -o %A_%a.220701_NCBI_Download.out
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 1
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 02:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

module load edirect/12.2

## Retrieve Bombyx mori protein and mRNA

esearch -db protein -query "Bombyx mori [ORGN]" | efetch -format fasta > B_mori_protein.fasta
```

Check how many proteins are in this file: 
```
grep ">" B_mori_protein.fasta | wc -l
```

35408 proteins are present. This file will now be used as protein evidence for BRAKER2


## (b) Run [ProtHint](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) to create protein gff file.

Use the following script execution on the both the masked genome and protein fasta as input files:

```
sbatch -J Al_ProtHint prothint.sh /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/masked_genome.fasta /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/B_mori_protein.fasta
```

```
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

dates;hostname;pwd

genome=${1}
protein=${2}

module load prothint/2.6.0
module load genemark_es/4.69

prothint.py --threads ${SLURM_CPUS_ON_NODE:-1} ${genome} ${protein}
```

## (c) Run BRAKER2 with protein evidence

Now that the ProtHint .gff file is completed, we can use this as protein evidence for runnning BRAKER2. Use the following command to execute the code below it.

```
sbatch -J Al_prot_braker2 Al_braker2_protein.sh /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/masked_genome.fasta /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/prothint_augustus.gff Actias_luna
```
```
#!/bin/bash
#SBATCH --job-name=%j_Al_braker2_prot
#SBATCH --output=%j_Al_braker2_prot.log
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1}
protein_gff=${2}
species=${3}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/rkeating.godfrey/Hyles_lineata_genome/Hl_busco/Augustus/config \
--genome=${genome} --species ${species} --hints=${protein_gff} --softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```

## Step 2 of BRAKER2 Feature Annotation – Running with RNAseq data

Directories for running BRAKER with RNAseq data:
```
### raw RNAseq data files
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/trimmomatic/RAW_rnaseq 

### trimmed RNAseq data files
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/trimmomatic/trimmed 

### parent folder for scripts and outputs related to braker
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2 
```

List of RNAseq files and important metadata:
```
2017RNALibPool02-10_S21_L002_R1_001.fastq.gz    2017RNALibPool02-10_S21_L002_R2_001.fastq.gz        luna_10_head		Pupa
2017RNALibPool02-11_S22_L002_R1_001.fastq.gz    2017RNALibPool02-11_S22_L002_R2_001.fastq.gz        luna_10_thorax		Pupa
2017RNALibPool03-1_S23_L003_R1_001.fastq.gz     2017RNALibPool03-1_S23_L003_R2_001.fastq.gz         luna_11_head		Pupa
2017RNALibPool03-2_S24_L003_R1_001.fastq.gz     2017RNALibPool03-2_S24_L003_R2_001.fastq.gz         luna_11_thorax		Pupa
2017RNALibPool03-3_S25_L003_R1_001.fastq.gz     2017RNALibPool03-3_S25_L003_R2_001.fastq.gz         luna_11_abdomen		Pupa
2017RNALibPool03-4_S26_L003_R1_001.fastq.gz     2017RNALibPool03-4_S26_L003_R2_001.fastq.gz         luna_12_head		Pupa
2017RNALibPool03-5_S27_L003_R1_001.fastq.gz     2017RNALibPool03-5_S27_L003_R2_001.fastq.gz         luna_12_thorax		Pupa
2017RNALibPool03-6_S28_L003_R1_001.fastq.gz     2017RNALibPool03-6_S28_L003_R2_001.fastq.gz         luna_12_abdomen		Pupa
2017RNALibPool03-7_S29_L003_R1_001.fastq.gz     2017RNALibPool03-7_S29_L003_R2_001.fastq.gz         luna_14_head		Adult
2017RNALibPool03-8_S30_L003_R1_001.fastq.gz     2017RNALibPool03-8_S30_L003_R2_001.fastq.gz         luna_14_thorax		Adult
2017RNALibPool03-9_S31_L003_R1_001.fastq.gz     2017RNALibPool03-9_S31_L003_R2_001.fastq.gz         luna_14_abdomen		Adult
2017RNALibPool03-10_S32_L003_R1_001.fastq.gz    2017RNALibPool03-10_S32_L003_R2_001.fastq.gz        luna_15_head		Adult
2017RNALibPool03-11_S33_L003_R1_001.fastq.gz    2017RNALibPool03-11_S33_L003_R2_001.fastq.gz        luna_15_thorax		Adult
2017RNALibPool03-12_S34_L003_R1_001.fastq.gz    2017RNALibPool03-12_S34_L003_R2_001.fastq.gz        luna_15_abdomen		Adult
2017RNALibPool04-1_S35_L004_R1_001.fastq.gz     2017RNALibPool04-1_S35_L004_R2_001.fastq.gz         luna_16_wholebody	        1st instar
2017RNALibPool01-1_S1_L001_R1_001.fastq.gz      2017RNALibPool01-1_S1_L001_R2_001.fastq.gz          luna_01_head		4th instar
2017RNALibPool01-2_S2_L001_R1_001.fastq.gz      2017RNALibPool01-2_S2_L001_R2_001.fastq.gz          luna_01_thorax		4th instar
2017RNALibPool01-3_S3_L001_R1_001.fastq.gz      2017RNALibPool01-3_S3_L001_R2_001.fastq.gz          luna_01_abdomen		4th instar
2017RNALibPool01-4_S4_L001_R1_001.fastq.gz      2017RNALibPool01-4_S4_L001_R2_001.fastq.gz          luna_02_head		4th instar
2017RNALibPool01-5_S5_L001_R1_001.fastq.gz      2017RNALibPool01-5_S5_L001_R2_001.fastq.gz          luna_02_thorax		4th instar
2017RNALibPool07-1_S68_L007_R1_001.fastq.gz     2017RNALibPool07-1_S68_L007_R2_001.fastq.gz         luna_02_abdomen		4th instar
2017RNALibPool01-6_S6_L001_R1_001.fastq.gz      2017RNALibPool01-6_S6_L001_R2_001.fastq.gz          luna_03_wholebody           6-eggs
2017RNALibPool01-7_S7_L001_R1_001.fastq.gz      2017RNALibPool01-7_S7_L001_R2_001.fastq.gz          luna_04_wholebody           1st instar
2017RNALibPool01-8_S8_L001_R1_001.fastq.gz      2017RNALibPool01-8_S8_L001_R2_001.fastq.gz          luna_05_wholebody           1st instar
2017RNALibPool01-10_S10_L001_R1_001.fastq.gz    2017RNALibPool01-10_S10_L001_R2_001.fastq.gz        luna_06_thorax		4th instar
2017RNALibPool01-11_S11_L001_R1_001.fastq.gz    2017RNALibPool01-11_S11_L001_R2_001.fastq.gz        luna_06_abdomen		4th instar
2017RNALibPool01-9_S9_L001_R1_001.fastq.gz      2017RNALibPool01-9_S9_L001_R2_001.fastq.gz          luna_06_head		4th instar
2017RNALibPool02-1_S12_L002_R1_001.fastq.gz     2017RNALibPool02-1_S12_L002_R2_001.fastq.gz         luna_07_head		Pre-pupa
2017RNALibPool02-2_S13_L002_R1_001.fastq.gz     2017RNALibPool02-2_S13_L002_R2_001.fastq.gz         luna_07_thorax		Pre-pupa
2017RNALibPool02-3_S14_L002_R1_001.fastq.gz     2017RNALibPool02-3_S14_L002_R2_001.fastq.gz         luna_07_abdomen		Pre-pupa
2017RNALibPool02-4_S15_L002_R1_001.fastq.gz     2017RNALibPool02-4_S15_L002_R2_001.fastq.gz         luna_08_head		Pre-pupa
2017RNALibPool02-5_S16_L002_R1_001.fastq.gz     2017RNALibPool02-5_S16_L002_R2_001.fastq.gz         luna_08_thorax		Pre-pupa
2017RNALibPool02-6_S17_L002_R1_001.fastq.gz     2017RNALibPool02-6_S17_L002_R2_001.fastq.gz         luna_08_abdomen		Pre-pupa
2017RNALibPool02-7_S18_L002_R1_001.fastq.gz     2017RNALibPool02-7_S18_L002_R2_001.fastq.gz         luna_09_head		Pre-pupa
2017RNALibPool02-8_S19_L002_R1_001.fastq.gz     2017RNALibPool02-8_S19_L002_R2_001.fastq.gz         luna_09_thorax		Pre-pupa
2017RNALibPool02-9_S20_L002_R1_001.fastq.gz     2017RNALibPool02-9_S20_L002_R2_001.fastq.gz         luna_09_abdomen		Pre-pupa
2017RNALibPool05-8_S53_L005_R1_001.fastq.gz     2017RNALibPool05-8_S53_L005_R2_001.fastq.gz         luna_13_thorax		Adult
2017RNALibPool07-10_S77_L007_R1_001.fastq.gz    2017RNALibPool07-10_S77_L007_R2_001.fastq.gz        luna_13_abdomen		Adult
    	
```

## (a) Trim RNAseq reads

To trim the reads, I use the following Trimmomatic script for all of the .fastq files in the list above.

```
#!/bin/bash
#SBATCH --job-name=trimmomatic_Al
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load trimmomatic/0.39

for sample in $(ls *fastq.gz | cut -d "_" -f 1,2,3 | sort | uniq)
do
    fq1=$(ls ${sample}_R1*)
    fq2=$(ls ${sample}_R2*)
    trimmomatic PE -threads 16 \
    ${fq1} ${fq2} \
    ${sample}_R1_clean.fastq.gz ${sample}_R1_unpaired.fastq.gz \
    ${sample}_R2_clean.fastq.gz ${sample}_R2_unpaired.fastq.gz LEADING:3 TRAILING:3 MINLEN:36
done
```

## (b) Combine tissues of same samples per life stage

Due to the nature of pooling our RNAseq data, there are reads for different tissues of the same individual. For this annotation, I use one sample from the following life stages: 4th instar (luna_01), eggs (luna_03), 1st instar (luna_05), prepupa (luna_07), pupa (luna_12), and adult (luna_14).

To combine samples of the same tissue, I use the following cat commands:
```
Pre-pupa (luna_07)
### combine R1 files 
```cat 2017RNALibPool02-1_S12_L002_R1_clean.fastq.gz 2017RNALibPool02-2_S13_L002_R1_clean.fastq.gz 2017RNALibPool02-3_S14_L002_R1_clean.fastq.gz > luna_07_prepupa_R1.fastq.gz```
### combine R2 files
```cat 2017RNALibPool02-1_S12_L002_R2_clean.fastq.gz 2017RNALibPool02-2_S13_L002_R2_clean.fastq.gz 2017RNALibPool02-3_S14_L002_R2_clean.fastq.gz > luna_07_prepupa_R2.fastq.gz```


Pupa (luna_12)
```cat 2017RNALibPool03-4_S26_L003_R1_clean.fastq.gz 2017RNALibPool03-5_S27_L003_R1_clean.fastq.gz 2017RNALibPool03-6_S28_L003_R1_clean.fastq.gz > luna_12_pupa_R1.fastq.gz```
```cat 2017RNALibPool03-4_S26_L003_R2_clean.fastq.gz 2017RNALibPool03-5_S27_L003_R2_clean.fastq.gz 2017RNALibPool03-6_S28_L003_R2_clean.fastq.gz > luna_12_pupa_R2.fastq.gz```


4th Instar (luna_01)
```cat 2017RNALibPool01-1_S1_L001_R1_clean.fastq.gz 2017RNALibPool01-2_S2_L001_R1_clean.fastq.gz 2017RNALibPool01-3_S3_L001_R1_clean.fastq.gz > luna_01_instar4_R1.fastq.gz```
```cat 2017RNALibPool01-1_S1_L001_R2_clean.fastq.gz 2017RNALibPool01-2_S2_L001_R2_clean.fastq.gz 2017RNALibPool01-3_S3_L001_R2_clean.fastq.gz > luna_01_instar4_R2.fastq.gz```


Adult (luna_14)
```cat 2017RNALibPool03-7_S29_L003_R1_clean.fastq.gz 2017RNALibPool03-8_S30_L003_R1_clean.fastq.gz 2017RNALibPool03-9_S31_L003_R1_clean.fastq.gz > luna_14_adult_R1.fastq.gz```
```cat 2017RNALibPool03-7_S29_L003_R2_clean.fastq.gz 2017RNALibPool03-8_S30_L003_R2_clean.fastq.gz 2017RNALibPool03-9_S31_L003_R2_clean.fastq.gz > luna_14_adult_R2.fastq.gz```

Eggs (luna_03)
```mv 2017RNALibPool01-6_S6_L001_R1_clean.fastq.gz luna_03_egg_R1.fastq.gz```
```mv 2017RNALibPool01-6_S6_L001_R2_clean.fastq.gz luna_03_egg_R2.fastq.gz```


1st Instar RNAseq (luna_05)
```mv 2017RNALibPool01-8_S8_L001_R1_clean.fastq.gz luna_05_instar1_R1.fastq.gz```
```mv 2017RNALibPool01-8_S8_L001_R2_clean.fastq.gz luna_05_instar1_R2.fastq.gz```
```

I created a braker folder for the RNA evidence run and moved these files there:

```
cd /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna
mv /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/trimmomatic/RAW_rnaseq/TRIMMED/luna_* /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna
```

## (c) Run hisat2 to align transcripts to the genome assembly

I use the following script, with the combined cat files above as my transcript input files:

```
#!/bin/bash
#SBATCH --job-name=Al_hisat2
#SBATCH --output=hisat2_%j.log
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

module load hisat2
module load samtools

hisat2-build /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/masked_genome.fasta /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna/Al_index

for sample in $(ls /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna/*fastq.gz | cut -d "_" -f 1,2,3,4,5,6,7 | sort | uniq)
do
    fq1=$(ls ${sample}_R1*)
    fq2=$(ls ${sample}_R2*)
    name=$(echo ${sample} | cut -d "_" -f 1,2,3)
    hisat2 -p 32 \
    -x /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna/Al_index \
    -1 ${fq1} -2 ${fq2} --phred33 --rna-strandness FR | samtools sort -@ 10 -O BAM -o Al_Al_aln.bam
done
```

The .log output file will contain mapping prercentages. For these files the mapping rates were:
- 4th instar (luna_01): 90.82%
- eggs (luna_03): 89.12%
- 1st instar (luna_05): 90.12%
- prepupa (luna_07): 90.31%
- pupa (luna_12): 88.07%
- adult (luna_14): 89.33%

## (d) Run BRAKER2 with aln.bam file

Lastly, I ran BRAKER with the output from hisat2 using the following command and script:

```sbatch -J Al_braker_RNAseq /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/scripts/braker2_RNAseq.sh /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/masked_genome.fasta aluna```

```bash
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1}
species=${2}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/config \
--genome=${genome} --species ${species} \
--bam=/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Al.Al.aln.bam \
--softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```


## Step 3 of BRAKER2 Feature Annotation – Running with IsoSeq (long-read) data

For the third run of BRAKER, I am using a new [modified pipeline](https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md#braker2) specifically for PacBio HiFi long reads in annotation. I use this modified pipeline up until the mapping step. For mapping and collapsing, I use [PacBio bulk IsoSeq workflow](https://isoseq.how/classification/workflow.html) using the pbmm and collapse funcitons to ensure script compatibility.


## (a) Merge all subreads into a single FASTQ file

Note: Files come back from ICBR in the form of .bam files per each barcoded instar. To convert them to FASTQ files, I used samtools bam2fq function:
```
samtools bam2fq SAMPLE.bam > SAMPLE.fastq
```

This can also be done using bedtools as follows:
```
srundev
bedtools bamtofastq -i instar1.bam -fq instar1.fq
bedtools bamtofastq -i instar2.bam -fq instar2.fq
bedtools bamtofastq -i instar3.bam -fq instar3.fq
bedtools bamtofastq -i instar4.bam -fq instar4.fq
bedtools bamtofastq -i instar5.bam -fq instar5.fq
```

Once .bam files have been converted to .fq, I used this code to concat each instar into all.subreads.fq file:
```
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_isoseq
cat instar1.fq instar2.fq instar3.fq instar4.fq instar5.fq > all.subreads.fq
```

Next, I ran pbmm (PacBio minimap) to map transcripts back to my A.luna masked genome using the following script:
```
#!/bin/bash
#SBATCH --job-name=%x_pbmm2_%j
#SBATCH -o %x_minimap_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 24:00:00
#SBATCH -c 32

module load pbmm2
module load isoseq3

pbmm2 align --preset ISOSEQ --sort /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_isoseq/all.subreads.bam \
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/masked_genome.fasta \
al_isoseq_mapped.bam
```

After mapping, I use the PacBio IsoSeq collapse funciton which serves the same purpose as Cupcake, and collapses redundant IsoForms:
```
#!/bin/bash
#SBATCH --job-name=%x_collapse_%j
#SBATCH -o %x_collapse_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 24:00:00
#SBATCH -c 32

module load isoseq3

isoseq3 collapse al_isoseq_mapped.bam al_isoseq_collapse.gff
```

After collapsing, I return to the modified BRAKER2 pipeline to conduct GeneMarkS-T predictions. I first run the Augustus stringtie2fa.py script within the BRAKER2 pipeline using the masked genome, and output collapsed .gff file from the previous isoseq collapse step
```
stringtie2fa.py -g genome.fa -f cupcake.collapsed.gff -o cupcake.fa
```

stringtie2fa.py script:
```
#!/usr/bin/env python3

# Author: Katharina J. Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Last modified on November 8th 2021
#
# This Python script extracts exon features from a GTF file, excises
# corresponding sequence windows from a genome FASTA file, stitches the
# codingseq parts together, makes reverse complement
# if required
# Output file is:
#    * file with mRNA squences in FASTA format
# Beware: the script assumes that the gtf input file is sorted by coordinates!
# This script is also compatible with cupcake gtf format

try:
    import argparse
except ImportError:
    raise ImportError(
        'Failed to import argparse. Try installing with \"pip3 install argparse\"')

try:
    import re
except ImportError:
    raise ImportError(
        'Failed to import argparse. Try installing with \"pip3 install re\"')

try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError(
        'Failed to import biophython modules. Try installing with \"pip3 install biopython\"')


parser = argparse.ArgumentParser(
    description='Generate *.codingseq and *.aa FASTA-format files from genes \
                 in a GTF-file produced by AUGUSTUS auxprogs tool joingenes \
                 and a corresponding genomic FASTA-file.')
parser.add_argument('-g', '--genome', required=True,
                    type=str, help='genome sequence file (FASTA-format)')
parser.add_argument('-o', '--out', required=True, type=str,
                    help="name stem pf output file with coding sequences and \
                    protein sequences (FASTA-format); will be extended by \
                    .codingseq/.aa")
parser.add_argument('-p', '--print_format_examples', required=False, action='store_true',
                    help="Print gtf input format examples, do not perform analysis")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-f', '--gtf',
                    type=str, help='file with CDS coordinates (GTF-format)')
args = parser.parse_args()


if args.print_format_examples:
    print('This script requires an annotation of a transcript with exon ' +
        'features in GTF format. We here provide a compatible ' +
        'example. ' +
        'This script will only process the exon lines. The input data may contain other feature lines that ' +
        'are ignored.')
    print('\nGTF format example:\n')
    print('1\tStringTie\ttranscript\t1555206\t1572018\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; cov "5.737374"; FPKM "5.261884"; TPM "18.775906";\n' +
         '1\tStringTie\texon\t1555206\t1555441\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "1"; cov "6.080509";\n' +
         '1\tStringTie\texon\t1565008\t1565346\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "2"; cov "5.917404";\n' +
         '1\tStringTie\texon\t1571901\t1572018\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "3"; cov "4.533898";\n')
    print('\nThis script has successfully been tested with GTF format produced by Stringtie2.')
    exit(0)

# output file names:
mrnaFile = args.out + ".mrna"

# Read GTF file exon entries for transcripts
tx2seq = {}
tx2str = {}
mrna = {}

if args.gtf:
    try:
        with open(args.gtf, "r") as gtf_handle:
            for line in gtf_handle:
                if re.match(
                        r"\S+\t[^\t]+\texon\t\d+\t\d+\t\S+\t\S\t\.\t.*transcript_id \"(\S+)\"", line):
                    #print("I attempt to store exon info")
                    seq_id, st, en, stx, tx_id = re.match(
                        r"(\S+)\t[^\t]+\texon\t(\d+)\t(\d+)\t\S+\t(\S)\t\.\t.*transcript_id \"(\S+)\"", line).groups()
                    if seq_id not in mrna:
                        mrna[seq_id] = {}
                    if tx_id not in mrna[seq_id]:
                        mrna[seq_id][tx_id] = []
                    mrna[seq_id][tx_id].append(
                        {'start': int(st), 'end': int(en), 'strand': stx})
                    if not tx_id in tx2seq:
                        tx2seq[tx_id] = seq_id
                        tx2str[tx_id] = stx
    except IOError:
        print("Error: Failed to open file " + args.gtf + "!")
        exit(1)
else:
    print("Error: No annotation file in GTF format was provided!")
    exit(1)

# Read genome file (single FASTA entries are held in memory, only), extract
# exon sequence windows
seq_len = {}
mrnaseq = {}
try:
    with open(args.genome, "r") as genome_handle:
        for record in SeqIO.parse(genome_handle, "fasta"):
            seq_len[record.id] = len(record.seq)
            if record.id in mrna:
                for tx in mrna[record.id]:
                    #print("I do something for tx")
                    if tx not in mrnaseq:
                        if mrna[record.id][tx][0]['strand'] == '.':
                            descr = tx + ' strand_unknown'
                        else:
                            descr = tx
                        mrnaseq[tx] = SeqRecord(Seq(""), id=tx, description=descr)
                    nExons = len(mrna[record.id][tx])
                    for i in range(0, nExons):
                        mrna_line = mrna[record.id][tx][i]
                        mrnaseq[tx].seq += record.seq[mrna_line['start'] - 1:mrna_line['end']]
                        if i == (nExons - 1) and mrna_line['strand'] == '-':
                            mrnaseq[tx].seq = mrnaseq[tx].seq.reverse_complement()
except IOError:
    print("Error: Failed to open file " + args.genome + "!")
    exit(1)

# Print mRNA sequences to file
try:
    with open(mrnaFile, "w") as mrna_handle:
        for tx_id, seq_rec in mrnaseq.items():
            SeqIO.write(seq_rec, mrna_handle, "fasta")
except IOError:
    print("Error: Failed to open file " + mrnaFile + "!")
    exit(1)
```

After using stringtie2fa.py, I ran GeneMarkST to use the IsoSeq data to train the gene model:
```
module load genemark_s/t-3.10.001
gmst.pl --strand direct cupcake.fa.mrna --output gmst.out --format GFF
```

Lastly, I use the GeneMarkS-T coordinates and the long-read transcripts to create a gene set in GTF format. Note, the gmst2globalCoords.py script is only in the long_read BRAKER documentation, so you must export this script following the [installation instructions](https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md). The output file should be gmst.global.gtf:
```
gmst2globalCoords.py -t al_isoseq_collapse.gff -p gmst.out -o gmst.global.gtf -g /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/masked_genome.fasta
```

I want to see how the gene model will look by combining all three outputs (protein evidence, RNA-seq evidence, and IsoSeq long read evidence), as well as the IsoSeq long read evidence on it's own. To assess combined model quality, I use TSEBRA to combine the gene models and then run BUSCO on the resulting .gtf file (converted to .aa) to assess quality. Please see the "Evaluate gene models" section for these results.


## Evaluate gene models produced by Braker2 

To evaluate the success of our gene models (of all 3 BRAKER2 runs), I used the BUSCO Lepidoptera ortholog database (odb10_lepidoptera).

## (a) From Bombyx mori protein evidence 

```
sbatch Al_Bm_prot_model_busco.sh
```

```
#!/bin/bash
#SBATCH --job-name=Al_lep_all_genemodel_busco
#SBATCH -o Al_lep_all_genemodel_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="/blue/kawahara/amanda.markee/Aluna_genome/aluna_assembly/hifiasm/BUSCO/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/config"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_prot/augustus.hints.aa \
 -o ./Al_Bm_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12
```
```
Results:
	- C:95.1%[S:88.3%,D:6.8%],F:1.2%,M:3.7%,n:5286	   
	- 5028	Complete BUSCOs (C)			   
	- 4669	Complete and single-copy BUSCOs (S)	   
	- 359	Complete and duplicated BUSCOs (D)	   
	- 66	Fragmented BUSCOs (F)			   
	- 192	Missing BUSCOs (M)			   
	- 5286	Total BUSCO groups searched
```

To determine how many genes were predicted using B.mori protein evidence, I used this command:
```
grep ">" augustus.hints.aa | wc -l
20801 genes
```

## (b) From Actias luna RNAseq evidence

```
sbatch Al_RNAseq_model_busco.sh
```

```
#!/bin/bash
#SBATCH --job-name=Al_lep_RNA_genemodel_busco
#SBATCH -o Al_lep_RNA_genemodel_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="/blue/kawahara/amanda.markee/Aluna_genome/aluna_assembly/hifiasm/BUSCO/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/config"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna/braker_rna_out/augustus.hints.aa \
 -o ./Al_RNA_genemod_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12 
```
```
Results:
	- C:97.2%[S:84.0%,D:13.2%],F:1.0%,M:1.8%,n:5286	   
	- 5140	Complete BUSCOs (C)			   
	- 4440	Complete and single-copy BUSCOs (S)	   
	- 700	Complete and duplicated BUSCOs (D)	   
	- 51	Fragmented BUSCOs (F)			   
	- 95	Missing BUSCOs (M)			   
	- 5286	Total BUSCO groups searched
```

To determine how many genes were predicted using A.luna RNAseq evidence, I used this command:
```
grep ">" augustus.hints.aa | wc -l
21773 genes
```

## (c) From Actias luna long-read IsoSeq evidence

```
sbatch Al_IsoSeq_model_busco.sh
```

```
#!/bin/bash
#SBATCH --job-name=Al_lep_IsoSeq_genemodel_busco
#SBATCH -o l_lep_IsoSeq_genemodel_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="/blue/kawahara/amanda.markee/Aluna_genome/aluna_assembly/hifiasm/BUSCO/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/config"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_isoseq/augustus.hints.aa \
-o ./Al_IsoSeq_genemod_busco_out \
-l /data/reference/busco/v5/lineages/endopterygota_odb10 \
-m protein -c 12
```



## Genome Annotation: Combine BRAKER2 outputs in TSEBRA, and assess combined gene models

## (a) Combine gene models from protein and transcriptome evidence

[TSEBRA Github](https://github.com/Gaius-Augustus/TSEBRA) is a transcript selector that allows you to combine the Braker outputs from different evidence into a single list of transcripts / predicted amino acid sequences.

First I cloned the TSEBRA directory into my working annotation directory on the cluster:
```
cd /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/tsebra
git clone https://github.com/Gaius-Augustus/TSEBRA
```

Next, I gather the necessary files as we will use them in the TSEBRA script below:

- augustus.hints.gtf and hintsfile.gff from BRAKER2 protein evidence
- augustus.hints.gff and hintsfile.gff from BRAKER2 RNA evidence
- gmst.global.gtf file from BRAKER2 long-read IsoSeq evidence
- A configuration (.cfg) file as descirbed [here](https://github.com/Gaius-Augustus/TSEBRA#configuration-file) 
_Note: You can use the default file located in /TSEBRA/config/default.ctg as a start, but the default is quite strict_

Code for combining all evidence (augustus.hints.gtf and hintsfile.gff, then appending the long_read protocol gmst.global.gtf). Output file should be tsebra_longread.gtf:
```
cd /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/tsebra/TSEBRA
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/tsebra/TSEBRA/bin/tsebra.py -g /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna/braker_rna_out/augustus.hints.gtf,/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_prot/augustus.hints.gtf -e /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna/braker_rna_out/hintsfile.gff,/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_prot/hintsfile.gff -l /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_isoseq/gmst.global.gtf -c /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/tsebra/TSEBRA/config/long_reads.cfg -o tsebra_longread.gtf
```

Once I have the combined TSEBRA model containing the long read evidence (.gtf), I have to convert this to .aa format so it can be read by Augustus (BUSCO):
```
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/scripts/gtf2aa.pl \
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/masked_genome.fasta \
/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/tsebra/TSEBRA/tsebra_longread.gtf \
tsebra_longread_aa.fa
```

Check how many genes: 
```
grep ">" tsebra_longread_aa.fa | wc -l
20144
```


## (b) Run BUSCO on the combined gene model

```
sbatch all_model_busco.sh
```

```
#!/bin/bash
#SBATCH --job-name=all_model_busco
#SBATCH -o all_model_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="/blue/kawahara/amanda.markee/Aluna_genome/aluna_assembly/hifiasm/BUSCO/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/config"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/tsebra/TSEBRA/tsebra_longread_aa.fa \
-o ./all_model_busco_out \
-l /data/reference/busco/v5/lineages/endopterygota_odb10 \
-m protein -c 12
```

The results for running BUSCO 5.3.0 are here:
```
***** Results: *****

	C:95.1%[S:75.4%,D:19.7%],F:0.9%,M:4.0%,n:2124	   
	2020	Complete BUSCOs (C)			   
	1602	Complete and single-copy BUSCOs (S)	   
	418	Complete and duplicated BUSCOs (D)	   
	19	Fragmented BUSCOs (F)			   
	85	Missing BUSCOs (M)			   
	2124	Total BUSCO groups searched		   

```


