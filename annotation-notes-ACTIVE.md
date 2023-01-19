## _STOP: Please see the [assembly](https://github.com/amandamarkee/actias-luna-genome/blob/main/assembly-notes-ACTIVE.md) note for the first half of the instructions before proceeding_ ##

## **10/9/2022; Annotation – RepeatModeler2**
Now that my genome assembly is complete, we have a high-quality and accurate assembly to begin building our annotation from. The ultimate goal is to perform both feature and functional annotation, but we will begin with feature annotation. For my annotation, I will be following [Dr. YiMing Weng's](https://github.com/yimingweng/Kely_genome_project/blob/main/note.md#09072022) and [Dr. Daren Card's](https://darencard.net/blog/2022-07-09-genome-repeat-annotation/) workflow, which is a modified version of the [MAKER](https://darencard.net/blog/2017-05-16-maker-genome-annotation/) annotation protocol. 

First, we need to run [RepeatModeler2](https://www.pnas.org/doi/10.1073/pnas.1921046117) to identify and model for repeat elements in the genome.

Navigate to your working directory for annotation, and execute the following script on your genome assembly using the following code:
```
[amanda.markee@login1 aluna_annotation]$ pwd
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation
sbatch repeatmodeler.sh /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/aluna_final_assembly.fasta aluna_repeatmodeler
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
cd /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_masker
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
/blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/aluna_final_assembly.fasta &> aluna_repeatmasker_step1.out
```

## Step 2: Feature Annotation – Mask Repeats Based on Existing Databases

Once we mask the simple repeat elements, we will use the existing database [Repbase](https://www.girinst.org/repbase/) to continue building our gene model. First we navigate too our RepeatMasker directory, then run the following script in that directory.

```
cd /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_masker
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
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_masker/aluna_repeatmasker_step1/aluna_final_assembly.fasta.masked &> aluna_repeatmasker_step2.out

mv aluna_repeatmask_step2.out out_files
```

## Step 3: Feature Annotation – Mask Genome Based on RepeatModeler2 Output

Next, we use the output from RepeatModeler2 (from the beginning of the annotation workflow) to continue masking repeat regions. 

```
cd /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_masker
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
-lib /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_modeler/aluna-genome-families.prefix.fa.known \
-dir aluna_repeatmasker_step3 \
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_masker/aluna_repeatmasker_step2/aluna_final_assembly.fasta.masked.masked &> aluna_repeatmasker_step3.out

mv aluna_repeatmask_step3.out aluna_repeatmasker_step3
```

## Step 4: Feature Annotation – Calculate Masking Percentage
The final step is to calculate the masking percentage to determine the soft masking rate vs. the hard masking rate. I used the following script.

```
cd /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_masker/aluna_repeatmasker_step3
bash /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/repeat_masker/aluna_repeatmasker_step3/maskrate.sh aluna_final_assembly.fasta.masked.masked.masked
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
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/funannotate
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
13. 

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
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/funannotate/predict
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
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2
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
sbatch -J Al_ProtHint prothint.sh /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/masked_genome.fasta /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/B_mori_protein.fasta
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
sbatch -J Al_prot_braker2 Al_braker2_protein.sh /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/masked_genome.fasta /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/prothint_augustus.gff Actias_luna
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
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/trimmomatic/RAW_rnaseq 

### trimmed RNAseq data files
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/trimmomatic/trimmed 

### parent folder for scripts and outputs related to braker
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2 
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
cd /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_rna
mv /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/trimmomatic/RAW_rnaseq/TRIMMED/luna_* /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_rna
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

hisat2-build /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/masked_genome.fasta /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_rna/Al_index

for sample in $(ls /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_rna/*fastq.gz | cut -d "_" -f 1,2,3,4,5,6,7 | sort | uniq)
do
    fq1=$(ls ${sample}_R1*)
    fq2=$(ls ${sample}_R2*)
    name=$(echo ${sample} | cut -d "_" -f 1,2,3)
    hisat2 -p 32 \
    -x /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_rna/Al_index \
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

```sbatch -J Al_braker_RNAseq /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/scripts/braker2_RNAseq.sh /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/masked_genome.fasta aluna```

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
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/Augustus/config \
--genome=${genome} --species ${species} \
--bam=/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/Al.Al.aln.bam \
--softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```


## Step 3 of BRAKER2 Feature Annotation – Running with IsoSeq (long-read) data

For the third run of BRAKER, I am using a new [modified pipeline](https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md#braker2) specifically for PacBio HiFi long reads.

From this protocol, I will be picking up starting at the BRAKER2 step, since I have already run BRAKER on protein evidence following Dr. Godfrey's protocol in Step 1. 

## (a) Convert IsoSeq raw .bam reads to .fastq

All IsoSeq BRAKER steps will be run in the following directory:
```
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_isoseq
```

First I requested a development node, and ran the following commands using the bedtools module to convert all .bam files to .fastq:

```
srundev
bedtools bamtofastq -i instar1.bam -fq instar1.fq
bedtools bamtofastq -i instar2.bam -fq instar2.fq
bedtools bamtofastq -i instar3.bam -fq instar3.fq
bedtools bamtofastq -i instar4.bam -fq instar4.fq
bedtools bamtofastq -i instar5.bam -fq instar5.fq

mkdir isoseq_raw_bam 
mv *.bam isoseq_raw_bam/
```

Next, I concatenate all subread libraries (instar1-instar5) files into a single file:
```
cat instar1.fq instar2.fq instar3.fq instar4.fq instar5.fq > all.subreads.fastq
```

## (b) Use minimap2 to map transcripts to genome

Use Minimap2 to map the transcripts to the genome sequence, and then sort the result.The final output will be called Al_isoseq_final.bam:
```
cd /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_isoseq
sbatch minimap2.sh /blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2 all.subreads.fastq Al_isoseq_final
```

```
#!/bin/bash
#SBATCH --job-name=%x_minimap_%j
#SBATCH -o %x_minimap_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 24:00:00
#SBATCH -c 32

module load minimap/2.21
module load samtools/1.15

genome=${1}
read=${2}
prefix=${3}

minimap2 -t 8 -ax map-pb ${genome} ${read} --secondary=no | samtools sort -m 1G -o ${prefix}.bam -T tmp.ali
```

## (c) Collapse redundant isoforms in Cupcake


Next, we collapse redundant isoforms using a script from Cupcake. Make sure to module load cupcake to have access to the collapse isoforms script. Note: [the GitHub](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step#collapse-redundant-isoforms-has-genome) mentions implementing this as isoseq3? Look into this later.

```
module load cupcake
collapse_isoforms_by_sam.py --input all.subreads.fastq --fq  -s Al_isoseq_final.sam --dun-merge-5-shorter -o cupcake

```

Note: This code gives me the following error, which I will pick up on tomorrow.
```
File "/apps/cupcake/22.0.0/bin/collapse_isoforms_by_sam.py", line 245, in <module>
    main(args)
  File "/apps/cupcake/22.0.0/bin/collapse_isoforms_by_sam.py", line 169, in main
    check_ids_unique(args.input, is_fq=args.fq)
  File "/apps/cupcake/22.0.0/lib/python3.9/site-packages/cupcake/tofu/utils.py", line 13, in check_ids_unique
    raise Exception("Duplicate id {0} detected. Abort!".format(r.id))
Exception: Duplicate id m64219e_220708_202551/155/ccs detected. Abort!
```

## (d) Final organization

I organized my output files for each braker run in the following directories:
```
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_prot
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_rna
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/braker_isoseq
```

I also moved all of my scripts to this folder:
```
/blue/kawahara/amanda.markee/insect_genomics_2022/aluna_annotation/braker2/scripts
```
