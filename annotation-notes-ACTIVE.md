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
