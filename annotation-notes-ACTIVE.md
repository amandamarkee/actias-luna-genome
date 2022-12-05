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

In this case, I do have RNAseq data, and will use the annotation workflow from [FunAnnotate](https://funannotate.readthedocs.io/en/latest/index.html#) to include long-read RNAseq data, and information from the UniProt/SwissProt protein database. 

Funannotate uses Evidence Modeler to combine ab initio gene model predictions with evidence (transcripts or proteins) aligned to the genome. Therefore, the evidence that you supply at runtime for --transcript_evidence and --protein_evidence are important. By default, funannotate will use the UniProtKb/SwissProt curated protein database for protein evidence. However, you can specify other forms of protein evidence, perhaps from a well-annotated closely related species, using the --protein_evidence option. Multiple files can be passed to both --transcript_evidence or --protein_evidence by separating the files by spaces, for example:

```
funannotate predict --input genome.fa --species "Awesome species" --transcript_evidence --pacbio_isoseq isoseq.fasta myESTs.fa \
    -o output --protein_evidence closely_related.fasta $FUNANNOTATE_DB/uniprot_sprot.fasta
```


