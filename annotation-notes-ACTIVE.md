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
#SBATCH -t 08:00:00
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
