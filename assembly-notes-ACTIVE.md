
## **09/19/2022**; General Information & Setup

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
<br />

## **09/19/2022; Data Accumulation and File Explanations**

- Samples were sequenced at the University of Florida's Interdisciplinary Center for Biotechnology Research (ICBR)
- Once ICBR has completed HiFi sequencing, the core emailed me regarding my run reports (see run-reports directory)
- To initiate data delivery, I coordinated with ICBR's bioinformatics team following [these instructions](https://biotech.ufl.edu/wp-content/uploads/2022/04/Configuring-Globus-with-ICBR-as-a-HiPerGator-Account-Holder.pdf)
- All raw data for the genomic sequence data (NS2497) can be found in the following HiPerGator directory

```
ssh amanda.markee@hpg.rc.ufl.edu
## type password here
cd /blue/kawahara/amanda.markee/ICBR/NS-2497
```
- ICBR uses the Globus file transfer system to upload the following directories of data related to your sequences:
  - 0000000200 = folder that contains full completed hifi reads w/metadata
  - 0000000205 = folder that contains css / .json  
  - 1_A05 = folder that contains raw data, and ccs reports
<br />

## **09/19/2022; Raw Read Quality Assessment with FastQC**
-[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality assessment tool for next generation sequencing, often used to assess raw reed quality and highlight problem areas in the form of visualizations (see results below).

First, to run FastQC, create a directory for the output data and script.
```
mkdir fastqc_raw
```

I then copied the following script into the new directory that contained my raw hifi reads:
```
#!/bin/bash
#SBATCH --job-name=A_luna_fastqc
#SBATCH --output=A_luna_raw_reads_fastqc-%j.out
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=02:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

## This script runs a fast qc on the fastq file output from PacBio HiFi Sequel IIe

## load fastqc module
module load fastqc

## run fastqc on the data 
fastqc m64219e_220210_175238.hifi_reads.fastq.gz
```

The FastQC output files include an .html file, which will contain visualizations for the following results:
- Basic Statistics
- Per base sequence quality
- Per sequence quality scores (Phred scores)
- Per base sequence content
- Per sequence GC content
- Per base N content
- Sequence Length Distribution
- Sequence Duplication Levels
- Overrepresented sequences
- Adapter content


![Screen Shot 2022-10-03 at 12 19 27 PM](https://user-images.githubusercontent.com/56971761/193628513-bdf8ac54-fd40-4d6d-bd21-3fc35f2a1766.png)

![Screen Shot 2022-10-03 at 12 19 37 PM](https://user-images.githubusercontent.com/56971761/193628527-02fc686c-ed19-4043-b758-22dbdc1e1fb0.png)
<br />

## **09/19/2022; Genome Assembly with hifiasm**

 - [hifiasm](https://hifiasm.readthedocs.io/en/latest/) is a fast and easy haplotype-resolved de novo assembly software for PacBio HiFi reads
 - hifiasm documentation explaining input parameters: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html
 - hifiasm documentation explaining output files: https://hifiasm.readthedocs.io/en/latest/interpreting-output.html
 - Script for hifiasm with aggressive duplicate purging (option -l 2)

```
#!/bin/bash
#SBATCH --job-name=markee_hifi_assembly
#SBATCH -o %A_%a.2022-06-14.hifiasm_assembly.out
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 32
#SBATCH --mem-per-cpu=10gb
#SBATCH -t 30:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara-b

date;hostname;pwd

module load ufrc
module load hifiasm

hifiasm -o /blue/kawahara/amanda.markee/ICBR/hifiasm/aclu_hifi_assembly_06-14-2022.asm 
-l 2 -t 32 /blue/kawahara/amanda.markee/ICBR/NS-2497/0000000200/outputs/m64219e_220210_175238.hifi_reads.fastq.gz

# output files moved to following directory
/blue/kawahara/amanda.markee/ICBR/hifiasm/hifiasm_output_files
```
<br />

## **10/03/2022; Genome Assembly Quality Assessment with assemblystats.py**

- After assembly with hifiasm, we can assess assembly quality using the [assemblystats.py script](https://github.com/MikeTrizna/assembly_stats/tree/0.1.4) created by Mike Trizna.
- The version of assemblystats.py used here was modified by Paul Frandsen (Brigham Young University).

First, I copied this script into my working directory, and called it assemblystats.py

```
#!/usr/bin/env python

import numpy as np
from itertools import groupby
import json
import sys


def fasta_iter(fasta_file):
    """Takes a FASTA file, and produces a generator of Header and Sequences.
    This is a memory-efficient way of analyzing a FASTA files -- without
    reading the entire file into memory.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    header: str
        The string contained in the header portion of the sequence record
        (everything after the '>')
    seq: str
        The sequence portion of the sequence record
    """

    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq


def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])

        stats['L' + str(level)] = l_level
        stats['N' + str(level)] = n_level
    return stats


if __name__ == "__main__":
    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    print(json.dumps(stat_output, indent=2, sort_keys=True))
```

Next, I changed permissions as follows to allow execution permissions.
```
chmod +x assemblystats.py
```

Then, I produced a FASTA file from the initial GFA output file from the hifiasm assembly output. I used the primary contig file, as indicated with asm.bp.p_ctg.fa (p_ctg meaning primary contig)
```
awk '/^S/{print ">"$2;print $3}' aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.gfa > aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa 
```

Lastly, I ran the assemblystats.py script on the newly generated fasta file of the fga in the format of scriptfilepath/scirptname.py nameofassembly.fa
```
/blue/kawahara/amanda.markee/ICBR/hifiasm assemblystats.py /blue/kawahara/amanda.markee/ICBR/hifiasm/hifiasm_output_files/aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa
```
Save results as a text file as shown.
```
./assemblystats.py aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa >> aclu_hifi_assembly_06-14-2022.txt
```
The results will look like the following table:

![Screen Shot 2022-10-03 at 11 54 34 AM](https://user-images.githubusercontent.com/56971761/193622529-568bc8f8-f936-4995-91c8-989f8da21759.png)
<br />

## **10/03/2022; Genome Completeness with BUSCO**

- As another quality assessment step, I used [BUSCO](https://busco.ezlab.org/busco_userguide.html) (Benchmarking Universal Single-Copy Orthologs) to determine genome completeness. BUSCO looks for the presence or absence of highly conserved genes in an assembly.
- Per Pacific Bioscience's [documentation](https://busco.ezlab.org/busco_userguide.html) of quality assessment post-sequencing, generally a BUSCO score of over 95% is goood, indicating a high percentage of genes identified in the assembly.


First, make a BUSCO directory in your working directory, and copy the [UF HiPerGator config.ini file](https://help.rc.ufl.edu/doc/BUSCO) into this directory.
```
mkdir BUSCO
cp $HPC_BUSCO_CONF/config.ini .
```

Then, edit the SLURM submission script to ensure that the paths to your assembly file and BUSCO config file are correct.
```
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH -o Luna_busco5_endo_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 4:00:00
#SBATCH -c 6

export BUSCO_CONFIG_FILE=/blue/kawahara/amanda.markee/ICBR/BUSCO/config.ini
export AUGUSTUS_CONFIG_PATH=/blue/kawahara/amanda.markee/ICBR/BUSCO/

echo $BUSCO_CONFIG_FILE

module load busco/5.2.0

busco -f -i /blue/kawahara/amanda.markee/ICBR/hifiasm/aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa \
 -o BUSCO_Luna_endopterygotap -l /data/reference/busco/v5/lineages/endopterygota_odb10   \
 -m genome -c 6

#	endopterygota_odb10
#	insecta_odb10
```

Once BUSCO has finished running, a short summary output table will show you your...
- C: completeness [S:single-copy, D:duplicated]
- F: fragmented
- M: missing
- n: number of genes used

```
***** A.luna BUSCO Results: *****

	C:99.4%[S:99.0%,D:0.4%],F:0.2%,M:0.4%,n:2124	   
	2111	Complete BUSCOs (C)			   
	2102	Complete and single-copy BUSCOs (S)	   
	9	Complete and duplicated BUSCOs (D)	   
	5	Fragmented BUSCOs (F)			   
	8	Missing BUSCOs (M)			   
	2124	Total BUSCO groups searched
```
<br />

## **10/03/2022; Genome Size Estimation (kmer) with KMC and GenomeScope**

- To estimate the size, heterozygosity and repetitiveness of the genome, I used a kmer distribution-based approach following the methods described in [Kawahara et al. 2022 publcation](https://gigabytejournal.com/articles/64).

- Following the short and simple instructions on the [GenomeScope documentation page](http://qb.cshl.edu/genomescope/genomescope2.0/), and the script below, I generated kmer count estimations and a histogram of kmer distribution using GenomeScope (v.2.0)
- Note, make sure to use GenomeScope version 2.0, shown in the documentation link above. KMC will count and estimate genome size based on kmer distribution, while GenomeScope visializes the results.

For running KMC in a development node, paste the following code directly into command line after requesting resources (ensuring that the command formatis input path, output file name, output path:
```
kmc -k21 -t10 -m64 -ci1 -cs10000 /blue/kawahara/amanda.markee/ICBR/NS-2497/0000000200/outputs/m64219e_220210_175238.hifi_reads.fastq.gz kmer_count_aluna_genome.kmc /blue/kawahara/amanda.markee/insect_genomics_2022/kmc_temp_dir_dev/kmer_count_aluna_genome.tmp
```

For running KMC using a SLURM submission scirpt, use the following script:
```
#!/bin/bash
#SBATCH --job-name=Aluna_kmc
#SBATCH -o Aluna_kmc_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH -c 3
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 00:30:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara-b

module load kmc/3.2.1

# create directory for kmc temporary files
mkdir kmc_tmp
 
kmc -k21 /blue/kawahara/amanda.markee/ICBR/NS-2497/0000000200/outputs/m64219e_220210_175238.hifi_reads.fastq.gz 21mers kmc_tmp

# Having the k-mers counted it is possible to dump KMC binary database to textual form with kmc_tools.

kmc_tools transform 21mers dump 21mers.txt

kmc_tools transform 21mers histogram 21mer_reads.histo
```

Lastly, once your KMC script is completed (either via submission or development node), the output file ending with .histo can be uploaded directly into the [GenomeScope GUI](http://qb.cshl.edu/genomescope/genomescope2.0/) to visualize results.

A link to the results I am going to display below can be found [here](http://genomescope.org/genomescope2.0/analysis.php?code=yZmFQCw5YGmWjOazjOE7)


![Screen Shot 2022-10-03 at 11 42 28 AM](https://user-images.githubusercontent.com/56971761/193621855-c7eecdc6-1c9a-4876-aebc-74101c825927.png)

![Screen Shot 2022-10-03 at 11 43 32 AM](https://user-images.githubusercontent.com/56971761/193621902-49f8ec6d-4499-48a8-ad3b-e666f2568c12.png)

![Screen Shot 2022-10-03 at 11 43 49 AM](https://user-images.githubusercontent.com/56971761/193622067-762f4454-57a1-4f71-a6a5-4c56707e117e.png)

![Screen Shot 2022-10-03 at 11 43 56 AM](https://user-images.githubusercontent.com/56971761/193622077-77435740-5c58-460b-9739-1136b698468d.png)

_GenomeScope Results Summarized:_
- The estimated genome size is about 420 Mbp.
- The heterozygosity is fairly low, about 1.3%
- The mean coverage is about 19.6X, which is in line with what is expected from the model (1.959e+01).
<br />

## **10/04/2022; Contaminaiton Filtering with BlobTools**

- We can use Blob Tools to assess genome contamination (version 1.0)

![Screen Shot 2022-10-04 at 9 40 26 AM](https://user-images.githubusercontent.com/56971761/193835309-9ddf6dfc-7cd4-4f33-a94d-ad6ab778ca45.png)

- To run blob tools we need three files:
	- (1) nodes.dmp and names.dmp files (from NCBI tax dump) 
	- (2) .bam alignment file (from minimap2 & samtools script) 
	- (3) .nt blast match file (from NCBI megablast)

(1) To download NCBI taxdump and create nodes.dmp and names.dmp (copy and paste this command in terminal, it will create data directory and create files)

```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
. blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
```

(2) To generate mapping files (BAM), I used the following script which utilizes minimap and samtools:
```
#!/bin/sh
#SBATCH --job-name=Al_minimap2
#SBATCH --output=Al_minimap2_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60gb
#SBATCH --time=08:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

pwd; hostname; date

module load minimap2

minimap2 -ax map-pb /blue/kawahara/amanda.markee/ICBR/hifiasm/hifiasm_output_files/aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa /blue/kawahara/amanda.markee/ICBR/NS-2497/0000000200/outputs/m64219e_220210_175238.hifi_reads.fastq.gz > Al.aln.sam

module load samtools
#convert SAM file to BAM file
samtools view -S -b Al.aln.sam > Al.aln.bam

#Use samtools sort to convert the BAM file to a coordinate sorted BAM file
samtools sort Al.aln.bam > Al.aln.sorted.bam

#index a genome sorted bAM file for quick alignment
samtools index Al.aln.sorted.bam > Al_indexed_sorted_bam
```
Note: The output for the megablast script will contain multiple files. The one used for input for BlobTools is the XX.aln.sorted.bam file


(3) Lastly, to generate taxanomic hit files, I used the following script which utilizes megablast:
```
#!/bin/sh
#SBATCH --job-name=Al_megablast
#SBATCH --output=Al_megablast_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amarkee@floridamuseum.ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20gb
#SBATCH --time 8-00:00:00
#SBATCH --qos=kawahara
#SBATCH --account=kawahara

pwd; hostname; date

module load ncbi_blast
blastn -db nt -task megablast -query /blue/kawahara/amanda.markee/ICBR/hifiasm/hifiasm_output_files/aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa -out Aluna_megablast.nt -evalue 1e-5 -outfmt "6 qseqid staxids bitscore sgi sskingdoms sscinames" -max_target_seqs 1 -num_threads=16
```
Note: The output for the megablast script will contain a file with a ".nt" and ".txt" extension. The ".nt" extension is used as the hit file input for BlobTools.

Once all three input files are created, I ran the following BlobTools script to generate the contamination visualization (BlobPlot):
```
#SBATCH --job-name=Al_blob
#SBATCH --output=Al_blob_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb
#SBATCH --time 04:00:00
#SBATCH --qos=kawahara
#SBATCH --account=kawahara

pwd; hostname; date

## Run blob tools
 
module load blobtools/1.0

blobtools create -i /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa \
-b /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/Al.aln.sorted.bam \
-t /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/Aluna_megablast.nt \
--nodes /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/nodes.dmp \
--names /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/names.dmp \
-o Al_blob_result

## You can then view and plot
#blobtools view -i result.blobDB.json
#blobtools plot -i result.blobDB.json
```

![Al_blob_result blobDB json bestsum phylum p7 span 100 blobplot bam0](https://user-images.githubusercontent.com/56971761/193838591-0404d830-aff1-4b7c-bbe2-23b6420c56ad.png)

![Al_blob_result blobDB json bestsum phylum p7 span 100 blobplot read_cov bam0](https://user-images.githubusercontent.com/56971761/193838619-2fa4dc30-3de9-4237-a1e5-0a369d08b213.png)

Interestingly, while most of the DNA came back in line with arthropoda (72.26%), there was a percentage that came back as microsporidia (2.05%) which could be indicative of a fungal infection in the organism I sequenced. Evidence of this is shown in [some literature for Saturniidae](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6842381/) and I will be looking further into this.

<br />

## **10/9/2022; Contaminaiton Removal and BUSCO Re-Run**

Blobtools revealed three non-target contigs in the assembly (ptg000026l, ptg000035l, and ptg000043l) belonging to Microsporidia, and eight non-target contigs in the assembly (ptg000047l, ptg000085l, ptg000090l, ptg000096l, ptg000097l, ptg000098l, ptg000106l, ptg000123l) beloging to Streptophyta. 

1) I conducted a quick check on the original BUSCO output result to see if any contribution was made from these contigs,
- If the contig showed up in the original output table, I removed them manually, and redid BUSCO.
- After checking the first contig, it did show up in the table, so I manually removed all contigs and re-ran BUSCO.

```
[amanda.markee@login1 run_endopterygota_odb10]$ pwd
/blue/kawahara/amanda.markee/ICBR/hifiasm/BUSCO/BUSCO_Luna_endopterygotap/run_endopterygota_odb10
[amanda.markee@login1 run_endopterygota_odb10]$ cat full_table.tsv | grep "ptg000026l"
# this contig came up multiple times in the BUSCO output
```

2) To manually make the final assembly, remove non-target contigs from the genome using a simple grep/sed command. I could probably edit this into a prettier for-loop, but for the sake of time, and due to there only being a few contaminant contigs, I wrote the script this way:
```
#!/bin/sh
#SBATCH --job-name=Al_final_asm
#SBATCH --output=Al_final_asm_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb
#SBATCH --time 04:00:00
#SBATCH --qos=kawahara
#SBATCH --account=kawahara

cd /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools

sed -e '/ptg000026l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa > temp1_aluna_final_assembly.fasta

sed -e '/ptg000035l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp1_aluna_final_assembly.fasta > temp2_aluna_final_assembly.fasta

sed -e '/ptg000043l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp2_aluna_final_assembly.fasta > temp3_aluna_final_assembly.fasta

sed -e '/ptg000047l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp3_aluna_final_assembly.fasta > temp4_aluna_final_assembly.fasta

sed -e '/ptg000085l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp4_aluna_final_assembly.fasta > temp5_aluna_final_assembly.fasta

sed -e '/ptg000090l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp5_aluna_final_assembly.fasta > temp6_aluna_final_assembly.fasta

sed -e '/ptg000096l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp6_aluna_final_assembly.fasta > temp7_aluna_final_assembly.fasta

sed -e '/ptg000097l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp7_aluna_final_assembly.fasta > temp8_aluna_final_assembly.fasta

sed -e '/ptg000098l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp8_aluna_final_assembly.fasta > temp9_aluna_final_assembly.fasta

sed -e '/ptg000106l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp9_aluna_final_assembly.fasta > temp10_aluna_final_assembly.fasta

sed -e '/ptg000123l/,+1d' /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/temp10_aluna_final_assembly.fasta > aluna_final_assembly.fasta
```

To confirm that these removals worked, we can take the original wc of the original assembly, and subtract 2n (where n = number of contigs removed, since two lines should have been removed per non-target contig):
- 310 original lines - 2(3 Microsporidia + 8 Streptophyta) = 288 lines remaining
```
[amanda.markee@c0709a-s30 blobtools]$ wc aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa 
      310       310 533454976 aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa # this is the original file, containing 310 lines
      
[amanda.markee@c0709a-s30 blobtools]$ wc aluna_final_assembly.fasta 
      288       288 518035841 aluna_final_assembly.fasta # this is the new and final assembly file, containing 288 lines, aka what we want
```

We can also confirm the removals worked by using grep to look for the original sequence names we used in the removal script. If they are no longer in the final assembly file, then we removed the seqs successfully:
```
[amanda.markee@c0709a-s30 blobtools]$ grep ptg000026l aluna_final_assembly.fasta # no results for this example from Microsporidia
[amanda.markee@c0709a-s30 blobtools]$ grep ptg000096l aluna_final_assembly.fasta # no results for this example from Streptophyta
```

3) I reran busco on the final verison of the assembly to ensure that these contaminated contigs didn't effect duplication.

```
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH -o Luna_busco5_endo_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 4:00:00
#SBATCH -c 6

export BUSCO_CONFIG_FILE=/blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/BUSCO_RERUN/config.ini
export AUGUSTUS_CONFIG_PATH=/blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/BUSCO_RERUN

echo $BUSCO_CONFIG_FILE

module load busco/5.2.0

busco -f -i /blue/kawahara/amanda.markee/insect_genomics_2022/blobtools/BUSCO_RERUN/aluna_final_assembly.fasta \
 -o BUSCO_Luna_endopterygotap -l /data/reference/busco/v5/lineages/endopterygota_odb10   \
 -m genome -c 6

#       endopterygota_odb10
#       insecta_odb10
```
<br />
