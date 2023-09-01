## Genome Annotation: BRAKER3 (YiMing Weng Pipeline) for Protien and RNAseq (8/29/23)
For the final version of the genome, Dr. YiMing Weng helped run BRAKER3 using a combined protein evidence from _B.mori_ and RNAseq data from 36 samples of varying tissue types (head, thorax and abdomen). I have uploaded his verbatim workflow to this repository, and explain in depth below. 

Note: We chose to use BRAKER3 instead of BRAKER2 with TSEBRA (of all evidence: _B.mori_ protien, RNAseq for each instar and long-read IsoSeq data of silk glands for each instar) because the long-read IsoSeq pipeline is not well developed. Isoforms were not properly collapsed with the pipeilne used above, and hidnered the model's duplication rate. Thus, we use BRAKER3 with just _B.mori_ protein, and RNAseq for all samples.

## (a) File Setup (PICK UP HERE)
```
# soft masked genome: /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/aluna_masked_genome.fa
# protien fasta: /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_prot/B_mori_protein.fasta
```
[yimingweng@login5 Actias_luna]$ pwd
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/ETP/Actias_luna

# get the softmasked genome assembly
cp /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/aluna_masked_genome.fa ./

# get the protein
cp /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_prot/B_mori_protein.fasta ./

# trim the sequence name to prevent the potential error 
cat B_mori_protein.fasta  | cut -d " " -f 1 >> B_mori_protein_trk.fasta
rm B_mori_protein.fasta

# get the RNA sequences
mkdir RNA
cd RNA
cp /orange/kawahara/amanda.markee/prelim_luna_dge/RAW_rnaseq/*001.fastq.gz ./

# rename RNA fastq files so that BRAKER can read them as pairs
ls | cut -d "_" -f 1 | sort -u >> list
IFS=$'\n'
for id in $(cat list)
do
    mv ${id}_*_R1_* ${id}_1.fastq.gz
    mv ${id}_*_R2_* ${id}_2.fastq.gz
done
rm list 
#### note: I used the mapped reads in bam file, provided by Amanda
cp /blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/braker_rna/hisat/Al_Al_aln.bam ./

# braker3 ETP pipeline
sbatch -J Alun_brk3 /blue/kawahara/yimingweng/universal_scripts/braker3_ETP_bam.slurm \
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/ETP/Actias_luna/aluna_masked_genome.fa \
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/ETP/Actias_luna/B_mori_protein_trk.fasta \
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/ETP/Actias_luna/Al_Al_aln.bam

# finish running after ~6 hrs, evaluate the completeness for braker.aa and augustus.hints.aa
[yimingweng@login5 braker]$ pwd
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/ETP/Actias_luna/braker
sbatch -J Alun_brk /blue/kawahara/yimingweng/universal_scripts/busco_gene_model.slurm braker.aa
        --------------------------------------------------
        |Results from dataset lepidoptera_odb10           |
        --------------------------------------------------
        |C:90.5%[S:75.1%,D:15.4%],F:0.6%,M:8.9%,n:5286    |
        |4784   Complete BUSCOs (C)                       |
        |3968   Complete and single-copy BUSCOs (S)       |
        |816    Complete and duplicated BUSCOs (D)        |
        |31     Fragmented BUSCOs (F)                     |
        |471    Missing BUSCOs (M)                        |
        |5286   Total BUSCO groups searched               |
        --------------------------------------------------


[yimingweng@login5 Augustus]$ pwd
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/ETP/Actias_luna/braker/Augustus
sbatch -J Alun_brk /blue/kawahara/yimingweng/universal_scripts/busco_gene_model.slurm augustus.hints.aa
cat Alun_gene_model_busco_3681843.log
        --------------------------------------------------
        |Results from dataset lepidoptera_odb10           |
        --------------------------------------------------
        |C:98.9%[S:86.5%,D:12.4%],F:0.4%,M:0.7%,n:5286    |
        |5230   Complete BUSCOs (C)                       |
        |4573   Complete and single-copy BUSCOs (S)       |
        |657    Complete and duplicated BUSCOs (D)        |
        |19     Fragmented BUSCOs (F)                     |
        |37     Missing BUSCOs (M)                        |
        |5286   Total BUSCO groups searched               |
        --------------------------------------------------

# get primary isoform model
[yimingweng@login5 Augustus]$ pwd
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/ETP/Actias_luna/braker/Augustus
echo "Actias_luna.aa" >> list
cat augustus.hints.aa >> Actias_luna.aa
sbatch -J Alun /blue/kawahara/yimingweng/universal_scripts/brakeraa4primary.slurm  list
python3 /blue/kawahara/yimingweng/universal_scripts/primary_transcript.py Alun.fa
cat ./primary_transcripts/Alun.fa >> Alun_primary.fa
rm -r primary_transcripts
sbatch -J Alun_primary /blue/kawahara/yimingweng/universal_scripts/busco_gene_model.slurm Alun_primary.fa

cat Alun_primary_gene_model_busco_3683604.log
        --------------------------------------------------
        |Results from dataset lepidoptera_odb10           |
        --------------------------------------------------
        |C:98.9%[S:98.1%,D:0.8%],F:0.4%,M:0.7%,n:5286     |
        |5229   Complete BUSCOs (C)                       |
        |5186   Complete and single-copy BUSCOs (S)       |
        |43     Complete and duplicated BUSCOs (D)        |
        |19     Fragmented BUSCOs (F)                     |
        |38     Missing BUSCOs (M)                        |
        |5286   Total BUSCO groups searched               |
        --------------------------------------------------

