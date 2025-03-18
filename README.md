# Project: Chromosome Y extended

## Reference sequences - notes

### T2Tv2 / chm13 + HG002-Y

1. whole-genome assembly: [GCF_009914755.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)
    - source: [Rhie et al., "T2T Y paper" DOI:10.1038/s41586-023-06457-y](https://pmc.ncbi.nlm.nih.gov/articles/PMC10752217/)
2. HG002-Y part of T2T/chm13 v2 release: [Genbank CP086569.2](https://www.ncbi.nlm.nih.gov/nuccore/CP086569.2) and RefSeq NC_060948.1
    - sequence download: [Genbank](https://www.ncbi.nlm.nih.gov/nuccore/CP086569.2?report=fasta&to=62460029)
    - reference sequence length: 62460029 bp
3. HG002-Y subregion annotation: provided by Pille, copied from suppl. tables 21 and 22 in [Rhie et al., "T2T Y paper" DOI:10.1038/s41586-023-06457-y](https://pmc.ncbi.nlm.nih.gov/articles/PMC10752217/)

### GRCh38

The Y sequence class annotation includes N-gaps that were skipped when dumping the 'region db' file:

```
Pandas(Index=73, chrom='chrY', start=26644163, end=27078488, name='other2', fasta_header='chrY_hg38_other2') N gap: 93.3%
Pandas(Index=74, chrom='chrY', start=27078488, end=56887902, name='HET', fasta_header='chrY_hg38_HET') N gap: 99.4%
Pandas(Index=0, chrom='chrY', start=0, end=10000, name='01_unlabeled', fasta_header='chrY_hg38_01_unlabeled', size=10000) N gap: 100.0%
Pandas(Index=83, chrom='chrY', start=26644163, end=27078488, name='84_other2', fasta_header='chrY_hg38_84_other2', size=434325) N gap: 93.3%
Pandas(Index=84, chrom='chrY', start=27078488, end=56887902, name='85_HET', fasta_header='chrY_hg38_85_HET', size=29809414) N gap: 99.4%
Pandas(Index=86, chrom='chrY', start=57217416, end=57227415, name='87_unlabeled', fasta_header='chrY_hg38_87_unlabeled', size=9999) N gap: 100.0%
```

## Verkko assemblies - known issues

1. sample `NA21093`: contig `haplotype2-0000204` was identified as Chromosome Y and assigned to the haplotype 1 FASTA w/o renaming:

```
$ zgrep -F haplotype2-0000204 NA21093.assembly.refOriented.haplotype1.fasta.gz
>chrY_haplotype2-0000204
```

2. sample `NA19700`: several contigs identified as Chromosome Y in haplotype 2 w/o assigning them to the haplotype 1 FASTA:

```
haplotype1-0000009 0 54962076 chrY_haplotype1-0000009 104.563

haplotype2-0000217 0 562410 chrY_haplotype2-0000217 126.671
haplotype2-0000218 357424 0 chrY_haplotype2-0000218 118.864
haplotype2-0000219 367032 0 chrY_haplotype2-0000219 105.449
```

```
$ zgrep -F haplotype2-0000217 NA19700.assembly.refOriented.haplotype2.fasta.gz
>chrY_haplotype2-0000217
$ zgrep -F haplotype2-0000218 NA19700.assembly.refOriented.haplotype2.fasta.gz
>chrY_haplotype2-0000218
$ zgrep -F haplotype2-0000219 NA19700.assembly.refOriented.haplotype2.fasta.gz
>chrY_haplotype2-0000219
```

3. sample HG03270 likely exhibits an assembly error in (one of) the PAR region; sample is female, but best hit in minimap alignments is to chrY
    - see Snakefile `extract-sex-chromosomes::workflow::Snakefile::reassign_sequences_by_chrom`

```
sample  chromosome      source_name     name    name_id seq_length      num_A   num_C   num_G   num_T   num_N   top_hit orientation     top_matching_bp top_matching_pct
HG03270 chrX    chrX_haplotype2-0000104 chrX|HG03270|XX|UNK|hap2|SQN:0000104|FRG        3C94FF6F        245047  67079   56689   54635   66644   0       chrY    1       201156  82.09
```

4. sample HG00423 likely exhibits an assembly error in (one of) the PAR region; sample is female, but best hit in minimap alignments is to chrY
    - see Snakefile `extract-sex-chromosomes::workflow::Snakefile::reassign_sequences_by_chrom`

```
sample  chromosome      source_name     name    name_id seq_length      num_A   num_C   num_G   num_T   num_N   top_hit orientation     top_matching_bp top_matching_pct
HG00423 chrX    chrX_haplotype1-0000026 chrX|HG00423|XX|UNK|hap1|SQN:0000026|FRG        D2ACC43B        211785  52863   53550   51223   54149   0       chrY    1       170003  80.27
```
