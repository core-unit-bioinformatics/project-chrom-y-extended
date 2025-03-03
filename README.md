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
