
SUB_WD = WD.joinpath("70-annotate-regions")


localrules: dump_genome_seq_sizes
rule dump_genome_seq_sizes:
    input:
        bed = rules.dump_input_regions_bed.output.bed
    output:
        tsv = SUB_WD.joinpath("suppl", "genome_sizes", "{sample}.seq-sizes.tsv")
    run:
        with open(input.bed) as regions:
            with open(output.tsv, "w") as sizes:
                for line in regions:
                    columns = line.strip().split()
                    sizes.write(f"{columns[0]}\t{columns[2]}\n")
    # END OF RUN BLOCK


rule create_gap_track:
    input:
        sizes = rules.dump_genome_seq_sizes.output.tsv,
        regions = rules.merge_hmmer_hits_into_aln_kmer_regions.output.bed
    output:
        bed = SUB_WD.joinpath("suppl", "gap_tracks", "{sample}_gaps.bed")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    shell:
        "bedtools complement -i {input.regions} -g {input.sizes} > {output}"


localrules: merge_gaps_into_seqclass_labels
rule merge_gaps_into_seqclass_labels:
    input:
        labels = rules.merge_hmmer_hits_into_aln_kmer_regions.output.bed,
        gaps = rules.create_gap_track.output.bed
    output:
        bed = SUB_WD.joinpath("suppl", "seqclass_gaps", "{sample}.{ref}.chrY-regions.gaps.tsv")
    run:
        raise


rule intersect_labels_and_qc:
    input:
        qc_win = rules.merge_qc_track_intersections.output.tsv,
        labels = rules.merge_gaps_into_seqclass_labels.output.bed
    output:
        isect = SUB_WD.joinpath("suppl", "label_qc_isect", "{sample}.{ref}.chrY-regions.qc-win.isect.tsv")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bedtools intersect -wao -a {input.labels} -b {input.qc_win} > {output}"


rule run_all_annotate_regions:
    input:
        isect = expand(
            rules.intersect_labels_and_qc.output.isect,
            sample=SAMPLES,
            ref=list(MODULE_REF_GENOMES.keys())
        ),
