
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
        bed = SUB_WD.joinpath("suppl", "gap_tracks", "{sample}.{ref}.gaps.bed")
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
        tsv = SUB_WD.joinpath("suppl", "seqclass_gaps", "{sample}.{ref}.chrY-regions.gaps.tsv")
    run:
        import pandas as pd
        gaps = pd.read_csv(input.gaps, sep="\t", header=None, names=["#seq", "start", "end"])

        gaps["name"] = "GAP"
        gaps["score"] = 0
        gaps["strand"] = "."
        gaps["thickStart"] = gaps["start"]
        gaps["thickEnd"] = gaps["end"]
        gaps["assign_method"] = "complement"
        gaps["second_best_guess"] = "GAP"
        gaps["kmer_top_enrich"] = 0.
        gaps["other_support"] = "none"
        gaps["other_orientation"] = "."
        gaps["cluster_id"] = -1

        labels = pd.read_csv(input.labels, sep="\t", header=0)
        labels = pd.concat([labels, gaps], axis=0, ignore_index=False)
        labels.sort_values(["#seq", "start", "end"], inplace=True)

        assert not pd.isnull(labels).any(axis=0).any()

        labels.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule intersect_labels_and_qc:
    input:
        qc_win = rules.merge_qc_track_intersections.output.tsv,
        labels = rules.merge_gaps_into_seqclass_labels.output.tsv
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
