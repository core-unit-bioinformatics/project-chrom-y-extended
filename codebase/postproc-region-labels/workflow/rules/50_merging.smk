
SUB_WD = WD.joinpath("50-merging")

rule merge_aln_and_kmer_annotation:
    """The merge largely refers to merging the two
    different k-mer based annotation files. Apart from
    that, the "join" between alignment-based and k-mer-based
    sequence labels (essentially, for the Y repeat blocks)
    implements a priority for the k-mer-derived labels.
    No information is really discarded, just condensed.
    """
    input:
        aln_bed = rules.check_realign_precision.output.bed,
        strand_kmer_tsv = rules.normalize_kmer_based_color_annotation.output.tsv_dir,
        uniq_color_kmer_tsv = rules.normalize_color_unique_kmer_counts.output.tsv_dir
    output:
        bed = SUB_WD.joinpath(
            "00_aln_kmer",
            "{sample}.{ref}.chrY-regions.mrg.aln-kmer.bed"
        )
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    params:
        script=PROJECT_REPO_ROOT.joinpath(
            "codebase", "postproc-region-labels", "workflow",
            "scripts", "merge_kmer_annotation.py"
        ).resolve(strict=True),
        # 20-postprocess-labels/20-color-blocks-compare/kmer_annotation/norm_tsv/
        # NA21093.chrY-colorblock-kmers.bed
        kmer_strand_info = lambda wildcards, input: pathlib.Path(
            input.strand_kmer_tsv,
            f"{wildcards.sample}.chrY-colorblock-kmers.bed"
        ),
        # 20-postprocess-labels/30-color-blocks-uniqk/kmer_annotation/norm_tsv/
        # NA21093_chrY-regions_uniq-kmer-counts.tsv
        kmer_color_info = lambda wildcards, input: pathlib.Path(
            input.uniq_color_kmer_tsv,
            f"{wildcards.sample}_chrY-regions_uniq-kmer-counts.tsv"
        ),
    shell:
        "{params.script} --aln-region-labels {input.aln_bed} "
        "--kmer-strands {params.kmer_strand_info} "
        "--kmer-region-labels {params.kmer_color_info} "
        "--output {output.bed}"


localrules: merge_hmmer_hits_into_aln_kmer_regions
rule merge_hmmer_hits_into_aln_kmer_regions:
    """This merge simply replaces all regions
    that represent a motif also identified (placed)
    by HMMER with the respective annotation, which
    has been vetted by Mark Loftus.
    """
    input:
        regions = rules.merge_aln_and_kmer_annotation.output.bed,
        hmmer_hits = expand(
            WD.joinpath("40-hmmer-motifs").joinpath(
                "hmmer_filtered_hits", "{motif}",
                "{sample}.{motif}.chrY-hmmer-filtered-hits.bed"
            ),
            motif=["tspy", "dyz19", "yq12"],
            allow_missing=True
        )
    output:
        bed = SUB_WD.joinpath(
            "10_aln_kmer_hmmer",
            "{sample}.{ref}.chrY-regions.mrg.aln-kmer-hmmer.bed"
        )
    run:
        import pandas as pd
        import re

        hmmer = []
        for bed_file in input.hmmer_hits:
            motif = pd.read_csv(bed_file, sep="\t", header=0)
            hmmer.append(motif)
        hmmer = pd.concat(hmmer, axis=0, ignore_index=False)

        regions = pd.read_csv(input.regions, sep="\t", header=0)

        tag_motifs = build_redundant_motif_filter(set(
            ["tspy", "dyz19", "yq12"] + hmmer["name"].tolist()
        ))

        regions["is_replaced"] = regions["name"].apply(tag_motifs) | regions["second_best_guess"].apply(tag_motifs)
        regions = regions.loc[~regions["is_replaced"], :].copy()

        regions.drop("is_replaced", axis=1, inplace=True)

        hmmer["thickStart"] = hmmer["start"]
        hmmer["thickEnd"] = hmmer["end"]
        hmmer["assign_method"] = "hmmer"
        hmmer["second_best_guess"] = hmmer["name"]
        hmmer["kmer_top_enrich"] = 0
        hmmer["other_support"] = "unknown"
        hmmer["other_orientation"] = "."
        hmmer["cluster_id"] = -1

        regions = pd.concat([regions, hmmer], axis=0, ignore_index=False)
        regions.sort_values(["#seq", "start", "end"], inplace=True)

        assert not pd.isnull(regions).any(axis=0).any()

        regions.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_merging:
    input:
        mrg_aln_kmer = expand(
            rules.merge_aln_and_kmer_annotation.output.bed,
            sample=[sample for sample in SAMPLES if sample not in ["RFGRC38-R1", "RFCHM13-J1"]],
            ref=list(MODULE_REF_GENOMES.keys())
        ),
        mrg_aln_kmer_hmmer = expand(
            rules.merge_hmmer_hits_into_aln_kmer_regions.output.bed,
            sample=[sample for sample in SAMPLES if sample not in ["RFGRC38-R1", "RFCHM13-J1"]],
            ref=list(MODULE_REF_GENOMES.keys())
        )
