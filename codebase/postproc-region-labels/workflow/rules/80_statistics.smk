
import pathlib

SUB_WD = WD.joinpath("80-statistics")


localrules: compare_to_previous
rule compare_to_previous:
    """This exists essentially for debugging

    The || true at the end is necessary because non-empty diffs
    yield an exit of 1, which will be turned into a failed rule
    by Snakemake
    """
    input:
        prev_labels = pathlib.Path(
            "/gpfs/project/projects/medbioinf/data/00_RESTRUCTURE/shares/globus/outgoing/hgsvc/sig_chry/v2/verkko-v2.2.1/annotations/seq_classes/2026-01_final",
            "{ref}", "{sample}.{ref}.chrY-regions.ngaps.err-struct-base.bed"
        ),
        curr_labels = rules.add_gap_fillers_to_annotation.output.bed
    output:
        diff = SUB_WD.joinpath(
            "suppl", "diff_prev_curr",
            "{ref}", "{sample}.{ref}.delta.txt"
        )
    shell:
        "diff --suppress-common-lines {input.prev_labels} {input.curr_labels} > {output.diff} || true"


rule run_all_diff_previous:
    input:
        txt = expand(
            rules.compare_to_previous.output.diff,
            sample=SAMPLES,
            ref=list(MODULE_REF_GENOMES.keys())
        )


### BELOW: actual statistics for manuscript


localrules: compute_label_dist_stats
rule compute_label_dist_stats:
    input:
        check = rules.check_all_bases_covered.output.check,
        labels = rules.add_gap_fillers_to_annotation.output.bed,
        gsize = rules.dump_genome_seq_sizes.output.tsv
    output:
        tsv = SUB_WD.joinpath(
            "suppl", "label_dist_stats",
            "{sample}.{ref}.label-dist-stats.tsv"
        )
    run:
        import pandas as pd
        pd.set_option('future.no_silent_downcasting', True)

        seq_sizes = dict(
            (row.seq, int(row.length)) for row in
            (pd.read_csv(
                input.gsize, header=None, sep="\t", names=["seq", "length"]
            )).itertuples()
        )
        labels = pd.read_csv(input.labels, header=0, sep="\t")
        labels["length"] = labels["end"] - labels["start"]
        total_length = sum(seq_sizes.values())

        # aggregate on per-seq level
        agg = labels.groupby(["#seq", "name"])["length"].sum().reset_index(drop=False, inplace=False)
        agg["seq_length"] = agg["#seq"].replace(seq_sizes, inplace=False).astype(int)
        agg["seq_pct_cov"] = (agg["length"] / agg["seq_length"] * 100).round(5)

        # aggregate on all-seq level
        agg2 = labels.groupby("name")["length"].sum().reset_index(drop=False, inplace=False)
        agg2["total_pct_cov"] = (agg2["length"] / total_length * 100).round(5)
        agg2 = agg2.set_index("name", inplace=False).to_dict()["total_pct_cov"]

        agg["total_pct_cov"] = agg["name"].replace(agg2, inplace=False).astype(float)

        agg.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: merge_label_dist_stats
rule merge_label_dist_stats:
    input:
        stats_tables = expand(
            rules.compute_label_dist_stats.output.tsv,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        tsv = SUB_WD.joinpath(
            "results", "ref_merged_dist_stats",
            "{ref}.label-dist-stats.tsv"
        )
    run:
        import pathlib as pl
        import pandas as pd

        def set_seq_type(seq_name):
            if seq_name.endswith("_chrY"):
                return "main"
            else:
                assert "random" in seq_name
                return "rand"

        merged = []
        for table_file in input.stats_tables:
            sample = pl.Path(table_file).name.split(".")[0]
            df = pd.read_csv(table_file, sep="\t", header=0)
            df["sample"] = sample
            df["seq_type"] = df["#seq"].apply(set_seq_type)
            merged.append(df)
        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged.sort_values(["sample", "#seq", "name"], inplace=True)
        merged.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_label_dist_stats:
    input:
        tsv = expand(
            rules.merge_label_dist_stats.output.tsv,
            ref=list(MODULE_REF_GENOMES.keys())
        )


### self-overlap statistics

rule compute_region_self_overlap:
    input:
        regions = rules.add_gap_fillers_to_annotation.output.bed
    output:
        isect = SUB_WD.joinpath(
            "suppl", "self_overlap", "{sample}.{ref}.annot-isect.tsv"
        )
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bedtools intersect -wo -a {input.regions} -b {input.regions} > {output.isect}"


localrules: aggregate_self_overlap_table
rule aggregate_self_overlap_table:
    input:
        tsv = rules.compute_region_self_overlap.output.isect
    output:
        tsv = SUB_WD.joinpath(
            "suppl", "agg_self_ovl", "{sample}.{ref}.agg-isect.tsv"
        )
    run:
        import pandas as pd

        plain_header = ["seq", "start", "end", "name", "score", "strand"]
        header1 = [f"{hd}1" for hd in plain_header]
        header2 = [f"{hd}2" for hd in plain_header]
        header = header1 + header2 + ["overlap_bp"]

        df = pd.read_csv(input.tsv, sep="\t", header=None, names=header)
        # drop self-overlap
        df = df.loc[df["name1"] != df["name2"], :].copy()

        grouping = header1 + ["name2"]

        agg = df.groupby(grouping)["overlap_bp"].sum()
        agg = agg.reset_index(drop=False, inplace=False)
        agg["length"] = agg["end1"] - agg["start1"]
        agg["overlap_pct"] = (agg["overlap_bp"] / agg["length"] * 100).round(3)

        reheader = [c.strip("1") for c in agg.columns]
        agg.columns = reheader

        agg.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: merge_self_overlap_stats
rule merge_self_overlap_stats:
    input:
        ovl_stats = expand(
            rules.aggregate_self_overlap_table.output.tsv,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        tsv = SUB_WD.joinpath(
            "results", "ref_merged_ovl_stats",
            "{ref}.self-ovl-stats.tsv"
        )
    run:
        import pathlib as pl
        import pandas as pd

        def set_seq_type(seq_name):
            if seq_name.endswith("_chrY"):
                return "main"
            else:
                assert "random" in seq_name
                return "rand"

        merged = []
        for table_file in input.ovl_stats:
            sample = pl.Path(table_file).name.split(".")[0]
            df = pd.read_csv(table_file, sep="\t", header=0)
            df["sample"] = sample
            df["seq_type"] = df["seq"].apply(set_seq_type)
            merged.append(df)

        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged.sort_values(["sample", "seq", "start"], inplace=True)
        merged.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


# for convenience only
# intersect the final regions with the early variant
# after adding error annotations to check which region
# type exhibits most errors and N gaps

rule intersect_final_regions_with_draft:
    input:
        draft = rules.add_kmer_blocks.output.bed,
        final = rules.add_gap_fillers_to_annotation.output.bed
    output:
        isect = SUB_WD.joinpath(
            "suppl", "isect_final_draft", "{sample}.{ref}.final-draft-isect.tsv"
        )
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bedtools intersect -wo -a {input.draft} -b {input.final} > {output.isect}"


localrules: compute_assign_umbrella_errors
rule compute_assign_umbrella_errors:
    input:
        isect = rules.intersect_final_regions_with_draft.output.isect,
        json = rules.determine_umbrella_label_matchings.output.rename_smp
    output:
        tsv = SUB_WD.joinpath(
            "suppl", "umbrella_errors", "{sample}.{ref}.umbrella-errors.tsv"
        )
    params:
        script=PROJECT_REPO_ROOT.joinpath(
            "codebase", "postproc-region-labels", "workflow",
            "scripts", "match_umbrella_labels.py"
        ).resolve(strict=True)
    shell:
        "{params.script} --intersect {input.isect} --umbrella-labels {input.json} --output {output.tsv}"


localrules: aggregate_umbrella_errors
rule aggregate_umbrella_errors:
    input:
        tsv = expand(
            rules.compute_assign_umbrella_errors.output.tsv,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        tsv = SUB_WD.joinpath(
            "results", "umbrella_errors", "{ref}.umbrella-errors.tsv"
        )
    run:
        import pathlib as pl
        import pandas as pd

        concat = []
        for table in input.tsv:
            sample = pl.Path(table).name.split(".")[0]
            assert sample in SAMPLES
            df = pd.read_csv(table, sep="\t", header=0)
            df["sample"] = sample
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["sample", "seq", "start"], inplace=True)
        concat.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: aggregate_final_draft_overlap_table
rule aggregate_final_draft_overlap_table:
    input:
        tsv = rules.intersect_final_regions_with_draft.output.isect
    output:
        tsv = SUB_WD.joinpath(
            "suppl", "isect_final_draft", "{sample}.{ref}.agg-final-draft-isect.tsv"
        )
    run:
        import pandas as pd

        plain_header = ["seq", "start", "end", "name", "score", "strand"]
        header1 = [f"{hd}1" for hd in plain_header]
        header2 = [f"{hd}2" for hd in plain_header]
        header = header1 + header2 + ["overlap_bp"]

        df = pd.read_csv(input.tsv, sep="\t", header=None, names=header)
        # drop self-overlap
        df = df.loc[df["name1"] != df["name2"], :].copy()

        grouping = header1 + ["name2"]

        agg = df.groupby(grouping)["overlap_bp"].sum()
        agg = agg.reset_index(drop=False, inplace=False)
        agg["length"] = agg["end1"] - agg["start1"]
        agg["overlap_pct"] = (agg["overlap_bp"] / agg["length"] * 100).round(3)

        reheader = [c.strip("1") for c in agg.columns]
        agg.columns = reheader

        agg.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: merge_final_draft_overlap_stats
rule merge_final_draft_overlap_stats:
    input:
        ovl_stats = expand(
            rules.aggregate_final_draft_overlap_table.output.tsv,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        tsv = SUB_WD.joinpath(
            "results", "ref_merged_final_draft_ovl",
            "{ref}.final-draft-ovl-stats.tsv"
        )
    run:
        import pathlib as pl
        import pandas as pd

        def set_seq_type(seq_name):
            if seq_name.endswith("_chrY"):
                return "main"
            else:
                assert "random" in seq_name
                return "rand"

        merged = []
        for table_file in input.ovl_stats:
            sample = pl.Path(table_file).name.split(".")[0]
            df = pd.read_csv(table_file, sep="\t", header=0)
            df["sample"] = sample
            df["seq_type"] = df["seq"].apply(set_seq_type)
            merged.append(df)

        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged.sort_values(["sample", "seq", "start"], inplace=True)
        merged.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_self_overlaps:
    input:
        tsv_self = expand(
            rules.merge_self_overlap_stats.output.tsv,
            ref=list(MODULE_REF_GENOMES.keys())
        ),
        tsv_draft = expand(
            rules.merge_final_draft_overlap_stats.output.tsv,
            ref=list(MODULE_REF_GENOMES.keys())
        ),
        tsv_umbrella = expand(
            rules.aggregate_umbrella_errors.output.tsv,
            ref=list(MODULE_REF_GENOMES.keys())
        )
