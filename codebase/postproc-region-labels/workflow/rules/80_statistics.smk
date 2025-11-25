
SUB_WD = WD.joinpath("80-statistics")


rule compute_label_dist_stats:
    input:
        labels = rules.label_and_merge_windows.output.bed,
        gsize = rules.dump_genome_seq_sizes.output.tsv
    output:
        tsv = SUB_WD.joinpath(
            "results", "label_dist_stats",
            "{sample}.{ref}.label-dist-stats.tsv"
        )
    run:
        import pandas as pd

        seq_sizes = {
            (row.seq, int(row.length)) for row in
            (pd.read_csv(
                input.gsize, header=None, sep="\t", names=["seq", "length"]
            )).itertuples()
        }
        labels = pd.read_csv(input.labels, header=0, sep="\t")
        labels["length"] = labels["end"] - labels["start"]

        agg = labels.groupby(["#seq", "name"])["length"].sum().reset_index(drop=False, inplace=False)
        agg["seq_length"] = agg["#seq"].replace(seq_sizes, inplace=False).astype(int)
        agg["pct_cov"] = (agg["length"] / agg["seq_length"] * 100).round(3)

        agg.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_label_dist_stats:
    input:
        tsv = expand(
            rules.compute_label_dist_stats.output.tsv,
            sample=SAMPLES,
            ref=list(MODULE_REF_GENOMES.keys())
        )
