
SUB_WD = WD.joinpath("60-qc-prep")


localrules: filter_sequences_to_sex_chrom
rule filter_sequences_to_sex_chrom:
    input:
        fasta = lambda wildcards: get_sample_file(SAMPLE_SHEET, wildcards.sample, "input_path"),
        qc_bed = lambda wildcards: get_sample_file(SAMPLE_SHEET, wildcards.sample, wildcards.qc_track)
    output:
        bed = SUB_WD.joinpath("{qc_track}", "{sample}_{qc_track}.bed")
    run:
        import dnaio
        import pandas as pd

        known_seqs = set()
        with dnaio.open(input.fasta) as fasta:
            for record in fasta:
                known_seqs.add(record.name)

        qc_header = [
            "seq", "start", "end", "label",
            "score", "strand",
            "thick_start", "thick_end", "color"
        ]
        qc_flagged_regions = pd.read_csv(input.qc_bed, sep="\t", header=None, names=qc_header)
        selector = qc_flagged_regions["seq"].isin(known_seqs)
        qc_flagged_regions = qc_flagged_regions.loc[selector, :].copy()
        assert not qc_flagged_regions.empty, f"No seqs selected: {known_seqs}"
        qc_flagged_regions.sort_values(["seq", "start", "end"], inplace=True)
        qc_flagged_regions.rename({"seq": "#seq"}, axis=1, inplace=True)
        qc_flagged_regions.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: dump_input_regions_bed
rule dump_input_regions_bed:
    input:
        fasta = lambda wildcards: get_sample_file(SAMPLE_SHEET, wildcards.sample, "input_path")
    output:
        bed = SUB_WD.joinpath("suppl", "regions", "{sample}.regions.bed")
    run:
        import dnaio

        regions = []
        with dnaio.open(input.fasta) as fasta:
            for record in fasta:
                regions.append(
                    "\t".join(
                        [
                            record.name,
                            "0",
                            str(len(record.sequence)),
                            record.name
                        ]
                    )
                )

        with open(output.bed, "w") as dump:
            _ = dump.write("\n".join(regions) + "\n")
    # END OF RUN BLOCK


rule make_input_regions_windows:
    input:
        bed = rules.dump_input_regions_bed.output.bed
    output:
        bed = SUB_WD.joinpath("suppl", "windows", "{sample}.win-1k.bed")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    shell:
        "bedtools makewindows -b {input.bed} -w 1000 -i srcwinnum > {output.bed}"


rule intersect_windows_qc_regions:
    input:
        windows = rules.make_input_regions_windows.output.bed,
        qc_track = rules.filter_sequences_to_sex_chrom.output.bed
    output:
        tsv = SUB_WD.joinpath("suppl", "isect", "{sample}_{qc_track}.win-1k.tsv")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bedtools intersect -wao -a {input.windows} -b {input.qc_track} > {output.tsv}"


localrules: normalize_qc_track_intersections
rule normalize_qc_track_intersections:
    input:
        tsv = rules.intersect_windows_qc_regions.output.tsv
    output:
        tsv = SUB_WD.joinpath("suppl", "norm_isect", "{sample}_{qc_track}.win-1k.norm.tsv")
    run:
        import pandas as pd

        def compute_window_ranking(df):
            """Compute rank bins for windows to summarize qc tracks
            by relative sequence position
            """
            rankings = []
            for seq, windows in df.groupby("seq"):
                ranks = (windows["end"].rank(pct=True) * 100).round(0).astype(int)
                rankings.append(ranks)

            df["rank_bin"] = -1
            for ranks in rankings:
                df.loc[ranks.index, "rank_bin"] = ranks.values

            return df

        def simplify_labeling(df):

            if wildcards.qc_track == "flagger":
                df["flagger_is_clean"] = 0
                df.loc[df["flagger_label"] == "Hap", "flagger_is_clean"] = 1
                # hold because of 1 kbp binning
                assert df.shape[0] == df["window"].nunique()
            elif wildcards.qc_track == "nucflag":
                df["nucflag_label"] = df["nucflag_label"].replace({".": "Hap"}, inplace=False)
                df["nucflag_is_clean"] = 0
                df.loc[df["nucflag_label"] == "Hap", "nucflag_is_clean"] = 1
                df.drop_duplicates(
                    subset=["seq", "window", "nucflag_label", "nucflag_is_clean"],
                    keep="first", inplace=True
                )

                set_new_labels = []
                drop_indices = []
                for window, regions in df.loc[df.duplicated(subset="window", keep=False), :].groupby("window"):
                    merged_labels = "|".join(sorted(set(regions["nucflag_label"].values)))
                    keep_index = regions.index.min()
                    drop_indices.extend(idx for idx in regions.index if idx != keep_index)
                    set_new_labels.append((keep_index, merged_labels))
                for idx, label in set_new_labels:
                    df.loc[idx, "nucflag_label"] = label
                df.drop(drop_indices, axis=0, inplace=True)
                df.reset_index(inplace=True, drop=True)
                assert df.shape[0] == df["window"].nunique()
            else:
                raise
            return df

        intersect_header = [
            "seq", "start", "end", "window",
            "seq2", "start2", "end2",
            f"{wildcards.qc_track}_label",
            "score", "strand", "thick_start", "thick_end",
            "color", "overlap_bp"
        ]
        usecols = [
            "seq", "start", "end", "window",
            f"{wildcards.qc_track}_label"
        ]

        qc_regions = pd.read_csv(
            input.tsv, sep="\t", header=None,
            names=intersect_header, usecols=usecols
        )
        qc_regions = simplify_labeling(qc_regions)
        qc_regions = compute_window_ranking(qc_regions)
        qc_regions.sort_values(["seq", "start", "end"], inplace=True)
        qc_regions.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: merge_qc_track_intersections
rule merge_qc_track_intersections:
    input:
        tables = expand(
            rules.normalize_qc_track_intersections.output.tsv,
            qc_track=["flagger", "nucflag"],
            allow_missing=True
        )
    output:
        tsv = SUB_WD.joinpath("suppl", "merge_qc", "{sample}_qclabels.win-1k.tsv"),
    run:
        import pandas as pd
        assert len(input.tables) == 2

        use_index = ["seq", "start", "end", "window", "rank_bin"]
        qc1 = pd.read_csv(input.tables[0], sep="\t", header=0, index_col=use_index)
        qc2 = pd.read_csv(input.tables[1], sep="\t", header=0, index_col=use_index)

        merge = qc1.join(qc2, how="outer")
        assert merge.shape[0] == qc1.shape[0] == qc2.shape[0]

        merge.to_csv(output.tsv, sep="\t", header=True, index=True)
    # END OF RUN BLOCK


rule run_all_prep_qc:
    input:
        qc_beds = expand(
            rules.filter_sequences_to_sex_chrom.output.bed,
            qc_track=["flagger", "nucflag"],
            sample=SAMPLES
        ),
        qc_labels = expand(
            rules.merge_qc_track_intersections.output.tsv
            sample=SAMPLES
        )
