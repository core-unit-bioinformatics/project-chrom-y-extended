
SUB_WD = WD.joinpath("70-annotate-regions")


localrules: dump_genome_seq_sizes
rule dump_genome_seq_sizes:
    input:
        bed = rules.dump_input_regions_bed.output.bed
    output:
        tsv = SUB_WD.joinpath("suppl", "genome_sizes", "{sample}.seq-sizes.tsv")
    run:
        seq_sizes = []
        with open(input.bed) as regions:
            for line in regions:
                columns = line.strip().split()
                seq_sizes.append((columns[0], columns[2]))
        with open(output.tsv, "w") as sizes:
            for name, size in sorted(seq_sizes):
                sizes.write(f"{name}\t{size}\n")
    # END OF RUN BLOCK


rule create_gap_track:
    input:
        sizes = rules.dump_genome_seq_sizes.output.tsv,
        regions = rules.merge_centromere_into_seqclasses.output.bed
    output:
        bed = SUB_WD.joinpath("suppl", "gap_tracks", "{sample}.{ref}.gaps.bed")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    shell:
        "bedtools complement -i {input.regions} -g {input.sizes} > {output}"


localrules: merge_gaps_into_seqclass_labels
rule merge_gaps_into_seqclass_labels:
    input:
        labels = rules.merge_centromere_into_seqclasses.output.bed,
        gaps = rules.create_gap_track.output.bed
    output:
        tsv = SUB_WD.joinpath("suppl", "seqclass_gaps", "{sample}.{ref}.chrY-regions.gaps.tsv"),
        header = SUB_WD.joinpath("suppl", "seqclass_gaps", "{sample}.{ref}.chrY-regions.gaps.header"),
    run:
        import pandas as pd
        gaps = pd.read_csv(input.gaps, sep="\t", header=None, names=["#seq", "start", "end"])

        gaps["name"] = "UNASSIGNED"
        gaps["score"] = 0
        gaps["strand"] = "."
        gaps["thickStart"] = gaps["start"]
        gaps["thickEnd"] = gaps["end"]
        gaps["assign_method"] = "complement"
        gaps["second_best_guess"] = "UNASSIGNED"
        gaps["kmer_top_enrich"] = 0.
        gaps["other_support"] = "none"
        gaps["other_orientation"] = "."
        gaps["cluster_id"] = -1

        labels = pd.read_csv(input.labels, sep="\t", header=0)
        labels = pd.concat([labels, gaps], axis=0, ignore_index=False)
        labels.sort_values(["#seq", "start", "end"], inplace=True)

        assert not pd.isnull(labels).any(axis=0).any()

        labels.to_csv(output.tsv, sep="\t", header=True, index=False)

        # in prep for rule set_error_windows
        labels.rename({"#seq": "seq"}, axis=1, inplace=True)
        with open(output.header, "w") as dump:
            _ = dump.write(",".join(labels.columns) + "\n")
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


localrules: set_error_windows
rule set_error_windows:
    input:
        isect = rules.intersect_labels_and_qc.output.isect,
        qc_header = rules.merge_qc_track_intersections.output.header,
        seqclass_header = rules.merge_gaps_into_seqclass_labels.output.header,
    output:
        tsv = SUB_WD.joinpath("suppl", "add_err_windows", "{sample}.{ref}.chrY-regions.qc-win.err-struct-basewin.tsv")
    run:
        import pandas as pd

        def load_header(fp):
            with open(fp) as hd:
                columns = hd.readline().strip().split(",")
            return columns

        columns = load_header(input.seqclass_header) + load_header(input.qc_header) + ["overlap_bp"]

        df = pd.read_csv(input.isect, sep="\t", header=None, names=columns)
        df.drop(["cluster_id"], axis=1, inplace=True)

        select_flagger_hifi_dirty = df["flagger_hifi_is_clean"] == 0  # False / not clean
        select_flagger_ont_dirty = df["flagger_ont_is_clean"] == 0  # False / not clean
        select_nucflag_hifi_dirty = df["nucflag_hifi_is_clean"] == 0  # False / not clean
        select_nucflag_ont_dirty = df["nucflag_ont_is_clean"] == 0  # False / not clean

        select_kmer_dirty = df["kmer_errors_is_clean"] == 0  # False / not clean

        # 2026-01-06
        # decision: flagger/nucflag errors are labeled as "structural errors";
        # bases flagged by erroneous k-mers are labeled as "base errors"
        # strict: require error flag from all tools
        select_struct_errors = (
            select_flagger_hifi_dirty
            &
            select_flagger_ont_dirty
            &
            select_nucflag_hifi_dirty
            &
            select_nucflag_ont_dirty
        )
        df["error_struct"] = 0
        df.loc[select_struct_errors, "error_struct"] = 1

        df["error_base_win"] = 0
        df.loc[select_kmer_dirty, "error_base_win"] = 1

        add_error_regions = []
        drop_rows = set()
        for row in df.itertuples(index=True):
            # Decision from chrY call on Nov. 11.
            # only consider windows flagged by both
            # tools (strict criterion) as true errors
            if row.error_struct == 1:
                err_label = "ERRSTRUCT"
                assign_method = "qcflag"
                new_row = row._asdict()
                del new_row["Index"]
                # change seq coord to window coord
                new_row["start"] = row.win_start
                new_row["end"] = row.win_end
                new_row["name"] = err_label
                new_row["score"] = 0
                new_row["second_best_guess"] = err_label
                new_row["assign_method"] = assign_method
                add_error_regions.append(new_row)
                if row.name == "UNASSIGNED":
                    # unassigned blocks that are labeled as errors can
                    # just be dropped from the list of regions
                    drop_rows.add(row.Index)
            if row.error_base_win == 1:
                err_label = "ERRBASEWIN"
                assign_method = "kmer"
                new_row = row._asdict()
                del new_row["Index"]
                # change seq coord to window coord
                new_row["start"] = row.win_start
                new_row["end"] = row.win_end
                new_row["name"] = err_label
                new_row["score"] = 0
                new_row["second_best_guess"] = err_label
                new_row["assign_method"] = assign_method
                add_error_regions.append(new_row)
                if row.name == "UNASSIGNED":
                    # unassigned blocks that are labeled as errors can
                    # just be dropped from the list of regions
                    drop_rows.add(row.Index)

        df.drop(drop_rows, axis=0, inplace=True)

        add_error_regions = pd.DataFrame.from_records(add_error_regions)
        # the mix-in of the curated hmmer calls can lead to duplicates;
        # drop them and arbitrarily keep the first error windows
        add_error_regions.drop_duplicates(["seq", "start", "end", "name"], keep="first", inplace=True)
        df = pd.concat([df, add_error_regions], axis=0, ignore_index=False)
        df.sort_values(["seq", "start", "end"], axis=0, inplace=True)
        df.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: add_kmer_high_res_blocks
rule add_kmer_high_res_blocks:
    input:
        kmer_track = expand(
            rules.filter_sequences_to_sex_chrom.output.bed,
            qc_track="kmer_errors",
            allow_missing=True
        ),
        window_track = rules.set_error_windows.output.tsv
    output:
        tsv = SUB_WD.joinpath("suppl", "add_high_res_kmer", "{sample}.{ref}.chrY-regions.qc-win.err-struct-base-win.tsv")
    run:
        import pandas as pd
        # NB: due to the expand (?), input.kmer_track
        # is a Snakemake Namedlist
        kmers = pd.read_csv(
            input.kmer_track[0], sep="\t", header=0,
            usecols=[
                "#seq", "start", "end", "strand",
            ]
        )
        kmers.rename({"#seq": "seq"}, axis=1, inplace=True)
        kmers["seq2"] = kmers["seq"]
        kmers["thickStart"] = kmers["start"]
        kmers["thickEnd"] = kmers["end"]
        kmers["name"] = "ERRBASE"
        kmers["score"] = 0
        kmers["second_best_guess"] = "ERRBASE"
        kmers["assign_method"] = "kmer"
        kmers["error_base"] = 1
        kmers["error_base_win"] = -1
        kmers["error_struct"] = -1
        kmers["kmer_top_enrich"] = 0.
        kmers["other_support"] = "none"
        kmers["other_orientation"] = "."
        kmers["win_pctile"] = -1
        kmers["overlap_bp"] = 0
        kmers["win_start"] = -1
        kmers["win_end"] = -1
        kmers["win_name"] = "UNK"

        err_win = pd.read_csv(input.window_track, sep="\t", header=0)
        concat = pd.concat([err_win, kmers], axis=0, ignore_index=False)
        concat.sort_values(["seq", "start", "end"], inplace=True)
        na_cols = pd.isna(concat).any(axis=0)
        if na_cols.any():
            column_names = concat.columns[na_cols]
            for cn in column_names:
                if "is_clean" in cn:
                    concat[cn] = concat[cn].fillna(-1)
                if "_label" in cn:
                    concat[cn] = concat[cn].fillna("UNK")
            print(column_names)
            raise ValueError(f"missing values: {column_names}")
        concat.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule label_and_merge_windows:
    input:
        tsv = rules.set_error_windows.output.tsv
    output:
        bed = SUB_WD.joinpath(
            "results", "seq_class", "{sample}.{ref}.chrY-regions.err-strict.bed"
        ),
        tsv = SUB_WD.joinpath(
            "suppl", "win_merge_debug", "{sample}.{ref}.chrY-regions.err-strict.debug.tsv"
        )
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=PROJECT_REPO_ROOT.joinpath(
            "codebase", "postproc-region-labels", "workflow",
            "scripts", "finalize_labels.py"
        ).resolve(strict=True)
    shell:
        "{params.script} --isect-table {input.tsv} --output {output.bed} --debug-out {output.tsv}"


rule run_all_annotate_regions:
    input:
        label_beds = expand(
            rules.add_kmer_high_res_blocks.output.tsv,
            sample=SAMPLES,
            ref=list(MODULE_REF_GENOMES.keys())
        ),

        # label_beds = expand(
        #     rules.label_and_merge_windows.output.bed,
        #     sample=SAMPLES,
        #     ref=list(MODULE_REF_GENOMES.keys())
        # ),
