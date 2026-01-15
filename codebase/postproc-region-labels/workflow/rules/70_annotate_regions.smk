
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

        # this is kept here for reference, but the erroneous k-mers are
        # added later with their original coordinates to keep the resolution
        select_kmer_dirty = df["kmer_errors_is_clean"] == 0  # False / not clean
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
            # see comment above about k-mers
            if False and row.error_base_win == 1:
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


rule label_and_merge_windows:
    input:
        tsv = rules.set_error_windows.output.tsv
    output:
        bed = SUB_WD.joinpath(
            "suppl", "seq_class_draft", "{sample}.{ref}.chrY-regions.err-struct.bed"
        ),
        tsv = SUB_WD.joinpath(
            "suppl", "win_merge_debug", "{sample}.{ref}.chrY-regions.err-struct.debug.tsv"
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


localrules: add_kmer_blocks
rule add_kmer_blocks:
    input:
        regions = rules.label_and_merge_windows.output.bed,
        kmers = expand(
            rules.filter_sequences_to_sex_chrom.output.bed,
            qc_track="kmer_errors",
            allow_missing=True
        )
    output:
        bed = SUB_WD.joinpath(
            "suppl", "seq_class_kmers", "{sample}.{ref}.chrY-regions.err-struct-base.bed"
        )
    run:
        import pandas as pd

        regions = pd.read_csv(input.regions, sep="\t", header=0)
        kmers = pd.read_csv(input.kmers[0], sep="\t", header=0, usecols=["#seq", "start", "end"])
        kmers["name"] = "ERRBASE"
        kmers["score"] = 0
        kmers["strand"] = "+"

        concat = pd.concat([regions, kmers], axis=0, ignore_index=False)
        concat.sort_values(["#seq", "start", "end"], inplace=True)
        na_cols = pd.isnull(concat).any(axis=0)
        if na_cols.any():
            col_names = concat.columns[na_cols]
            print(col_names)
            assert not pd.isnull(concat).any(axis=0).any()

        concat.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule fill_remaining_gaps:
    """There can still be some small gaps
    in the annotation following the way the
    error windows are defined (fixed 1 kbp
    boundaries).
    """
    input:
        sizes = rules.dump_genome_seq_sizes.output.tsv,
        regions = rules.add_kmer_blocks.output.bed
    output:
        bed = SUB_WD.joinpath("suppl", "fill_gaps", "{sample}.{ref}.fillers.bed")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    shell:
        "bedtools complement -i {input.regions} -g {input.sizes} > {output}"


localrules: add_gap_fillers_to_annotation
rule add_gap_fillers_to_annotation:
    input:
        gaps = rules.fill_remaining_gaps.output.bed,
        regions = rules.add_kmer_blocks.output.bed
    output:
        bed = SUB_WD.joinpath(
            "results", "seq_annotation",
            "{sample}.{ref}.chrY-regions.err-struct-base.bed"
        )
    run:
        import pandas as pd

        regions = pd.read_csv(input.regions, sep="\t", header=0)
        gaps = pd.read_csv(input.gaps, sep="\t", header=None, names=["#seq", "start", "end"])
        if gaps.empty:
            regions.to_csv(output.bed, sep="\t", header=True, index=False)
        else:
            gaps["score"] = 500
            gaps["name"] = "UNASSIGNED"
            gaps["strand"] = "+"
            regions = pd.concat([regions, gaps], axis=0, ignore_index=False)
            assert not pd.isnull(regions).any(axis=0).any()
            regions.sort_values(["#seq", "start", "end"], inplace=True)
            regions.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule check_all_bases_covered:
    input:
        sizes = rules.dump_genome_seq_sizes.output.tsv,
        regions = rules.add_gap_fillers_to_annotation.output.bed
    output:
        check = SUB_WD.joinpath("suppl", "sanity_check", "{sample}.{ref}.chrY-regions.check.ok")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import numpy as np

        regions = pd.read_csv(input.regions, sep="\t", header=0)
        seq_sizes = dict()
        with open(input.sizes) as listing:
            for line in listing:
                name, length = line.strip().split()
                seq_sizes[name] = int(length)

        for seq, seq_regions in regions.groupby("seq"):
            indicator = np.zeros(seq_sizes[seq], dtype=bool)
            for row in seq_regions.itertuples():
                indicator[row.start:row.end] |= True
            total_covered = indicator.sum()
            if total_covered != indicator.size:
                err_msg = (
                    f"Gaps remaining: {wildcards.sample} / {wildcards.ref}: "
                    f"Should: {indicator.size} - Is: {total_covered}"
                )
                raise ValueError(err_msg)
        with open(output.check):
            pass
    # END OF RUN BLOCK


rule run_all_annotate_regions:
    input:
        check = expand(
            rules.check_all_bases_covered.output.check,
            sample=SAMPLES,
            ref=list(MODULE_REF_GENOMES.keys())
        ),
