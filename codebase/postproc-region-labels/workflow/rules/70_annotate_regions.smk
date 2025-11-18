
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


localrules: set_error_windows
rule set_error_windows:
    input:
        isect = rules.intersect_labels_and_qc.output.isect
    output:
        tsv = SUB_WD.joinpath("suppl", "add_err_windows", "{sample}.{ref}.chrY-regions.qc-win.err-strict.tsv")
    run:
        import pandas as pd
        # reheader intersection
        columns = [
            "seq", "start", "end", "name", "score", "strand",
            "thickStart", "thickEnd", "assign_method",
            "second_best_guess", "kmer_top_enrich", "other_support",
            "other_orientation", "cluster_id",
            "seq2", "win_start", "win_end", "win_name", "win_pctile",
            "flagger_label", "flagger_is_clean",
            "nucflag_label", "nucflag_is_clean",
            "overlap_bp"
        ]
        df = pd.read_csv(input.isect, sep="\t", header=None, names=columns)
        df.drop(["cluster_id"], axis=1, inplace=True)

        select_flagger_dirty = df["flagger_is_clean"] == 0  # False / not clean
        select_nucflag_dirty = df["nucflag_is_clean"] == 0  # False / not clean

        # strict: require error flag from both tools
        select_strict_dirty = select_flagger_dirty & select_nucflag_dirty
        df["error_strict"] = 0
        df.loc[select_strict_dirty, "error_strict"] = 1
        # lenient: require only one error flag
        select_lenient_dirty = select_flagger_dirty | select_nucflag_dirty
        df["error_lenient"] = 0
        df.loc[select_lenient_dirty, "error_lenient"] = 1

        add_error_regions = []
        drop_rows = []
        for row in df.itertuples(index=True):
            # Decision from chrY call on Nov. 11.
            # only consider windows flagged by both
            # tools (strict criterion) as true errors
            if row.error_strict == 1:
                new_row = row._asdict()
                del new_row["Index"]
                # change seq coord to window coord
                new_row["start"] = row.win_start
                new_row["end"] = row.win_end
                new_row["name"] = "ERR"
                new_row["score"] = 0
                new_row["second_best_guess"] = "ERR"
                new_row["assign_method"] = "qcflag"
                add_error_regions.append(new_row)
                if row.name == "GAP":
                    # GAPs that are labeled as errors can
                    # just be dropped from the list of regions
                    drop_rows.append(row.Index)
        df.drop(drop_rows, axis=0, inplace=True)

        add_error_regions = pd.DataFrame.from_records(add_error_regions)
        # the mix-in of the curated hmmer calls can lead to duplicates;
        # drop them and arbitrarily keep the first error windows
        add_error_regions.drop_duplicates(["seq", "start", "end"], keep="first", inplace=True)
        df = pd.concat([df, add_error_regions], axis=0, ignore_index=False)
        df.sort_values(["seq", "start", "end"], axis=0, inplace=True)
        df.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_annotate_regions:
    input:
        err_win = expand(
            rules.set_error_windows.output.tsv,
            sample=SAMPLES,
            ref=list(MODULE_REF_GENOMES.keys())
        ),
