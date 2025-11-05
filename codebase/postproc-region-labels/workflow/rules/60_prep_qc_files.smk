
SUB_WD = WD.joinpath("60-qc-prep")


rule filter_sequences_to_sex_chrom:
    input:
        bed = rules.merge_aln_and_kmer_annotation.output.bed,
        qc_bed = lambda wildcards: get_sample_file(SAMPLE_SHEET, wildcards.sample, wildcards.qc_track)
    output:
        bed = SUB_WD.joinpath("{qc_track}", "{sample}_{qc_track}.bed")
    run:
        import pandas as pd

        known_seqs = pd.read_csv(input.bed, sep="\t", header=0)
        known_seqs = set(known_seqs["#seq"].values)
        qc_header = [
            "seq", "start", "end", "label",
            "score", "strand",
            "thick_start", "thick_end", "color"
        ]
        qc_flagged_regions = pd.read_csv(input.qc_bed, sep="\t", header=None, names=qc_header)
        selector = qc_flagged_regions["seq"].isin(known_seqs)
        qc_flagged_regions = qc_flagged_regions.loc[selector, :].copy()
        qc_flagged_regions.sort_values(["seq", "start", "end"], inplace=True)
        qc_flagged_regions.rename({"seq": "#seq"}, axis=1, inplace=True)
        qc_flagged_regions.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_prep_qc:
    input:
        qc_beds = expand(
            rules.filter_sequences_to_sex_chrom.output.bed,
            qc_track=["flagger", "nucflag"],
            sample=SAMPLES
        )
