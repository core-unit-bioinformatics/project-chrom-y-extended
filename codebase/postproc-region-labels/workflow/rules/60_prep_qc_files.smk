
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


rule run_all_prep_qc:
    input:
        qc_beds = expand(
            rules.filter_sequences_to_sex_chrom.output.bed,
            qc_track=["flagger", "nucflag"],
            sample=SAMPLES
        )
