"""The data flow of this module is a bit involved:
    HMMER motif annotations are produced at HHU, sent off to
    JAX/Mark Loftus, who filters the motif hits for the
    highest-quality ones and sends those lists back to
    HHU to have the remaining matches integrated into
    the Y sequence class annotations.
"""

SUB_WD = WD.joinpath("40-hmmer-motifs")

localrules: normalize_filtered_hmmer_hits
rule normalize_filtered_hmmer_hits:
    """This rule is highly application-specific because
    it depends on externally filtered high-quality HMMER
    hit (provided by Mark Loftus). These data come in slightly
    non-uniform format and must be normalized first.
    """
    input:
        hmmer_filtered = [
            pathlib.Path(PROJECT_CONFIG["remote_hilbert_prefix"]).joinpath(PROJECT_CONFIG["hmmer_filtered_dyz19"]),
            pathlib.Path(PROJECT_CONFIG["remote_hilbert_prefix"]).joinpath(PROJECT_CONFIG["hmmer_filtered_tspy"]),
            pathlib.Path(PROJECT_CONFIG["remote_hilbert_prefix"]).joinpath(PROJECT_CONFIG["hmmer_filtered_yq12"])
        ]
    output:
        hmmer_motifs = expand(
            SUB_WD.joinpath(
                "hmmer_filtered_hits", "{motif}",
                "{sample}.{motif}.chrY-hmmer-filtered-hits.bed"
            ),
            sample=SAMPLES,
            motif=["tspy", "dyz19", "yq12"]
        )
    run:

        for folder in input.hmmer_filtered:
            folder = pathlib.Path(folder).resolve(strict=True)
            motif_name = folder.name.replace("hmmerFilter", "")
            assert motif_name in ["tspy", "dyz19", "yq12"]
            for csv_file in folder.glob("**/*.csv"):
                if "old" in str(csv_file.parent):
                    continue
                sample = csv_file.name.split(".")[0]
                bed_df = load_ml_motif_hits(csv_file, motif_name)
                out_file = SUB_WD.joinpath(
                    "hmmer_filtered_hits", f"{motif_name}",
                    f"{sample}.{motif_name}.chrY-hmmer-filtered-hits.bed"
                )
                out_file.parent.mkdir(exist_ok=True, parents=True)
                bed_df.to_csv(out_file, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_norm_filtered_hmmer:
    input:
        tables = rules.normalize_filtered_hmmer_hits.output.hmmer_motifs
