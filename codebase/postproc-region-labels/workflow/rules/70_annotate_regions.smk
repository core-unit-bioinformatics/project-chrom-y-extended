
SUB_WD = WD.joinpath("70-annotate-regions")


rule intersect_labels_and_qc:
    input:
        qc_win = rules.merge_qc_track_intersections.output.tsv,
        labels = rules.merge_hmmer_hits_into_aln_kmer_regions.output.bed
    output:
        isect = SUB_WD.joinpath("label_qc_isect", "{sample}.{ref}.chrY-regions.qc-win.isect.tsv")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bedtools intersect -wao -a {input.labels} -b {input.qc_win} > {output}"


rule run_all_annotate_regions:
    input:
        mrg_aln_kmer = expand(
            rules.intersect_labels_and_qc.output.tsv,
            sample=SAMPLES,
            ref=list(MODULE_REF_GENOMES.keys())
        ),
