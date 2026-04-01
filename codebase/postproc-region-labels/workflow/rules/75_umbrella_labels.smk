
import pathlib

SUB_WD = WD.joinpath("75-umbrella-labels")

localrules: determine_umbrella_label_matchings
rule determine_umbrella_label_matchings:
    """This rule creates a mapping for the
    umbrella labels between reference and
    samples in an attempt to create some form
    of unified labeling that makes evaluation
    easier.
    Only for umbrella terms!
    """
    input:
        ref_labels = lambda wildcards: lambda wildcards: WD.joinpath("references", f"{MODULE_REF_LABELINGS[wildcards.ref]}"),
        smp_labels = rules.add_gap_fillers_to_annotation.output.bed,
        chk_ok = rules.check_all_bases_covered.output.check
    output:
        rename_ref = SUB_WD.joinpath(
            "results", "unified_umbrella",
            "{ref}", "{ref}.{sample}.unified-umbrella.json"
        ),
        rename_smp = SUB_WD.joinpath(
            "results", "unified_umbrella",
            "{ref}", "{sample}.{ref}.chrY-unified-umbrella.json"
        )
    params:
        script=PROJECT_REPO_ROOT.joinpath(
            "codebase", "postproc-region-labels", "workflow",
            "scripts", "match_umbrella_labels.py"
        ).resolve(strict=True)
    shell:
        "{params.script} --ref-regions {input.ref_labels} --sample-regions {input.smp_labels} "
        "--ref-out {output.rename_ref} --output {output.rename_smp}"


rule intersect_reference_labelings:
    """Decision chrY call 2026-03-31
    Merge both reference labelings for the umbrella
    terms to (ideally) close some of the unassigned
    regions by simple union.
    """
    input:
        ref_hg38 = expand(
            rules.add_gap_fillers_to_annotation.output.bed,
            ref="hg38",
            allow_missing=True
        ),
        ref_t2t = expand(
            rules.add_gap_fillers_to_annotation.output.bed,
            ref="t2tv2",
            allow_missing=True
        )
    output:
        isect = SUB_WD.joinpath("suppl", "ref_label_isect", "{sample}.uniref.chrY-regions.isect.tsv")
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    shell:
        "bedtools intersect -wao -a {input.ref_t2t} -b {input.ref_hg38} > {output.isect}"


localrules: simplify_umbrella_region_labels
rule simplify_umbrella_region_labels:
    input:
        rename_smp = rules.determine_umbrella_label_matchings.output.rename_smp,
        labels = rules.add_gap_fillers_to_annotation.output.bed
    output:
        disjoined = SUB_WD.joinpath(
            "results", "seq_annotation", "{ref}",
            "{sample}.{ref}.chrY-umbrella.ngaps.err-struct-base.disjoined.bed"
        ),
        stitched = SUB_WD.joinpath(
            "results", "seq_annotation", "{ref}",
            "{sample}.{ref}.chrY-umbrella.ngaps.err-struct-base.stitched.bed"
        ),
    params:
        script = PROJECT_REPO_ROOT.joinpath(
            "codebase", "postproc-region-labels", "workflow",
            "scripts", "simplify_umbrella_labels.py"
        ).resolve(strict=True)
    shell:
        "{params.script} --region-labels {input.labels} --umbrella-labels {input.rename_smp} "
        "--disjoined {output.disjoined} --stitched {output.stitched}"


rule run_all_umbrella_computations:
    input:
        renamer = expand(
            rules.determine_umbrella_label_matchings.output,
            ref=list(MODULE_REF_GENOMES.keys()),
            sample=SAMPLES
        ),
        merged = expand(
            rules.simplify_umbrella_region_labels.output,
            ref=list(MODULE_REF_GENOMES.keys()),
            sample=SAMPLES
        ),
        uniref = expand(
            rules.intersect_reference_labelings.output,
            sample=SAMPLES
        )
