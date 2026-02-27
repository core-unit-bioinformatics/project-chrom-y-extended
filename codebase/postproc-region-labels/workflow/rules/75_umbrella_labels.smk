
import pathlib

SUB_WD = WD.joinpath("80-statistics")

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

