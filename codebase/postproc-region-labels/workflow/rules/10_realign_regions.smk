"""
This module:
check quality of labeled regions by realigning them
back to the reference genome and assessing the precision
of that alignment.
"""


rule merge_region_by_coord_and_name:
    """The script first simplifies the region labels
    and then merges them by overlap (default: 1 bp) and name.

    The script of this rule is inherently "dumb". The label
    simplification corrects for the naming differences
    between the "region DB" and the "labeled ref" approach.
    Apart from that, it just uses pyranges to merge regions
    by position and name.

    """
    input:
        sample_region_labels=lambda wildcards: WF_SEQANN_WD.joinpath(
            "results", "annotations", "combined"
        ).glob(f"{wildcards.sample}.*cmb-{wildcards.ref}.concat.bed.gz")
    output:
        bed_file = WD.joinpath(
            "region_labels", "merged",
            "{sample}.{ref}.chrY-regions-merged.bed.gz"
        )
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    params:
        script=PROJECT_REPO_ROOT.joinpath(
            "codebase", "postproc-region-labels", "workflow",
            "scripts", "merge_by_ovl_and_name.py"
        ).resolve(strict=True)
    shell:
        "{params.script} --bed-in {input.sample_region_labels} "
        "--merge-dist 1 --bed-out {output.bed_file}"


rule extract_labeled_sequences:
    """Based on the merged region windows (previous rule), we extract the
    labeled sequences again from the de novo assemblies to align them
    back to the references. If the sequences precisely realign to the
    respective reference region, we can be quite certain that the assembled
    sequence represents that chrY sequence class (irrespective of any
    errors that are still in the sequence).
    """
    input:
        bed_file = rules.merge_region_by_coord_and_name.output.bed_file,
        fasta = lambda wildcards: get_sample_file(SAMPLE_SHEET, wildcards.sample, "input_path")
    output:
        fasta = WD.joinpath(
            "region_seqs", "{sample}.{ref}.chrY-regions.fasta"
        ),
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed_file} -name -fo {output.fasta}"


rule ref_align_extracted_sequences:
    input:
        fasta = rules.extract_labeled_sequences.output.fasta,
        ref = lambda wildcards: WD.joinpath("references", f"{MODULE_REF_GENOMES[wildcards.ref]}")
    output:
        paf = WD.joinpath(
            "region_alignments",
            "{sample}.{ref}.chrY-regions-realigned.paf.gz"
        )
    conda:
        GLOBAL_CONDA_ENVS.joinpath("align_tools.yaml")
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt * attempt,
    shell:
        "minimap2 -x asm20 --eqx -c -t {threads} --secondary=no "
        "{input.ref} {input.fasta} | gzip > {output.paf}"


rule normalize_realigned_region_seqs:
    input:
        paf = rules.ref_align_extracted_sequences.output.paf
    output:
        tsv = WD.joinpath(
            "region_alignments",
            "{sample}.{ref}.chrY-regions-realigned.norm-paf.tsv.gz"
        )
    conda:
        GLOBAL_CONDA_ENVS.joinpath("align_tools.yaml")
    params:
        script=GLOBAL_SCRIPTS.joinpath("normalizer", "normalize_paf.py")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule check_realign_precision:
    """TODO
    this script is a bit involved - simplify
    """
    input:
        tsv = rules.normalize_realigned_region_seqs.output.tsv,
        labels = lambda wildcards: lambda wildcards: WD.joinpath("references", f"{MODULE_REF_LABELINGS[wildcards.ref]}")
    output:
        tmp = temp(WD.joinpath(
            "region_alignments",
            "{sample}.{ref}.chrY-regions-realigned.tmp.bed"
        )),
        bed = WD.joinpath(
            "region_alignments",
            "{sample}.{ref}.chrY-regions-realigned.bed.gz"
        ),
    conda:
        GLOBAL_CONDA_ENVS.joinpath("seqtools.yaml")
    params:
        script=PROJECT_REPO_ROOT.joinpath(
            "codebase", "postproc-region-labels", "workflow",
            "scripts", "check_realign_precision.py"
        ).resolve(strict=True)
    shell:
        "{params.script} --input-aln {input.tsv} --input-regions {input.labels} "
        "--output-regions {output.tmp}"
            " && "
        "bgzip -c {output.tmp} > {output.bed}"
            " && "
        "tabix -p bed {output.bed}"


rule run_all_realign_regions:
    input:
        bed = expand(
            rules.check_realign_precision.output.bed,
            sample=sample=[sample for sample in SAMPLES if sample not in ["RFGRC38-R1", "RFCHM13-J1"]],
            ref=list(MODULE_REF_GENOMES.keys())
        )
