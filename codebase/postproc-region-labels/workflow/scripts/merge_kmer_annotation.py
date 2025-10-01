#!/usr/bin/env python3

import argparse as argp
import collections as col
import difflib as diffl
import pathlib as pl
import re

import pandas as pd
import pyranges as pr


LABEL_ORDER_PREFIX = re.compile(r"^[0-9]{2}(n|u|s)_")

UNCERTAIN_SCRAMBLED = "".join(sorted("uncertain"))

def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--aln-region-labels", "-aln", "-a",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_aln",
        required=True,
        help=(
            "The labeled ref. / alignment-based inferred region "
            "labels. In the Snakefile, this is the output of the "
            "rule: check_realign_precision"
        )
    )
    parser.add_argument(
        "--kmer-region-labels", "-kmer", "-k",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_kmer",
        required=True,
        help=(
            "These files contain the k-mer based repeat block "
            "annotation/structrue as determined by Mark Loftus. "
            "(without orientation info). In addition, these files "
            "list unique k-mers per region/repeat block. "
            "In the Snakefile, this is the output of the agg rule: "
            "run_all_colorblocks_uniqk"
        )
    )

    parser.add_argument(
        "--kmer-strands", "-ks",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_kmer_strands",
        required=True,
        help=(
            "These files contain full orientation info and represent the output "
            "of Mark Loftus' comparison of his k-mer based repeat annotation "
            "to the annotation that Arang Rhie produces (typically, these two "
            "do agree). In the Snakefile, this is the output of the agg "
            "rule: run_all_colorblocks_compare"
        )
    )

    parser.add_argument(
        "--output-regions", "-out", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output_reg",
        required=True
    )

    args = parser.parse_args()

    return args


def plainify_label(label):
    """
    Convert a label to a plain version by removing special characters and converting to lowercase.
    """
    if LABEL_ORDER_PREFIX.match(label):
        label = label.split("_", 1)[-1]  # Remove order prefix
    plainifed = re.findall(r"[a-z]", label.lower())
    return "".join(sorted(plainifed))


def check_label_identity(label_a, label_b):

    if label_a == label_b:
        return True
    # more involved - ignore special characters
    label_a_plain = plainify_label(label_a)
    label_b_plain = plainify_label(label_b)
    return label_a_plain == label_b_plain


def collapse_label_groups(sub):

    # outer_boundary
    outer_start = sub["start"].min()
    outer_end = sub["end"].max()

    # inner boundary
    inner_start = sub["start"].max()
    inner_end = sub["end"].min()

    assigned_is_unique = sub.loc[sub["source"] == "uniq_kmer", "assigned_label"].nunique()
    assert assigned_is_unique, "Expected exactly one assigned label from kmer source"
    assigned_label = sub.loc[sub["source"] == "uniq_kmer", "assigned_label"].iloc[0]
    if assigned_label == "uncertain":
        # can happen for some k-mer labels
        kmer_strand = "+"
        top_enrich = 0
        score = 0
    else:
        kmer_strand = sub.loc[(sub["source"] == "kmer_strand"), "strand"].iloc[0]
        top_enrich = sub.loc[(sub["source"] == "uniq_kmer") & (sub["kmerColor"] == assigned_label), "top_enrich"].iloc[0]
        score = 1000

    if "aln" in sub["source"].values:
        # do we have a full match?
        aln_labels = sub.loc[sub["source"] == "aln", "name"].str.contains(assigned_label, regex=False)
        if aln_labels.any():
            aln_support = "match"
            aln_orientation = sub.loc[aln_labels.index[aln_labels], "strand"].iloc[0]
        else:
            aln_support = "approx"
            aln_orientation = "."
    else:
        aln_support = "none"
        aln_orientation = "."

    seq = sub["seq"].iloc[0]

    collapsed_region = (
        seq, outer_start, outer_end,
        assigned_label, score, kmer_strand,
        inner_start, inner_end,
        "kmer", top_enrich,
        aln_support, aln_orientation,
        sub["Cluster"].iloc[0]
    )

    return collapsed_region


def extract_second_best_guess(asm_cutout_window):

    if asm_cutout_window == "none":
        return "none"
    else:
        # looks like this
        # 01n_PAR1::NA21093_chrY:757-148195
        parts = asm_cutout_window.split("::")[0]
        label = parts.split("_", 1)[-1]  # remove order prefix
        return label


def cluster_kmer_annotations(concat):
    """This function encapsulates the merge operation
    that combines the two different k-mer annotation inputs:
    one giving the strand/orientation of the repeat block
    and the other giving some proxy for confidence in the
    assignment by stating the number of unique k-mers in
    the respective region.
    These two annotations tend to exhibit some wiggle room
    (off-by-one errors in the coordinates?), which is implicitly
    corrected for here by clustering with a slack of 3 bp.

    This function returns the merged k-mer regions concatenated
    to the other (alignment-derived) regions.
    """

    intervals = pr.from_dict(
        {
            "Chromosome": concat["seq"].values,
            "Start": concat["start"].values,
            "End": concat["end"].values,
            "Strand": concat["strand"].values,
            "kmer_label": concat["is_kmer_label"].values,
            "plain_label": concat["plain_label"].values,
            "pd_idx": concat.index.values,
        }
    )
    # note here: slack -3 to avoid overlaps that look like off-by-one errors
    # in the various k-mer annotations
    intervals = intervals.cluster(by=["kmer_label", "plain_label"], count=False, strand=False, slack=-3).df
    intervals.set_index("pd_idx", inplace=True)
    intervals.drop(["Chromosome", "Start", "End", "Strand", "kmer_label", "plain_label"], axis=1, inplace=True)
    # join the cluster info back into concat
    concat = concat.join(intervals, how="outer", rsuffix="_cluster")

    kmer_regions = []
    collapsed_indices = set()
    for cluster_id, group in concat.groupby("Cluster"):
        if group.shape[0] == 1:
            # single label, nothing to do
            continue
        if group["source"].nunique() == 1:
            # not a mix of kmer and aln labels, nothing to do
            assert group["source"].iloc[0] == "aln"
            continue
        # not recorded - is this helpful/informative?
        # source_mix = group["source"].value_counts()
        sub = group.sort_values(["start", "end"], inplace=False)
        collapsed_region = collapse_label_groups(sub)
        assert collapsed_region[-1] == cluster_id
        kmer_regions.append(collapsed_region)
        collapsed_indices.update(set(sub.index))

    kmer_regions = pd.DataFrame.from_records(
        kmer_regions,
        columns=[
            "seq", "start", "end", "name", "score", "strand",
            "thickStart", "thickEnd",
            "assign_method",
            "kmer_top_enrich",
            "other_support", "other_orientation", "cluster_id"
        ]
    )
    kmer_regions.sort_values(by=["seq", "start", "end"], inplace=True)
    kmer_regions["asm_cutout"] = "none"

    # what is left in concat now must have been assigned by the
    # alignment-based labelings
    concat = concat.loc[~concat.index.isin(collapsed_indices), :].copy()

    # add columns that are used in the merge operation
    # of the k-mer regions above
    concat["assign_method"] = "aln"
    concat["kmer_top_enrich"] = concat["top_enrich"]
    concat["cluster_id"] = concat["Cluster"]
    concat["other_support"] = "none"
    concat["other_orientation"] = "."
    concat["thickStart"] = concat["start"]
    concat["thickEnd"] = concat["end"]
    concat.drop(
        [
            "Cluster", "source", "is_kmer_label", "plain_label", "top_enrich",
            "ProportionOfTotalKmers", "TotalKmersOfRegion", "totalUniqueKmerHits",
            "kmerColor", "assigned_label"
        ], axis=1, inplace=True
    )

    merged_regions = pd.concat([concat, kmer_regions], axis=0, ignore_index=False)
    merged_regions.sort_values(by=["seq", "start", "end"], inplace=True)
    merged_regions.reset_index(drop=True, inplace=True)

    return merged_regions


def main():

    args = parse_command_line()

    # nothing special: read the alignment-based labeling
    print(f"Processing alignment file: {args.input_aln.name}")
    aln_labels = pd.read_csv(
        args.input_aln,
        sep="\t", header=0
    )
    aln_labels.rename({"#chrom": "seq"}, axis=1, inplace=True)
    # in prep for pyranges clustering
    aln_labels["strand"] = aln_labels["strand"].replace({1: "+", -1: "-"}, inplace=False).astype(str)
    aln_labels["plain_label"] = aln_labels["name"].apply(plainify_label)
    aln_labels["source"] = "aln"

    # collect labels also assigned via k-mer approach
    # NB: not all labels / Y sequence classes are derived
    # with the alternative approaches, e.g., PAR1 or XDR/XTR
    # regions do not exist in these annotations. Only repeat
    # (color) blocks such as blue, green, red etc.
    kmer_based_labels = set()

    # example of the input read here (after normalization - *chrY-regions_uniq-kmer-counts.tsv):
    # #Contig Start   End     kmerColor       totalUniqueKmerHits     TotalKmersOfRegion      ProportionOfTotalKmers  assigned_label  top_enrich
    # NA21093_chrY    22871428        23039139        blue1   14334   167612  0.0855189365916521      blue1   11.95
    print(f"Processing k-mer label/block file: {args.input_kmer.name}")
    kmer_labels = pd.read_csv(
        args.input_kmer,
        sep="\t", header=0
    )
    kmer_labels.rename({"#Contig": "seq", "Start": "start", "End": "end"}, axis=1, inplace=True)
    kmer_labels["plain_label"] = kmer_labels["kmerColor"].apply(plainify_label)
    kmer_labels["name"] = kmer_labels["kmerColor"]
    kmer_labels["asm_cutout"] = "n/a"
    kmer_labels["source"] = "uniq_kmer"
    kmer_labels["strand"] = "+"
    kmer_labels["score"] = 1000

    kmer_based_labels.update(set(kmer_labels["plain_label"].values))

    # example of the input read here (after normalization - *chrY-colorblock-kmers.bed):
    # #seq_name       start   end     name    score   strand
    # NA21093_chrY    22871429        23039136        blue1   1000    -
    # NA21093_chrY    23039135        23154782        teal2   1000    -
    print(f"Processing k-mer strand/compare file: {args.input_kmer_strands.name}")
    kmer_strands = pd.read_csv(
        args.input_kmer_strands,
        sep="\t", header=0
    )
    kmer_strands.rename({"#seq_name": "seq"}, axis=1, inplace=True)
    kmer_strands["plain_label"] = kmer_strands["name"].apply(plainify_label)
    kmer_strands["source"] = "kmer_strand"
    kmer_strands["asm_cutout"] = "n/a"

    kmer_based_labels.update(set(kmer_strands["plain_label"].values))

    concat = pd.concat([aln_labels, kmer_labels, kmer_strands], axis=0, ignore_index=False)
    concat.sort_values(by=["seq", "start", "end"], inplace=True)
    concat.reset_index(drop=True, inplace=True)

    # fill disjoint fields / columns with reasonable defaults
    concat["top_enrich"] = concat["top_enrich"].fillna(0, inplace=False).astype(float)
    concat["assigned_label"] = concat["assigned_label"].fillna("n/a", inplace=False).astype(str)
    concat["ProportionOfTotalKmers"] = concat["ProportionOfTotalKmers"].fillna(0, inplace=False).astype(float)
    concat["TotalKmersOfRegion"] = concat["TotalKmersOfRegion"].fillna(0, inplace=False).astype(int)
    concat["totalUniqueKmerHits"] = concat["totalUniqueKmerHits"].fillna(0, inplace=False).astype(int)
    concat["kmerColor"] = concat["kmerColor"].fillna("n/a", inplace=False).astype(str)

    concat["is_kmer_label"] = concat["plain_label"].apply(lambda n: n in kmer_based_labels or n == UNCERTAIN_SCRAMBLED)

    merged_regions = cluster_kmer_annotations(concat)

    # fallback option - this label just assumes that the alignment-based
    # label is correct and extracts that info from the FASTA header of the
    # sequence window that was cut out of the assembly
    merged_regions["second_best_guess"] = merged_regions["asm_cutout"].apply(
        extract_second_best_guess
    )
    merged_regions.drop("asm_cutout", axis=1, inplace=True)

    merged_regions["name"] = merged_regions["name"].str.replace(
        LABEL_ORDER_PREFIX, "", regex=True
    )

    args.output_reg.parent.mkdir(parents=True, exist_ok=True)
    merged_regions = merged_regions[
        [
            "seq", "start", "end", "name", "score", "strand",
            "thickStart", "thickEnd",
            "assign_method", "second_best_guess",
            "kmer_top_enrich", "other_support", "other_orientation",
            "cluster_id"
        ]
    ]
    with open(args.output_reg, "w") as out_fh:
        out_fh.write("#")
        merged_regions.to_csv(out_fh, sep="\t", index=False, header=True)

    return 0


if __name__ == "__main__":
    main()
