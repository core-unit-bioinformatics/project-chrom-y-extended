#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl
import re
import sys

import pandas as pd
import pyranges as pr


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input-aln", "-aln", "-a",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_aln",
        required=False,
        default="NA20905.hg38.chrY-regions-realigned.norm-paf.tsv.gz"
    )

    hg38 = "GRCh38-Y_regions_repeat-details.bed"
    t2t = "HG002-Y_T2Tv2_regions_repeat-details.bed"

    parser.add_argument(
        "--input-regions", "-reg", "-r",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_reg",
        required=False,
        default=hg38
    )

    parser.add_argument(
        "--output-regions", "-out", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output_reg",
        required=False,
        default="out.NA20905.hg38.bed"
    )

    args = parser.parse_args()

    return args


def join_region_labels(alignments, regions):
    """This function joins (= overlaps) the input
    regions, i.e. the regions of interest labeled
    in the reference, with the resulting alignment
    records of the inferred labeled regions in the
    de novo assemblies (the query in the PAF) aligned
    back to the labeled reference.
    This process is supposed to check if the inferred
    labels in the de novo assembly / query do
    actually align to the correct region as annotated
    in the reference.

    Args:
        alignments (pandas.DataFrame): _description_
        regions (pandas.DataFrame): _description_

    Returns:
        pandas.DataFrame: _description_
    """

    # as usual for Pyranges interval operations, carry
    # additional information from the PAF
    # The 'query' below refers to the regions in the
    # de novo assemblies with the inferred label.
    iv_align = pr.from_dict(
        {
            "Chromosome": alignments["target_name"].values,
            "Start": alignments["target_start"].values,
            "End": alignments["target_end"].values,
            "pd_idx_aln": alignments.index.values,
            "asm_seq_name": alignments["query_name"].values,
            "asm_seq_start": alignments["query_start"].values,
            "asm_seq_end": alignments["query_end"].values,
            "aln_strand": alignments["align_orient"].values,
            "mapq": alignments["mapq"].values
        }
    )

    # this is just the reference annotation
    # set of regions (Y seq. classes)
    iv_regions = pr.from_dict(
        {
            "Chromosome": regions["chrom"].values,
            "Start": regions["start"].values,
            "End": regions["end"].values,
            "Name": regions["name"].values
        }
    )

    joined = iv_align.join(iv_regions, how="outer", report_overlap=True, suffix="_region")

    return joined.df


def compute_asm_seq_length(fasta_header):
    """compute_asm_seq_length _summary_

    Args:
        fasta_header (str): fasta seq. header like this "01n_PAR1::HG02040_chrY:0-40573"

    Raises:
        err: invalid fasta header

    Returns:
        int: length of the seq., e.g., (40573 - 0) = 40753
    """
    try:
        seq_coord = fasta_header.split(":")[-1]
        start, end = seq_coord.split("-")
        length = int(end) - int(start)
    except Exception as err:
        err.add_note("Cannot process FASTA header: {}".format(fasta_header))
        raise err
    return length


def add_asm_seq_offset(fasta_header):
    """add_asm_seq_offset _summary_

    Args:
        fasta_header (str): fasta seq. header like this "01n_PAR1::HG02040_chrY:0-40573"

    Returns:
        int: offset, i.e., the start position of the sequence fragment cut out from
            the entire assembled sequence after inferring the window of the
            labeled region, e.g. HG02040_chrY:0-40573 = 0
    """

    seq_coord = fasta_header.split(":")[-1]
    start, _ = seq_coord.split("-")
    offset = int(start)
    return offset


def add_asm_seq_start(fasta_header):
    return add_asm_seq_offset(fasta_header)


def add_asm_seq_end(fasta_header):
    """add_asm_seq_end _summary_

    Args:
        fasta_header (str): fasta seq. header like this "01n_PAR1::HG02040_chrY:0-40573"

    Returns:
        int: end, i.e., HG02040_chrY:0-40573 = 40573
    """

    seq_coord = fasta_header.split(":")[-1]
    _, end = seq_coord.split("-")
    end = int(end)
    return end


def add_asm_region_label(fasta_header):
    """add_asm_region_label _summary_

    Args:
        fasta_header (str): fasta seq. header like this "01n_PAR1::HG02040_chrY:0-40573"

    Returns:
        str: region label, e.g., PAR1
    """

    label = fasta_header.split("::")[0].split("_", 1)[-1]
    return label


def add_asm_seq(fasta_header):
    """add_asm_seq _summary_

    Args:
        fasta_header (str): fasta seq. header like this "01n_PAR1::HG02040_chrY:0-40573"

    Returns:
        str: sequence name, e.g., HG02040_chrY
    """

    seq_name = fasta_header.split(":")[2]
    return seq_name


def add_asm_seq_length(joined):

    joined["asm_seq_length"] = joined["asm_seq_name"].apply(compute_asm_seq_length)
    # note that the overlap_pct can be larger than 100
    # due to breaks/skips in the alignment, but will typically
    # exceed 100 by only a few points
    joined["overlap_pct"] = (joined["Overlap"] / joined["asm_seq_length"] * 100).round(2)

    return joined


def check_same_color_family(region_label, matched_names):

    possible_colors = ["green", "blue", "yellow", "red", "teal", "gray"]
    region_color = None
    for c in possible_colors:
        if c not in region_label:
            continue
        assert region_color is None, "Multiple colors found in region label"
        region_color = c
    if region_color is None:
        # that region is simply not a 'color region'
        return []

    family_match = []
    for idx, name in matched_names.items():
        if region_color not in name:
            continue
        family_match.append(idx)
    return family_match


def merge_aligned_regions(joined):

    out_regions = []
    for seq_name, alignments in joined.groupby("asm_seq_name"):
        # this iterates over aggregates like this:
        # 01n_PAR1::HG02040_chrY:0-40573
        # i.e., it's labeled sequence fragments
        offset = alignments["asm_seq_offset"].iloc[0]

        align_target_is_unique_name = alignments["Name"].nunique() == 1
        align_target_is_not_unique_name = not align_target_is_unique_name
        query_name_in_target = alignments["asm_region_label"].iloc[0] in alignments["Name"].values
        query_name_not_in_target = not query_name_in_target

        if align_target_is_not_unique_name or query_name_not_in_target:
            if query_name_not_in_target:
                # no match at all
                # two options: if it's a 'color family' match, it's likely just
                # an alignment artifact (alignment-based matching is not sensitive enough)
                # otherwise, no clue what the problem is / where the mismatch comes from
                region_label_name = alignments["asm_region_label"].iloc[0]
                matched_names = alignments["Name"]
                matches_by_family = check_same_color_family(region_label_name, matched_names)
                if matches_by_family:
                    # use only 'color' matches, e.g., accept that the aligner
                    # matches green2 to green3 etc.
                    subset = alignments.loc[matches_by_family, :].copy()
                    orientation = subset.groupby("aln_strand")["Overlap"].sum()
                    orientation = orientation.index[orientation.argmax()]
                    region = {
                        "seq": subset["asm_seq"].iloc[0],
                        "start": subset["asm_seq_start"].min() + offset,
                        "end": subset["asm_seq_end"].max() + offset,
                        "name": region_label_name, # label
                        "score": 500,
                        "strand": orientation,
                        "asm_cutout": subset["asm_seq_name"].iloc[0]
                    }
                else:
                    # uncertain - might be an alignment artifact
                    orientation = alignments.groupby("aln_strand")["Overlap"].sum()
                    orientation = orientation.index[orientation.argmax()]
                    region = {
                        "seq": alignments["asm_seq"].iloc[0],
                        "start": alignments["asm_seq_start"].min() + offset,
                        "end": alignments["asm_seq_end"].max() + offset,
                        "name": "uncertain", # label
                        "score": 0,
                        "strand": orientation,
                        "asm_cutout": alignments["asm_seq_name"].iloc[0]
                    }
            else:
                # only look at the subset with matching label;
                # this is a heuristic that ignores alignment artifacts
                # and effectively changes/reduces the size of the matched
                # sequence part
                selector = alignments["Name"] == alignments["asm_region_label"]
                subset = alignments.loc[selector, :].copy()
                orientation = subset.groupby("aln_strand")["Overlap"].sum()
                orientation = orientation.index[orientation.argmax()]
                region = {
                    "seq": subset["asm_seq"].iloc[0],
                    "start": subset["asm_seq_start"].min() + offset,
                    "end": subset["asm_seq_end"].max() + offset,
                    "name": subset["Name"].iloc[0], # label
                    "score": 750,
                    "strand": orientation,
                    "asm_cutout": subset["asm_seq_name"].iloc[0]
                }
        else:
            # easy case, simple 1:1 label mapping
            orientation = alignments.groupby("aln_strand")["Overlap"].sum()
            orientation = orientation.index[orientation.argmax()]
            region = {
                "seq": alignments["asm_seq"].iloc[0],
                "start": alignments["asm_seq_start"].min() + offset,
                "end": alignments["asm_seq_end"].max() + offset,
                "name": alignments["Name"].iloc[0], # label
                "score": 1000,
                "strand": orientation,
                "asm_cutout": alignments["asm_seq_name"].iloc[0]
            }
        out_regions.append(region)

    df = pd.DataFrame.from_records(out_regions)
    df.sort_values(["seq", "start"], inplace=True)

    return df


def get_un_labels(regions):

    match_un_label = re.compile("^[0-9]+[un]{1}")

    un_labels = []
    for row in regions.itertuples():
        if not match_un_label.match(row.seqclass):
            continue
        un_labels.append(row.seqclass)
        un_labels.append(row.name)

    return un_labels


def discard_enclosed_un_alignments(align, un_labels):

    global_labels = col.defaultdict(list)

    for row in align.itertuples():
        if not any(label in row.query_name for label in un_labels):
            continue
        idx = row.Index
        seq_name = add_asm_seq(row.query_name)
        seq_start = add_asm_seq_start(row.query_name)
        seq_end = add_asm_seq_end(row.query_name)
        seq_label = add_asm_region_label(row.query_name)
        seq_length = seq_end - seq_start
        global_labels["Chromosome"].append(seq_name)
        global_labels["Start"].append(seq_start)
        global_labels["End"].append(seq_end)
        global_labels["Name"].append(seq_label)
        global_labels["length"].append(seq_length)
        global_labels["pd_idx"].append(idx)

    iv = pr.from_dict(global_labels)
    iv = iv.cluster(slack=-1, count=True).df
    drop_alignments = col.defaultdict(set)
    for cluster, aln in iv.groupby("Cluster"):
        # annoying that PyRanges seems not to be able
        # to cluster and compute the fractional overlap
        # at the same time... anyway...
        for a in aln.itertuples():
            for b in aln.itertuples():
                if a.Name == b.Name:
                    continue
                ovl = min(a.End, b.End) - max(a.Start, b.Start)
                dropped = 0
                frac_a = ovl/a.length
                frac_b = ovl/b.length
                #if frac_a >= 1:  # deactivated for heuristic below
                if frac_a > 0.999:
                    dropped += 1
                    drop_alignments[a.Name].add(a.pd_idx)
                #if frac_b >= 1:  # deactivate for heuristic below
                if frac_b > 0.999:
                    dropped += 1
                    drop_alignments[b.Name].add(b.pd_idx)
                if dropped > 1:
                    raise ValueError(f"double dropping: {a} - {b}")

    # Heuristic: check if the above would lead to a u/n label being
    # removed entirely for this sample, which we do not want.
    # In this case, retain the single alignment representing the
    # largest match - heuristic of last resort...
    for label, label_alns in drop_alignments.items():
        # check if dropping those alignments would lead to
        # removing an un label entirely
        check_drop = align.loc[~align.index.isin(label_alns), align.columns[:10]]
        if check_drop["query_name"].str.contains(label).any():
            # label still contained, we can drop those alignments
            align.drop(label_alns, inplace=True)
        else:
            # dropping all those alignments would remove an 'un'
            # label from the sample, which we do not want
            # -> keep the alignment w/ most matching bases
            keep_idx = align.loc[sorted(label_alns), "align_matching"].idxmax()
            drop_others = label_alns - set([keep_idx])
            align.drop(drop_others, inplace=True)

    align = align.reset_index(drop=True, inplace=False)

    return align


def check_missing_un_labels(sample_regions, un_labels):

    sub_select = sample_regions["name"] == "uncertain"
    fixes = []
    for row in sample_regions.loc[sub_select, :].itertuples():
        if not any(label in row.asm_cutout for label in un_labels):
            continue
        label = add_asm_region_label(row.asm_cutout)
        if label in sample_regions["name"].values:
            # label already exists, do not foce another
            # one into the annotation
            continue
        fixes.append((row.Index, label))

    for idx, label in fixes:
        sample_regions.loc[idx, "name"] = label
        sample_regions.loc[idx, "score"] = 250

    return sample_regions


def main():

    args = parse_command_line()

    regions = pd.read_csv(args.input_reg, sep="\t", header=0)
    # to simplify matching, we use the sequence class label as
    # region name, i.e., switch from 01n_PAR1 to just PAR1
    regions[["name", "seqclass"]] = regions[["seqclass", "name"]]
    regions.rename({"#chrom": "chrom"}, axis=1, inplace=True)

    align = pd.read_csv(args.input_aln, sep="\t", header=0)

    # 2026-02-06 // add heuristics
    # add external knowledge to the process and filter as follows:
    # (1) umbrella (u) or default (n) seq. class labels that have
    # alignments that are fully enclosed by other u or n seq. class
    # label alignments are discarded (the respective alignments).
    # (2) [final step below]: "uncertain" regions that represent
    # an umbrella (u) or default (n) label are added back into
    # the region set; this fixes the missing AMPL1 label in
    # the T2Tv2-based labeling because of a reference bias;
    # In the T2Tv2 labeling, the AMPL1_IR3 (hg38) is split into
    # AMPL1 and IR3 and the AMPL1 candidate alignments
    # (example: sample HG04228) have their best matches to
    # AMPL2 (end) and IR3 and there is thus no label match
    # from AMPL1 to AMPL1

    # apply heuristic number 1: discard fully enclosed
    # umbrella/global seq. class alignments; this is NOT
    # done for the "subsumed" (s) seq. class terms.
    un_labels = get_un_labels(regions)

    align = discard_enclosed_un_alignments(align, un_labels)

    joined = join_region_labels(align, regions)

    # if any region (labels) do not exist in the assembly,
    # they will be dropped here with only an error message
    missing_labels = joined["asm_seq_name"] == "-1"
    if missing_labels.any():
        # note here: at this stage, it is not possible with
        # certainty to say whether or not the missing label
        # is an actual assembly error or just a recalcitrant
        # alignment artifact
        label_names = sorted(joined.loc[missing_labels, "Name"].unique())
        sys.stderr.write(
            f"\nWarning: the following region labels are missing in the assembly: {label_names}\n"
        )
        joined = joined.loc[~missing_labels, :].copy()
    joined = add_asm_seq_length(joined)
    joined["asm_region_label"] = joined["asm_seq_name"].apply(add_asm_region_label)
    joined["asm_seq"] = joined["asm_seq_name"].apply(add_asm_seq)
    joined["asm_seq_offset"] = joined["asm_seq_name"].apply(add_asm_seq_offset)

    df = merge_aligned_regions(joined)

    # apply second heuristic
    df = check_missing_un_labels(df, un_labels)

    df.rename({"seq": "#chrom"}, axis=1, inplace=True)

    args.output_reg.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(args.output_reg, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
