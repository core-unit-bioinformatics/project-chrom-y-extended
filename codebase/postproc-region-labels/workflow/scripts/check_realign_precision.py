#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import pandas as pd
import pyranges as pr


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input-aln", "-aln", "-a",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_aln",
        required=True
    )
    parser.add_argument(
        "--input-regions", "-reg", "-r",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_reg",
        required=True
    )

    parser.add_argument(
        "--output-regions", "-out", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output_reg",
        required=True
    )

    args = parser.parse_args()

    return args


def join_region_labels(alignments, regions):

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

    try:
        seq_coord = fasta_header.split(":")[-1]
        start, end = seq_coord.split("-")
        length = int(end) - int(start)
    except Exception as err:
        err.add_note("Cannot process FASTA header: {}".format(fasta_header))
        raise err
    return length


def add_asm_seq_offset(fasta_header):

    seq_coord = fasta_header.split(":")[-1]
    start, _ = seq_coord.split("-")
    offset = int(start)
    return offset


def add_asm_region_label(fasta_header):

    label = fasta_header.split("::")[0]
    return label


def add_asm_seq(fasta_header):

    seq_name = fasta_header.split(":")[2]
    return seq_name


def add_asm_seq_length(joined):

    joined["asm_seq_length"] = joined["asm_seq_name"].apply(compute_asm_seq_length)
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
        offset = alignments["asm_seq_offset"].iloc[0]
        if alignments["Name"].nunique() != 1:
            if alignments["asm_region_label"].iloc[0] not in alignments["Name"].values:
                # two options: if it's a 'color family' match, it's likely just
                # an alignment artifact (alignment-based matching is not sensitive enough)
                # otherwise, no clue what the problem is / where the mismatch comes from
                region_label_name = alignments["asm_region_label"].iloc[0]
                matched_names = alignments["Name"]
                matches_by_family = check_same_color_family(region_label_name, matched_names)
                if matches_by_family:
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
                selector = alignments["Name"] == alignments["asm_region_label"]
                subset = alignments.loc[selector, :].copy()
                orientation = subset.groupby("aln_strand")["Overlap"].sum()
                orientation = orientation.index[orientation.argmax()]
                region = {
                    "seq": subset["asm_seq"].iloc[0],
                    "start": subset["asm_seq_start"].min() + offset,
                    "end": subset["asm_seq_end"].max() + offset,
                    "name": subset["Name"].iloc[0], # label, # label
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


def main():

    args = parse_command_line()

    regions = pd.read_csv(args.input_reg, sep="\t", header=0)
    regions.rename({"#chrom": "chrom"}, axis=1, inplace=True)

    align = pd.read_csv(args.input_aln, sep="\t", header=0)

    joined = join_region_labels(align, regions)
    joined = add_asm_seq_length(joined)
    joined["asm_region_label"] = joined["asm_seq_name"].apply(add_asm_region_label)
    joined["asm_seq"] = joined["asm_seq_name"].apply(add_asm_seq)
    joined["asm_seq_offset"] = joined["asm_seq_name"].apply(add_asm_seq_offset)

    df = merge_aligned_regions(joined)
    df.rename({"seq": "#chrom"}, axis=1, inplace=True)

    args.output_reg.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(args.output_reg, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
