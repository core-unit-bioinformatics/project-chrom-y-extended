#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import pandas as pd
import pyranges as pr


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--bed-in", "-i",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_bed",
        required=True
    )

    parser.add_argument(
        "--merge-dist", "-d",
        type=int,
        default=1,
        dest="merge_dist"
    )

    parser.add_argument(
        "--bed-out", "-o",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="output_bed",
        required=True
    )

    args = parser.parse_args()

    return args


def simplify_label(label):

    prefix = label.split("[")[0]
    if prefix.startswith("chrY_"):
        new_label = prefix.split("_", 2)[-1]
    else:
        new_label = prefix
    return new_label


def main():

    args = parse_command_line()

    bed_regions = pd.read_csv(
        args.input_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "long_name"]
    )
    bed_regions["name"] = bed_regions["long_name"].apply(simplify_label)

    intervals = pr.from_dict(
        {
            "Chromosome": bed_regions["chrom"],
            "Start": bed_regions["start"],
            "End": bed_regions["end"],
            "Name": bed_regions["name"],
            "pd_idx": bed_regions.index.values
        }
    )
    intervals = intervals.merge(by="Name", slack=args.merge_dist, count=True, count_col="n_merged")

    intervals = intervals.df
    intervals.rename(
        dict((col, new_col) for col, new_col in
        zip(intervals.columns, ["#chrom", "start", "end", "name", "n_merged"])),
        axis=1, inplace=True
    )
    intervals.sort_values(["#chrom", "start"], inplace=True)
    intervals.to_csv(args.output_bed, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
