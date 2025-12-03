#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import pandas as pd
import pyranges as pyr


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--label-file", "-l",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="label_file"
    )

    parser.add_argument(
        "--centromere", "-c",
        type=str,
        dest="centromere"
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output"
    )

    args = parser.parse_args()

    return args


def main():

    args = parse_command_line()

    seq_classes = pd.read_csv(args.label_file, sep="\t", header=0)

    CEN_LABELS = ["CEN", "CEN1", "CEN2"]

    cen_specs = seq_classes.loc[seq_classes["name"].isin(CEN_LABELS), :].copy()
    cen_specs["name"] = "periCEN"
    cen_specs["second_best_guess"] = "periCEN"
    column_sort_order = cen_specs.columns

    # kept for re-merging
    seq_classes.drop(cen_specs.index, axis=0, inplace=True)

    # parse cen coordinate string
    # HG01243_chrY:10359455-10679866
    # NA19443_chrY:10313852-11285628

    seq, window = args.centromere.split(":")
    start, end = window.split("-")

    cen_window = pyr.from_dict(
        {
            "Chromosome": [seq],
            "Start": [int(start)],
            "End": [int(end)]
        }
    )

    cen_regions = pyr.from_dict(
        {
            "Chromosome": cen_specs["#seq"].values,
            "Start": cen_specs["start"].values,
            "End": cen_specs["end"].values,
            "pd_idx": cen_specs.index
        }
    )

    cen_regions = cen_regions.subtract(cen_window).df.rename(
        {"Start": "start", "End": "end", "Chromosome": "#seq"},
        axis=1, inplace=False
    )
    if not cen_regions.empty:
        select_cen_to_keep = cen_specs.index.isin(cen_regions["pd_idx"])
        # some initial CEN regions might have been replaced/subsumed by
        # the external CEN annotation as a whole, reduce original
        # set to what is left over
        # NB: because cen_regions was not empty, there must be
        # something left here
        cen_specs = cen_specs.loc[select_cen_to_keep, :].copy()
        cen_specs.drop(["start", "end", "#seq"], axis=1, inplace=True)

        cen_regions.set_index("pd_idx", inplace=True)
        cen_specs = cen_specs.join(cen_regions)
        assert cen_specs.shape[0] == cen_regions.shape[0]

        cen_specs = cen_specs[column_sort_order]

    # now add the external CEN window
    cen_window = cen_window.df.rename(
        {"Start": "start", "End": "end", "Chromosome": "#seq"},
        axis=1, inplace=False
    )
    cen_window["name"] = "CEN"
    cen_window["score"] = 1000
    cen_window["strand"] = "+"
    cen_window["thickStart"] = cen_window["start"]
    cen_window["thickEnd"] = cen_window["end"]
    cen_window["assign_method"] = "exp"
    cen_window["second_best_guess"] = cen_window["name"]
    cen_window["kmer_top_enrich"] = 0.
    cen_window["other_support"] = "none"
    cen_window["other_orientation"] = "."
    cen_window["cluster_id"] = -1

    if not cen_regions.empty:
        cen_specs = pd.concat([cen_specs, cen_window], axis=0, ignore_index=False)
    else:
        # empty case: what's left to add is just the CEN window from
        # the external expert annotation
        cen_specs = cen_window

    seq_classes = pd.concat([seq_classes, cen_specs], axis=0, ignore_index=False)
    seq_classes.sort_values(["#seq", "start", "end"], inplace=True)
    seq_classes.reset_index(drop=True, inplace=True)
    assert not pd.isnull(seq_classes).any(axis=1).any()

    args.output.parent.mkdir(exist_ok=True, parents=True)
    seq_classes.to_csv(args.output, sep="\t", header=True, index=False)

    return 0




if __name__ == "__main__":
    main()
