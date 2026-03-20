#!/usr/bin/env python3

import argparse as argp
import collections as col
import json
import pathlib as pl

import pandas as pd


ISECT_COL_NAMES = list(
    (pos, name) for pos, name in
    enumerate([
        "seq", "start", "end", "label",  # 0-3
        "score", "strand",
        "src_start", "src_end",
        "seqtype", "group_label", # 8-9
        "seq2", "start2", "end2", "name2",  # 13
        "score2", "strand2",
        "overlap_bp" # 16
    ], start=0)
)

ISECT_HEADER = [hd for _, hd in ISECT_COL_NAMES]
SELECT_COLS = [0,1,2,3,8,9,13,16]
USE_HEADER = [hd for pos,hd in ISECT_COL_NAMES if pos in SELECT_COLS]

OVERLAP_ISSUES = ["ERRBASE", "ERRSTRUCT", "UNASSIGNED"]
ALL_ISSUES = OVERLAP_ISSUES + ["NGAP"]


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "-i", "--input",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input"
    )
    parser.add_argument(
        "-u", "--umbrella",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="umbrella"
    )
    parser.add_argument(
        "-o", "--output",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output"
    )

    args = parser.parse_args()

    return args


def load_isect_table(isect_table, umbrella_labels):

    df = pd.read_csv(isect_table, sep="\t", header=None, names=ISECT_HEADER, usecols=USE_HEADER)
    df = df.loc[~df["label"].isin(OVERLAP_ISSUES), :].copy()
    df["length"] = df["end"] - df["start"]
    df = df.sort_values(["seq", "start"], inplace=False).reset_index(inplace=False, drop=True)

    with open(umbrella_labels) as dump:
        umbrellas = json.load(dump)


    def unify_label_name(name, umbrellas):
        try:
            uni_label = umbrellas[name]["unified"]
            if "," in uni_label:
                uni_label = name
            if uni_label.startswith("other"):
                uni_label = "OTHER"
        except KeyError:
            uni_label = name
        return uni_label

    df["label"] = df["label"].apply(unify_label_name, args=(umbrellas,))

    return df, umbrellas


def assign_ngaps(isect, seqtype_groups, grouping_labels):

    ngaps = col.Counter()
    for st_group, st_label in zip(seqtype_groups, grouping_labels):
        sub = isect.loc[isect["seqtype"].isin(st_group), :].copy()

        ngap_indices = sub.loc[sub["label"] == "NGAP", :].index
        if ngap_indices.empty:
            continue

        for ngap in ngap_indices:
            gap_row = sub.loc[ngap, :]
            gap_size = gap_row.end - gap_row.start

            row_before = sub.loc[ngap-1, :]
            row_after = sub.loc[ngap+1, :]
            # record - ?
            if row_before.label == row_after.label:
                assigned = "enclosed"
            else:
                assigned = "left"
            ngaps[(st_label, row_before.label, "ovl_num")] += 1
            ngaps[(st_label, row_before.label, "ovl_bp")] += gap_size
    return ngaps


def main():

    args = parse_command_line()
    isect, umbrellas = load_isect_table(args.input, args.umbrella)

    sample, ref = args.input.name.split(".")[:2]

    seqtype_groups = [("main",), ("main", "rand")]
    grouping_labels = ["main", "all"]

    ngaps = assign_ngaps(isect, seqtype_groups, grouping_labels)

    merged = []
    for st_group, st_label in zip(seqtype_groups, grouping_labels):
        sub = isect.loc[isect["seqtype"].isin(st_group), :].copy()

        region_lengths = sub.groupby("label")["length"].unique().apply(sum)

        agg_ovl = sub.groupby(["label", "name2"])["overlap_bp"].sum()
        agg_num = sub.groupby(["label", "name2"])["overlap_bp"].size()

        for label, length in region_lengths.items():
            if label == "NGAP":
                continue
            base_key = [ref, sample, st_label, label]
            merged.append(
                (tuple(base_key + ["SIZE", "bp"]), length)
            )

            ngap_num = ngaps[(st_label, label, "ovl_num")]
            ngap_bp = ngaps[(st_label, label, "ovl_bp")]
            merged.append(
                (tuple(base_key + ["NGAP", "bp"]), ngap_bp)
            )
            merged.append(
                (tuple(base_key + ["NGAP", "num"]), ngap_num)
            )
            pct = round(ngap_bp / length * 100, 4)
            merged.append(
                (tuple(base_key + ["NGAP", "pct"]), pct)
            )
            for issue in OVERLAP_ISSUES:
                try:
                    ovl_num = agg_num.xs((label, issue))
                    ovl_bp = agg_ovl.xs((label, issue))
                    merged.append(
                        (tuple(base_key + [issue, "bp"]), ovl_bp)
                    )
                    merged.append(
                        (tuple(base_key + [issue, "num"]), ovl_num)
                    )
                    pct = round(ovl_bp / length * 100, 4)
                    merged.append(
                        (tuple(base_key + [issue, "pct"]), pct)
                    )
                except KeyError:
                    merged.append(
                        (tuple(base_key + [issue, "bp"]), 0)
                    )
                    merged.append(
                        (tuple(base_key + [issue, "num"]), 0)
                    )
                    merged.append(
                        (tuple(base_key + [issue, "pct"]), 0)
                    )

    midx = pd.MultiIndex.from_tuples([t[0] for t in merged], names=["ref", "sample", "seqtype", "label", "other", "statistic"])
    merged = pd.DataFrame.from_records(merged, index=midx, columns=["tidx", "value"])
    merged.drop("tidx", axis=1, inplace=True)

    # sanity checking
    data_labels = set(merged.index.get_level_values("label"))
    umbrella_labels = set(d["unified"] for d in umbrellas.values())
    umbrella_labels.add("TELO")
    unknowns = data_labels - umbrella_labels
    assert len(unknowns) == 0, unknowns

    args.output.parent.mkdir(exist_ok=True, parents=True)
    merged.to_csv(args.output, sep="\t", index=True, header=True)

    return 0


if __name__ == "__main__":
    main()

