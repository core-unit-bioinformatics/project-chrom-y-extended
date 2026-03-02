#!/usr/bin/env python3

import argparse as argp
import collections as col
import json
import pathlib as pl

import pandas as pd

DEBUG = False

HEADER = [
    "seq", "start", "end", "draft_label", "score", "strand",
    "seq2", "start2", "end2", "final_label", "score2", "strand2",
    "overlap_bp"
]
BED_HEADER = ["seq", "start", "end", "name", "score", "strand"]
ISSUE_LABELS = ["UNASSIGNED", "ERRBASE", "ERRSTRUCT", "NGAP"]


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "-i", "--intersect",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="intersect",
        required=True,
    )
    parser.add_argument(
        "-l", "--umbrella-labels",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="umbrella_labels",
        required=True,
    )
    parser.add_argument(
        "-o", "--output",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="output",
        required=True
    )
    args = parser.parse_args()

    return args


def load_umbrella_labels(fp):
    umbrella_labels = json.load(open(fp))
    return umbrella_labels


def find_previous(regions, seq, start_idx):
    # index selection relies on sort-order of regions
    select_index = regions.index < start_idx
    select_umbrella = regions["is_umbrella"]
    select_seq = regions["seq2"] == seq
    selector = select_index & select_umbrella & select_seq
    return regions.loc[selector, :].iloc[-1]


def find_next(regions, seq, start_idx):
    # index selection relies on sort-order of regions
    select_index = regions.index > start_idx
    select_umbrella = regions["is_umbrella"]
    select_seq = regions["seq2"] == seq
    selector = select_index & select_umbrella & select_seq
    try:
        next_umbrella = regions.loc[selector, :].iloc[0]
    except IndexError:
        next_umbrella = None
    return next_umbrella


def check_overlap(issue_start, issue_end, umbrella_region):
    if umbrella_region is None:
        return -1
    ovl = min(issue_end, umbrella_region.end2) - max(issue_start, umbrella_region.start2)
    return ovl


def extract_issue_information(overlaps, self=None):

    issues = col.Counter()
    for row in overlaps.itertuples():
        if row.draft_label not in ISSUE_LABELS:
            continue
        issue = row.draft_label
        if self is not None and issue == self:
            continue
        issues[f"{issue}_n"] += 1
        issues[f"{issue}_bp"] += row.overlap_bp
    return issues


def assign_label(regions, umbrellas, overlaps, issue_name):

    if DEBUG:
        print(overlaps)
    if issue_name == "UNASSIGNED":
        # by construction, cannot overlap
        sub = overlaps.loc[~overlaps["final_label"].isin(ISSUE_LABELS), :].copy()
    elif issue_name == "NGAP":
        sub = overlaps.loc[~overlaps["draft_label"].isin(ISSUE_LABELS), :].copy()
    else:
        raise
    if sub.empty:
        min_index = overlaps.index.min()
        max_index = overlaps.index.max()
        seq_name = overlaps["seq2"].iloc[0]
        start = overlaps["start2"].iloc[0]
        end = overlaps["end2"].iloc[0]
        previous_umbrella = find_previous(regions, seq_name, min_index)
        prev_ovl = check_overlap(start, end, previous_umbrella)

        if DEBUG:
            print("prev")
            print(previous_umbrella)
            print(prev_ovl)

        next_umbrella = find_next(regions, seq_name, max_index)
        next_ovl = check_overlap(start, end, next_umbrella)

        if DEBUG:
            print("next")
            print(next_umbrella)
            print(next_ovl)

        if next_umbrella is not None and (previous_umbrella.final_label == next_umbrella.final_label):
            # assign-enclosed
            assign = previous_umbrella.final_label
            assign_rule = "enclosed"
        elif prev_ovl == 0 and next_ovl < 0:
            # assign-prev
            assign = previous_umbrella.final_label
            assign_rule = "touchleft"
        elif prev_ovl < 0 and next_ovl == 0:
            # assign-next
            assert next_umbrella is not None
            assign = next_umbrella.final_label
            assign_rule = "touchright"
        elif prev_ovl == 0 and next_ovl == 0:
            # assign-left rule
            assign = previous_umbrella.final_label
            assign_rule = "assignleft"
        else:
            raise
    else:
        sub.sort_values("overlap_bp", ascending=True, inplace=True)
        assign = sub["draft_label"].iloc[0]
        assign_rule = "overlap"
    label_group = umbrellas[assign]["group"]
    if DEBUG:
        print(assign)
        print(label_group)
        print(assign_rule)
    assert assign in umbrellas
    return assign, label_group, assign_rule


def load_intersection_table(fp, umbrellas):

    df = pd.read_csv(fp, sep="\t", header=None, names=HEADER)

    umbrella_labels = list(umbrellas.keys())

    all_labels = ISSUE_LABELS + umbrella_labels

    select_draft = df["draft_label"].isin(all_labels)
    select_final = df["final_label"].isin(all_labels)
    select_rows = select_draft & select_final

    df = df.loc[select_rows, :].copy()
    df["is_umbrella"] = df["final_label"].apply(lambda l: l in umbrellas)

    # relevant sorting for finding prev/next umbrella terms
    df = df.sort_values(["seq2", "start2", "end2"], inplace=False).reset_index(drop=True, inplace=False)

    return df


def main():

    args = parse_command_line()

    umbrellas = load_umbrella_labels(args.umbrella_labels)

    df = load_intersection_table(args.intersect, umbrellas)

    regions = []

    for region, overlaps in df.groupby(["seq2", "start2", "end2", "final_label", "score2", "strand2"]):
        region_infos = dict(
            (h, v) for h, v in zip(BED_HEADER, list(region))
        )
        start_match = True  # region_infos["start"] == 18475829
        end_match = True  # region_infos["end"] == 18483329
        name_match = False  # region_infos["name"] == "PAR1"
        DEBUG = start_match & end_match & name_match
        if region_infos["name"] in ["ERRBASE", "ERRSTRUCT"]:
            label_group = "ERROR"
            region_infos["label_group"] = label_group
            if region_infos["name"] == "ERRSTRUCT":
                region_issues = extract_issue_information(overlaps, region_infos["name"])
                region_infos.update(region_issues)
            region_infos["assigned"] = region_infos["name"]
            region_infos["assign_method"] = "identity"
            regions.append(region_infos)
            continue
        try:
            label_group = umbrellas[region_infos["name"]]["group"]
            region_infos["label_group"] = label_group
            region_infos["assigned"] = region_infos["name"]
            region_infos["assign_method"] = "identity"
            region_issues = extract_issue_information(overlaps)
        except KeyError:
            assert region_infos["name"] in ["NGAP", "UNASSIGNED"]
            region_issues = extract_issue_information(overlaps, region_infos["name"])
            assign, label_group, assign_method = assign_label(df, umbrellas, overlaps, region_infos["name"])
            region_infos["label_group"] = label_group
            region_infos["assigned"] = assign
            region_infos["assign_method"] = assign_method
        region_infos.update(region_issues)
        regions.append(region_infos)


    regions = pd.DataFrame.from_records(
        regions
    )

    for issue_label in ISSUE_LABELS:
        col_name = f"{issue_label}_n"
        try:
            regions[col_name] = regions[col_name].fillna(0, inplace=False).astype(int)
        except KeyError:
            regions[col_name] = 0
        col_name = f"{issue_label}_bp"
        try:
            regions[col_name] = regions[col_name].fillna(0, inplace=False).astype(int)
        except KeyError:
            regions[col_name] = 0

    regions.sort_values(["seq", "start", "end"], inplace=True)
    args.output.parent.mkdir(exist_ok=True, parents=True)
    regions.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
