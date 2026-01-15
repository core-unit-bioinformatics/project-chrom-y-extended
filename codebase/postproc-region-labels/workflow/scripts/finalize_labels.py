#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl
import re

import pandas as pd
import pyranges as pyr


DebugRegion = col.namedtuple(
    "DebugRegion",
    [
        "seq", "start", "end", "name", "score", "strand",
        "option", "alt", "region_start", "region_end"
    ]
)

MergedRegion = col.namedtuple(
    "MergedRegion",
    (
        "seq", "start", "end", "name", "score", "strand"
    )
)


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--isect-table", "-i",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="isect_table"
    )

    # parser.add_argument(
    #     "--region-annotation", "-r",
    #     type=lambda fp: pl.Path(fp).resolve(strict=True),
    #     dest="region_annotation"
    # )

    parser.add_argument(
        "--output", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output"
    )

    parser.add_argument(
        "--debug-out", "-d",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="debug_out"
    )

    args = parser.parse_args()

    return args


# def contain_splits(row1, row2):
#     """Generate containment splits:

#     s1 ================ e1
#           s2 ==== e2

#     s1 > s2
#     s2 > e2
#     e2 > e1

#     Args:
#         row1 (_type_): _description_
#         row2 (_type_): _description_

#     Returns:
#         _type_: _description_
#     """
#     splits = set(
#         [
#             (row1.start, row2.start, (row1.Index, row1.Index)),
#             (row2.start, row2.end, (row1.Index, row2.Index)),
#             (row2.end, row1.end, (row1.Index, row1.Index))

#         ]
#     )
#     return splits


# def overlap_splits(row1, row2):
#     """Generate overlap splits:

#     s1 ================ e1
#             s2 ============== e2

#     or

#     s1 ================ e1
#             s2 ======== e2

#     and so on...

#     s1 > s2
#     s2 > e1
#     e1 > e2

#     Args:
#         row1 (_type_): _description_
#         row2 (_type_): _description_

#     Returns:
#         _type_: _description_
#     """
#     splits = set(
#         [
#             (row1.start, row2.start, (row1.Index, row1.Index)),
#             (row2.start, row1.end, (row1.Index, row2.Index)),
#             (row1.end, row2.end, (row2.Index, row2.Index))
#         ]
#     )
#     return splits


# def bookend_splits(row1, row2):
#     splits = set(
#         [
#             (row1.start, row1.end, (row1.Index, row1.Index)),
#             (row2.start, row2.end, (row2.Index, row2.Index))
#         ]
#     )
#     return splits


# def identity_splits(row1, row2):
#     splits = set(
#         [
#             (row1.start, row1.end, (row1.Index, row2.Index)),
#         ]
#     )
#     return splits


# def apply_heuristics(row1, row2):

#     if row1["assign_method"] == "kmer":
#         return row1
#     if row2["assign_method"] == "kmer":
#         return row2
#     if row1["score"] > row2["score"]:
#         return row1
#     if row1["score"] < row2["score"]:
#         return row2
#     if row1["overlap_bp"] > row2["overlap_bp"]:
#         return row1
#     if row1["overlap_bp"] < row2["overlap_bp"]:
#         return row2

#     return None


# def determine_dominant_label(row1, row2, start, end, ref_regions):

#     option1 = row1["name"]
#     alt1 = row1["second_best_guess"]

#     option2 = row2["name"]
#     alt2 = row2["second_best_guess"]

#     if option1 == "ERR":
#         dbr = [window_to_debug_region(start, end, row1)]
#     elif option2 == "ERR":
#         dbr = [window_to_debug_region(start, end, row2)]
#     elif option1 == option2:
#         dbr = [window_to_debug_region(start, end, row1)]
#     elif is_proper_label(option1) and not is_proper_label(option2):
#         dbr = [window_to_debug_region(start, end, row1)]
#     elif is_proper_label(option2) and not is_proper_label(option1):
#         dbr = [window_to_debug_region(start, end, row2)]
#     else:
#         # both names must be proper labels or not
#         if is_proper_label(option1) and is_proper_label(option2):
#             if overlaps(option1, option2, ref_regions):
#                 # accept both labels
#                 dbr1 = window_to_debug_region(start, end, row1)
#                 dbr2 = window_to_debug_region(start, end, row2)
#                 dbr = [dbr1, dbr2]
#             else:
#                 selected = apply_heuristics(row1, row2)
#                 if selected is None:
#                     raise RuntimeError(f"{start} -> {end}: ", row1, row2)
#                 dbr = [window_to_debug_region(start, end, selected)]

#         elif is_proper_label(alt1) and is_proper_label(alt2):
#             if overlaps(alt1, alt2, ref_regions):
#                 # accept both labels
#                 dbr1 = window_to_debug_region(start, end, row1)
#                 dbr2 = window_to_debug_region(start, end, row2)
#                 dbr = [dbr1, dbr2]
#             else:
#                 selected = apply_heuristics(row1, row2)
#                 if selected is None:
#                     raise RuntimeError(f"{start} -> {end}: ", row1, row2)
#                 dbr = [window_to_debug_region(start, end, selected)]

#         else:
#             raise RuntimeError(f"{start} -> {end}: ", row1, row2)

#     return dbr


# def overlaps(label1, label2, ref_regions):

#     try:
#         start1, end1 = ref_regions.loc[ref_regions["seqclass"] == label1, ["start", "end"]].values[0]
#         start2, end2 = ref_regions.loc[ref_regions["seqclass"] == label2, ["start", "end"]].values[0]
#     except IndexError:
#         return False

#     overlap = min(end1, end2) - max(start1, start2)
#     # negative = distance / no overlap
#     # zero = book-ended
#     return overlap > 0


# def determine_window_labels(isect_table, ref_regions):

#     debug_output = set()
#     for cluster_id, regions in isect_table.groupby("cluster_id"):
#         primary_is_unique = regions["name"].nunique() == 1
#         secondary_is_unique = regions["second_best_guess"].nunique() == 1
#         # the 'is_label' criterion is only checked in case
#         # of a unique region, hence .iloc[0] is informative
#         primary_is_label = is_proper_label(regions["name"].iloc[0])
#         secondary_is_label = is_proper_label(regions["second_best_guess"].iloc[0])

#         if primary_is_unique and primary_is_label:
#             debug_output.add(merge_to_debug_region(regions))

#         elif secondary_is_unique and secondary_is_label:
#             debug_output.add(merge_to_debug_region(regions, False))

#         else:
#             # more involved case: the cluster consists of overlapping
#             # regions with different labels; this may be easy to resolve
#             # in case of strict labels such as ERR or problematic in case
#             # of umbrella/subregion categories
#             # Important to recall here that 'regions' is sorted
#             all_splits = set()
#             for row1, row2 in itt.pairwise(regions.itertuples()):
#                 if row2.start < row1.end and row2.end < row1.end:
#                     # row2 region is fully contained in row1 region
#                     all_splits.update(contain_splits(row1, row2))
#                 elif row2.start < row1.end and row2.end >= row1.end:
#                     # overhang situation
#                     all_splits.update(overlap_splits(row1, row2))
#                 elif row2.start == row1.end:
#                     # bookended regions
#                     all_splits.update(bookend_splits(row1, row2))
#                 elif row1.start == row2.start and row1.end == row2.end:
#                     all_splits.update(identity_splits(row1, row2))
#                 elif row1.end < row2.start:
#                     # no more overlap
#                     continue
#                 else:
#                     raise RuntimeError(row1, row2)
#             # splits can be empty, e.g., in the case of identical
#             # start coordinates; those cases can just be skipped
#             iter_nonempty_splits = filter(lambda split: split[0] != split[1], all_splits)
#             last_label = None
#             for (start, end, (idx1, idx2)) in sorted(iter_nonempty_splits):
#                 if idx1 == idx2:
#                     dbg_region = window_to_debug_region(start, end, isect_table.loc[idx1, :])
#                     debug_output.add(dbg_region)
#                     last_label = dbg_region.name
#                 else:
#                     row1 = isect_table.loc[idx1, :]
#                     row2 = isect_table.loc[idx2, :]
#                     dbg_regions = determine_dominant_label(row1, row2, start, end, ref_regions)
#                     # dbg_region = DebugRegion(
#                     #     row1["seq"], start, end, "unset",
#                     #     row1["name"], row1["second_best_guess"],
#                     #     row2["name"], row2["second_best_guess"]
#                     # )
#                     debug_output.update(dbg_regions)
#                     last_label = "unset"

#     debug_output = pd.DataFrame.from_records(sorted(debug_output), columns=DebugRegion._fields)
#     debug_output.sort_values(["seq", "start"], inplace=True)

#     return debug_output


def cluster_regions(windowed_regions):
    """
    Args:
        windowed_regions (pandas.DataFrame): _description_

    Returns:
        _type_: _description_
    """

    ranges = pyr.from_dict(
        {
            "Chromosome": windowed_regions["seq"].values,
            "Start": windowed_regions["start"].values,
            "End": windowed_regions["end"].values,
            "Strand": windowed_regions["strand"].values,
            "Name": windowed_regions["name"].values,
            "pd_idx": windowed_regions.index.values
        }
    )

    clustered = ranges.cluster(strand="same", by="Name", slack=2)
    # CAUTION
    # PyRanges has the very annoying habit of reorganizing
    # the underlying data by chromosome (and strand), which
    # means a simple assignment like this:
    # isect_table["cluster_id"] = clustered.Cluster.values
    # will NOT assign the correct cluster IDs to the original
    # entries in the DataFrame. Need to build a new df with
    # the matching index
    tmp_cluster_df = clustered.df.set_index("pd_idx", inplace=False)
    tmp_cluster_df.sort_index(inplace=True)
    windowed_regions["cluster_id"] = tmp_cluster_df.Cluster
    windowed_regions.sort_values(["seq", "start", "end"], axis=0, inplace=True)
    windowed_regions.reset_index(drop=True, inplace=True)

    return windowed_regions


def is_proper_label(label):
    return label not in ["uncertain", "none", "unknown"]


def merge_to_debug_regions(windows):

    assert windows.shape[0] > 1

    # Decision on 2026-01-14 / PH
    # the structural errors should be represented in the same way
    # as the k-mer (base-level) errors on top of the other labels;
    # added if False here to skip the code path that would ignore
    # all other labels and just generate an ERRSTRUCT window
    # for the respective region
    if False and windows["name"].isin(["ERRSTRUCT"]).any():
        dbr = window_to_debug_region(
            windows.loc[windows["name"].isin(["ERRSTRUCT"]), :].iloc[0].to_dict()
        )
        dbrs = [dbr]
    else:
        dbrs = [
            window_to_debug_region(row._asdict()) for
            row in windows.itertuples()
        ]
    return dbrs


def window_to_debug_region(window):

    assert isinstance(window, dict), window

    seq = window["seq"]
    start = max(window["start"], window["win_start"])
    end = min(window["end"], window["win_end"])
    score = 1000
    strand = window["strand"]
    if strand == ".":
        strand = "+"

    primary_label = window["name"]
    secondary_label = window["second_best_guess"]

    if primary_label == "ERRSTRUCT":
        score = 0
        strand = "+"
    if primary_label == "UNASSIGNED":
        score = 500
        strand = "+"

    if is_proper_label(primary_label):
        select_label = primary_label
    else:
        select_label = secondary_label

    dbg_region = DebugRegion(
        seq, start, end, select_label, score, strand,
        primary_label, secondary_label, window["start"], window["end"]
    )
    return dbg_region


def assign_label_per_window(isect_table):

    debug_output = []
    for win_name, labels in isect_table.groupby("win_name"):
        if labels.shape[0] > 1:
            dbrs = merge_to_debug_regions(labels)
            debug_output.extend(dbrs)
        else:
            dbr = window_to_debug_region(labels.iloc[0].to_dict())
            debug_output.append(dbr)

    debug_output = pd.DataFrame.from_records(sorted(debug_output), columns=DebugRegion._fields)
    debug_output.sort_values(["seq", "start"], inplace=True)

    return debug_output


def merge_clustered_windows(clustered_windows):

    final_regions = []
    for cluster_id, windows in clustered_windows.groupby("cluster_id"):
        reg = MergedRegion(
            windows["seq"].iloc[0], windows["start"].min(),
            windows["end"].max(), windows["name"].iloc[0],
            windows["strand"].iloc[0], windows["score"].max()
        )
        final_regions.append(reg)
    final_regions = pd.DataFrame.from_records(final_regions, columns=MergedRegion._fields)
    final_regions.sort_values(["seq", "start"], inplace=True)
    return final_regions


def simplify_color_block_labels(label):
    colors = "(" + "|".join(["blue", "red", "green", "yellow", "gray"]) +")"
    mobj = re.search(colors, label)
    if mobj is None:
        return label
    s,e = mobj.span()
    color = label[s:e]
    qualifier = "(plus|IR[0-9])"
    mobj = re.search(qualifier, label)
    if mobj is not None:
        s, e = mobj.span()
        qual = label[s:e]
        new_label = f"{color}-{qual}"
    else:
        new_label = color
    return new_label


def main():

    args = parse_command_line()

    isect_table = pd.read_csv(args.isect_table, sep="\t", header=0)
    #ref_regions = pd.read_csv(args.region_annotation, sep="\t", header=0)

    # simplify: drop score zero / check w/ Pille
    drop_names = ["uncertain", "TSPY-small", "TSPY-large"]
    isect_table = isect_table.loc[~isect_table["name"].isin(drop_names), :].copy()
    isect_table["name"] = isect_table["name"].apply(simplify_color_block_labels)
    isect_table["second_best_guess"] = isect_table["second_best_guess"].apply(simplify_color_block_labels)

    debug_output = assign_label_per_window(isect_table)
    args.debug_out.parent.mkdir(exist_ok=True, parents=True)
    debug_output.to_csv(args.debug_out, sep="\t", header=True, index=False)

    debug_output = cluster_regions(debug_output)
    final_regions = merge_clustered_windows(debug_output)
    args.output.parent.mkdir(exist_ok=True, parents=True)

    # make BED-like header
    final_regions.rename({"seq": "#seq"}, axis=1, inplace=True)
    final_regions.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
