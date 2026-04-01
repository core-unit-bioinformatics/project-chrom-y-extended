#!/usr/bin/env python3

import argparse as argp
import collections as col
import json
import pathlib as pl

import pandas as pd


ISECT_COLUMNS = [
    "seq", "start", "end", "label", "score", "strand"
]

HEADER = [f"{c}_ref" for c in ISECT_COLUMNS]
HEADER += [f"{c}_patch" for c in ISECT_COLUMNS]
HEADER += ["overlap_bp"]

ISSUES = ["ERRBASE", "ERRSTRUCT", "NGAP", "UNASSIGNED"]

EXTEND_THRESHOLD = 1000


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--isect", "-i",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="isect"
    )
    parser.add_argument(
        "--seq-sizes", "-s",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="seq_sizes"
    )
    parser.add_argument(
        "--reference-umbrellas", "-r",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="ref_umbrellas"
    )
    parser.add_argument(
        "--patch-umbrellas", "-p",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="patch_umbrellas"
    )
    parser.add_argument(
        "--output", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output"
    )
    args = parser.parse_args()

    return args


def load_umbrella_map(file_path):

    with open(file_path) as dump:
        umbrellas = json.load(dump)
    unified_terms = set()
    for term, alts in umbrellas.items():
        if alts["unified"] == "PAR1,PAR2":
            unified_terms.add("TELO")
            unified_terms.add("TELOp")
            unified_terms.add("TELOq")
        else:
            unified_terms.add(alts["unified"])
    return umbrellas, unified_terms


def load_seq_sizes(file_path):

    sizes = dict()
    with open(file_path) as listing:
        for line in listing:
            name, size = line.strip().split()
            sizes[name] = int(size)
    return sizes


def subset_intersect_table(file_path, ref_umbrellas, patch_umbrellas):

    df = pd.read_csv(file_path, sep="\t", header=None, names=HEADER)

    # 1 - select all rows where both labels are umbrellas
    df["umbrella_ref"] = df["label_ref"].apply(lambda l: 1 if l in ref_umbrellas else 0)
    df["umbrella_patch"] = df["label_patch"].apply(lambda l: 1 if l in patch_umbrellas else 0)
    select_umbrella = (df["umbrella_ref"] + df["umbrella_patch"]) > 1

    # 2 - select issue on issue (i.e. NGAP)
    df["issue_ref"] = df["label_ref"].apply(lambda l: 1 if l in ISSUES else 0)
    df["issue_patch"] = df["label_patch"].apply(lambda l: 1 if l in ISSUES else 0)

    select_issue = (df["issue_ref"] + df["issue_patch"]) > 1

    # 3 - select mixed rows (umbrella or issue)
    select_mixed = ((df["umbrella_ref"] + df["umbrella_patch"]) > 0) & ((df["issue_ref"] + df["issue_patch"]) > 0)

    select_rows = select_umbrella | select_issue | select_mixed

    sub = df.loc[select_rows, :].copy()
    return sub


def extract_overlapping_region(isect_row):

    max_start = max(isect_row.start_ref, isect_row.start_patch)
    min_end = min(isect_row.end_ref, isect_row.end_patch)
    assert max_start < min_end
    return max_start, min_end


def patch_regions(regions, ref_umbrellas, ref_unified, patch_umbrellas):

    label_regions = set()
    issue_regions = set()

    for row in regions.itertuples():
        if row.label_ref == "UNASSIGNED" and row.label_patch not in ISSUES:
            # patching
            start, end = extract_overlapping_region(row)
            umbrellas = patch_umbrellas
            label = row.label_patch
        elif row.label_ref == "UNASSIGNED" and row.label_patch == "UNASSIGNED":
            # unclear is necessary for any sample
            start, end = extract_overlapping_region(row)
            label = "UNASSIGNED"
            umbrellas = dict()
        else:
            # just use the ref label
            start, end = row.start_ref, row.end_ref
            label = row.label_ref
            umbrellas = ref_umbrellas
        if label in ISSUES:
            region = (row.seq_ref, start, end, label, 0, "+", "ISSUE", "NOLABEL", "NOGROUP", "KEEP")
            issue_regions.add(region)
            continue

        unified = umbrellas[label]["unified"]
        group = umbrellas[label]["group"]
        if unified == "PAR1,PAR2":
            if start < int(5e6):
                unified = "TELOp"
                group = "TEL"
            else:
                unified = "TELOq"
                group = "TEL"
        if unified == "periCEN":
            # 2026-04-01
            # change requested by Arang Rhie
            group = "SAT"
        try:
            assert unified in ref_unified, unified
        except AssertionError:
            if unified == "AMPL1" and "AMPL1_IR3d" in ref_unified:
                unified = "AMPL1_IR3d"
            elif unified == "AMPL1_IR3d" and "AMPL1" in ref_unified:
                unified = "AMPL1"
            elif unified == "OTHER" and row.start_ref > int(20e6):
                unified = "other2"
                assert unified in ref_unified, "remap failed"
            else:
                raise

        region = (
            row.seq_ref, start, end, unified,
            1000, row.strand_ref,
            "LABEL", unified, group, "KEEP"
        )
        label_regions.add(region)

    return label_regions, issue_regions


def disjoin_regions(current_region, next_region):

    overlap = current_region[2] - next_region[1]
    # AMPL* labels do not have priority
    if current_region[3].startswith("AMPL"):
        # subtract all from AMPL
        current_region[2] -= overlap
        assert current_region[1] < current_region[2]
        assert current_region[2] == next_region[1]
    elif next_region[3].startswith("AMPL"):
        next_region[1] += overlap
        assert next_region[1] < next_region[2]
        assert current_region[2] == next_region[1]
    else:
        if overlap % 2 == 0:
            half_rev, half_fwd = overlap // 2, overlap // 2
        else:
            half_rev, half_fwd = overlap // 2, (overlap // 2) + 1
        current_region[2] -= half_rev
        next_region[1] += half_fwd
        assert current_region[1] < current_region[2]
        assert next_region[1] < next_region[2]
        assert current_region[2] == next_region[1]
    return tuple(current_region), tuple(next_region)


def extend_regions(current_region, next_region):

    if current_region[3] == "TELOp":
        assert next_region[3] == "PAR1"
        # extend PAR1 backwards to TELOp
        next_region[1] = current_region[2]
    elif current_region[3] == "PAR2":
        assert next_region == "TELOq"
        # extend PAR2 forward to TELOq
        current_region[2] = next_region[1]
    else:
        pass
    return tuple(current_region), tuple(next_region)


def stitch_up_label_regions(label_regions, seq_sizes):

    label_regions = col.deque(sorted(label_regions, key=lambda r: (r[0], r[1])))
    label_regions.append(None)

    final_regions = []
    current_region = label_regions.popleft()
    # 0 - seq
    # 1 - start /  2 - end
    # 3 - label(name)
    while 1:
        next_region = label_regions.popleft()
        if next_region is None:
            final_regions.append(current_region)
            break
        if current_region[0] != next_region[0]:
            # seq switch
            # check proper begin/end
            current_region = list(current_region)
            if current_region[1] < EXTEND_THRESHOLD:  #type: ignore
                current_region[1] = 0  #type: ignore
            if seq_sizes[current_region[0]] - current_region[2] < EXTEND_THRESHOLD:
                current_region[2] = seq_sizes[current_region[0]]
            if current_region[3] == "PAR2":
                # end w/o TELOq
                # extend up to any distance
                # from seq end
                current_region[2] = seq_sizes[current_region[0]]
            current_region = tuple(current_region)

            next_region = list(next_region)
            if next_region[3] == "PAR1":
                # begin w/o TELOp
                next_region[1] = 0
            next_region = tuple(next_region)

            final_regions.append(current_region)
            current_region = next_region
            continue
        # same label - stitch together
        if current_region[3] == next_region[3]:
            new_current = list(current_region)
            new_current[2] = next_region[2]
            assert new_current[1] < new_current[2]
            current_region = tuple(new_current)
            continue
        # bookended - just iterate through
        if current_region[2] == next_region[1]:
            final_regions.append(current_region)
            current_region = next_region
            continue
        if current_region[2] > next_region[1]:
            current_region, next_region = disjoin_regions(list(current_region), list(next_region))
            final_regions.append(current_region)
            current_region = next_region
        else:
            current_region, next_region = extend_regions(list(current_region), list(next_region))
            final_regions.append(current_region)
            current_region = next_region
    return final_regions


def merge_labels_and_issues(label_regions, issue_regions, seq_sizes):

    merged = col.deque(sorted(label_regions, key=lambda r: (r[0], r[1])))
    merged.append(None)

    issue_regions = [issue for issue in issue_regions if issue[3] != "UNASSIGNED"]

    current_region = merged.popleft()
    region_per_seq = col.Counter()
    while 1:
        next_region = merged.popleft()
        if next_region is None:
            merged.append(current_region)
            break
        if current_region[0] != next_region[0]:
            if region_per_seq[current_region[0]] == 0:
                # first and only region for sequence
                if current_region[1] != 0:
                    unassigned_start = 0
                    unassigned_end = current_region[1]
                    assert unassigned_start < unassigned_end
                    issue_regions.append(
                        (
                            current_region[0], unassigned_start, unassigned_end,
                            "UNASSIGNED", 0, "+", "LABEL", "NOLABEL", "UNASSIGNED", "KEEP"
                        )
                    )
            if current_region[2] != seq_sizes[current_region[0]]:
                unassigned_start = current_region[2]
                unassigned_end = seq_sizes[current_region[0]]
                assert unassigned_start < unassigned_end
                issue_regions.append(
                    (
                        current_region[0], unassigned_start, unassigned_end,
                        "UNASSIGNED", 0, "+", "LABEL", "NOLABEL", "UNASSIGNED", "KEEP"
                    )
                )
            merged.append(current_region)
            region_per_seq[current_region[0]] += 1
            current_region = next_region
            continue
        if current_region[2] < next_region[1]:
            unassigned_start = current_region[2]
            unassigned_end = next_region[1]
            assert unassigned_start < unassigned_end
            issue_regions.append(
                (
                    current_region[0], unassigned_start, unassigned_end,
                    "UNASSIGNED", 0, "+", "LABEL", "NOLABEL", "UNASSIGNED", "KEEP"
                )
            )
        merged.append(current_region)
        region_per_seq[current_region[0]] += 1
        current_region = next_region

    merged = sorted(merged)
    merged += sorted(issue_regions, key=lambda r: (r[0], r[1]))

    merged = pd.DataFrame.from_records(
        merged,
        columns=[
            "seq", "start", "end", "label",
            "score", "strand",
            "category", "assign_label", "assign_group", "filter"
        ]
    )
    merged = merged.sort_values(["seq", "start"], inplace=False).reset_index(drop=True, inplace=False)

    set_values = []
    for (q, s, e, l), region in merged.groupby(["seq", "start", "end", "label"]):
        if l in ISSUES:
            continue
        select_seq = merged["seq"] == q
        select_start = merged["start"] >= s
        select_end = merged["end"] <= e
        select_sub = select_seq & select_start & select_end
        if select_sub.sum() == 1:
            continue
        subsumed = merged.loc[select_sub, :]
        for row in subsumed.itertuples():
            if row.category == "LABEL":
                continue
            assert row.label in ISSUES
            set_values.append((row.Index, "assign_label", region.iloc[0].assign_label))
            set_values.append((row.Index, "assign_group", region.iloc[0].assign_group))
            if row.label == "UNASSIGNED":
                set_values.append((row.Index, "filter", "DROP"))

    for idx, column, value in set_values:
        merged.loc[idx, column] = value

    return merged


def main():

    args = parse_command_line()

    seq_sizes = load_seq_sizes(args.seq_sizes)

    ref_umbrellas, ref_unified = load_umbrella_map(args.ref_umbrellas)
    patch_umbrellas, _ = load_umbrella_map(args.patch_umbrellas)

    df = subset_intersect_table(args.isect, ref_umbrellas, patch_umbrellas)

    label_regions, issue_regions = patch_regions(df, ref_umbrellas, ref_unified, patch_umbrellas)

    label_regions = stitch_up_label_regions(label_regions, seq_sizes)

    merged = merge_labels_and_issues(label_regions, issue_regions, seq_sizes)

    merged.drop("filter", axis=1, inplace=True)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    merged.rename({"seq": "#seq"}, axis=1, inplace=True)
    merged.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
