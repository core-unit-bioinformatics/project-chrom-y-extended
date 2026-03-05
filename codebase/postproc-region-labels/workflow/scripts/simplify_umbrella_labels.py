#!/usr/bin/env python3

import argparse as argp
import collections as col
import itertools as itt
import json
import pathlib as pl

import pandas as pd


ISSUE_LABELS = ["UNASSIGNED", "ERRBASE", "ERRSTRUCT", "NGAP"]


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "-r", "--region-labels",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="region_labels",
        required=True,
    )
    parser.add_argument(
        "-l", "--umbrella-labels",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="umbrella_labels",
        required=True,
    )
    parser.add_argument(
        "-d", "--disjoined",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="disjoined",
        required=True
    )
    parser.add_argument(
        "-s", "--stitched",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="stitched",
        required=True
    )

    args = parser.parse_args()

    return args


def load_umbrella_labels(fp):
    umbrella_labels = json.load(open(fp))
    return umbrella_labels


def load_sample_regions(fp, umbrellas):

    regions = pd.read_csv(fp, sep="\t", header=0)
    regions.rename({"#seq": "seq"}, axis=1, inplace=True)
    regions.sort_values(["seq", "start"], inplace=True)

    regions["skip"] = regions["name"].apply(lambda n: n not in umbrellas and n not in ISSUE_LABELS)
    regions = regions.loc[~regions["skip"], :].copy()
    regions.drop("skip", axis=1, inplace=True)

    return regions


def to_dict(region):
    d = region._asdict()
    d["source_start"] = d["start"]
    d["source_end"] = d["end"]
    if region.seq.endswith("_chrY"):
        d["seq_type"] = "main"
    else:
        d["seq_type"] = "rand"
    return d


def overlap(region1, region2):
    if region1["seq"] != region2["seq"]:
        return -1
    ovl = min(region1["end"], region2["end"]) - max(region1["start"], region2["start"])
    return ovl


def check_region(region):

    l = region["end"] - region["start"]
    if l < 2:
        raise ValueError(f"invalid region: {region}")
    assert region["end"] <= region["source_end"]
    assert region["start"] >= region["source_start"]
    return


def idname(r1, r2):
    return r1["seq"] == r2["seq"] and r1["name"] == r2["name"]


def is_enclosed(r1, r2):
    r1_len = r1["end"] - r1["start"]
    r2_len = r2["end"] - r2["start"]
    len_ratio = r1_len / r2_len
    start_beyond = r2["start"] <= r1["start"]
    end_before = r1["end"] <= r2["end"]
    group_match = r1["group"] == r2["group"]
    mismap = len_ratio < 0.05
    group_enclosed = start_beyond & end_before & group_match
    mismap_enclosed = start_beyond & end_before & mismap
    return group_enclosed or mismap_enclosed


def adjust_buffered_regions(buffer):

    buffer.append(None)
    adjusted = False
    enclosed = 0
    while 1:
        r1 = buffer.popleft()
        r2 = buffer.popleft()
        if r2 is None:
            buffer.append(r1)
            buffer.append(None)
            if not adjusted:
                break
            adjusted = False
            # next cycle
            continue
        ovl = overlap(r1, r2)
        if ovl < 0:
            raise RuntimeError(f"logic error: {ovl} / {r1} / {r2}")
        elif ovl == 0:
            buffer.append(r1)
            buffer.appendleft(r2)
        else:
            # special check for full enclosure
            if is_enclosed(r1, r2):
                # r1 contained in r2
                # keep only r2
                r2["disjoined"] = True
                buffer.appendleft(r2)
                adjusted = True
                enclosed += 1
                continue
            elif is_enclosed(r2, r1):
                # r2 contained in r1
                r1["disjoined"] = True
                buffer.appendleft(r1)
                adjusted = True
                enclosed += 1
                continue
            else:
                pass

            try:
                check_region(r1)
                check_region(r2)
            except ValueError:
                print(r1)
                print(r2)
                print(ovl)
                raise ValueError("pre-adjust error")

            half = ovl // 2
            r1["end"] = r1["end"] - half
            r2["start"] = r1["end"]
            try:
                check_region(r1)
                check_region(r2)
            except ValueError:
                print(r1)
                print(r2)
                print(ovl)
                raise ValueError("post-adjust error")
            r1["disjoined"] = True
            r2["disjoined"] = True
            buffer.append(r1)
            buffer.appendleft(r2)
            adjusted = True

    for r1, r2 in itt.pairwise(buffer):
        if r2 is None:
            break
        ovl = overlap(r1, r2)
        if ovl != 0:
            raise RuntimeError(f"adjustment failed: {ovl} / {r1} / {r2}")
    return buffer, enclosed


def disjoin_regions(regions, umbrellas):

    # processed regions
    final = col.deque()
    # active regions to check for overlap
    active = col.deque()
    # buffer overlapping regions
    buffer = col.deque()
    last_issue = -1
    row_n = -1
    total_rows = regions.shape[0]
    total_enclosed = 0
    regions_per_seq = regions.agg("seq").value_counts().to_dict()
    for seq, seq_regions in regions.groupby("seq"):
        for row_n, region in enumerate(seq_regions.itertuples(index=False), start=1):
            rd = to_dict(region)
            rd["disjoined"] = False
            if region.name in ISSUE_LABELS:
                rd["group"] = "ISSUE"
                final.append(rd)
                last_issue = row_n
                continue
            label_infos = umbrellas[rd["name"]]
            rd["group"] = label_infos["group"]
            if rd["name"] != "TELO":
                if "," in label_infos["unified"]:
                    rd["name"] = label_infos["group"]
            while 1:
                try:
                    other = active.popleft()
                    print("other ", other)
                    if overlap(rd, other) > 0:
                        buffer.append(other)
                    else:
                        final.append(other)
                except IndexError:
                    # active is empty now
                    buffer.append(rd)
                    break
            if len(buffer) > 0:
                adj_regions, enclosed = adjust_buffered_regions(buffer)
                total_enclosed += enclosed
                while 1:
                    r1 = adj_regions.popleft()
                    r2 = adj_regions.popleft()
                    if r2 is None:
                        # implies that r1 is 'region' above
                        assert len(active) == 0
                        if row_n == regions_per_seq[seq]:
                            # last region of this seq does not need
                            # to go back into the cycle
                            final.append(r1)
                        else:
                            active.append(r1)
                        buffer = col.deque()
                        break
                    else:
                        final.append(r1)
                        adj_regions.appendleft(r2)
    if len(active) == 1 and row_n in [last_issue + 1, last_issue]:
        final.append(active.pop())
    elif len(active) == 0:
        pass
    else:
        raise RuntimeError(f"active not empty: {active.pop()} / {row_n}")
    assert len(buffer) == 0

    assert len(final) == (total_rows - total_enclosed)

    disjoined = pd.DataFrame.from_records(
        sorted(final, key=lambda d: (d["seq"],d["start"]))
    )
    return disjoined


def stitch_regions(disjoined):

    # exclude N-gap, not stitching over that
    to_stitch = col.deque(
        to_dict(row) for row in
        disjoined.loc[~disjoined["name"].isin(ISSUE_LABELS[:-1]), :].itertuples(index=False)
    )
    to_stitch.append(None)

    stitched = col.deque()
    active = None

    while 1:
        r1 = to_stitch.popleft()
        r2 = to_stitch.popleft()
        if r2 is None:
            if active is None:
                stitched.append(r1)
                break
            elif idname(active, r1):
                active["end"] = r1["end"]
                active["strand"] = "+"
                stitched.append(active)
                break
            else:
                stitched.append(active)
                stitched.append(r2)
                break
        if idname(r1, r2):
            if active is not None:
                active["end"] = r2["end"]
                active["source_end"] = r2["source_end"]
                active["strand"] = "+"
                to_stitch.appendleft(r2)
                continue
            else:
                active = dict(r2)
                active["source_start"] = r1["source_start"]
                active["start"] = r1["start"]
                active["strand"] = "+"
                to_stitch.appendleft(r2)
                continue
        else:
            if active is None:
                if r2["name"] == "NGAP":
                    r1["end"] = r2["start"]
                if r1["name"] == "NGAP":
                    r2["start"] = r1["end"]
                stitched.append(r1)
                to_stitch.appendleft(r2)
                continue
            elif idname(active, r1):
                # should be identity op
                active["end"] = r1["end"]
                if r2["name"] == "NGAP":
                    active["end"] = r2["start"]
                active["strand"] = "+"
                active["source_end"] = r1["source_end"]
                stitched.append(active)
                active = None
                to_stitch.appendleft(r2)
                continue
            else:
                stitched.append(active)
                active = None
                stitched.append(r1)
                to_stitch.appendleft(r2)
                continue

    if active is not None:
        stitched.append(active)
    stitched = pd.DataFrame.from_records(
        sorted(stitched, key=lambda d: (d["seq"],d["start"]))
    )
    stitched.drop("disjoined", axis=1, inplace=True)

    return stitched


def main():

    args = parse_command_line()

    umbrellas = load_umbrella_labels(args.umbrella_labels)
    regions = load_sample_regions(args.region_labels, umbrellas)

    disjoined = disjoin_regions(regions, umbrellas)
    stitched = stitch_regions(disjoined)
    stitched = pd.concat(
        [stitched, disjoined.loc[disjoined["name"].isin(ISSUE_LABELS[:-1]), :].copy()],
        axis=0, ignore_index=False
    ).sort_values(["seq", "start"], inplace=False)
    stitched.drop("disjoined", axis=1, inplace=True)

    args.disjoined.parent.mkdir(exist_ok=True, parents=True)
    args.stitched.parent.mkdir(exist_ok=True, parents=True)
    disjoined.to_csv(args.disjoined, sep="\t", header=True, index=False)
    stitched.to_csv(args.stitched, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
