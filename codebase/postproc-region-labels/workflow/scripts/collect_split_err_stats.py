#!/usr/bin/env python3

import argparse as argp
import collections as col
import json
import pathlib as pl

import pandas as pd


DEBUG = True

DBG_ISECT_TABLE = pl.Path("/home/ebertp/work/projects/chrom-y-extended/wf-data/region-labels/split_qc/NA19043.t2tv2.chrY-regions-split-errors.rhd.tsv")
DBG_UMBRELLA_MAP = pl.Path("/home/ebertp/work/projects/chrom-y-extended/wf-data/region-labels/umbrellas/t2tv2/NA19043.t2tv2.chrY-unified-umbrella.json")


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "-i", "--isect-table",
        type=lambda fp: pl.Path(fp).resolve(strict=not DEBUG),
        dest="isect_table",
        required=not DEBUG
    )

    parser.add_argument(
        "-u", "--umbrella-map",
        type=lambda fp: pl.Path(fp).resolve(strict=not DEBUG),
        dest="umbrella_map",
        required=not DEBUG
    )

    parser.add_argument(
        "-o", "--output",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="output",
        required=not DEBUG
    )

    args = parser.parse_args()

    if DEBUG:
        setattr(args, "isect_table", DBG_ISECT_TABLE)
        setattr(args, "umbrella_map", DBG_UMBRELLA_MAP)

    return args


def load_umbrella_mapping(fp):

    with open(fp) as dump:
        umbrella_labels = json.load(dump)
    return umbrella_labels


def load_isect_table(fp):

    df = pd.read_csv(fp, sep="\t", header=0)
    filename = fp.name
    sample, ref = filename.split(".")[:2]
    return df, sample, ref


def main():

    ISSUE_LABELS = ["UNASSIGNED", "ERRSTRUCT", "ERRBASE", "NGAP"]

    args = parse_command_line()
    umbrellas = load_umbrella_mapping(args.umbrella_map)
    isect_table, sample, ref = load_isect_table(args.isect_table)

    stats = col.defaultdict(col.Counter)
    for (seq, start, end, name), overlaps in isect_table.groupby(["seq", "start", "end", "name"]):
        if name in ISSUE_LABELS:
            stat_label = name
        elif name in umbrellas:
            stat_label = umbrellas[name]["unified"]
            if "," in stat_label:
                if name == "TELO" and start < int(2.5e6):
                    stat_label = "TELOp"
                elif name == "TELO" and start > int(4e7):
                    stat_label = "TELOq"
                else:
                    raise ValueError(seq, name, stat_label)
        else:
            continue
        size = end - start
        stats[stat_label][(ref, sample, "size")] += size

        flagger_is_dirty = (overlaps["flagger_hifi_is_clean"] < 1) & (overlaps["flagger_ont_is_clean"] < 1)
        nucflag_is_dirty = (overlaps["nucflag_hifi_is_clean"] < 1) & (overlaps["nucflag_ont_is_clean"] < 1)
        both_hit = flagger_is_dirty & nucflag_is_dirty

        hit_counts = {
            (ref, sample, "flagger_err"): overlaps.loc[flagger_is_dirty, "overlap_bp"].sum(),
            (ref, sample, "nucflag_err"): overlaps.loc[nucflag_is_dirty, "overlap_bp"].sum(),
            (ref, sample, "both_err"): overlaps.loc[both_hit, "overlap_bp"].sum()
        }
        stats[stat_label].update(hit_counts)

    stats = pd.DataFrame(stats)
    args.output.parent.mkdir(exist_ok=True, parents=True)
    stats.to_csv(args.output, sep="\t", header=True, index=True, index_label=["ref", "sample", "statistic"])

    return 0


if __name__ == "__main__":
    main()
