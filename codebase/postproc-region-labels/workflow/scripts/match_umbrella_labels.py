#!/usr/bin/env python3

import argparse as argp
import collections as col
import json
import pathlib as pl
import re

import pandas as pd


DEV_REF1 = None
DEV_REF2 = None

DEV_SMP1 = None
DEV_SMP2 = None


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "-r", "--ref-regions",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="ref_regions",
        required=True,
    )
    parser.add_argument(
        "-s", "--sample-regions",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="sample_regions",
        required=True,
    )
    parser.add_argument(
        "-n", "--ref-name",
        type=str,
        required=True,
        dest="ref_name"
    )
    parser.add_argument(
        "-o", "--output",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="output",
        required=True,
    )
    parser.add_argument(
        "-t", "--test",
        type=int,
        choices=[1,2],
        dest="test",
        default=None
    )

    args = parser.parse_args()

    if False:
        if args.ref_regions is None:
            if args.test is None:
                raise
            if args.test == 1:
                setattr(args, "ref_regions", DEV_REF1)
                setattr(args, "sample_regions", DEV_SMP1)
            else:
                setattr(args, "ref_regions", DEV_REF2)
                setattr(args, "sample_regions", DEV_SMP2)

    return args


def make_alphanum_label(label):
    reduced = "".join(re.findall("[a-z0-9]+", label, flags=re.IGNORECASE)).upper()
    return reduced


def make_char_label(label):
    reduced = "".join(re.findall("[a-z]+", label, flags=re.IGNORECASE)).upper()
    return reduced


def is_umbrella_label(label):
    return re.match("^[0-9]+[un]_", label) is not None


def get_unified_label(label):
    """This implements the wishes of the project lead.

    Args:
        label (_type_): _description_

    Returns:
        _type_: _description_
    """
    if "DYZ3" in label:
        return "CEN-DYZ3"
    if "DYZ17" in label:
        return "periCEN"
    if "DYZ19" in label:
        return "DYZ19"
    if "SAT" in label:
        return "periCEN"
    if "CEN1" in label or "CEN2" in label:
        return "periCEN"
    if "SAT" in label:
        return "periCEN"
    if "CEN" == label:
        return "CEN-DYZ3"
    return label


def get_prefix_group(label):
    if any(part in label for part in ["CEN", "DYZ"]):
        return label
    if "AMP" in label:
        return "AMPL"
    return label[:3].upper()


def read_regions(file_path):
    df = pd.read_csv(file_path, sep="\t", header=0)
    return df


def read_reference_regions(file_path):

    df = read_regions(file_path)

    df["is_umbrella"] = df["name"].apply(is_umbrella_label)
    color_blocks = ["red", "green", "yellow", "blue", "grey", "gray", "teal"]

    non_umbrella_labels = set(df.loc[~df["is_umbrella"], "seqclass"].values)

    color_only = "|".join(sorted(make_char_label(label) for label in color_blocks))
    alnum = "|".join(sorted(make_alphanum_label(label) for label in non_umbrella_labels))
    exclude_match = re.compile(
        f"^({color_only}|{alnum})$",
        flags=re.IGNORECASE
    )

    df = df.loc[df["is_umbrella"], :].copy()
    df["check_label"] = df["seqclass"].apply(make_alphanum_label)
    df["unified_label"] = df["seqclass"].apply(get_unified_label)
    df["prefix_group"] = df["unified_label"].apply(get_prefix_group)

    label_mult = col.Counter(df["check_label"].values)
    seen = col.Counter()
    ref_labels = col.defaultdict(list)
    for row in df.itertuples():
        if label_mult[row.check_label] > 1:
            seen[row.check_label] += 1
            enum_suffix = str(seen[row.check_label])
        else:
            enum_suffix = ""
        ref_labels[row.check_label].append(
            (
                row.name, row.seqclass, row.unified_label, row.prefix_group, f"{row.seqclass}{enum_suffix}"
            )
        )
    return dict(ref_labels), exclude_match


def reduce_match_info(matchings):

    labels = set()
    unified = set()
    prefixes = set()
    for m in matchings:
        try:
            _, sq, uni, pfx, sq2 = m
        except ValueError:
            print(matchings)
            raise
        labels.add(sq)
        labels.add(sq2)
        unified.add(uni)
        prefixes.add(pfx)
    labels = ",".join(sorted(labels))
    unified = ",".join(sorted(unified))
    prefixes = ",".join(sorted(prefixes))
    return labels, unified, prefixes


def find_complex_match(check_label, ref_labels):
    matchings = None
    if check_label == "TELO":
        matchings = [ref_labels["PAR1"][0], ref_labels["PAR2"][0]]
    unified_label = get_unified_label(check_label)
    for check_label, label_infos in ref_labels.items():
        for label_info in label_infos:
            _, _, uni, _, _ = label_info
            if unified_label.lower() == uni.lower():
                if matchings is None:
                    matchings = [label_info]
                else:
                    matchings.append(label_info)
    return matchings


def find_matching_reference_label(sample_regions, ref_labels, match_non_umbrella):

    issue_labels = ["UNASSIGNED", "ERRBASE", "ERRSTRUCT", "NGAP"]
    # HMMER motifs were added as part of the process and can thus not
    # be recognized given the reference labels
    # --- this here is the explicit = safe way of dealing with them
    motif_labels = ["TSPY", "YQHET3K1BP", "YQHET2K7BP", "DYZ18YQ", "DYZ1YQ", "DYZ2CON"]
    # umbrella motif labels are:
    # (CEN-)DYZ3, DYZ19 --- hg38
    # DYZ17, DYZ19 --- t2t

    matched_labels = []
    unified_labels = []
    matched_prefixes = []
    for row in sample_regions.itertuples():
        check_label = make_alphanum_label(row.name)

        is_non_umbrella = match_non_umbrella.search(check_label) is not None
        is_issue_label = check_label in issue_labels
        is_motif_label = check_label in motif_labels

        if is_non_umbrella or is_issue_label or is_motif_label:
            if is_non_umbrella:
                print("skipping ", check_label)
            matched_labels.append("non-umbrella")
            unified_labels.append("non-umbrella")
            matched_prefixes.append("non-umbrella")
            continue
        try:
            matchings = ref_labels[check_label]
        except KeyError:
            matchings = find_complex_match(check_label, ref_labels)
            if matchings is None:
                print("no match ", check_label)
                matched_labels.append("non-umbrella")
                unified_labels.append("non-umbrella")
                matched_prefixes.append("non-umbrella")
                continue
        labels, unified, prefixes = reduce_match_info(matchings)
        matched_labels.append(labels)
        unified_labels.append(unified)
        matched_prefixes.append(prefixes)
    sample_regions["ref_umbrella_label"] = matched_labels
    sample_regions["unified_label"] = unified_labels
    sample_regions["umbrella_prefix"] = matched_prefixes
    return sample_regions


def write_ref_renamer(output_path, ref_name, ref_labels):

    # row.name, row.seqclass, row.unified_label, row.prefix_group, f"{row.seqclass}{enum_suffix}"
    renamer = dict()
    for _, label_info in ref_labels.items():
        full_name, seqclass, unified, group = label_info[0][:4]
        renamer[full_name] = {
            "seqclass": seqclass,
            "unified": unified,
            "group": group
        }
    output_path.parent.mkdir(exist_ok=True, parents=True)
    out_json = output_path.parent.joinpath(f"{ref_name}.unified-umbrella.json")
    if out_json.is_file():
        return None
    with open(out_json, "w") as dump:
        json.dump(renamer, dump, indent=2, ensure_ascii=True)
    return None


def write_sample_renamer(output_path, sample_regions):

    renamer = dict()
    for row in sample_regions.itertuples():
        if row.unified_label == "non-umbrella":
            continue
        renamer[row.name] = {
            "seqclass": row.ref_umbrella_label,
            "unified": row.unified_label,
            "group": row.umbrella_prefix
        }

    with open(output_path, "w") as dump:
        json.dump(renamer, dump, indent=2, ensure_ascii=True)
    return None


def main():

    args = parse_command_line()
    ref_labels, match_non_umbrella = read_reference_regions(args.ref_regions)
    sample_regions = read_regions(args.sample_regions)
    sample_regions = find_matching_reference_label(sample_regions, ref_labels, match_non_umbrella)

    args.output.parent.mkdir(exist_ok=True, parents=True)

    sample_regions.to_csv(args.output, header=True, index=False, sep="\t")

    write_ref_renamer(args.output, args.ref_name, ref_labels)
    write_sample_renamer(args.output, sample_regions)

    return 0



if __name__ == "__main__":
    main()
