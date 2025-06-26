import pandas
import collections

_POSTPROC_REGION_LABELS_PYUTILS_CACHE = None


def _init_label_cache():

    labels_lut = collections.defaultdict(set)
    for reference, ref_label_file in MODULE_REF_LABELINGS.items():
        ref_label_file = WD.joinpath("references", ref_label_file).resolve(strict=True)
        labels = set(pandas.read_csv(
            ref_label_file, sep="\t",header=0
        )["seqclass"])
        labels_lut[reference].update(labels)
        labels_lut["all"].update(labels)
    global _POSTPROC_REGION_LABELS_PYUTILS_CACHE
    _POSTPROC_REGION_LABELS_PYUTILS_CACHE = labels_lut
    return None


def normalize_label_name(label_name, reference="all"):
    """Normalize label name to a standard form.

    Args:
        label_name (str): The label name to normalize.
        reference (str): The reference genome to use for normalization.

    Returns:
        str: The normalized label name.
    """
    global _POSTPROC_REGION_LABELS_PYUTILS_CACHE
    if _POSTPROC_REGION_LABELS_PYUTILS_CACHE is None:
        _init_label_cache()

    reference_labels = _POSTPROC_REGION_LABELS_PYUTILS_CACHE[reference]

    if label_name in reference_labels:
        norm_label = label_name
    elif "bIR" in label_name:
        norm_label = label_name.replace("bIR", "IR")
        parts = norm_label.split("-", 1)
        norm_label = parts[0] + "_" + parts[1]
    elif "gIR" in label_name:
        norm_label = label_name.replace("gIR", "IR")
        parts = norm_label.split("-", 1)
        norm_label = parts[0] + "_" + parts[1]
    elif "spacer" in label_name:
        spacer_num = label_name[-1]
        assert int(spacer_num) in [1,2,3,4,5,6,7,8], "Invalid spacer number"
        norm_label = f"P{spacer_num}-spacer"
    else:
        pass

    if norm_label not in reference_labels:
        raise ValueError(f"Normalization failed: from {label_name} to {norm_label} (ref: {reference})")

    return norm_label


def load_ml_motif_hits(file_path, motif_name):

    # these are the names that show up in the files
    # processed here --- see annotation/norm/motifs
    known_motif_names = [
        'DYZ18_Yq', 'DYZ19_Yq', 'DYZ1_Yq','DYZ2_Con',
        'DYZ3-prim_Ycentro', 'DYZ3-sec_Ycentro',
        'TSPY', 'Yqhet_2k7bp', 'Yqhet_3k1bp'
    ]

    if motif_name in ["tspy"]:
        header = ["#seq", "name", "start", "end", "strand"]
        columns = [1, 2, 5, 6, 7]
    elif motif_name in ["dyz19"]:
        header = ["#seq", "name", "strand", "start", "end"]
        columns = [0, 1, 2, 3, 4]
    elif motif_name in ["yq12"]:
        header = ["#seq", "name", "strand", "start", "end"]
        columns = [1, 2, 3, 4, 5]
    else:
        raise ValueError(f"Unknown motif: {motif_name}")

    df = pandas.read_csv(
        file_path, sep=",", header=None,
        names=header, usecols=columns,
        skiprows=1
    )

    # manual fixes for ... stuff
    df["name"] = df["name"].replace({"TSPY2": "TSPY"}, inplace=False)
    # my brain hurts...
    df["strand"] = df["strand"].replace({"PossibleInversion": ".", "Possible_Inversion": "."}, inplace=False)

    df["score"] = 1000
    df = df[["#seq", "start", "end", "name", "score", "strand"]]
    assert (df["strand"].isin(["+", "-", "."])).all(), f"malformed file / strand: {file_path.name} / {df}"
    assert (df["name"].isin(known_motif_names)).all(), f"malformed file / name: {file_path.name} / {df}"
    assert (df["start"] < df["end"]).all(), f"malformed file / coord: {file_path.name} / {df}"
    df.sort_values(["#seq", "start", "end"], inplace=True)
    return df
