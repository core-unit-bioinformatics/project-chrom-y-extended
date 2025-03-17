import pathlib
import xopen
import datetime as dt


def read_set_listing(file_path):
    with xopen.xopen(file_path) as listing:
        items = set(listing.read().strip().split())
    return items


def get_sample_file(sample_sheet, sample, file_key):

    if "-" in sample:
        sample, hg = sample.split("-")
    file_path = sample_sheet.loc[sample, file_key]
    _ = pathlib.Path(file_path).resolve(strict=True)
    return file_path


def get_timestamp():
    ts = dt.datetime.today().strftime("%Y%m%dT%H%M")
    return ts


_REVCOMP_CHAR_MAP = {
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "n": "n",
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N"
}

_REVCOMP_TABLE = str.maketrans(_REVCOMP_CHAR_MAP)


def revcomp(sequence):
    return sequence.translate(_REVCOMP_TABLE)[::-1]

assert revcomp("ACTTG") == "CAAGT"
