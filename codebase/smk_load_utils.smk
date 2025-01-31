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
