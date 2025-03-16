import pathlib
import json
import subprocess

import pandas


def find_project_repository_root(starting_point=None):
    """Temp/workaround until it is clear how/where in
    Snakemake 8 basic information such as the location
    of the Snakefile is accessible via the new API
    """
    if starting_point is None:
        # unlikely to work, but we can try
        starting_point = pathlib.Path(".").resolve()

    git_cmd = ["git", "rev-parse", "--show-toplevel"]
    output = subprocess.check_output(git_cmd, cwd=starting_point)
    output = output.decode("utf-8").strip()
    repo_root = pathlib.Path(output).resolve(strict=True)

    return repo_root


def load_sample_sheet(sample_sheet_name, repo_root):

    if not sample_sheet_name.endswith(".tsv"):
        sample_sheet_name += ".tsv"

    sample_sheet_path = repo_root.joinpath(
        "samples", sample_sheet_name
    ).resolve(strict=True)

    sample_sheet = pandas.read_csv(
        sample_sheet_path, sep="\t",
        comment="#", header=0
    )

    sample_sheet.set_index("sample", inplace=True)

    return sample_sheet



PROJECT_REPO_ROOT = find_project_repository_root(config.get("cwd", None))

GLOBAL_CONDA_ENVS = PROJECT_REPO_ROOT.joinpath("codebase", "global_envs").resolve(strict=True)
GLOBAL_SCRIPTS = PROJECT_REPO_ROOT.joinpath("codebase", "global_scripts").resolve(strict=True)

SAMPLE_SHEETS = {
    "vrk_assm": load_sample_sheet("verkko_assemblies", PROJECT_REPO_ROOT),
    "vrk_chrom": load_sample_sheet("verkko_sex_chrom", PROJECT_REPO_ROOT)
}

SUBFOLDERS = {
    "extract-sex-chromosomes": "00-extract",
    "process-region-annotation": "10-process-regions"
}

PROJECT_CONFIG_JSON = PROJECT_REPO_ROOT.joinpath("codebase", "project-config.json").resolve(strict=True)

PROJECT_CONFIG = json.load(open(PROJECT_CONFIG_JSON, "r"))["GENERIC"]
