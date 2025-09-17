

import pathlib


SUB_WD = WD.joinpath("30-color-blocks-uniqk")

GLOBUS_SHARE_BED_INPUT_DIR = pathlib.Path(
    PROJECT_CONFIG["remote_hilbert_prefix"],
    PROJECT_CONFIG["kmer_color_blocks_uniqk"]
).resolve(strict=True)


localrules: unzip_color_unique_kmer_counts
rule unzip_color_unique_kmer_counts:
    """2025-09-17 - DEPRECATED
    update after dataset freeze: these annotation data
    are now shared via Globus directly - read BED input
    in the subsequent rules from the share location on
    the file system
    """
    input:
        zip = PROJECT_REPO_ROOT.joinpath(
            "annotation", "raw", "20250707_Peter_Bed_Files_ColorBlocks.ML.zip",
        ).resolve(strict=True)
    output:
        bed_dir = directory(WD.joinpath("region_kmer_annot", "color_unique_kmers"))
    shell:
        "exit 1"
        "unzip -o -j -d {output.bed_dir} {input.zip} "


localrules: normalize_color_unique_kmer_counts
rule normalize_color_unique_kmer_counts:
    input:
        #bed_dir = rules.unzip_color_unique_kmer_counts.output.bed_dir
        bed_dir = GLOBUS_SHARE_BED_INPUT_DIR
    output:
        tsv_dir = directory(SUB_WD.joinpath("kmer_annotation", "norm_tsv"))
    run:
        import hashlib as hl
        import numpy as np
        import pandas as pd
        import pathlib as pl
        import sys
        np.seterr(divide = 'raise')

        for bed_file in pl.Path(input.bed_dir).glob("*.bed"):
            # the format is of course idiosyncratic
            file_name = bed_file.name
            sample = file_name.split("_")[0]
            if sample not in SAMPLES:
                sys.stderr.write(f"\nWarning: skipping file - likely not a sample file: {file_name}\n")
                continue  # skip files that do not match any sample name
            df = pd.read_csv(
                bed_file,
                sep="\t",
                header=0,
            )
            df["region_id"] = df.apply(
                lambda row: hl.md5(
                    f"{row['#Contig']}{row['Start']}{row['End']}".encode("utf-8")
                ).hexdigest(),
                axis=1
            )
            df["kmerColor"] = df["kmerColor"].apply(normalize_label_name)

            # following: rank the possible labels by
            # their unique k-mer count and compute some
            # proxy for the confidence in the assignment
            # (enrichment = km(best label) / km(second-best label)
            ext_df = []
            for region, labels in df.groupby("region_id"):
                if labels.shape[0] == 1:
                    # only one label to choose from,
                    # take it as best hit, set enrichment
                    # to factor 100
                    ext_df.append(
                        (labels.index[0], labels["kmerColor"].iloc[0], 100)
                    )
                    continue
                # note 'ascending=false': rank 1 is the best hit, rank 2 second and so on
                ranks = labels["totalUniqueKmerHits"].rank(method="dense", ascending=False)
                top_rank = ranks.argmin()
                last_rank = ranks.argmax()
                if top_rank == last_rank:
                    # Several potential labels, but all have the same
                    # unique k-mer count and thus, no best match can be
                    # assigned here. Note that ML solves these cases by
                    # looking at all shared k-mers and then assigns the
                    # label that maximizes the number of shared k-mers.
                    # But we don't have that information available here.
                    ext_df.extend(
                        [(idx, "uncertain", 0) for idx in labels.index]
                    )
                else:
                    # at least two different ranks implies there is a best match
                    top_kmer = labels.loc[ranks.idxmin(), "totalUniqueKmerHits"]
                    assigned_label = labels.loc[ranks.idxmin(), "kmerColor"]
                    # note here: several labels could have been ranked as two,
                    # but all of them must have the same unique k-mer count,
                    # hence arbitrarily taking the first does not affect the
                    # computation
                    second_rank = ranks.index[ranks == 2][0]
                    second_kmer = labels.loc[second_rank, "totalUniqueKmerHits"]
                    try:
                        top_enrich = min(100, round(top_kmer / second_kmer, 2))
                    except FloatingPointError:
                        top_enrich = 100
                    ext_df.extend(
                        [(idx, assigned_label, top_enrich) for idx in labels.index]
                    )
            ext_df = pd.DataFrame.from_records(ext_df, columns=["index", "assigned_label", "top_enrich"])
            ext_df.set_index("index", inplace=True)

            df = df.join(ext_df)
            df.drop("region_id", axis=1, inplace=True)
            df.sort_values(["#Contig", "Start", "totalUniqueKmerHits"], inplace=True, ascending=[True, True, False])

            out_filename = f"{sample}_chrY-regions_uniq-kmer-counts.tsv"
            out_filepath = pl.Path(output.tsv_dir).joinpath(out_filename)
            out_filepath.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(out_filepath, sep="\t", index=False, header=True)
    # END OF RUN BLOCK


rule run_all_colorblocks_uniqk:
    input:
        norm_dir = rules.normalize_color_unique_kmer_counts.output.tsv_dir
