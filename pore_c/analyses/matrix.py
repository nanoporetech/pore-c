from logging import getLogger
from pathlib import Path

import pandas as pd
from cooler import Cooler, annotate


logger = getLogger(__name__)


def correlate(x_cool: Cooler, y_cool: Cooler, xy_path: Path, coefficients_path: Path, resolution: int):
    metadata = {}

    for key, cool in (("x", x_cool), ("y", y_cool)):
        i = cool.info
        metadata[key] = {"contacts": int(i["sum"]), "non_zero_pixels": int(i["nnz"]), "num_bins": int(i["nbins"])}

    chrom_dtype = pd.CategoricalDtype(list(set(x_cool.chromnames).union(set(x_cool.chromnames))))
    xy_df = join_count(x_cool, y_cool).astype({"chrom1": chrom_dtype, "chrom2": chrom_dtype})
    logger.debug("Found {} non-zero overlapping bins".format(len(xy_df)))
    xy_df.to_parquet(xy_path)

    subsets = dict(
        [
            ("all", None),
            ("diagonal", "distance == 0"),
            ("trans", "distance == -1"),
            ("cis_off_diagonal", "distance > 0"),
            ("cis", "distance >= 0"),
        ]
    )
    res = []
    for label, query in subsets.items():
        if query:
            _df = xy_df.query(query)
        else:
            _df = xy_df
        num_bins = len(_df)
        for method in ["pearson", "spearman"]:
            logger.debug(f"Calculating {method} correlation coefficent for {label} contacts in {num_bins} bins")
            corr = _df[["x", "y"]].corr(method).loc["x", "y"]
            logger.info(f"{method} correlation coefficent for {label} contacts: {corr}")
            res.append({"subset": label, "corr_method": method, "coefficient": corr})
    coeff_df = (
        pd.DataFrame(res).assign(resolution=resolution).set_index(["subset", "corr_method"])["coefficient"].unstack()
    )
    logger.debug("Writing correlation coefficients to {}".format(str(coefficients_path)))
    coeff_df.to_csv(coefficients_path)
    metadata["pearson"] = {key: float(val) for key, val in coeff_df["pearson"].to_dict().items()}
    return metadata


def get_count_df(cool_mat):
    logger.debug(f"Creating counts dataframe for {cool_mat}")
    df = (
        annotate(cool_mat.pixels()[:], cool_mat.bins()[:]["chrom"])
        .eval("is_cis = (chrom1 == chrom2)")
        .pipe(get_distance, cool_mat.binsize)
        .set_index(["bin1_id", "bin2_id"])
        .sort_index()
    )
    num_bins = len(df)
    logger.info(f"Read {num_bins} bins from {cool_mat}")
    return df


def join_count(x_cool, y_cool):
    x_df = get_count_df(x_cool).rename(columns={"count": "x"})
    y_df = get_count_df(y_cool).rename(columns={"count": "y"})[["y"]]
    res = x_df.join(y_df, how="inner")
    return res


def get_distance(df, resolution):
    distance = df.eval(f"abs(bin1_id - bin2_id) * {resolution}")
    distance[~df.is_cis] = -1
    df["distance"] = distance
    return df
