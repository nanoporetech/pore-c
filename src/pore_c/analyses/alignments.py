import logging
from pore_c.datasources import NameSortedBamSource
import pandas as pd
import networkx as nx
import numpy as np
from tqdm import tqdm
import pyarrow as pa
from pyarrow import parquet as pq
from pysam import AlignmentFile
import intake
import yaml
from pathlib import Path
from intake.catalog.local import YAMLFileCatalog


logger = logging.getLogger(__name__)

FILTER_REASON_DTYPE = pd.CategoricalDtype(
    ['pass', 'unmapped', 'singleton', 'low_mq', 'overlap_on_read', "not_on_shortest_path"],
    ordered=True
)

class TableWriter(object):
    def __init__(self, path):
        self.path = path
        self._writer = None

    def write(self, df):
        table = pa.Table.from_pandas(df, preserve_index=False)
        if self._writer is None:
            self._writer = pq.ParquetWriter(self.path, schema=table.schema)
        self._writer.write_table(table)

    def close(self):
        self._writer.close()


def filter_alignments(input_bam: str, pass_bam: str = None, fail_bam: str = None, align_table: str = None, read_table: str = None, catalog_file: str = None, chunksize=50000):
    """Filter alignments to keep only alignments that contribute to contacts

    Parameters
    ----------

    input_bam : str
                Path to a namesorted bam with unfiltered alignments
    chunksize: int
                The alignments are batched for processing, this controls the batch size

    """
    source_aligns = NameSortedBamSource(input_bam, metadata={})

    pbar = tqdm(total=None, unit=" reads")
    writers = dict(
        pass_align = AlignmentFile(pass_bam, 'wb', template=AlignmentFile(input_bam)) if pass_bam else None,
        fail_align = AlignmentFile(fail_bam, 'wb', template=AlignmentFile(input_bam)) if fail_bam else None,
        align_table = TableWriter(align_table) if align_table else None,
        read_table = TableWriter(read_table) if read_table else None
    )

    running_total = None
    for chunk_idx, (aligns, align_df) in enumerate(source_aligns.read_chunked(chunksize=chunksize, yield_aligns=True)):

        # apply filters to the alignments
        res = filter_read_alignments(align_df)

        # calculate pre-read statistics on the result of the filtering
        read_stats = calculate_read_stats(res)

        if writers['pass_align']:
            [writers['pass_align'].write(aligns[x]) for x in res.index[res.pass_filter]]

        if writers['fail_align']:
            [writers['fail_align'].write(aligns[x]) for x in res.index[res.pass_filter == False].values]

        if writers['align_table']:
            writers['align_table'].write(res)

        if writers['read_table']:
            writers['read_table'].write(read_stats)

        pass_stats = res['reason'].value_counts(dropna=False)
        if chunk_idx == 0:
            running_total = pass_stats
        else:
            running_total += pass_stats
        counts = running_total.sort_index().to_frame().T
        percents = (100.0 * counts).div(counts.sum(axis=1), axis=0)
        pbar.set_postfix(percents.loc['reason', :].to_dict())
        pbar.update(chunksize)

        if chunk_idx == 1:
            break

    cat = create_catalog_yaml(writers, catalog_file)
    return cat


def create_catalog_yaml(writers, catalog_file):
    catalog = {
        "name": "pore_c_filtered_files",
        "description": "Output files of pore-c tools alignment filtering",
        "sources": {}
    }
    for file_key, writer in writers.items():
        if not writer:
            pass
        entry = {
            "args": {}
        }
        writer.close()
        if isinstance(writer, TableWriter):
            entry['args']['urlpath'] = "{{ CATALOG_DIR }}/" + str(writer.path.name)
            entry['driver'] = 'parquet'
        elif isinstance(writer, AlignmentFile):
            entry['args']['urlpath'] = "{{ CATALOG_DIR }}/" + Path(writer.filename.decode('utf8')).name
            entry['driver'] = 'pore_c.datasources.NameSortedBamSource'
        catalog['sources'][file_key] = entry

    with open(catalog_file, "w") as fh:
        fh.write(yaml.dump(catalog))
    return YAMLFileCatalog(str(catalog_file))


def calculate_read_stats(read_df):

    # total number of reads and how many aligments per read
    num_aligns = (
        read_df
        .groupby(["read_name"])[['read_length']]
        .max()
        .join(
            read_df
            .query("chrom != 'NULL'")
            .groupby(["read_name"])
            .size()
            .rename("num_aligns")
            .to_frame()
        )
        .fillna(0)
        .astype({"num_aligns": np.uint8})
    )

    # subset of alignments that passed filters
    pass_aligns = (
        read_df
        .query("pass_filter == True")
        .eval("perc_read_aligned = 100.0 * (read_end - read_start) / read_length")
    )

    # how many alignments pass filtering and how many unique chromosomes
    # do they hit
    read_stats = (
        pass_aligns
        .groupby(['read_name'])['chrom']
        .agg(['nunique', 'size'])
        .rename(columns={"nunique": "num_chroms", "size": "num_pass_aligns"})
    )

    # figure out what the "main" chromosome is and how much of the read hits there
    main_chrom = (
        pass_aligns
        .groupby(['read_name', 'chrom'])
        .agg({"perc_read_aligned": "sum", "read_end": "size"})
        .sort_values("perc_read_aligned", ascending=False)
        .groupby("read_name")
        .head(1)
        .reset_index()
        .rename(columns={
            "chrom": "main_chrom",
            "perc_read_aligned": "main_chrom_perc_read_aligned",
            "read_end": "main_chrom_num_aligns"}
        )
        .set_index(['read_name'])
    )

    # create a merged dataframe with one row per read
    res = (
        num_aligns.join(read_stats.join(main_chrom), how="outer")
    )

    unmapped = res['main_chrom'].isnull()
    res.loc[unmapped, 'main_chrom'] = 'NULL'
    res.loc[unmapped, 'num_chroms'] = 0
    res.loc[unmapped, 'num_pass_aligns'] = 0
    res.loc[unmapped, 'main_chrom_num_aligns'] = 0
    res.loc[unmapped, 'main_chrom_perc_read_aligned'] = 0.0
    res = (
        res
        .astype({
            "num_aligns": np.uint8,
            "num_chroms": np.uint8,
            "num_pass_aligns": np.uint8,
            "main_chrom_num_aligns": np.uint8
        })
    )
    return res


def filter_read_alignments(df, mapping_quality_cutoff=1):
    res = (
        df.assign(
            pass_filter=True,
            reason="pass"
        )
        .astype({'reason': FILTER_REASON_DTYPE})
    )
    # apply alignment-level filters
    unmapped_mask = res.mapping_type == 'unmapped'
    res.loc[unmapped_mask, "pass_filter"] = False
    res.loc[unmapped_mask, "reason"] = 'unmapped'

    # if not unmapped, but mapping quality below cutoff then fail
    fail_mq_mask = (~unmapped_mask & (df.mapping_quality <= mapping_quality_cutoff))
    res.loc[fail_mq_mask, "pass_filter"] = False
    res.loc[fail_mq_mask, "reason"] = 'low_mq'

    # no need to do other checks if nothing left
    if not res['pass_filter'].any():
        return res

    # for the remaining alignments filtering happens on a per-read basis
    by_read_res = (
        res.query("pass_filter == True")
        .groupby("read_name", sort=False, as_index=False)
        .apply(
            apply_per_read_filters
        )
    )
    res.update(by_read_res[['pass_filter', 'reason']])
    return res.astype({'reason': FILTER_REASON_DTYPE})



def apply_per_read_filters(read_df):
    return (read_df
        .pipe(filter_singleton)
        .pipe(filter_overlap_on_query)
        .pipe(filter_shortest_path)
    )
    return read_df


def filter_singleton(read_df):
    if len(read_df) == 1: # if you have a single alignment at this point you fail
        read_df.loc[:, 'pass_filter'] = False
        read_df.loc[:, 'reason'] = 'singleton'
    return read_df


def filter_overlap_on_query(read_df):
    overlaps_on_read = read_df.duplicated(subset = ['read_start', 'read_end'], keep=False)
    if overlaps_on_read.any():
        raise ValueError(read_df.loc[overlaps_on_read, :])
    return read_df

    query_endpoints = defaultdict(list)
    for index, key in enumerate([(a.query_alignment_start, a.query_alignment_end) for a in aligns]):
        query_endpoints[key].append(index)
    results = np.zeros(len(aligns), dtype=bool)
    for endpoints, indices in query_endpoints.items():
        if len(indices) > 1:
            tie_break = sorted([(aligns[i].get_tag("AS"), i) for i in indices], reverse=True)[0]
            results[tie_break[1]] = True
        else:
            results[indices[0]] = True
    return results


def minimap_gapscore(length, o1=4, o2=24, e1=2, e2=1):
    return min([o1 + int(length) * e1, o2 + int(length) * e2])


def bwa_gapscore(length, O=5, E=2):
    # O=5 E=2 is default for bwa bwasw
    # O=6 E=1 is default for bwa mem
    return O + length * E


def create_align_graph(aligns, gap_fn):
    # we'll visit the alignments in order of increasing endpoint on the read, need to keep
    # the ids as the index in the original list of aligns for filtering later
    aligns = (
        aligns[['read_start', 'read_end', 'read_length', 'score']]
        .copy()
        .sort_values(['read_end'])
    )
    node_ids = list(aligns.index)
    graph = nx.DiGraph()
    # initialise graph with root and sink node, and one for each alignment
    # edges in the graph will represent transitions from one alignment segment
    # to the next
    graph.add_nodes_from(['ROOT', 'SINK'] + node_ids)
    for align in aligns.itertuples():
        align_idx = align.Index
        align_score = align.score
        gap_penalty_start = gap_fn(align.read_start)
        graph.add_edge("ROOT", align_idx, weight=gap_penalty_start - align_score)
        gap_penalty_end = gap_fn(int(align.read_length - align.read_end))
        graph.add_edge(align_idx, "SINK", weight=gap_penalty_end)

    # for each pair of aligned segments add an edge
    for idx_a, align_idx_a in enumerate(node_ids[:-1]):
        align_a_end = aligns.at[align_idx_a, 'read_end']
        for align_idx_b in node_ids[idx_a+1:]:
            align_b_score = aligns.at[align_idx_b, 'score']
            align_b_read_start = aligns.at[align_idx_b, 'read_start']
            gap_penalty = gap_fn(abs(int(align_b_read_start) - int(align_a_end)))
            graph.add_edge(align_idx_a, align_idx_b, weight=gap_penalty - align_b_score)

    return graph



def filter_shortest_path(read_df, aligner="minimap2"):
    aligns = read_df[read_df.pass_filter]
    num_aligns = len(aligns)
    if num_aligns < 2:
        # can't build a graph, so nothing is filtered by this method
        return read_df

    if aligner == "minimap2":
        gap_fn = minimap_gapscore
    elif aligner == 'bwa':
        gap_fn = bwa_gapscore
    else:
        raise ValueError(f"Unrecognised aligner: {aligner}")
    graph = create_align_graph(aligns, gap_fn)
    distance, shortest_path = nx.single_source_bellman_ford(graph, "ROOT", "SINK")
    for idx  in aligns.index:
        if idx not in shortest_path:
            read_df.at[idx, 'pass_filter'] = False
            read_df.at[idx, 'reason'] = 'not_on_shortest_path'
    return read_df
