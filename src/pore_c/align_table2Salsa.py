#!/usr/bin/env python

from pathlib import Path
import pandas as pd
from itertools import combinations
import sys
import numpy as np


def to_salsa(df):
    colname = {'chrom': str,
               'start':int,
               'end':int,
               'readid':str,
               'align_score':int,
               'strand':str}
    
    res = []
    keep_segments = (
        df.query("pass_filter == True").replace({"strand": {True: "+", False: "-"}})
        # convert strand to +/- and chrom to string for lexographical sorting to get upper-triangle
        .astype({"strand": colname["strand"], "chrom": colname['chrom']})
    )
    for x, (read_idx, read_df) in enumerate(keep_segments.groupby("read_idx", as_index=False)):
        if len(read_df) <= 1:
            continue
        rows = list(
            read_df.sort_values(["read_start"], ascending=True)
            .itertuples()
        )

        for combi,(pos1, pos2) in enumerate(combinations(rows, 2)):
            if pos1.align_idx == pos2.align_idx:
                raise ValueError
            for rec, paired_id in zip([pos1, pos2], [1,2]):
                new_read_id = f'{rec.read_name}_{combi}/{paired_id}'
                yield(to_salsa_bed_tuple(rec,new_read_id))


def to_salsa_bed_tuple(pos,read_id):
    return (
        pos.chrom,
        pos.start,
        pos.end,
        read_id,
        pos.mapping_quality,
        pos.strand
    )



def main():
    #fn = '/data/20190830_porec_anaylsis/align_table/GM12878_NlaIII_native_20190813_972fa_guppy_3.0.5.hac_new123878assm_bwa_0.7.17_sensitive.at.alignment.parquet'
    fn = sys.argv[1]
    fout = sys.argv[2]
    columns=['align_idx',"chrom", 'read_start','start','end', "read_name",'mapping_quality','strand',"pass_filter", 'read_idx']

    df = pd.read_parquet(fn,columns=columns)
    salsa = to_salsa(df)
    with open(fout, 'w') as bedout:
        for n,line in enumerate(salsa):
            print(f'{n:,}', end='\r', flush=True)
            bedout.write('{}\n'.format('\t'.join(map(str,line))))

if "__main__" == __name__:
    main()
