import pysam

file_in = pysam.AlignmentFile(snakemake.input[0])
file_out = open(snakemake.output[0],'w')

template = "{ch}\t{t_st}\t{t_en}\t{read_id}\t{q_st}\t{q_en}\t{mapq}\t{strand}\n"

for align in file_in:
    file_out.write(template.format(ch=align.reference_name,t_st = align.reference_start, t_en = align.reference_end,
                                   read_id=align.query_name, q_st = align.query_alignment_start, q_en = align.query_alignment_end,
                                   mapq = align.mapq, strand = 16 if align.is_reverse else 0))

file_out.close()
