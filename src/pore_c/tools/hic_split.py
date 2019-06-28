def split_hic_data(input_hictxt: str, output_hictxt_prefix: str, output_inter_hictxt: str) -> None:
    inter_out = open(output_inter_hictxt, "w")
    seen = set()
    # because the files are sort-ordered by chromosome, they can mostly be traversed
    current_chr = False
    for entry in open(input_hictxt):
        l = entry.strip().split()
        if l[2] != l[6]:
            inter_out.write(entry)

        elif l[2] == l[6] and current_chr == False:
            current_chr = l[2]
            seen.add(l[2])
            intra_out = open("{filebase}.{}.hic.txt", "w")
            intra_out.write(entry)
        elif l[2] == l[6] and current_chr == l[2]:
            intra_out.write(entry)
        elif l[2] == l[6] and current_chr != l[2] and l[2] in seen:
            raise ValueError("hictxt file {} is not properly ordered.".format(input_hictxt))
        elif l[2] == l[6] and current_chr != l[2] and l[2] not in seen:
            intra_out.close()
            intra_out = open("{filebase}.{}.hic.txt", "w")
    inter_out.close()
