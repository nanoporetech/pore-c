from collections import defaultdict, Counter

import gzip
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def plot_contact_distances(EC_matrix_file_in: str, graph_file_out: str, ref_bin_file: str) -> None:

    chr_sizes =  {}

    bin_size = False
    for entry in gzip.open(ref_bin_file):
        l = entry.strip().split()
        if not bin_size:
            bin_size = int(l[2]) - int(l[1])

        if chr_sizes[l[0]] < l[3]:
            chr_sizes[l[0]] = l[3]

        min_max_size = len(min(chr_sizes.values()))

    bin_data = {}
    for entry in gzip.open(ref_bin_file):
        l = entry.strip().split()
        l[1] = int(l[1])
        l[2] = int(l[2])
        l[3] = int(l[3])
        bin_data[int(l[3])] = l

    data = Counter()

    for entry in open(EC_matrix_file_in):
        l = list(map(int,entry.strip().split()))
        bin1 = bin_data[l[0]]
        bin2 = bin_data[l[1]]
        if bin1[0] == bin2[0]:
            #intrachromosomal
            if bin1[3] - bin2[3] <= max_bin_size:
                data[bin_size * (bin1[3] - bin2[3])] += l[3]

    fig, ax = plt.subplots(1,figsize=(12,6))
    plt.subtitle("{} Contact Distance Distribution".format(EC_matrix_file_in))
    distances, counts = zip(*sorted(data.items(),key = lambda x: x[0]))
    distances_Mb = np.array(np.array(distances) / 10**6,dtype = int)
    ax.plot(distances_Mb,counts)
    ax.xlabel("distance (Mbp)", fontsize = "x-small")
    ax.ylabel("Contacts", fontsize = "x-small")
    ax.set_xlim(0,min_max_size * bin_size)
    fig.savefig(graph_file_out)


construction= """
def plot_corrected_contact_map(EC_matrix_file_in: str, heat_map_file_out: str, ref_bin_file: str) -> None:

    names = []
    markers = []
    lastChr = False
    size = 0
    for idx, entry in enumerate(gzip.open(ref_bin_file)):
        l = entry.strip().split()
        if not lastChr:
            lastChr = l[0]
        if lastChr != l[0]:
            markers.append(idx)
            names.append(lastChr)
        size = idx
    
    matrix = np.zeros((size,size))
    for entry in open(EC_matrix_file_in):
        l = entry.strip().split()
        l[0] = int(l[0])
        l[1] = int(l[1])

        matrix[


    fig, ax 
"""
