from typing import Optional

from collections import defaultdict, Counter

import numpy as np
import gzip
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.stats import pearsonr
from dataclasses import dataclass

@dataclass
class Matrix_Entry:
    bin1: int
    bin2: int
    raw_counts: float

    contact_probability: float = -1.0
    corrected_counts: float = -1.0
    E1: float = 0.0
    E2: float = 0.0

    @classmethod
    def from_string(cls, entry):
        l = entry.strip().split()
        for x in range(2):
            l[x] = int(l[x])
        if len(l) > 3:
            for x in range(2,len(l)):
                l[x] = float(l[x])
        return cls(*l)

    def to_string(self):
        if self.contact_probability != 0.0 and self.corrected_counts != 0.0 and self.E1 != 0.0 and self.E2 != 0.0:
            return "{bin1} {bin2} {raw_counts} {contact_probability} {corrected_counts} {E1} {E2}\n".format(bin1 = self.bin1,
                                                                                                          bin2 = self.bin2,
                                                                                                          raw_counts = self.raw_counts,
                                                                                                          contact_probability = self.contact_probability,
                                                                                                          corrected_counts = self.corrected_counts,
                                                                                                          E1 = self.E1,
                                                                                                          E2 = self.E2)
        else:
            return "{bin1} {bin2} {raw_counts}\n".format(bin1 = self.bin1,
                                                       bin2 = self.bin2,
                                                       raw_counts = self.raw_counts)

def plot_contact_distances(EC_matrix_file_in: str,ref_bin_file: str,  graph_file_out: str) -> None:

    chr_sizes =  {}

    bin_size = -1
    for entry in gzip.open(ref_bin_file,'rt'):
        l = entry.strip().split()
        if l in ["M"]:
            continue
        l[3] = int(l[3])
        if bin_size == -1:
            bin_size = int(l[2]) - int(l[1])
        if l[0] not in chr_sizes:
            chr_sizes[l[0]] = l[3]
        elif chr_sizes[l[0]] < l[3]:
            chr_sizes[l[0]] = l[3]

    min_max_size = min(chr_sizes.values())
    print('distance cap: {}'.format(min_max_size))

    bin_data = {}

    for entry in gzip.open(ref_bin_file):
        l = entry.strip().split()
        l[1] = int(l[1])
        l[2] = int(l[2])
        l[3] = int(l[3])
        bin_data[int(l[3])] = l

    data = Counter()

    for entry in map(Matrix_Entry.from_string, open(EC_matrix_file_in)):
        bin1 = bin_data[entry.bin1]
        bin2 = bin_data[entry.bin2]
        if bin1[1] > bin2[1]:
            continue #forces removal of redundant entries
        if bin1[0] == bin2[0]:
            #intrachromosomal
            if bin2[3] - bin1[3] <= min_max_size:
                data[bin_size * (bin2[3] - bin1[3])] += l[3]

    fig, ax = plt.subplots(1,figsize=(12,6))
    plt.title("{} Contact Distance Distribution".format(EC_matrix_file_in))
    distances, counts = zip(*sorted(data.items(),key = lambda x: x[0]))
    print(distances[:100])
    print(counts[:100])

    distances_Mb = np.array(np.array(distances) / 10**6,dtype = int)
    print(distances_Mb[:100])    
    ax.plot(distances_Mb,counts, 'k:')
    ax.set_xlabel("distance (Mbp)", fontsize = "x-small")
    ax.set_ylabel("Contacts", fontsize = "x-small")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(0,min_max_size)
    ax.set_ylim(100,max(counts))
    fig.savefig(graph_file_out)


def plot_contact_map(matrix_file_in: str,ref_bin_file: str, heat_map_file_out: str, matrix_type: Optional[str] = "raw") -> None:

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
    
    matrix = np.zeros((size+1,size+1))
    for entry in map(Matrix_Entry.from_string, open(matrix_file_in)):
#        try:
#            assert entry.contact_probability != -1.0
#        except:
#            raise ValueError ("This matrix file has not been balanced.")
#        if entry.bin1 == entry.bin2:
#            #the diagonal is never informative and only serves to scale down the rest of the data in the colorspace
            continue 
        if matrix_type == "corrected":
            matrix[entry.bin1,entry.bin2] = entry.corrected_counts
            matrix[entry.bin2,entry.bin1] = entry.corrected_counts
        elif matrix_type == "raw":
            matrix[entry.bin1,entry.bin2] = entry.raw_counts
            matrix[entry.bin2,entry.bin1] = entry.raw_counts
        elif matrix_type == "compare":
            matrix[entry.bin1,entry.bin2] = entry.raw_counts
            matrix[entry.bin2,entry.bin1] = entry.corrected_counts

    fig, ax = plt.subplots(1,figsize= (12,6))

    plt.imshow(matrix,norm=colors.LogNorm(vmin=1, vmax=matrix.max()), cmap="viridis")

    #TODO: chromosome names halfway between the major ticks
    ax.set_yticks(markers, names)
    ax.set_xticks(markers, names)
    if matrix_type == "compare":
        ax.set_xlabel("corrected counts")
        ax.set_ylabel("raw counts")

    ax.yaxis.set_label_position("right")
    plt.savefig(heat_map_file_out)


#this is done on corrected values
def cis_trans_analysis(EC_matrix_file_in: str, ref_bin_file: str, data_file_out:str, results_file_out: str, scatter_map_file_out: str ) -> None:

    #coordinate the bins with their chromosomes based on the 
    chrs = {}
    for entry in gzip.open(ref_bin_file,'rt'):
        l = entry.strip().split()
        chrs[l[3]] = l[0]

    #key these on matrix bin with values equal to count
    intra = Counter()
    inter = Counter()
    for entry in open(EC_matrix_file_in):
        l = entry.strip().split()
        l[0] = int(l[0])
        l[1] = int(l[1])
        if l[0] > l[1]:
            continue
        c_count = float(l[2]) # raw counts
        #c_count = float(l[2]) #corrected counts
        
        if chrs[l[0]] == chrs[l[1]]:
            intra[l[0]] += c_count
        else:
            inter[l[0]] += c_count


    shared = sorted(list(set(intra.keys()).intersection(set(inter.keys()))))
    print("data sizes:")
    print("intra:",len(intra))
    print("inter:",len(inter))
    print("shared data coordinates:",len(shared))
    fig, ax = plt.subplots(1,figsize= (12,6))
    intra_values = np.array([intra[x] for x in shared])
    inter_values = np.array([inter[x] for x in shared])

    print("max intra:",max(intra_values))
    print("max inter:",max(inter_values))
    ratios = intra_values / inter_values
    plt.hexbin(intra_values,inter_values,gridsize = (200,200))

    ax.set_xlabel("cis counts")
    ax.set_ylabel("trans counts")
#    ax.set_yscale('log')
#    ax.set_xscale('log')

    plt.savefig(scatter_map_file_out)

    f_out = open(data_file_out,'w')
    for idx, val in enumerate( shared):
        f_out.write("{}\t{}\n".format(val, ratios[idx]))

    f_out.close()

    f_out = open(results_file_out,'w')
    f_out.write("intra,inter,ratio\n")
    x = np.sum(intra_values)
    y = np.sum(inter_values)

    ratio = x / y
    print(x,y, ratio)

    f_out.write("{},{},{}".format(x, y, ratio))
    f_out.close()
    

####
#plots a log-log scatter plot of point-matched raw count values and calculates the pearson correlation coefficient for that distribution.
def matrix_correlation(matrix1_file_in: str, matrix2_file_in: str, plot_out: str, result_out:str) -> None:
    matrix1_data = {}
    matrix2_data = {}
    for entry in map(Matrix_Entry.from_string, open(matrix1_file_in)):
        matrix1_data[(entry.bin1,entry.bin2)] = entry.raw_counts

    for entry in map(Matrix_Entry.from_string, open(matrix2_file_in)):
        matrix2_data[(entry.bin1,entry.bin2)] = entry.raw_counts

    matrix1_nonzero = set(matrix1_data.keys())
    shared_nonzero = list(matrix1_nonzero.intersection(set(matrix2_data.keys())))

    shared_matrix1 = [matrix1_data[x] for x in shared_nonzero]
    shared_matrix2 = [matrix2_data[x] for x in shared_nonzero]


    r, p = pearsonr(shared_matrix1, shared_matrix2)
    fig, ax = plt.subplots(1,figsize= (12,6))

    plt.plot(shared_matrix1,shared_matrix2, 'b,')
    ax.set_yscale('log')
    ax.set_xscale('log')

    plt.savefig(plot_out)

    f_out = open(result_out,'w')

    header = "matrix1,matrix2,shared_points,pearson_coeff,p_val\n"
    form = "{matrix1},{matrix2},{shared_points},{pearson_coeff},{p_val}\n"
    f_out.write(header)
    f_out.write(form.format(matrix1= matrix1_file_in,
                            matrix2 = matrix2_file_in,
                            shared_points = len(shared_nonzero),
                            pearson_coeff = r,
                            p_val = p)
            )
    f_out.close()


