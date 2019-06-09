from typing import Optional

from collections import defaultdict, Counter

import numpy as np
import gzip
import matplotlib
matplotlib.use('agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

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
        l[2] = float(l[2])
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



def plot_contact_distances(ref_bin_file: str,  graph_file_out: str, *EC_matrix_files_in) -> None:

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

    data_set = {}
    
    for fn in EC_matrix_files_in:
        data = Counter()
        for entry in map(Matrix_Entry.from_string, open( fn )):
            bin1 = bin_data[entry.bin1]
            bin2 = bin_data[entry.bin2]
            if bin1[1] > bin2[1]:
                continue #forces removal of redundant entries
            if bin1[0] == bin2[0]:
            #intrachromosomal
                if bin2[3] - bin1[3] <= min_max_size:
                    data[bin_size * (bin2[3] - bin1[3])] += l[3]
        data_set[fn] = data

    colors = cm.rainbow(np.linspace(0, 1, len(data_set)))

    fig, ax = plt.subplots(1,figsize=(12,6))
    plt.title("Contact Distance Distribution")
    for idx, d in enumerate(data_set.items()):
        name, data = d
        distances, counts = zip(*sorted(data.items(),key = lambda x: x[0]))
        distances_Mb = np.array(np.array(distances) / 10**6,dtype = int)
        ax.plot(distances_Mb,counts, c = colors[idx])

    ax.set_xlabel("distance (Mbp)", fontsize = "x-small")
    ax.set_ylabel("Corrected counts", fontsize = "x-small")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1,min_max_size)
    ax.set_ylim(1000,max(counts))
    fig.savefig(graph_file_out)


def plot_contact_map(matrix_file_in: str,ref_bin_file: str, heat_map_file_out: str, matrix_type: Optional[str] = "raw") -> None:

    names = []
    markers = []
    lastChr = False
    for idx, entry in enumerate(gzip.open(ref_bin_file,'rt')):
        l = entry.strip().split()
        if not lastChr:
            lastChr = l[0]
        if lastChr != l[0]:
            markers.append(idx - 1 - .5)
            names.append(lastChr)

        lastChr = l[0]

    #tail entry
    size = idx
    markers.append(idx - 0.5)
    _markers = [0] + markers
    minor_markers = [ (x+y) / 2 for x,y in zip(_markers[:-1],_markers[1:])]
    names.append(l[0])


    matrix = np.zeros((size+1,size+1))
    for entry in map(Matrix_Entry.from_string, open(matrix_file_in)):
        return "{bin1} {bin2} {raw_counts}\n".format(bin1 = self.bin1,
                                                     bin2 = self.bin2,
                                                     raw_counts = self.raw_counts)

        if matrix_type == "corrected":
            matrix[entry.bin1,entry.bin2] = entry.corrected_counts
            matrix[entry.bin2,entry.bin1] = entry.corrected_counts
        elif matrix_type == "raw":
            matrix[entry.bin1,entry.bin2] = entry.raw_counts
            matrix[entry.bin2,entry.bin1] = entry.raw_counts
        elif matrix_type == "compare":
            matrix[entry.bin1,entry.bin2] = entry.raw_counts
            matrix[entry.bin2,entry.bin1] = entry.corrected_counts
        elif matrix_type == "contactprobability":
             matrix[entry.bin1,entry.bin2] = entry.contact_probability
             matrix[entry.bin2,entry.bin1] = entry.contact_probability

    fig, ax = plt.subplots(1,figsize= (12,6), dpi = 500)

    matrix[matrix < 2] = .1

#    plt.imshow(matrix,norm=colors.LogNorm(vmin=.1, vmax=matrix.max()), cmap="gist_heat_r")
    plt.imshow(matrix,norm=colors.LogNorm(vmin=.1, vmax=matrix.max()), cmap="viridis")


    null_markers = [""] * len(markers)
    ax.set_yticks(markers)
    ax.set_yticks(minor_markers, minor = True)
    ax.set_yticklabels(null_markers)
    ax.set_yticklabels(names, minor = True)
    ax.set_xticks(markers)
    ax.set_xticklabels(null_markers)
    ax.set_xticks(minor_markers, minor = True)
    ax.set_xticklabels(names, minor = True,rotation=90)

    ax.tick_params( axis="both", which="minor",labelsize= 'xx-small',length=0)
    ax.tick_params( axis="both", which="major",labelsize= 'xx-small',length=3)

    #TODO: chromosome names halfway between the major ticks

#    print("markers:",markers)
#    print("minor markers:",minor_markers)
#    print("names:",names)
#    print("size:",size)

    ax.vlines(markers,0,size, linestyle = ":", linewidth = .5, alpha=0.4, color = '#357BA1')
    ax.hlines(markers,0,size, linestyle = ":", linewidth = .5, alpha=0.4, color = '#357BA1')

    if matrix_type == "compare":
        ax.set_xlabel("corrected counts")
        ax.set_ylabel("raw counts")
        ax.yaxis.set_label_position("right")

    plt.colorbar()
    plt.savefig(heat_map_file_out)


def comparison_contact_map(matrix1_file_in: str,matrix2_file_in: str,ref_bin_file: str, heat_map_file_out: str, matrix_type: Optional[str] = "raw", normalise: Optional[bool] = True, chr_file: Optional[str] = "None" ) -> None:

    if chr_file != "None":

        chr_strands = {}
        names = []
        for entry in open(chr_file):
            #first col is contig, second is strandedness as either + or -
            l = entry.strip().split()
            names.append(l[0])
            chr_strands[l[0]] = l[1]

        sizes = {}
        intrachr_binranges = {}
        absolute_binranges = {}
        lastChr = False
        for idx, entry in enumerate(gzip.open(ref_bin_file,'rt')):
            l = entry.strip().split()
            if not lastChr:
                lastChr = l[0]
                intrachr_idx = 0
                continue
            if lastChr != l[0]:
                sizes[lastChr] = max( 1, intrachr_idx)
                intrachr_binranges[lastChr] = (0,intrachr_idx )
                absolute_binranges[lastChr] = (idx - intrachr_idx ,idx )
                intrachr_idx = 0
            intrachr_idx += 1
            lastChr = l[0]

        #tail chromosome info
        sizes[lastChr] = intrachr_idx
        intrachr_binranges[lastChr] = (0,intrachr_idx)
        absolute_binranges[lastChr] = (idx - intrachr_idx ,idx )

        #keys are absolute bin position in genome
        #values are the position in the selection given the ordering and orientating of the selected chromosomes
        bin_mappings = {}
        markers = []
        pos = 0
        for entry in names:
            markers.append(pos + sizes[entry] - 0.5)
            start,stop = intrachr_binranges[entry]
            abs_start, abs_stop = absolute_binranges[entry]
            print(pos +start,pos+stop,abs_start,abs_stop)

            if chr_strands[entry] == "+":
                for x1, x2 in zip(  list(range(*absolute_binranges[entry])) , list(range(pos + start, pos + stop))):
                    bin_mappings[x1] = x2
            else:
                for x1, x2 in zip(  list(range(*absolute_binranges[entry])) , list(range(pos + start, pos + stop))[::-1]):
                    bin_mappings[x1] = x2
            pos += sizes[entry]

        size = pos + 1

        _markers = [0] + markers
        minor_markers = [ (x+y) / 2 for x,y in zip(_markers[:-1],_markers[1:])]


    else:

        names = []
        markers = []
        lastChr = False
        for idx, entry in enumerate(gzip.open(ref_bin_file,'rt')):
            l = entry.strip().split()
            if not lastChr:
                lastChr = l[0]
            if lastChr != l[0]:
                markers.append(idx - 1 - .5)
                names.append(lastChr)
            lastChr = l[0]

        markers.append(idx - 0.5)
        _markers = [-.5] + markers
        minor_markers = [ (x+y) / 2 for x,y in zip(_markers[:-1],_markers[1:])]
        names.append(l[0])
        size = idx + 1

######construct contact matrix
    matrix = np.zeros((size,size )) 

    upper_sum = 0
    for entry in map(Matrix_Entry.from_string, open(matrix1_file_in)):
        if chr_file != "None":
            if entry.bin1 in bin_mappings and entry.bin2 in bin_mappings:
                x = bin_mappings[entry.bin1]
                y = bin_mappings[entry.bin2]
                if matrix_type == "corrected":
                    matrix[x,y] = entry.corrected_counts
                    upper_sum += entry.corrected_counts
                elif matrix_type == "raw":
                    matrix[x,y] = entry.raw_counts
                    upper_sum += entry.raw_counts
            else:
                continue

        else:
            if matrix_type == "corrected":
                matrix[entry.bin1,entry.bin2] = entry.corrected_counts
                upper_sum += entry.corrected_counts
            elif matrix_type == "raw":
                matrix[entry.bin1,entry.bin2] = entry.raw_counts
                upper_sum += entry.raw_counts


    lower_sum = 0
    for entry in map(Matrix_Entry.from_string, open(matrix2_file_in)):
        if chr_file != "None":
            if entry.bin1 in bin_mappings and entry.bin2 in bin_mappings:
                x = bin_mappings[entry.bin1]
                y = bin_mappings[entry.bin2]
                if matrix_type == "corrected":
                    matrix[y,x] = entry.corrected_counts
                    lower_sum += entry.corrected_counts
                elif matrix_type == "raw":
                    matrix[y,x] = entry.raw_counts
                    lower_sum += entry.raw_counts
            else:
                continue
        else:
            if matrix_type == "corrected":
                matrix[entry.bin2,entry.bin1] = entry.corrected_counts
                lower_sum += entry.corrected_counts
            elif matrix_type == "raw":
                matrix[entry.bin2,entry.bin1] = entry.raw_counts
                lower_sum += entry.raw_counts

    if normalise:

        print('normalising two datasets.')
        print('lower factor: X  / {}'.format(lower_sum))
        print('upper factor: X  / {}'.format(upper_sum))

        if upper_sum > lower_sum:
            factor = lower_sum / upper_sum
            flag = True
        else:
            factor = upper_sum / lower_sum
            flag = False
        for x in range(size):
            for y in range(size):
                if x  == y:
                    continue
                elif x > y: 
                    if flag:
                        matrix[x,y] = matrix[x,y] / factor
                else:
                    if not flag:
                        matrix[x,y] = matrix[x,y] / factor

    #save numpy matrix before setting diagonal to zero in order to preserve
    # self contacts of small molecules
    np.save(heat_map_file_out.replace('.png','.npy'),matrix)                        
    
    fig, ax = plt.subplots(1,figsize= (12,12), dpi = 1000)

#    plt.imshow(matrix,norm=colors.LogNorm(vmin=1, vmax=matrix.max()), cmap="viridis")
    plt.imshow(matrix,norm=colors.LogNorm(vmin=1, vmax=matrix.max()), cmap="gist_heat_r")
#    plt.imshow(matrix, cmap="gist_heat_r")


    null_markers = [""] * len(markers)
    ax.set_yticks(markers)
    ax.set_yticks(minor_markers, minor = True)
    ax.set_yticklabels(null_markers)
    ax.set_yticklabels(names, minor = True)
    ax.set_xticks(markers)
    ax.set_xticklabels(null_markers)
    ax.set_xticks(minor_markers, minor = True)
    ax.set_xticklabels(names, minor = True,rotation=90)

    ax.set_ylabel(matrix1_file_in)
    ax.set_xlabel(matrix2_file_in)
    ax.yaxis.set_label_position("right")    

    ax.tick_params( axis="both", which="minor",labelsize= 'xx-small',length=0)
    ax.tick_params( axis="both", which="major",labelsize= 'xx-small',length=3)

    ax.vlines(markers,-0.5, size - 1, linestyle = ":", linewidth = .5, alpha=0.4, color = '#357BA1')
    ax.hlines(markers,-0.5, size - 1, linestyle = ":", linewidth = .5, alpha=0.4, color = '#357BA1')


    ax.set_xlim(right=size-1.5, left = -.5)
    ax.set_ylim(bottom=size-1.5, top =-.5)

#    plt.colorbar()
    plt.savefig(heat_map_file_out)



#this is done on corrected values
def cis_trans_analysis(EC_matrix_file_in: str, ref_bin_file: str, data_file_out:str, results_file_out: str, scatter_map_file_out: str ) -> None:

    #coordinate the bins with their chromosomes based on the 
    chrs = {}
    for entry in gzip.open(ref_bin_file,'rt'):
        l = entry.strip().split()
        chrs[int(l[3])] = l[0]

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
        
        if chrs[l[0]] == chrs[l[1]]:
            intra[l[0]] += c_count
        else:
            inter[l[0]] += c_count


    shared = sorted(list(set(intra.keys()).union(set(inter.keys()))))

    print("data sizes:")
    print("intra:",len(intra))
    print("inter:",len(inter))
    print("shared data coordinates:",len(shared))
    fig, ax = plt.subplots(1,figsize= (12,6))
    intra_values = np.array([intra[x] if x in intra else 0.0 for x in shared])
    inter_values = np.array([inter[x] if x in inter else 1.0 for x in shared])

    print("max intra:",max(intra_values))
    print("max inter:",max(inter_values))
    ratios = intra_values / inter_values
    plt.hexbin(intra_values,inter_values,gridsize = (100,100))

    ax.set_xlabel("cis counts")
    ax.set_ylabel("trans counts")
#    ax.set_yscale('log')
#    ax.set_xscale('log')
    ax.set_ylim(bottom = 0, top = 500)
    ax.set_xlim(left= 0, right = 500)

    plt.savefig(scatter_map_file_out)

    f_out = open(data_file_out,'w')
    for idx, val in enumerate( shared):
        f_out.write("{}\t{}\n".format(val, ratios[idx]))

    f_out.close()

    f_out = open(results_file_out,'w')
    f_out.write("intra,inter,ratio,pct_intra\n")
    x = np.sum(intra_values)
    y = np.sum(inter_values)

    ratio = x / y

    f_out.write("{},{},{},{}\n".format(x, y, ratio, x / (x+y)))
    f_out.close()
    

####
#plots a log-log scatter plot of point-matched raw count values and calculates the pearson correlation coefficient for that distribution.
def matrix_correlation(matrix1_file_in: str, matrix2_file_in: str, plot_out: str, result_out:str, contact_data_out: str,matrix_type: Optional[str] = "raw") -> None:
    matrix1_data = {}
    matrix2_data = {}
    for entry in map(Matrix_Entry.from_string, open(matrix1_file_in)):
        if matrix_type == "raw":
            matrix1_data[(entry.bin1,entry.bin2)] = entry.raw_counts
        if matrix_type == "corrected":
            matrix1_data[(entry.bin1,entry.bin2)] = entry.corrected_counts

    for entry in map(Matrix_Entry.from_string, open(matrix2_file_in)):
        if matrix_type == "raw":
            matrix2_data[(entry.bin1,entry.bin2)] = entry.raw_counts
        if matrix_type == "corrected":
            matrix2_data[(entry.bin1,entry.bin2)] = entry.corrected_counts

    matrix1_nonzero = set(matrix1_data.keys())
    shared_nonzero = list(matrix1_nonzero.intersection(set(matrix2_data.keys())))

    shared_matrix1 = [matrix1_data[x] for x in shared_nonzero]
    shared_matrix2 = [matrix2_data[x] for x in shared_nonzero]

    f_out = open(contact_data_out,'w')
    f_out.write("bin1,bin2,val1,val2\n")
    for entry in shared_nonzero:
        bin1, bin2 = entry
        f_out.write("{bin1},{bin2},{val1},{val2}\n".format(bin1 = bin1, bin2= bin2, val1 = matrix1_data[bin1,bin2],val2 = matrix2_data[bin1,bin2]))

    f_out.close()

    r, p = pearsonr(shared_matrix1, shared_matrix2)
    fig, ax = plt.subplots(1,figsize= (12,6))

    plt.plot(shared_matrix1,shared_matrix2, 'b,')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel("{} distances (bp)".format(matrix1_file_in))
    ax.set_ylabel("{} distances (bp)".format(matrix2_file_in))
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


