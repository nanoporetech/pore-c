from typing import Pattern, List, NamedTuple, Iterator, Union, Dict, Tuple, Optional

from dataclasses import dataclass

import scipy.stats as st
import numpy as np
import gzip

@dataclass
class HiCMap(object):
    def __init__(self, bin_ref):
        self.bin_count = 0

        for entry in gzip.open(bin_ref):
            self.bin_count += 1

        self.sparse_data = {}
#        self.matrix = np.zeros((self.bin_count,self.bin_count),dtype = float)
        self.matrix = np.ma.identity(self.bin_count,dtype = float)
        self.mask = None
        self.cP = None #np.copy(self.matrix)
        self.min_val = np.min(self.matrix)
        self.max_val = np.max(self.matrix)
        self.total_contacts = 0

    def add_datum(self,x,y, value):
        if x == y:
            #protect from losing total support
            if self.matrix[x,y] == 1:
                self.matrix[x,y] = value
            else:
                self.matrix[x,y] += value
        else:
            self.matrix[x,y] += value
            self.matrix[y,x] = self.matrix[x,y]

    def from_raw_matrix_file(self, matrix_file):
        for entry in open(matrix_file):
            l = entry.strip().split()
            l[0] = int(l[0])
            l[1] = int(l[1])
            l[2] = int(l[2])
            if l[0] > l[1]:
                self.add_datum(l[1], l[0], l[2])
            else:
                self.add_datum(*l)
            

        symmetry_test = set()

        for entry in open(matrix_file):
            l = entry.strip().split()
            try:
                assert len(l) == 3
            except:
                continue
            l[0] = int(l[0])
            l[1] = int(l[1])

            ts = tuple(sorted([l[0],l[1]]))
            if ts in symmetry_test:
                raise ValueError("This data set is not contact commutative. There are contacts such that there is an AB contact value that is not equal to BA contact value.")
            symmetry_test.add(ts)
            l[2] = float(l[2])
            self.sparse_data[(l[0],l[1])] = l
            if l[0] == l[1]:
                ###very clever cheating is happening here: force the matrix to have total support by including the notion that
                #  every bin on the genome is in contact with itself (i.e., add the constituted contact matrix with
                #  a diagonal matrix whose value is all 1). >:)
                self.matrix[l[0],l[1]] = np.max([self.matrix[l[0],l[1]],l[2]])
            else:
                self.matrix[l[0],l[1]] = self.matrix[l[1],l[0]] = l[2] 

            #test for matrix square-ness
 #           try:
 #               assert self.matrix.ndim == 2
 #               assert self.matrix.shape[0] == self.matrix.shape[1]
 #           except:
 #               raise ValueError("Contact matrix is not square.")

            #confirm that the resulting matrix is strictly positive matrix 
            # total support is tested in the set_thresholds step
            try:
                assert np.all(self.matrix >= 0)
            except:
                raise ValueError("Contact matrix is not strictly positive.")


    def set_thresholds(self, file_out, ci=0.9999, remove_zeros = True):
        """
        This protocol establishes a mask for the contact matrix
        such that values above the confidence interval are masked
        since they can underskew unsaturated data during correction.

        Additionally, low values are masked as their prominence 
        can be unduly amplified, and zero value
        """

        # mask zeros
        if remove_zeros:
            np.ma.masked_where(self.matrix == 0, self.matrix)
        
        raw_values_mean = np.median(np.diag(self.matrix))
        print(raw_values_mean)
        print( st.t.interval(ci,raw_values_mean))
        #interval(alpha,df,loc=0,scale=1)
        self.min_val, self.max_val = st.t.interval(ci,raw_values_mean)
        print (self.min_val, self.max_val)
        self.matrix[self.matrix < self.min_val] = np.ma.masked
        self.matrix[self.matrix > self.max_val] = np.ma.masked

        r = np.ones((self.bin_count,1))
        mdotr = self.matrix.dot(r)
        if not np.all(mdotr != 0):
            print(mdotr)
            raise ValueError("Raw contact matrix does not have total support, and may not converge.(1)")

        c = 1 / mdotr
        mdotc = self.matrix.dot(c)
        if not np.all(mdotc != 0):
            print(mdotc)
            raise ValueError("Raw contact matrix does not have total support, and may not converge.(2)")

        self.total_contacts = np.sum(self.matrix) #for calculated adjusted contact counts later

        print(self.matrix)

        import matplotlib
        matplotlib.use('agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1,figsize = (12,6))
        unique, counts = np.unique(self.matrix, return_counts=True)
        ax.plot(unique,counts, 'k.')
        ax.set_xlim(1,1000)
        ax.set_ylim(0, 100)

        ax.set_xlabel("total counts per bin")
        ax.set_ylabel("frqeuency")

        plt.savefig(file_out)



    def sinkhorn_knopp_correction(self, max_iter = 1000, epsilon = 10.**-6):
        """
        This protocol uses the Sinkhorn-knopp iterative correction algorithm
        to generate a double stochastic matrix from the raw contact count matrix.
        """

        #zeroeth iteration
        r = np.ones((self.bin_count,1))
 
        #keeping an original copy of the matrix may be unnecessary
#        self.cP = np.copy(self.matrix)
        self.cP = self.matrix

        e_min = 1 - epsilon
        e_max = 1 + epsilon
        
        iterations = 0
        while np.any(np.sum(self.cP,axis=1) < e_min) or \
              np.any(np.sum(self.cP, axis=1) > e_max) or \
              np.any(np.sum(self.cP, axis=0) < e_min) or \
              np.any(np.sum(self.cP, axis=0) > e_max):

            c = 1 / self.matrix.T.dot(r)
            r = 1 / self.matrix.dot(c)

            d1 = np.diag(np.squeeze(r))
            d2 = np.diag(np.squeeze(c))
            self.cP = d1.dot(self.matrix).dot(d2)

            iterations += 1

            if iterations > max_iter:
                break

        d1 = np.diag(np.squeeze(r))
        d2 = np.diag(np.squeeze(c))

        self.cP = d1.dot(self.matrix).dot(d2)

    def knight_ruiz_correction(self, remove_zeros = True, min = 1, max = 10**9, epsilon = 10.**-6):
        """
        This protocol uses the Knight-Ruiz iterative correction algorithm
        to generate a double stochastic matrix from the raw contact count matrix.

        For reference:
        http://msp.org/pjm/1967/21-2/pjm-v21-n2-p14-s.pdf
        """
        raise NotImplemented
        
    def calculate_eigenvectors(self):
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        self.pca = pca.fit_transform(self.cP)

    def write_out_sparse_probability_matrix(self, matrix_file_out):
        if self.cP is None:
            raise ValueError("Sparse matrix hasn't been calculated yet.")

#        template = "{row} {column} {raw_counts} {contact_probability} {corrected_counts}\n"
        template = "{row} {column} {raw_counts} {contact_probability} {corrected_counts} {E1} {E2}\n"

        f_out = open(matrix_file_out,'w')
        
        nonzero_coords = zip(*np.nonzero(self.cP))
        print("corrected:\n",  .01 * np.array( np.array(10000*self.cP, dtype = int),dtype=float))
        per_bin_contacts = self.total_contacts / float(self.matrix.shape[0])
        for x,y in list(nonzero_coords):
            if x > y:
                continue #makes the matrix non-redundant as the raw matrix was
            else:
                f_out.write(template.format(row = x,column = y, 
                                            raw_counts = self.matrix[x,y], 
                                            contact_probability = self.cP[x,y], 
                                            corrected_counts = self.cP[x,y] * per_bin_contacts,
                                            E1 = self.pca[x][0],
                                            E2 = self.pca[x][1]
                                    )
                        )
        f_out.close()

def compute_contact_probabilities( matrix_file_in: str, bin_ref:str, matrix_file_out: str, correction_method: str, ci: Optional[float] = 0.999, mask_zeros: Optional[bool] = True,  epsilon: Optional[float] = 10**-3 , max_iter: Optional[float] = 1000, eigen: Optional[bool] = False) -> None:

    #0 determine the size of the map from the reference file
    print("loading reference genome binning map...")
    contact_matrix = HiCMap(bin_ref)

    print("loading raw contact data...")
    #1. read in matrix from raw contact count file
    contact_matrix.from_raw_matrix_file(matrix_file_in)
    
    print("Thresholding data...")
    #2. mask values higher and lower than the confidence interval
    contact_matrix.set_thresholds(matrix_file_out.replace(".corrected.matrix",'.png'),ci = ci)
    
    print("performing correction...")
    #4. write out corrected matrix file in the same sparse format as raw contact matrix
    if correction_method == "SK":
        contact_matrix.sinkhorn_knopp_correction(max_iter = max_iter, epsilon = epsilon )

    if correction_method == "KR":
        raise NotImplementedError

    if eigen:
        print("Identifying principle components...")
        contact_matrix.calculate_eigenvectors()

    print("writing out corrected data...")
    contact_matrix.write_out_sparse_probability_matrix(matrix_file_out)


#takes in a set of contact matrix files and joins them together into a new .matrix file.
# if they've been balanced, the resulting matrix can be balanced if specified.
def join_contact_matrices( bin_ref_file_in, matrix_file_out, *matrix_filenames, correction=False):
    outputted_matrix = HiCMap(bin_ref_file_in)
    for filename in matrix_filenames:
        print("Adding in matrix {}".format(filename))
        #this works because raw matrix file loading is done using the add datum operation
        outputted_matrix.from_raw_matrix_file(filename)

    if correction:
        outputted_matrix.sinkhorn_knopp_correction()

    outputted_matrix.write_out_sparse_probability_matrix(matrix_file_out)
    
