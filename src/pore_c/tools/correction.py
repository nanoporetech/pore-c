from typing import Pattern, List, NamedTuple, Iterator, Union, Dict, Tuple, Optional

from dataclasses import dataclass

import scipy.stats as st
import numpy as np

@dataclass
class HiCMap(object):
    def __init__(self, bin_ref):
        self.bin_count = 0
        for entry in open(bin_ref):
            self.bin_count += 1

        self.matrix = np.ma.zeros((self.bin_count,self.bin_count),dtype = float)
        self.cP = None #np.copy(self.matrix)
        self.min_val = np.min(self.matrix)
        self.max_val = np.max(self.matrix)
        self.total_contacts = 0

    def populate_matrix(self, matrix_file):
        symmetry_test = set()

        for entry in open(matrix_file):
            l = entry.strip().split()
            try:
                assert len(l) == 3
            except:
                continue
            l[0] = int(l[0])
            l[1] = int(l[1])
            if tuple(sorted([l[0],l[1]])) in symmetry_test:
                raise ValueError("This data set is not contact commutative. There are contacts such that there is an AB contact value that is not equal to BA contact value.")
            symmetry_test.add(tuple(sorted([l[0],l[1]])))
            l[2] = float(l[2])
            if l[0] == l[1]:
                self.matrix[l[0],l[1]] = l[2]
            else:
                self.matrix[l[0],l[1]] = self.matrix[l[1],l[0]] = l[2] 

        #must mask any rows/columns that are zeros
        np.ma.masked_where(self.matrix == 0, self.matrix)
        self.mask = self.matrix.mask #store the mask to prevent overwriting information about blank row/columns

        #confirm that the resulting matrix is a square, strictly positive matrix with total support
        try:
            assert self.matrix.ndim == 2
            assert self.matrix.shape[0] == self.matrix.shape[1]
            assert np.all(self.matrix >= 0)
        except:
            raise ValueError("Contact matrix is not square and strictly positive.")

        size = self.matrix.shape[0]

        r = np.ones((size,1))
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

    def set_thresholds(self, ci=0.999, remove_zeros = True):
        """
        This protocol establishes a mask for the contact matrix
        such that values above the confidence interval are masked
        since they can skew the correction protocols.
        Additionally, low and zero values are masked as their prominence 
        can be unduly amplified. 
        """

        # mask zeros
        if remove_zeros:
            np.ma.masked_where(self.matrix == 0, self.matrix)
        
        raw_values_mean = np.mean(self.matrix)

        #interval(alpha,df,loc=0,scale=1)
        self.min_val, self.max_val = st.t.interval(ci,raw_values_mean)
        self.matrix[self.matrix < self.min_val] = np.ma.masked
        self.matrix[self.matrix > self.max_val] = np.ma.masked

    def sinkhorn_knopp_correction(self, max_iter = 1000, epsilon = 10.**-6):
        """
        This protocol uses the Sinkhorn-knopp iterative correction algorithm
        to generate a double stochastic matrix from the raw contact count matrix.
        """

        #zeroeth iteration
        size = self.matrix.shape[0]

        r = np.ones((size,1))
        
        self.cP = np.copy(self.matrix)

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

        print(iterations)

    def knight_ruiz_correction(self, remove_zeros = True, min = 1, max = 10**9, epsilon = 10.**-6):
        """
        This protocol uses the Knight-Ruiz iterative correction algorithm
        to generate a double stochastic matrix from the raw contact count matrix.

        For reference:
        http://msp.org/pjm/1967/21-2/pjm-v21-n2-p14-s.pdf
        """
        raise NotImplemented
        

    def write_out_sparse_probability_matrix(self, matrix_file_out):
        if self.cP is None:
            raise ValueError("Sparse matrix hasn't been calculated yet.")

        template = "{row} {column} {raw_counts} {contact_probability} {corrected_counts}\n"

        f_out = open(matrix_file_out,'w')
        
        nonzero_coords = zip(*np.nonzero(self.cP))
        print("corrected:\n",  .01 * np.array( np.array(10000*self.cP, dtype = int),dtype=float))
        for x,y in list(nonzero_coords):
            f_out.write(template.format(row = x,column = y, 
                                        raw_counts = self.matrix[x,y], 
                                        contact_probability = self.cP[x,y], 
                                        corrected_counts = int(self.cP[x,y] * self.total_contacts / len(self.cP[x]))
                                    )
                        )
        f_out.close()

def compute_contact_probabilities( matrix_file_in: str, bin_ref:str, corrected_matrix_file_out: str, output_eigenvector_data: str, correction_method: str, ci: Optional[float] = 0.999, mask_zeros: Optional[bool] = True,  epsilon: Optional[float] = 10**-3 , max_iter: Optional[float] = 1000) -> None:

    #0 determine the size of the map from the reference file
    contact_matrix = HiCMap(bin_ref)

    #1. read in matrix from raw contact count file
    contact_matrix.populate_matrix(matrix_file_in)
    
    #2. mask values higher and lower than the confidence interval
    contact_matrix.set_thresholds(ci = ci)
    
    #4. write out corrected matrix file in the same sparse format as raw contact matrix
    if correction_method == "SK":
        contact_matrix.sinkhorn_knopp_correction(max_iter = max_iter, epsilon = epsilon )
    if correction_method == "KR":
        raise NotImplementedError

    contact_matrix.write_out_sparse_probability_matrix(corrected_matrix_file_out)
