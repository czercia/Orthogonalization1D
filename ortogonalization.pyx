#orthogonalization module

import numpy as np
cimport numpy as np

def calculate_matrix(np.ndarray mpp, np.ndarray mpm):
    cdef np.ndarray[np.float64_t, ndim = 2] result
    m1 = np.concatenate((mpp, mpm), axis = 1)
    m2 = np.concatenate((mpm, mpp), axis = 1)
    result = np.concatenate((m1, m2))
    return result

def calculate_H(np.ndarray S, np.ndarray A, np.ndarray R, double d):
    cdef  np.ndarray[np.float64_t, ndim = 2] H
    #TODO macierz Mij =  Ej
    H =

