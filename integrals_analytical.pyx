from scipy import integrate
from scipy import special
import numpy as np
cimport numpy as np

cdef double f(double x, double ni):
    if x <=0:
        return special.pbdv(ni, -x)[0]
    else:
        return special.pbdv(ni, x)[0]

cdef double norm(double ni, double x_max):
    return 1./np.sqrt(integrate.quad(lambda x: abs(f(x, ni) * f(x, ni)) , -x_max, x_max)[0])

def spp_1d_integrate(np.ndarray matrix,  np.ndarray ni, double x_max):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * integrate.quad(lambda x: f(x, ni[i]) * f(x, ni[j]), -x_max, x_max)[0]
    return result

def spm_1d_integrate(np.ndarray matrix, np.ndarray ni, double x_max, double d):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * integrate.quad(lambda x: f(x - d, ni[i]) * f(x +d, ni[j]), -x_max, x_max)[0]
    return result