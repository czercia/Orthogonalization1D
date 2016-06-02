from scipy import integrate
from scipy import special
import numpy as np
cimport numpy as np

cdef double f(double x, double ni):
    if x <= 0:
        return special.pbdv(ni, -x)[0]
    else:
        return special.pbdv(ni, x)[0]

cdef double norm(double ni, double x_max):
    return 1. / np.sqrt(integrate.quad(lambda x: abs(f(x, ni) * f(x, ni)), -x_max, x_max)[0])

def spp_1d_integrate(np.ndarray matrix, np.ndarray ni, double x_max):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i <= j:
                result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * \
                           integrate.quad(lambda x: f(x, ni[i]) * f(x, ni[j]), -x_max, x_max)[0]
            else:
                result[i, j] = result[j, i]
    return result

def spm_1d_integrate(np.ndarray matrix, np.ndarray ni, double x_max, double d):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i <= j:
                result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * \
                           integrate.quad(lambda x: f(x - d, ni[i]) * f(x + d, ni[j]), -x_max, x_max)[0]
            else:
                result[i, j] = result[j, i]
    return result

def app_1d_integrate(np.ndarray matrix, np.ndarray ni, double x_max, double d):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i <= j:
                result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * \
                           integrate.quad(lambda x: f(x + d, ni[i]) * (x - d) * (x - d) * f(x + d, ni[j]), -x_max,
                                          x_max)[0]
            else:
                result[i, j] = result[j, i]
    return result

def apm_1d_integrate(np.ndarray matrix, np.ndarray ni, double x_max, double d):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i <= j:
                result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * \
                               integrate.quad(lambda x: f(x - d, ni[i]) * (x - d) * (x - d) * f(x + d, ni[j]), -x_max,
                                              x_max)[0]
            else:
                result[i, j] = result[j, i]
    return result

def rpp_1d_integrate(np.ndarray matrix, np.ndarray ni, double x_max, double d):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i <= j:
                result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * f(2 * d, ni[i]) * f(2 * d, ni[j])
            else:
                result[i, j] = result[j, i]
    return result

def rpm_1d_integrate(np.ndarray matrix, np.ndarray ni, double x_max, double d):
    cdef Py_ssize_t i, j
    cdef np.ndarray[np.float64_t, ndim = 2] result
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i <= j:
                result[i, j] = norm(ni[i], x_max) * norm(ni[j], x_max) * f(0, ni[i]) * f(2 * d, ni[j])
            else:
                result[i, j] = result[j, i]
    return result
