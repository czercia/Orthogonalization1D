#  integrals module
from scipy import integrate
from scipy import special
import numpy as np
cimport numpy as np
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)

ctypedef np.ndarray[np.float64_t, ndim = 2] arr_2D
ctypedef np.ndarray[np.float64_t, ndim = 1] arr_1D

cdef double f( double x, arr_1D x_vals, arr_1D y_vals):
    return np.interp(x, x_vals, y_vals)


def spp_1d_integrate(arr_2D matrix,  arr_1D rm, arr_1D x_vals, arr_1D y_vals):
    cdef double lim_low, lim_up
    cdef Py_ssize_t i, j
    result = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            lim_low = - min(rm[i], rm[j])
            lim_up = min(rm[i], rm[j])
            result[i, j] = integrate.quad(lambda x: f(i, x) * f(j, x), lim_low, lim_up)[0]
    return result


def spm_1d_integrate(i, j, rm, d):
    cdef double lim_low, lim_up
    lim_low = - rm[j] + d
    lim_up = rm[i] - d
    return integrate.quad(lambda x: f(i, x - d) * f(j, x + d), lim_low, lim_up)


def rpp_1d_integrate(i, j, rm, d):
    cdef double lim_low, lim_up
    lim_low = - min(rm[i], rm[j]) - d
    lim_up = min(rm[i], rm[j]) - d
    return integrate.quad(lambda x: f(i, x + d) * (1 / ((x - d) ** 2 + b ** 2) ** 2) * f(j, x + d), lim_low, lim_up)


def rpm_1d_integrate(i, j, rm, d):
    cdef double lim_low, lim_up
    lim_low = - rm[j] + d
    lim_up = rm[i] - d
    return integrate.quad(lambda x: f(i, x - d) * (1 / ((x - d) ** 2 + b ** 2) ** 2) * f(j, x + d), lim_low, lim_up)


def app_1d_integrate(i, j, rm, d):
    cdef double lim_low, lim_up
    lim_low = - min(rm[i], rm[j]) - d
    lim_up = min(rm[i], rm[j]) - d
    return integrate.quad(lambda x: f(i, x + d) * (x - d) ** 2 * f(j, x + d), lim_low, lim_up)


def apm_1d_integrate(i, j, rm, d):
    cdef double lim_low, lim_up
    lim_low = - rm[j] + d
    lim_up = rm[i] - d
    return integrate.quad(lambda x: f(i, x - d) * (x - d) ** 2 * f(j, x + d), lim_low, lim_up)


def integration(i, j, z_min, z_max, d_i, d_j, rm, d):
    return integrate.dblquad(lambda rho, zet: 2 * math.pi * rho * f(i, np.sqrt(rho ** 2 + (zet + d_i) ** 2))
                                              * f(j, np.sqrt(rho ** 2 + (zet + d_j) ** 2))
                                              * special.sph_harm(0, l[i], 0,
                                                                 ((zet + d_i) / np.sqrt(rho ** 2 + (zet + d_i) ** 2)))
                                              * special.sph_harm(0, l[j], 0,
                                                                 ((zet + d_j) / np.sqrt(rho ** 2 + (zet + d_j) ** 2))),
                             z_min, z_max,
                             lambda zet: 0,
                             lambda zet: np.sqrt((rm) ** 2 - (zet + d) ** 2),
                             epsabs=1.49e-03,
                             epsrel=1.49e-03)


def integration_r(i, j, z_min, z_max, d_i, d_j, rm, d):
    return integrate.dblquad(lambda rho, zet: 2 * math.pi * rho * f(i, np.sqrt(rho ** 2 + (zet + d_i) ** 2))
                                              * special.sph_harm(0, l[i], 0,
                                                                 ((zet + d_i) / np.sqrt(rho ** 2 + (zet + d_i) ** 2)))
                                              * special.sph_harm(0, l[j], 0,
                                                                 ((zet + d_j) / np.sqrt(rho ** 2 + (zet + d_j) ** 2)))
                                              * f(j, np.sqrt(rho ** 2 + (zet + d_j) ** 2)) / (
                                                  (np.sqrt((zet - d_j) ** 2 + rho ** 2) + b ** 2) ** 2),
                             z_min, z_max,
                             lambda zet: 0,
                             lambda zet: np.sqrt((rm) ** 2 - (zet + d) ** 2),
                             epsabs=1.49e-03,
                             epsrel=1.49e-03)


def integration_alpha(i, j, z_min, z_max, d_i, d_j, rm, d):
    return integrate.dblquad(
        lambda rho, zet: 0.5 * alpha * 2 * math.pi * rho * f(i, np.sqrt(rho ** 2 + (zet + d_i) ** 2))
                         * special.sph_harm(0, l[i], 0, ((zet + d_i) / np.sqrt(rho ** 2 + (zet + d_i) ** 2)))
                         * special.sph_harm(0, l[j], 0, ((zet + d_j) / np.sqrt(rho ** 2 + (zet + d_j) ** 2)))
                         * f(j, np.sqrt(rho ** 2 + (zet + d_j) ** 2)) * (rho ** 2 + (zet - d_j) ** 2),
        z_min, z_max,
        lambda zet: 0,
        lambda zet: np.sqrt((rm) ** 2 - (zet + d) ** 2),
        epsabs=1.49e-03,
        epsrel=1.49e-03)