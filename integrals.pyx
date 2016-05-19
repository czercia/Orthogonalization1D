#  integrals module
from scipy import integrate
from scipy import special


def spp_1d_integrate(i, j, rm):
    lim_low = - min(rm[i], rm[j])
    lim_up = min(rm[i], rm[j])
    return integrate.quad(lambda x: f(i, x) * f(j, x), lim_low, lim_up)


def spm_1d_integrate(i, j, rm, d):
    lim_low = - rm[j] + d
    lim_up = rm[i] - d
    return integrate.quad(lambda x: f(i, x - d) * f(j, x + d), lim_low, lim_up)


def rpp_1d_integrate(i, j, rm, d):
    lim_low = - min(rm[i], rm[j]) - d
    lim_up = min(rm[i], rm[j]) - d
    return integrate.quad(lambda x: f(i, x + d) * (1 / ((x - d) ** 2 + b ** 2) ** 2) * f(j, x + d), lim_low, lim_up)


def rpm_1d_integrate(i, j, rm, d):
    lim_low = - rm[j] + d
    lim_up = rm[i] - d
    return integrate.quad(lambda x: f(i, x - d) * (1 / ((x - d) ** 2 + b ** 2) ** 2) * f(j, x + d), lim_low, lim_up)


def app_1d_integrate(i, j, rm, d):
    lim_low = - min(rm[i], rm[j]) - d
    lim_up = min(rm[i], rm[j]) - d
    return integrate.quad(lambda x: f(i, x + d) * (x - d) ** 2 * f(j, x + d), lim_low, lim_up)


def apm_1d_integrate(i, j, rm, d):
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