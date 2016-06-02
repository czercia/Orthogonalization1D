import numpy as np
import integrals
import integrals_analytical
import multiprocessing
from time import time

import csv


#parameters
kappa = 1
d_min = 0
d_max = 6
delta_d = 0.1
x_max = 10
d = 1

#get the list of ni from file
ni = np.loadtxt("ni_list.dat", dtype=float)

matrix = np.zeros((5, 5))
start = time()
spp = integrals_analytical.spp_1d_integrate(matrix, ni, x_max)
spm = integrals_analytical.spm_1d_integrate(matrix, ni, x_max, d)
app = integrals_analytical.app_1d_integrate(matrix, ni, x_max, d)
apm = integrals_analytical.apm_1d_integrate(matrix, ni, x_max, d)
rpp = kappa * integrals_analytical.rpp_1d_integrate(matrix, ni, x_max, d)
rpm = kappa * integrals_analytical.rpm_1d_integrate(matrix, ni, x_max, d)

time =  time() - start
print ('calculations finished in ' + str(time) )

a = np.array([[1, 1], [1, 1]])
b = np.array([[2, 2], [2, 2]])
c = np.concatenate((a, b), axis = 1)
d = np.concatenate((b, a), axis=1)
print c
print d
e = np.concatenate((c, d))
print e

