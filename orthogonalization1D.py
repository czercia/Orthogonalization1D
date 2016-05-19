import math
import multiprocessing
import numpy as np
import warnings

import integrals

warnings.filterwarnings("ignore")
from contextlib import contextmanager


@contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


import csv

# ***********************************************************************************************
# parameters


d_init = 1
delta_d = 2
l_max = 11

alpha = 0.001
files_before = False
b = 0.1
plot_functions = False
ile_stanow = []
ile_stanow.append(10)

for i in range(1, l_max):
    ile_stanow.append(10)

# ************************************************************************************************
# read which functions to get


if not files_before:
    n_before = 0
if files_before:
    n_before = 0
    with open("spp_.txt") as data3:
        for line in csv.reader(data3, dialect="excel-tab"):
            n_before += 1

print ("n_before = " + str(n_before))

energy = []
n_states = 0
count = 1
l = []
n_l = []
n_l = ile_stanow

# *************************
for ll in range(l_max):
    with open("/home/martas/numerov_11_5/E_levels_l=" + str(ll) + "_b=" + str(b) + "_1.dat") as en:
        n = 0
        for line in csv.reader(en):
            if n < n_l[ll]:
                energy.append(float(line[0]))
		print (str(ll) + "   " + line[0])
                n += 1
n_states = np.sum(n_l)
# *************************
# with open("energy_levels.txt") as data3:
#     for line in csv.reader(data3):
#         energy.append(float(line[0]))
#         n_states += 1
#
# for i in range(n_states - 1):
#     if energy[i + 1] > energy[i]:
#         count += 1
#     else:
#         n_l.append(count)
#         count = 1
#     i += 1
# n_l.append(count)

j = 0
for value in n_l:
    print ('n_l[' + str(j) + '] = ' + str(n_l[j]))
    j += 1

for ll in range(l_max):
    for k in range(n_l[ll]):
        l.append(ll)

k = 0
for value in l:
    print ('l[' + str(k) + '] = ' + str(l[k]))
    k += 1

r_max_value = []
for i in range(n_states):
    if energy[i] <= 0:
        r_max_value.append(2)
    if 0 < energy[i] < 0.1:
        r_max_value.append(20)
    if 0.1 <= energy[i] < 0.3:
        r_max_value.append(40)
    if 0.3 <= energy[i] < 0.5:
        r_max_value.append(2.2 * np.sqrt(energy[i] / alpha))
    if 0.5 <= energy[i] < 4:
        r_max_value.append((2 - 0.1 * energy[i]) * np.sqrt(energy[i] / alpha))
    if 4 <= energy[i] < 10:
        r_max_value.append((2 - 0.05 * energy[i]) * np.sqrt(energy[i] / alpha))
    if energy[i] >= 10:
        r_max_value.append((2 - 0.01 * energy[i]) * np.sqrt(energy[i] / alpha) - 40)

print ("n_states = " + str(n_states))

# *************************************************************************************************
# get functions from file, normalize, interpolate
list_of_f_values = []
list_of_r_values = []

for j in range(n_states):
    r_values = []
    f_values = []
    with open("/home/martas/numerov_11_5_normalized_dr_const1D/psi_l=" + str(l[j]) + "_" + str(energy[j]) + '_normalized_dr_const1D.dat') as data:
        counter = 0
        for lin in csv.reader(data, dialect="excel-tab"):
            if float(lin[0]) != 0 and float(lin[0]) < r_max_value[j]:
                r_values.append(float(lin[0]))
                f_values.append(float(lin[1]) / math.fabs(float(lin[0])))
            counter += 1
        r_values.append(r_max_value[j])
        f_values.append(0)
list_of_r_values.append(r_values)
list_of_f_values.append(f_values)


def f(i, r):
    return np.interp(r, list_of_r_values[i], list_of_f_values[i])


# ************************************************************************************************************
# operations on files and matrices


def read_file(matrix, name, nb):
    if files_before:
        w = -1
        with open(str(name) + ".txt") as filename:
            for line in csv.reader(filename, dialect="excel-tab"):
                w += 1
                for j in range(nb):
                    matrix[w, j] = line[j]


def save_file(matrix, name, dim):
    plik = open(str(name) + '.txt', 'w')
    for i in range(dim * n_states):
        for j in range(dim * n_states):
            plik.write(str(matrix[i, j]))
            plik.write(' \t')
        plik.write('\n')
    plik.close()


def sym(matrix):
    for i in range(n_states):
        for j in range(n_states):
            if i > j:
                matrix[i, j] = matrix[j, i]


spp = np.zeros((n_states, n_states))
err_spp = np.zeros((n_states, n_states))
read_file(spp, 'spp_', n_before)
for i in range(n_states):
    for j in range(n_before, n_states):
        if i <= j:
            print("spp(i, j) = (" + str(i) + ", " + str(j) + ")")
            spp_values = integrals.spp_1d_integrate(i, j, r_max_value)
            spp[i, j] = spp_values[0]
            err_spp[i, j] = spp_values[1]
            print spp[i, j]
sym(spp)
save_file(spp, 'spp_', 1)
print ("spp:")
print spp


# def dd(p):
#     return d_init + p * delta_d
#
#
# def part(p_init_2, p_end_2):
#     for p in range(p_init_2, p_end_2):
#         d = dd(p)
#         print ("d = " + str(d))
#
#         spm = np.matrix(np.zeros(n_states))
#         err_spm = np.matrix(np.zeros(n_states))
#
#         n_before = 0
#         if files_before == True:
#             with open("spm_" + str(d) + ".txt") as data3:
#                 for line in csv.reader(data3, dialect="excel-tab"):
#                     n_before = n_before + 1
#
#         read_file(spm, 'spm_' + str(d), n_before)
#
#         for i in range(n_states):
#             for j in range(n_before, n_states):
#                 if i <= j:
#                     print("spm(" + str(i) + ", " + str(j) + "), d =" + str(d))
#                     spm_values = integrals.spm_1d_integrate(i, j, r_max_value, d)
#                     spm[i, j] = spm_values[0]
#                     err_spm[i, j] = spm_values[1]
#                     print spm[i, j]
#         sym(spm)
#         sym(err_spm)
#         save_file(spm, 'spm_' + str(d), 1)
#         save_file(err_spm, 'err_spm_' + str(d), 1)
#         print ("spm,  d = " + str(d))
#         print spm
#
#         s = np.matrix(np.identity(2 * n_states))
#
#         for i in range(2 * n_states):
#             for j in range(2 * n_states):
#                 if i < n_states:
#                     if j < n_states:
#                         s[i, j] = spp[i, j]
#                     if j >= n_states:
#                         k = j - n_states
#                         s[i, j] = spm[i, k]
#                 if i >= n_states:
#                     if j < n_states:
#                         s[i, j] = spm[i - n_states, j]  # mieszane
#                     if j >= n_states:
#                         s[i, j] = spp[i - n_states, j - n_states]  # minusy
#         save_file(s, 's_' + str(d), 2)
#
#         # # *********************************************************************************
#         ee = np.matrix(np.identity(2 * n_states))
#         i = -1
#         with open("spp_.txt") as data3:
#             for line in csv.reader(data3, dialect='excel-tab'):
#                 i += 1
#                 for j in range(n_states):
#                     ee[i, j] = (energy[j] - alpha * d ** 2) * float(line[j])
#         i = n_states - 1
#         with open("spp_.txt") as data3:
#             for line in csv.reader(data3, dialect='excel-tab'):
#                 i += 1
#                 for j in range(n_states, 2 * n_states):
#                     ee[i, j] = (energy[j - n_states] - alpha * d ** 2) * float(line[j - n_states])
#         i = -1
#         with open("spm_" + str(d) + ".txt") as data3:
#             for line in csv.reader(data3, dialect='excel-tab'):
#                 i += 1
#                 for j in range(n_states, 2 * n_states):
#                     ee[i, j] = (energy[j - n_states] - alpha * d ** 2) * float(line[j - n_states])
#         i = n_states - 1
#         with open("spm_" + str(d) + ".txt") as data3:
#             for line in csv.reader(data3, dialect='excel-tab'):
#                 i += 1
#                 for j in range(n_states):
#                     ee[i, j] = (energy[j] - alpha * d ** 2) * float(line[j])
#
#         rpp = np.matrix(np.identity(n_states))
#         err_rpp = np.matrix(np.identity(n_states))
#
#         n_before = 0
#         if files_before == True:
#             with open("rpp_" + str(d) + ".txt") as data3:
#                 for line in csv.reader(data3, dialect="excel-tab"):
#                     n_before = n_before + 1
#
#         for i in range(n_states):
#             for j in range(n_before, n_states):
#                 if i <= j:
#                     print("rpp(i, j) = (" + str(i) + ", " + str(j) + ")")
#                     rpp_values = integrals.rpp_1d_integrate(i, j, r_max_value, d)
#                     rpp[i, j] = rpp_values[0]
#                     err_rpp[i, j] = rpp_values[1]
#                     print rpp[i, j]
#                     print err_rpp[i, j]
#         sym(rpp)
#         sym(err_rpp)
#         save_file(rpp, 'rpp_' + str(d), 1)
#         save_file(err_rpp, 'err_rpp_' + str(d), 1)
#         print rpp
#
#         rpm = np.matrix(np.identity(n_states))
#         err_rpm = np.matrix(np.identity(n_states))
#
#         n_before = 0
#         if files_before == True:
#             with open("rpm_" + str(d) + ".txt") as data3:
#                 for line in csv.reader(data3, dialect="excel-tab"):
#                     n_before = n_before + 1
#
#         read_file(rpm, 'rpm_' + str(d), n_before)
#
#         for i in range(n_states):
#             for j in range(n_before, n_states):
#                 if i <= j:
#                     print("rpm(i, j) = (" + str(i) + ", " + str(j) + "), d =" + str(d))
#                     rpm_values = integrals.rpp_1d_integrate(i, j, r_max_value, d)
#                     rpm[i, j] = rpm_values[0]
#                     err_rpm[i, j] = rpm_values[1]
#                     print rpm[i, j]
#                     print err_rpm[i, j]
#         sym(rpm)
#         sym(err_rpm)
#         save_file(rpm, 'rpm_' + str(d), 1)
#         save_file(err_rpm, 'err_rpm' + str(d), 1)
#         print ("rpm (d = " + str(d) + ': ' + str(rpm))
#         app = np.matrix(np.identity(n_states))
#
#         err_app = np.matrix(np.identity(n_states))
#         n_before = 0
#         if files_before == True:
#             with open("app_" + str(d) + ".txt") as data3:
#                 for line in csv.reader(data3, dialect="excel-tab"):
#                     n_before = n_before + 1
#             read_file(app, 'app_' + str(d), n_before)
#
#         for i in range(n_states):
#             for j in range(n_before, n_states):
#                 if i <= j:
#                     print("app(i, j) = (" + str(i) + ", " + str(j) + ")")
#                     app_values = integrals.app_1d_integrate(i, j, r_max_value, d)
#                     app[i, j] = app_values[0]
#                     err_app[i, j] = app_values[1]
#                     print app[i, j]
#                     print err_app[i, j]
#         sym(app)
#         sym(err_app)
#         save_file(app, 'app_' + str(d), 1)
#         save_file(err_app, 'err_app' + str(d), 1)
#         print ("app (d = " + str(d) + ': ' + str(app))
#
#         apm = np.matrix(np.identity(n_states))
#         err_apm = np.matrix(np.identity(n_states))
#         n_before = 0
#         if files_before:
#             with open("apm_" + str(d) + ".txt") as data3:
#                 for line in csv.reader(data3, dialect="excel-tab"):
#                     n_before += 1
#             read_file(apm, 'apm_' + str(d), n_before)
#         for i in range(n_states):
#             for j in range(n_before, n_states):
#                 if i <= j:
#                     print("apm(i, j) = (" + str(i) + ", " + str(j) + ")")
#                     apm_values = integrals.apm_1d_integrate(i, j, r_max_value, d)
#                     apm[i, j] = apm_values[0]
#                     err_apm[i, j] = apm_values[1]
#                     print apm[i, j]
#                     print err_apm[i, j]
#         sym(apm)
#         sym(err_apm)
#         save_file(apm, 'apm_' + str(d), 1)
#         save_file(err_apm, 'err_apm' + str(d), 1)
#         print ("apm (d = " + str(d) + ': ' + str(apm))
#
#         h = np.matrix(np.identity(2 * n_states))
#
#         for i in range(2 * n_states):
#             for j in range(2 * n_states):
#                 if i < n_states:
#                     if j < n_states:
#                         h[i, j] = ee[i, j] - rpp[i, j] + 0.5 * alpha * app[i, j]
#                     if j >= n_states:
#                         k = j - n_states
#                         h[i, j] = - rpm[i, k] + ee[i, j] + 0.5 * alpha * apm[i, k]
#                 if i >= n_states:
#                     if j < n_states:
#                         h[i, j] = - rpm[i - n_states, j] + ee[i, j] + 0.5 * alpha * apm[i - n_states, j]  # mieszane
#                     if j >= n_states:
#                         h[i, j] = ee[i, j] - rpp[i - n_states, j - n_states] + 0.5 * alpha * app[
#                             i - n_states, j - n_states]  # minusy
#
#         save_file(h, 'H_' + str(d), 2)
#
#
# # print h
#
#
#
#
#
#
#
#
# # # ************************************************************************************
#
# num_workers = multiprocessing.cpu_count()
# print num_workers
# # n_integrals = n_states / (num_workers - 1)  # na kazdy proces ile calek, ostatni proces bierze to co zostalo
# n_integrals_last = n_states % (num_workers - 1)

# processes = []
#
# for x in range(num_workers - 1):
#     processes.append(multiprocessing.Process(target=part, args=(x, x + 1)))
# # processes.append(multiprocessing.Process(target=part, args=(
# #     (num_workers - 1) * n_integrals, num_workers * n_integrals + n_integrals_last)))
#
# if __name__ == '__main__':
#     for x in range(num_workers):
#         processes[x].start()