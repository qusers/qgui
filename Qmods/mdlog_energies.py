#!/usr/bin/env python3

# This file is part of Qgui.

# Qgui source file is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Qgui is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Qgui.  If not, see <http://www.gnu.org/licenses/>.

__description__ = \
    """
Functions to prepare extract energies from MD log files from Qdyn
"""

import numpy as np
import os

def is_md_log(file_to_test):
    """
    Check of a file is md log file from Qdyn.
    Returns True or False
    """
    mdlog = False
    logfile = open(file_to_test, 'r')
    for i in range(0, 20):
        if 'reading input from' in logfile.readline().lower():
            mdlog = True
            break

    return mdlog

def get_q_energies(logfiles=[], lambda_='1.00'):
    """
    Collects Q-energies at defined lambda from MD logfiles.
    """

    if len(logfiles) == 0:
        print('No MD logfiles specified! Aborting!')
        return

    qq_el = []
    qq_vdw = []
    qp_el = []
    qp_vdw = []
    qw_el = []
    qw_vdw = []
    qs_el = []
    qs_vdw = []

    qany_el = []
    qany_vdw = []
    qany_bnd = []
    qany_ang = []
    qany_tor = []
    qany_imp = []

    qtypes = ['Q-Q', 'Q-prot', 'Q-wat', 'Q-surr.', 'Q-any']
    qterms = [[qq_el, qq_vdw], [qp_el, qp_vdw], [qw_el, qw_vdw], [qs_el, qs_vdw],
              [qany_el, qany_vdw, qany_bnd, qany_ang, qany_tor, qany_imp]]

    #Go through logfiles and get Q energies:
    for logfile in logfiles:
        if os.path.isfile(logfile):
            with open(logfile) as mdlog:
                for line in mdlog:
                    try:
                        if line.split()[0] in qtypes:
                            if float(line.split()[2]) == float(lambda_):
                                #Decide qtype and get index:
                                for qindex in range(len(qtypes)):
                                    if qtypes[qindex] == line.split()[0]:
                                        insert_index = 0
                                        for eneterm in range(3, len(qterms[qindex]) + 3):
                                            qterms[qindex][insert_index].append(float(line.split()[eneterm]))
                                            insert_index += 1
                    except:
                        continue
        else:
            continue

    #Check if qterms are missing and if, add zero:
    for i in range(len(qterms)):
        for j in range(len(qterms[i])):
            if len(qterms[i][j]) == 0:
                qterms[i][j].append(0)

    qterms_ave = []
    qterms_stderr = []
    for interactions in range(len(qterms)):
        tmp_ave = []
        tmp_stderr = []
        for energy_type in range(len(qterms[interactions])):
            tmp_ave.append(np.average(qterms[interactions][energy_type]))
            tmp_stderr.append(estimate_error(qterms[interactions][energy_type]))
        qterms_ave.append(tmp_ave)
        qterms_stderr.append(tmp_stderr)

    return qterms, qterms_ave, qterms_stderr


def estimate_error(enelist):
    """
    Estimetaes the standard error of a sample by computing the statistical inefficiency
    """

    if np.average(enelist) == 0:
        return 0

    var_tot = (np.std(enelist)) ** 2
    ave_tot = np.average(enelist)
    s_list = list()

    t_max = len(enelist)

    #Estimate s as average of dividing sample in blocks of 4 - 10
    for b in range(5, 11):


        blocks = list()

        #Find block size
        t = int(t_max / float(b))

        for i in range(b):
            if i < b - 1:

                blocks.append(enelist[i * t: (i + 1) * t])

            elif i == b - 1:
                blocks.append(enelist[i * t:])

        block_ave = list()

        for j in blocks:
            block_ave.append(np.average(j))

        #bloc_var = (np.std(block_ave)) ** 2
        bloc_var = 0
        for ave in block_ave:
            bloc_var += (ave - ave_tot) ** 2

        bloc_var *= (1.0 / len(block_ave))

        s = float(t) * (bloc_var / var_tot)

        s_list.append(s)

    s_b = np.average(s_list)

    se = np.std(enelist) / np.sqrt(len(enelist) / s_b)

    return se
