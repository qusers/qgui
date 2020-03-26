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
Collection of functions commonly used by the different Qgui Classes.
"""

import numpy as np
import os
import re
from subprocess import call

def read_topology(topology, libfiles=list()):
    """
    :param topology: a topolgy file (*.top)
    :param libfiles: a list with absolute paths to library files in settings
    :return: atom_xyz ( dict {atom nr: [x,y,z]} ), res_nr ( dict {res nr: res name} ), libfiles (list with paths)
    """
    found_xyz = False
    found_res = False

    atomnr_xyz = dict()
    resnr_res = dict()

    #atom number
    atomnr = 0

    #residue number
    resnr = 0

    with open(topology, 'r') as top:
        for line in top:

            #Get library files:
            if 'LIB_FILES' in line:
                for lib in line.split()[1].split(';'):
                    if lib not in libfiles:
                        if os.path.isfile(lib):
                            libfiles.append(lib)
                        else:
                            print(('Warning: Did not find library file: %s' % lib))

            #Coordinates ends here:
            if 'No. of integer atom codes' in line:
                found_xyz = False

            #Residues ends here:
            if 'No. of separate molecules.' in line:
                found_res = False
                #for now I do not need any more info from the topology file, so stop reading it.
                break

            #Collect coordinates
            if found_xyz:
                if len(line.split()) > 2:
                    atomnr += 1
                    atomnr_xyz[atomnr] = list(map(float, line.split()[0:3]))

                if len(line.split()) == 6:
                    atomnr += 1
                    atomnr_xyz[atomnr] = list(map(float, line.split()[3:6]))

            #Collect residues
            if found_res:
                for residue in line.split():
                    resnr += 1
                    resnr_res[resnr] = residue

            #Coordinates starts here:
            if 'Coordinates: (2*3 per line)' in line:
                found_xyz = True

            #Residues starts here:
            if 'Sequence:' in line:
                found_res = True

        return atomnr_xyz, resnr_res, libfiles

def create_pdb_from_topology(topology, libfiles=list()):
    """
    :param topology:
    :param libfiles:
    :return: pdbfile (list)
    """
    pdbfile = list()

    #Get atoms coordinates, atom nr and residues from topology
    atommr_xyz, resnr_res, libs = read_topology(topology)

    for lib in libs:
        if lib not in libfiles:
            libfiles.append(lib)

    #Get a dictionary for all defined residues in library:
    lib = dict()
    for libfile in libfiles:
        lib.update(read_library(libfile))

    #go through residue list and create pdb file based on library and xyz
    atomnr = 0
    for resnr in sorted(resnr_res.keys()):
        residue = resnr_res[resnr]
        for nr in sorted(lib[residue].keys()):
            atomnr += 1
            atom = lib[residue][nr]['atom name']
            x,y,z = atommr_xyz[atomnr][:]

            pdbfile.append('ATOM  %5d  %4s%4s %4d    %8.3f%8.3f%8.3f\n' %
                           (atomnr, atom.ljust(4), residue.ljust(4), resnr, x, y, z))

    return pdbfile

def write_top_pdb(topology, pdbname, pdbpath, libfiles=list()):
    """
    Creates PDB file from topology file and writes it to pdbpath/pdbname.pdb
    :param topology: path to topology file
    :param pdbname: name of PDB file to write
    :param pdbpath: path to write PDB file
    :param libfiles: optional, libfiles (can be that libfiles defined in topology has been moved/deleted)
    :return: Nothing
    """

    pdb = create_pdb_from_topology(topology, libfiles)

    pdbout = open('%s/%s' % (pdbpath, pdbname.strip('/')), 'w')

    for line in pdb:
        pdbout.write(line)

    pdbout.close()
    print(('%s written to %s' % (pdbname, pdbpath)))


def get_pdb_resnr(pdbfile, res_nr):
    """
    :param pdbfile:
    :param resnr:
    :return: list of pdb lines for given residue number
    """
    found_res = False
    pdb_res = list()

    for line in pdbfile:
        if 'ATOM' in line:
            if found_res and int(line[22:26]) != int(res_nr):
                break
            if int(line[22:26]) == int(res_nr):
                found_res = True
                pdb_res.append(line)

    return pdb_res


def read_library(libfile):
    """
    :param libfile:
    :return: {res: {nr: {atom name, atom type, atom charge}}}
    """
    libdict = dict()
    res = None

    found_atoms = False
    with open(libfile, 'r') as lib:
        for line in lib:
            if re.search('{.*}', line):
                found_atoms = False
                res = re.search('{.*}', line).group(0).strip('{*}')
                libdict[res] = dict()
                atomnr = 0
            if found_atoms:
                if re.search('\[.*\]', line):
                    found_atoms = False

            if found_atoms:
                if len(line.split()) > 3 and not line.startswith('!'):
                    nr, name, atomtype, charge = line.split()[0:4]
                    nr = int(nr)
                    libdict[res][nr] = dict()

                    libdict[res][nr]['atom name'] = name
                    libdict[res][nr]['atom type'] = atomtype
                    libdict[res][nr]['atom charge'] = float(charge)

            if '[atoms]' in line:
                found_atoms = True

    return libdict


def get_fep_atoms(fepfile):
    """
    Simple function to get Q nr - atom nr - atom name from FEP file
    :param fepfile:
    :return: qnr_atomnr and qnr_atomname dicts
    """
    qnr_atomnr = dict()
    qnr_atomname = dict()

    found_atoms = False

    with open(fepfile, 'r') as fep:
        for line in fep:
            if found_atoms:
                if re.search('\[.*\]', line):
                    return qnr_atomnr, qnr_atomname
            if found_atoms and len(line.split()) > 1:
                q_nr = int(line.split()[0])
                atom_nr = int(line.split()[1])
                atom_name = line.split('!')[1].split()[0]

                qnr_atomnr[q_nr] = atom_nr
                qnr_atomname[q_nr] = atom_name

            if '[atoms]' in line:
                found_atoms = True

    return qnr_atomnr, qnr_atomname


def read_fep(fepfile, qoffset=0, atomoffset=None):
    """
    :param fepfile:
    :return: dictionary with all FEP info
    """
    #if atomsoffset is None, fill it with 0! {qnr: atomnumber offset}
    if not atomoffset:
        atomoffset=dict()
        for i in range(0,9999):
            atomoffset[i] = 0

    #TODO check that no terms are forgotten here!
    #Intitialize dictionary:
    fepdict = dict()

    #What part are we reading in the FEP file?
    section = None
    #Give ordered numbers to sections so that we can print in the same order as read!
    sectionnr = 0

    #Give ordered numbers to no_key sections so that we can print in the same order as read
    keynr = 0

    #sections without atomnumber(s)/Q number(s) etc.
    no_key = ['!info', '[angle_couplings]', '[torsion_couplings]', '[soft_pairs]']

    #How to handle keys if fepdict when modifying the: (nr of elements in key, type of modification to apply)
    key_type = {'[FEP]': (1, 'name'),
                '[atoms]': (1, qoffset),
                '[atom_types]': (1, 'name'),
                '[change_atoms]': (1, qoffset),
                '[change_charges]': (1, qoffset),
                '[bond_types]': (1, keynr),
                '[angle_types]': (1, keynr),
                '[torsion_types]': (1, keynr),
                '[improper_types]': (1, keynr),
                '[softcore]': (1, qoffset),
                '[el_scale]': (2, qoffset),
                '[change_bonds]': (2, atomoffset),
                '[changle_angles]': (3, atomoffset),
                '[change_torsions]': (4, atomoffset),
                '[change_impropers]': (4, atomoffset),
                '[excluded_pairs]': (2, atomoffset),
                '[soft_pairs]': (2, qoffset),
                '[torsion_couplings]': (2, keynr),
                '[angle_couplings]': (2, keynr),
                '[shake_constraints]': (2, atomoffset)}



    #Now read the FEP file and insert stuff to fepdict
    with open(fepfile, 'r') as fep:
        for line in fep:
            #New section?
            if re.search('\[.*]', line):
                sectionnr += 1
                section = line.split()[0]
                fepdict[sectionnr] = dict()
                fepdict[sectionnr][section] = dict()
                keynr = 0

            elif section and len(line.split()) >= 2 and not line.startswith('!'):
                keymod = key_type[section][1]

                if section in no_key:
                    keynr += 1
                    _key = keynr
                    _val = modify_fepline(keymod, line.split()[0:key_type[section][0]], section, key_type, qoffset)

                else:
                    _key = modify_fepline(keymod, line.split()[0:key_type[section][0]], section, key_type, qoffset)

                    if section == '[atoms]':
                        _val = [str(modify_fepline(atomoffset, line.split()[1:2], section, key_type, qoffset))] + \
                               line.split()[2:]
                    else:
                        _val = line.split()[key_type[section][0]:]

                fepdict[sectionnr][section][_key] = _val

            elif not section and len(line.strip('\n')) > 0:
                #Read info in top of file (if any):
                if not sectionnr in list(fepdict.keys()):
                    fepdict[sectionnr] = dict()
                if not '!info' in list(fepdict[sectionnr].keys()):
                    fepdict[sectionnr]['!info'] = dict()

                keynr += 1
                fepdict[sectionnr]['!info'][keynr] = [line.strip('\n')]

    return fepdict

def modify_fepline(keymod, modify, section, key_type, qoffset):
    """
    :param keymod: dictionary, integer or string
    :param modify: listelement from FEP file to modify
    :param section: section in FEP file f.ex [bonds]
    :param key_type: dictionary with info on FEP sections (see read_fep)
    :param qoffset: integer for qoffset
    :return: modify (modified line) as string or integer
    """
    #Make modefication to atomnumbers and Q numbers
    if isinstance(keymod, dict):

        if len(modify) == 1:
            modify = int(modify[0]) + keymod[int(modify[0])]
        else:
            modify = ' '.join(str(int(i) + keymod[int(i)]) for i in modify)

    elif isinstance(keymod, str):
        modify = ' '.join(modify)

    elif keymod == qoffset:
        if key_type[section][0] == 1:
            modify = int(modify[0]) + qoffset
        else:
            modify = ' '.join(str(int(i) + qoffset) for i in modify)

    return modify


def write_fepdict(fep, path=None, printfep=False):
    """
    Takes FEP dictionary and writes a regular FEP file to path
    :rtype : object
    :param fep: a dictionary in the format described in function read_fep in this file
    :param path: path to where the FEP file will be written
    :param printfep: print FEP=True/False. If True, file will not be written only printed.
    :return: fepprint: used for printing FEP file in Qgui main window
    """
    fepprint = list()

    if not printfep:
        fepout = open(path, 'w')

    #Some of the FEP sections are printed without the key (key is only for ordering in these cases)
    no_key = ['!info', '[angle_couplings]', '[torsion_couplings]', '[soft_pairs]']

    for nr in sorted(fep.keys()):
        #The nr is just for printing the sections in the "correct" order
        for section in sorted(fep[nr].keys()):
            if printfep:
                print(('\n%s' % section))
                fepprint.append('\n%s\n' % section)
            else:
                if section != '!info':
                    fepout.write('\n')
                fepout.write('%s\n' % section)

                #For some reason states needs to come first in [FEP]! adding this until we get a fix in Q...
                #TODO delete this if claus when fix for Q is available:
                if section == '[FEP]':
                    fepout.write('%6s %7s\n' % ('states'.ljust(6), fep[nr][section]['states'][0].ljust(7)))

            for _key in sorted(fep[nr][section].keys()):
                _val = fep[nr][section][_key]
                if section in no_key:
                    if printfep:
                        print((' '.join(['%7s' % w.ljust(7) for w in _val])))
                        fepprint.append(' '.join(['%7s' % w.ljust(7) for w in _val]))
                        fepprint.append('\n')
                    else:
                        fepout.write(' '.join(['%7s' % w.ljust(7) for w in _val]))
                        fepout.write('\n')
                else:
                    if printfep:
                        print(('%6s %s' % (str(_key).ljust(6), ' '.join(['%7s' % str(w).ljust(7) for w in _val]))))
                        fepprint.append('%6s %s' % (str(_key).ljust(6), ' '.join(['%7s' %
                                                                                  str(w).ljust(7) for w in _val])))
                        fepprint.append('\n')
                    else:
                        #TODO delete this if claus when fix for Q is available:
                        if _key != 'states':
                            fepout.write('%6s %s' %
                                         (str(_key).ljust(6), ' '.join(['%7s' % str(w).ljust(7) for w in _val])))
                            fepout.write('\n')
    if not printfep:
        fepout.close()

    else:
        return fepprint

def extend_fep_length(fep, nr_feps):
    """
    Takes a fepdict of lenth N, and extends the length to nr_feps by setting the end state (N) to all new steps.
    :param fep: fep dictionary
    :param nr_feps: how many steps should the FEP protocol be
    :return: fep
    """
    #sections to extend final state of final FEP to new FEP files (s1-->s2 in FEP6 becomes s2-->s2 in FEP7 f.ex)
    get_state = ['[change_atoms]', '[change_charges]', '[softcore]', '[el_scale]', '[change_bonds]', '[change_angles]',
                 '[change_torsions]', '[change_impropers]']

    #Last FEP step in existing dictionary:
    n = len(list(fep.keys()))
    #New number for new FEP protocols
    i = n

    #Find how many states FEP is:
    order_nr = get_fepdict_order_nr(fep[n], '[FEP]')
    s = int(fep[n][order_nr]['[FEP]']['states'][0])
    if s is None:
        s = 2

    while len(list(fep.keys())) < nr_feps:
        i += 1
        #Create new FEP file
        fep[i] = dict()

        for order_nr in sorted(fep[n].keys()):
            fep[i][order_nr] = dict()
            for section in list(fep[n][order_nr].keys()):
                fep[i][order_nr][section] = dict()
                for _key in sorted(fep[n][order_nr][section].keys()):
                    _val = fep[n][order_nr][section][_key]

                    #Append final state
                    if section in get_state:
                        final_val = fep[n][order_nr][section][_key][s - 1]
                        for alter in range(s):
                            _val[alter] = final_val

                    fep[i][order_nr][section][_key] = _val


    return fep

def get_offset(org, appending):
    """
    takes the section parts of two FEP dictinaries (feks [atoms]), and creates a translation dictionary for the
    appending FEP dictionary - for merging two FEPs.
    :param org: original FEP section
    :param appending: the FEP section that is being appended
    :return: offset: dictionary {old nr: new nr}
    """
    #{original nr: new nr}
    offset = dict()

    new_key = max(org.keys())

    for _key in sorted(appending.keys()):
        #Does value exist in org FEP section?
        if appending[_key] in list(org.values()):
            #OK, find original key for that value then
            for org_key in sorted(org.keys()):
                if org[org_key] == appending[_key]:
                    offset[_key] = org_key
                    break
        else:
            new_key += 1
            offset[_key] = new_key

    return offset

def get_fepdict_order_nr(fep, section):
    """
    Function that returns the order nr for a given section in fepdict
    :param fep1:
    :return: order_nr
    """
    order_nr = None

    for nr in sorted(fep.keys()):
        for sect in sorted(fep[nr].keys()):
            if sect == section:
                order_nr = nr
                break

    return order_nr


def merge_fep_dicts(fep1, fep2):
    """
    Will add fep2 to fep1 with modified Q nrs to depending on what is in fep1.
    :param fep1: dict with fepfiles {nr: stuff}
    :param fep2: dict with feofiles {nr: stuff}
    :return: fepdict
    merges two fep dictionaries so that nr of FEP files becomes equal!
    """
    #Check length of FEP protocols and make them equal in length
    if len(list(fep1.keys())) > len(list(fep2.keys())):
        fep2 = extend_fep_length(fep2, len(list(fep1.keys())))
    elif len(list(fep1.keys())) < len(list(fep2.keys())):
        fep1 = extend_fep_length(fep1, len(list(fep2.keys())))


    #How to handle keys if fepdict when modifying the: (nr of elements in key, type of modification to apply)
    key_type = {'[atoms]': ([0], ['q_offset']),
                '[change_atoms]': ([0], ['q_offset']),
                '[change_charges]': ([0], ['q_offset']),
                '[bond_types]': ([0], ['bonds_offset']),
                '[angle_types]': ([0], ['angles_offset']),
                '[torsion_types]': ([0], ['torsions_offset']),
                '[improper_types]': ([0], ['impropers_offset']),
                '[softcore]': ([0], ['q_offset']),
                '[el_scale]': ([0, 1], ['q_offset', 'q_offset']),
                '!info': ([0], ['info_offset'])
                }

    #How to handle values in fepdict when modifying ...
    val_type = {'[change_bonds]': ([0, 1], ['bonds_offset', 'bonds_offset']),
                '[changle_angles]': ([0, 1], ['angles_offset', 'angles_offset']),
                '[change_torsions]': ([0, 1], ['torsions_offset', 'torsions_offset']),
                '[change_impropers]': ([0, 1], ['impropers_offset', 'impropers_offset']),
                '[soft_pairs]': ([0, 1], ['q_offset', 'q_offset']),
                '[torsion_couplings]': ([0, 1], ['bonds_offset', 'torsions_offset']),
                '[angle_couplings]': ([0, 1], ['bonds_offset', 'angles_offset'])}

    for fep in sorted(fep2.keys()):
        #(re)-initialize translation dictionaries for FEP2
        offsets = {'q_offset': None,
               'bonds_offset': None,
               'angles_offset': None,
               'torsions_offset': None,
               'impropers_offset': None,
               'info_offset': None}

        for order_nr in sorted(fep2[fep].keys()):
            for section in sorted(fep2[fep][order_nr].keys()):
                #Ensure that section exists in fep1, it may not when merging different FEPs!
                insert_nr = get_fepdict_order_nr(fep1[fep], section)
                #Maybe this section did not exist in fep1
                if insert_nr is None:
                    if not order_nr in list(fep1[fep].keys()):
                        insert_nr = order_nr
                    else:
                        insert_nr = max(fep1[fep].keys()) + 1
                    fep1[fep][insert_nr] = dict()
                    fep1[fep][insert_nr][section] = dict()

                for _key in sorted(fep2[fep][order_nr][section].keys()):
                    _val = fep2[fep][order_nr][section][_key]

                    #info section?
                    #if section == '!info':
                    #    _key = len(fep1[fep][insert_nr][section].keys()) + 1

                    #Modify key?
                    if section in list(key_type.keys()):
                        new_key = list()
                        for i in range(len(key_type[section][0])):

                            #Do we have the correct offset for this section for this FEP file?
                            if not offsets[key_type[section][1][i]]:
                                appending = fep2[fep][order_nr][section]
                                #find section in in fep1
                                org = fep1[fep][insert_nr][section]

                                if len(list(org.keys())) < 1:
                                    org = appending

                                offsets[key_type[section][1][i]] = get_offset(org, appending)

                            if len(key_type[section][0]) > 1:
                                new_key.append(offsets[key_type[section][1][i]][_key.split()[key_type[section][0][i]]])
                            else:
                                new_key.append(offsets[key_type[section][1][i]][_key])

                        if len(new_key) > 1:
                            _key = ' '.join(str(int(i)) for i in new_key)
                        else:
                            _key = new_key[0]

                    #Modify value?
                    elif section in list(val_type.keys()):
                        new_val = list()
                        for i in range(len(val_type[section][0])):

                            #Do we have the correct offset for this section for this FEP file?
                            if not offsets[val_type[section][1][i]]:
                                appending = fep2[fep][order_nr][section]

                                #Find section in fep1
                                org = fep1[fep][insert_nr][section]

                                if len(list(org.keys())) < 1:
                                    org = appending

                                offsets[val_type[section][1][i]] = get_offset(org, appending)

                            new_val.append(offsets[val_type[section][1][i]][int(_val[val_type[section][0][i]])])

                        #Check if additional stuff exist in val part
                        if len(_val) > new_val:
                            for j in range(len(new_val), len(_val)):
                                new_val.append(_val[j])

                        #_val = ' '.join(str(i) for i in new_val)
                        _val = new_val

                    fep1[fep][insert_nr][section][_key] = _val

    return fep1

def get_md_settings(for_what='MD'):
    """
    :param: forwhat: keyword MD/FEP/EVB/resFEP/LIE used to invoke special settings
    :return:
    """
    #Molecular dynamics settings:
    md_settings = {'simtime': 0.01,
                   'stepsize': 1.0,
                   'inputfiles': 51,
                   'temperature': 'T_VAR',
                   'bath_coupling': 100,
                   'shake_solvent': 1,
                   'shake_solute': 0,
                   'shake_hydrogens': 0,
                   'lrf': 1,
                   'solute_solute_cut': 10,
                   'solute_solvent_cut': 10,
                   'solvent_solvent_cut': 10,
                   'q_atoms_cut': 99,
                   'lrf_cut': 99,
                   'shell_force': 10.0,
                   'shell_rad': 30.0,
                   'radial_force': 60.0,
                   'polarisation': 1,
                   'pol_force': 20.0,
                   'nonbond_list': 25,
                   'ene_summary': 5,
                   'ene_file': 10,
                   'trajectory': 1000,
                   'trajectory atoms': 'not excluded',
                   'seq_rest': [],
                   'atom_rest': [],
                   'dist_rest': [],
                   'wall_rest': [],
                   'fep_file': 'FEP_VAR'}

    if for_what == 'resFEP':
        md_settings['shake_hydrogens'] = 0

    return md_settings

def create_lambda_list(lambda_step, lambda_start=list(), lambda_end=list()):
    """
    creates a lambda list as they are seen in the Gui for FEP and EVB setup (nr l1 l2 l3 l4).
    Note that only one lambda-pair can be changed at the time when more than 2 states!
    :param lambda_step: lambda step size
    :param lambda_start: list with initial lambda values
    :param lambda_end: list with final lambda values
    :return: list of nr lambda values
    """
    if sum(lambda_start) != 1 or sum(lambda_end) != 1:
        return

    #Find out what values to change
    change_values = list()
    for i in range(len(lambda_start)):
        if lambda_start[i] != lambda_end[i]:
            change_values.append(i)

    if len(change_values) > 2:
        print('Can only vary two lambdas at the same time')
        print(('Will change lambda %d and lambda %d' % (change_values[0] + 1, change_values[1] + 1)))

    change_values = (change_values[0], change_values[1])

    #Decide sign on lambda step:
    if lambda_start[change_values[0]] > lambda_end[change_values[0]]:
        lambda_step = -lambda_step

    lambda_list = list()

    lambda1 = np.arange(lambda_start[change_values[0]], lambda_end[change_values[0]] + lambda_step,
                                 lambda_step)

    for i in range(len(lambda1)):
        tmp = lambda_start[:]
        tmp[change_values[0]] = abs(lambda1[i])
        tmp[change_values[1]] = 0
        tmp[change_values[1]] = abs(1. - sum(tmp))

        #it is list of strings: '%4d %s %s %s %s' % (i, l1, l2, l3, l4)
        # where %s comes from '%03.3f'
        lambda_list.append('%4d %s' % (i + 1, ' '.join(['%03.3f' % w for w in tmp])))

    return lambda_list


def write_md_inputfiles(inputfiles, md_settings, topology, lambda_list, qdyn, eq=None, submissionscript=False,
                        resFEP=False, pdbfile=None):
    """
    :param inputfiles: PATH to inputfiles
    :param md_settings: dictionary with MD file settings
    :param topology: path to topology
    :param fepname: Name of fepfile
    :param lambda_list: list of lambda values (floats)
    :param submissionscript: list with lines to put in submission script!
    :param qdyn: 'Qdyn5' or 'mpirun Qdyn5p' f. ex
    :param resFEP: True/False In resFEP special handling is needed for the equilibration
    :param pdbfile: pdb file in list format (optional)
    :return: some kind of message...
    """
    #TODO write fep file name into md_settings dict!

    #Make inputfiles directory and move to directory:
    if not os.path.exists(inputfiles):
        os.makedirs(inputfiles)

    #Open submission script and write head:
    submit = '%s/run.sh' % inputfiles
    submitfile = open(submit,'w')
    if not submissionscript:
        submissionscript = ['#!/bin/bash\n', '#Qdyn I/O\n']

    for line in submissionscript:
        if '#Qdyn I/O' in line:
            break
        else:
            submitfile.write(line)

    if not pdbfile:
        pdbfile = create_pdb_from_topology(topology)

    #Find atom nrs for solute and solvent:
    solute = []
    solvent = []
    all_atoms = []
    for line in pdbfile:
        if 'ATOM' in line:
            if 'HOH' in line or 'SPC' in line:
                solvent.append(line.split()[1])
            else:
                solute.append(line.split()[1])
            all_atoms.append(line.split()[1])

    #Get global settings
    on_off = {0: 'off', 1: 'on'}

    output_int = md_settings['ene_summary']
    trj_int = md_settings['trajectory']
    ene_int = md_settings['ene_file']
    non_bond_int = md_settings['nonbond_list']
    radial_force = md_settings['radial_force']
    polarisation = md_settings['pol_force']
    shell_force = md_settings['shell_force']
    shell_radius = md_settings['shell_rad']
    md_steps = int(round((float(md_settings['simtime']) * 1000000.00)/ float(md_settings['stepsize'])))
    md_stepsize = md_settings['stepsize']
    md_temp = md_settings['temperature']
    md_bath = md_settings['bath_coupling']
    shake_solvent = on_off[md_settings['shake_solvent']]
    shake_solute = on_off[md_settings['shake_solute']]
    shake_hydrogens = on_off[md_settings['shake_hydrogens']]
    lrf = on_off[md_settings['lrf']]
    lrf_cutoff = md_settings['lrf_cut']
    use_pol = on_off[md_settings['polarisation']]
    fepname = md_settings['fep_file']

    #Write default equilibration procedure
    #Leave base_name as md (names are long enough with the lambda values included)
    base_name = 'eq'
    random_seed = 'SEED_VAR'
    count = 0

    if eq:
        restart_file = None
        if resFEP:
            submitfile.write('if [ $index -lt 1 ]; then\n')
        for i in range(len(eq)):
            count += 1
            inputfile = base_name + '%d.inp' % count
            logfile = base_name + '%d.log' %count
            eq_file = open('%s/%s' % (inputfiles, inputfile), 'w')

            if eq[i][0] == 'End':
                temp = md_temp
            else:
                temp = eq[i][0]

            submitfile.write('%s %s > %s\n' % (qdyn, inputfile, logfile))

            eq_file.write('[MD]\n')
            eq_file.write('%25s %s\n' % ('steps'.ljust(25), eq[i][5]))
            eq_file.write('%25s %s\n' % ('stepsize'.ljust(25), eq[i][4]))
            eq_file.write('%25s %s\n' % ('temperature'.ljust(25), temp))
            eq_file.write('%25s %s\n' % ('bath_coupling'.ljust(25), eq[i][1]))
            if count == 1:
                eq_file.write('%25s %s\n' % ('random_seed'.ljust(25), random_seed))
                eq_file.write('%25s %s\n' % ('initial_temperature'.ljust(25), eq[i][0]))
                eq_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), 'on'))
            if count > 1:
                eq_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), shake_solvent))

            eq_file.write('%25s %s\n' % ('shake_hydrogens'.ljust(25), shake_hydrogens))
            eq_file.write('%25s %s\n' % ('shake_solute'.ljust(25), shake_solute))
            eq_file.write('%25s %s\n' % ('lrf'.ljust(25), lrf))

            eq_file.write('\n[cut-offs]\n')
            eq_file.write('%25s %s\n' % ('solute_solvent'.ljust(25), md_settings['solute_solvent_cut']))
            eq_file.write('%25s %s\n' % ('solute_solute'.ljust(25), md_settings['solute_solute_cut']))
            eq_file.write('%25s %s\n' % ('solvent_solvent'.ljust(25), md_settings['solvent_solvent_cut']))
            eq_file.write('%25s %s\n' % ('q_atom'.ljust(25), md_settings['q_atoms_cut']))
            if lrf == 'on':
                eq_file.write('%25s %s\n' % ('lrf'.ljust(25).ljust(25), lrf_cutoff))

            eq_file.write('\n[sphere]\n')
            eq_file.write('%25s %s\n' % ('shell_force'.ljust(25), shell_force))
            eq_file.write('%25s %s\n' % ('shell_radius'.ljust(25), shell_radius))

            eq_file.write(('\n[solvent]\n'))
            eq_file.write('%25s %s\n' % ('radial_force'.ljust(25), radial_force))
            eq_file.write(('%25s %s\n') % ('polarisation'.ljust(25), use_pol))
            if use_pol == 'on':
                eq_file.write('%25s %s\n' % ('polarisation_force'.ljust(25), polarisation))

            eq_file.write('\n[intervals]\n')
            eq_file.write('%25s %s\n' % ('output'.ljust(25), output_int))
            eq_file.write('%25s %s\n' % ('trajectory'.ljust(25), trj_int))
            eq_file.write('%25s %s\n' % ('non_bond'.ljust(25), non_bond_int))

            topology = topology.split('/')[-1]
            fepname = fepname.split('/')[-1]

            eq_file.write('\n[files]\n')
            eq_file.write('%25s %s\n' % ('topology'.ljust(25), topology))
            eq_file.write('%25s %s%d.dcd\n' % ('trajectory'.ljust(25), base_name, count))
            if count != 1:
                eq_file.write('%25s %s\n' % ('restart'.ljust(25), restart_file))
            eq_file.write('%25s %s%d.re\n' % ('final'.ljust(25), base_name, count))
            eq_file.write('%25s %s\n' % ('fep'.ljust(25), fepname))
            restart_file = logfile.split('.')[0]+'.re'

            eq_file.write('\n[trajectory_atoms]\n')
            eq_file.write('%s\n' % md_settings['trajectory atoms'])

            eq_file.write('\n[lambdas]\n')
            eq_file.write('%s\n' % ' '.join(lambda_list[0].split()[1:]))

            eq_file.write('\n[sequence_restraints]\n')
            force = float(eq[i][3])
            if eq[i][2] != 'None':
                if eq[i][2] == 'All':
                    atomlist = all_atoms
                    #eq_file.write(' not excluded  %4.1f 0  0\n' % force)
                elif eq[i][2] == 'Solute':
                    atomlist = solute
                elif eq[i][2] == 'Solvent':
                    atomlist = solvent
                try:
                    atom_i = atomlist[0]
                    atom_j = atomlist[-1]
                    eq_file.write('%6s %6s %4.1f 0  0\n' % (atom_i, atom_j, force))
                except:
                    continue

            if len(md_settings['seq_rest']) > 0:
                for restraint in md_settings['seq_rest']:
                    eq_file.write('%s\n' % restraint)

            if len(md_settings['dist_rest']) > 0:
                eq_file.write('\n[distance_restraints]\n')
                for restraint in md_settings['dist_rest']:
                    eq_file.write('%s\n' % restraint)

            if len(md_settings['atom_rest']) > 0:
                eq_file.write('\n[atom_restraints]\n')
                for restraint in md_settings['atom_rest']:
                    eq_file.write('%s\n' % restraint)

            if len(md_settings['wall_rest']) > 0:
                eq_file.write('\n[wall_restraints]\n')
                for restraint in md_settings['wall_rest']:
                    eq_file.write('%s\n' % restraint)

        if resFEP:
            submitfile.write('fi\n')

    #Write MD inputfiles for lambda steps:
    for lamda in lambda_list:
        lambda_value = ' '.join(lamda.split()[1:])
        base_name = 'md_'
        l = lambda_value.split()
        for i in range(0, len(l)):
            base_name += ''.join(l[i].split('.'))
            if i < len(l) - 1:
                base_name += '_'

        inputfile =  base_name + '.inp'
        logfile = base_name + '.log'
        md_file = open('%s/%s' % (inputfiles, inputfile), 'w')

        submitfile.write('%s %s > %s\n' % (qdyn, inputfile, logfile))

        md_file.write('[MD]\n')
        md_file.write('%25s %s\n' % ('steps'.ljust(25), md_steps))
        md_file.write('%25s %s\n' % ('stepsize'.ljust(25), md_stepsize))
        md_file.write('%25s %s\n' % ('temperature'.ljust(25), md_temp))
        md_file.write('%25s %s\n' % ('bath_coupling'.ljust(25), md_bath))

        md_file.write('%25s %s\n' % ('shake_hydrogens'.ljust(25), shake_hydrogens))
        md_file.write('%25s %s\n' % ('shake_solute'.ljust(25), shake_solute))
        md_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), shake_solvent))
        md_file.write('%25s %s\n' % ('lrf'.ljust(25), lrf))

        md_file.write('\n[cut-offs]\n')
        md_file.write('%25s %s\n' % ('solute_solvent'.ljust(25), md_settings['solute_solvent_cut']))
        md_file.write('%25s %s\n' % ('solute_solute'.ljust(25), md_settings['solute_solute_cut']))
        md_file.write('%25s %s\n' % ('solvent_solvent'.ljust(25), md_settings['solvent_solvent_cut']))
        md_file.write('%25s %s\n' % ('q_atom'.ljust(25), md_settings['q_atoms_cut']))
        if lrf == 'on':
            md_file.write('%25s %s\n' % ('lrf'.ljust(25).ljust(25), lrf_cutoff))

        md_file.write('\n[sphere]\n')
        md_file.write('%25s %s\n' % ('shell_force'.ljust(25), shell_force))
        md_file.write('%25s %s\n' % ('shell_radius'.ljust(25), shell_radius))

        md_file.write(('\n[solvent]\n'))
        md_file.write('%25s %s\n' % ('radial_force'.ljust(25), radial_force))
        md_file.write(('%25s %s\n') % ('polarisation'.ljust(25), use_pol))
        if use_pol == 'on':
            md_file.write('%25s %s\n' % ('polarisation_force'.ljust(25), polarisation))

        md_file.write('\n[intervals]\n')
        md_file.write('%25s %s\n' % ('output'.ljust(25), output_int))
        md_file.write('%25s %s\n' % ('energy'.ljust(25), ene_int))
        md_file.write('%25s %s\n' % ('trajectory'.ljust(25), trj_int))
        md_file.write('%25s %s\n' % ('non_bond'.ljust(25), non_bond_int))

        md_file.write('\n[files]\n')

        topology = topology.split('/')[-1]
        fepname = fepname.split('/')[-1]

        md_file.write('%25s %s\n' % ('topology'.ljust(25), topology))
        md_file.write('%25s %s.dcd\n' % ('trajectory'.ljust(25), base_name))
        md_file.write('%25s %s\n' % ('restart'.ljust(25), restart_file))
        md_file.write('%25s %s.en\n' % ('energy'.ljust(25), base_name))
        md_file.write('%25s %s.re\n' % ('final'.ljust(25), base_name))
        md_file.write('%25s %s\n' % ('fep'.ljust(25), fepname))
        restart_file = logfile.split('.')[0]+'.re'

        md_file.write('\n[trajectory_atoms]\n')
        md_file.write('%s\n' % md_settings['trajectory atoms'])

        md_file.write('\n[lambdas]\n')
        md_file.write('%s\n' % lambda_value)

        if len(md_settings['seq_rest']) > 0:
            md_file.write('\n[sequence_restraints]\n')
            for restraint in md_settings['seq_rest']:
                md_file.write('%s\n' % restraint)

        if len(md_settings['dist_rest']) > 0:
            md_file.write('\n[distance_restraints]\n')
            for restraint in md_settings['dist_rest']:
                md_file.write('%s\n' % restraint)

        if len(md_settings['atom_rest']) > 0:
            md_file.write('\n[atom_restraints]\n')
            for restraint in md_settings['atom_rest']:
                md_file.write('%s\n' % restraint)

        if len(md_settings['wall_rest']) > 0:
            md_file.write('\n[wall_restraints]\n')
            for restraint in md_settings['wall_rest']:
                md_file.write('%s\n' % restraint)
        md_file.close()

    #If use submission script, check for end statements (comes after #Qdyn I/O):
    write_end = False
    end_statements_start = None
    for k in range(len(submissionscript)):
        if '#Qdyn I/O' in submissionscript[k]:
            end_statements_start = k + 1
            write_end = True
    if write_end:
        for line in range(end_statements_start, len(submissionscript)):
            submitfile.write(submissionscript[line])

    submitfile.close()

    print(('info', 'FEP templates written to "/inputfiles"'))


def get_qnr_atomnr(fepdict):
    """
    Creates a dictionary {Q nr : atom nr} from fepdict
    :param fepdict:
    :return:
    """
    qnr_atomnr = dict()

    order_nr = get_fepdict_order_nr(fepdict, '[atoms]')

    for qnr in sorted(fepdict[order_nr]['[atoms]'].keys()):
        qnr_atomnr[qnr] = int(fepdict[order_nr]['[atoms]'][qnr][0])

    return qnr_atomnr


def get_ene_files(qpath):
    """
    :param qpath: direcory to look for .en files
    :return: enefiles
    """
    enefiles = list()

    for f in sorted(os.listdir(qpath)):
        if f.endswith('.en'):
            if f.startswith('md'):
                enefiles.append(f)

    print(('found %d energy files in in %s' % (len(enefiles), qpath)))
    return enefiles


def write_qfep_input(qfep_path, states=2, temp=298, linear_comb='1 0', skip=100, bins=50, minbinpts=50,
                     a=None, hij=None):
    """
    Write Qfep input file in qfep_path
    :param qfep_path: directory where qfep.inp should be written
    :param states: numbers of FEP states
    :param temp: Temperature
    :param linear_comb: linear combination of states --> '1 0' / '1 -1' / '-1 1'
    :param skip: numbers of bin points to skip
    :param bins: How many bins to use (energy gaps)
    :param minbinpts: minumum number of points in every bin
    :param a: list of alpha value(s) for EVB --> [alpha1, alpha2, ...]
    :param hij: list of off-diagonal value(s) for EVB --> ['1  2  Hij  0.0 0 0.000', ...]
    :return: True/False (something went right or wrong)
    """
    enefiles = get_ene_files(qfep_path)

    if len(enefiles) < 0:
        return False

    inputfile = open('%s/%s' % (qfep_path, 'qfep.inp'), 'w')

    inputfile.write('%d\n' % len(enefiles))
    inputfile.write('%d  0\n' % states)
    inputfile.write('%.4f  %d\n' % (0.001987209 * temp, skip))
    inputfile.write('%d\n' % bins)
    inputfile.write('%d\n' % minbinpts)

    #Write alpha values for state(s)
    for i in range(states - 1):
        if a:
            inputfile.write('%s\n' % str(a[i]))
        else:
            inputfile.write('0\n')

    #Number of off-diagonal elements:
    if hij:
        inputfile.write('%d\n' % (states - 1))
        for i in range(states - 1):
            inputfile.write(hij[i])
    else:
        inputfile.write('0\n')

    #Linear combination
    inputfile.write('%s\n' % linear_comb)

    #Write lambda (1,0) --> (0,1)
    enefiles.reverse()
    for ene in enefiles:
        inputfile.write('%s\n' % ene)

    inputfile.close()

    return True


def run_Qfep(qpath, qfep_inp='qfep.inp', qfep='Qfep5'):
    """
    Calls Qfep with qfep_inp and produces qfep.out
    :param qpath: path to run
    :param qfepin:
    :param qfep:
    :return: True/False
    """
    org_path = os.getcwd()

    if org_path != qpath:
        os.chdir(qpath)

    if not os.path.isfile(qfep_inp):
        print(('No Qfep inputfile (%s) in %s' % (qfep_inp, qpath)))
    else:
        tmpfile = open('qfep.out', 'w')

        call('%s <%s' % (qfep, qfep_inp), shell=True, stdout=tmpfile, stderr=tmpfile)

        tmpfile.close()

    if org_path != qpath:
        os.chdir(org_path)

def get_qfep_part1(qpath, qfep='qfep.out'):
    """
    Gets part 1 from qfep.out and returns a dictionary with data
    :param qpath: path to Qfep output file
    :param qfep: name of file generated by Qfep
    :return: part1: dict
    """
    part1 = {'lambda': list(), 'sum_dGf': list(), 'sum_dGr': list(), 'dG': list() }

    org_path = os.getcwd()

    if org_path != qpath:
        os.chdir(qpath)

    if not os.path.isfile(qfep):
        return False

    with open(qfep, 'r') as qfep_out:
        found_part1 = False
        for line in qfep_out:
            if line == 'Qfep5 terminated abnormally: Failed to read energies.\n':
                print(('Qfep failed to read energies in path: %s' % os.getcwd()))
                os.remove(qfep)
                try:
                    os.remove('qfep.inp')
                except:
                    continue
                return None
            elif 'ERROR:' in line:
                print('Unexpected end of .en file in first record.\nMaybe job is not done?')
                os.remove(qfep)
                try:
                    os.remove('qfep.inp')
                except:
                    continue
                return
            if found_part1:
                if len(line.split()) < 6:
                    break
                if not line.startswith('#'):
                    part1['lambda'].append(float(line.split()[0]))
                    part1['sum_dGf'].append(float(line.split()[2]))
                    part1['sum_dGr'].append(float(line.split()[4]))
                    part1['dG'].append(float(line.split()[5]))
            if '# lambda' in line and len(line.split()) > 5:
                found_part1 = True


    if org_path != qfep_out:
        os.chdir(org_path)

    return part1
