#!/usr/bin/env python

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
                            print 'Warning: Did not find library file: %s' % lib

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
                    atomnr_xyz[atomnr] = map(float, line.split()[0:3])

                if len(line.split()) == 6:
                    atomnr += 1
                    atomnr_xyz[atomnr] = map(float, line.split()[3:6])

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
                '[angle_couplings]': (2, keynr)}



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
                if not sectionnr in fepdict.keys():
                    fepdict[sectionnr] = dict()
                if not '!info' in fepdict[sectionnr].keys():
                    fepdict[sectionnr]['!info'] = dict()

                keynr += 1
                fepdict[sectionnr]['!info'][keynr] = [line]

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

    elif isinstance(keymod, basestring):
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
    :param fep: a dictionary in the format described in function read_fep in this file
    :param path: path to where the FEP file will be written
    :param printfep: print FEP=True/False. If True, file will not be written only printed.
    :return: True/False maybe?
    """

    if not printfep:
        fepout = open(path, 'w')

    #Some of the FEP sections are printed without the key (key is only for ordering in these cases)
    no_key = ['!info', '[angle_couplings]', '[torsion_couplings]', '[soft_pairs]']

    for nr in sorted(fep.keys()):
        #The nr is just for printing the sections in the "correct" order
        for section in sorted(fep[nr].keys()):
            if printfep:
                print('\n%s' % section)
                for _key in sorted(fep[nr][section].keys()):
                    _val = fep[nr][section][_key]
                    if section in no_key:
                        if printfep:
                            print ' '.join(['%7s' % w.ljust(7) for w in _val])
                        else:
                            fepout.write(' '.join(['%7s' % w.ljust(7) for w in _val]))
                    else:
                        if printfep:
                            print '%6s %s' % (str(_key).ljust(6), ' '.join(['%7s' % w.ljust(7) for w in _val]))
                        else:
                            fepout.write('%6s %s' % (str(_key).ljust(6), ' '.join(['%7s' % w.ljust(7) for w in _val])))

    fepout.close()

def merge_fep_dicts(fep1, fep2):
    """
    :param fep1: dict with fepfiles {nr: stuff}
    :param fep2: dict with feofiles {nr: stuff}
    :return: fepdict
    merges two fep dictionaries so that nr of FEP files becomes equal!
    """
    merged = dict()

    #Are FEP dicts of same length, just merge them:
    if len(fep1.keys()) == len(fep2.keys()):
        merged = fep1
        addon = fep2

    #If length is unequal, find which one is largest
    print fep1



