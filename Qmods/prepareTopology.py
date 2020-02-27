#!/usr/bin/env python2

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
Functions to prepare topolgy input file for Qprep5.
"""

import numpy as np
import sys
import re


def checkPdb(pdb):
    """
   Checks if pdb file is Q pdb or standard pdb format.
   Returns True if Q pdb or False if standard pdb.
   """
    print(pdb)
    try:
        pdbFile = open(pdb, 'r').readlines()
    except:
        print('Could not open file %s')
        pass
    qpdb = False
    if 'ATOM' in pdbFile[0].split()[0] or 'HETATM' in pdbFile[0].split()[0]:
        if len(pdbFile[0]) < 58:
            qpdb = True

    return qpdb


def convertPdb(pdb, workdir, atomnr=0, resnr=0):
    """
    Takes a standard PDB file format and converts it to Q PDB format.
    Automatically inserts GAPs where possible to interpret.
    :param pdb:
    :param atomnr:
    :param resnr:
    """

    pdbfile = open(pdb, 'r').readlines()
    residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TYR', 'TRP', 'MET', 'CYS', 'THR', 'GLU', 'ASP', 'LYS', 'ARG',
                'HIS', 'SER', 'PRO', 'ASN', 'GLN']

    newName = pdb.split('/')[-1].split('.')[0] + 'Q.pdb'

    newPdb = open(workdir + '/' + newName, 'w')
    #newPdb = open(newName, 'w')

    current_resnr = 0
    current_resletter = 'C'
    previousAtomtype = 'dummy'

    found_ter = False

    for line in pdbfile:
        try:
            if 'ATOM' in line.split()[0] or 'HETATM' in line.split()[0]:
                if atomnr == 'get':                          #Get atom nr from pdb file
                    atomnr = int(line[6:11])
                else:
                    atomnr += 1
                atomtype = line[11:17].strip().split()[0] #PDB atomtype
                if len(atomtype) == 4:
                    if atomtype[3] == 'A' or atomtype[3] == 'B':
                        atomtype = atomtype[0:3]
                resname = line[17:21].strip()             #Residue name
                if resnr == 'get':
                    resnr = int(line[22:26])               #Get residue nr from pdb file
                else:
                    if current_resnr != int(line[22:26]):
                        #Check if GAP shoud be inserted:
                        if current_resnr + 1 != int(line[22:26]) and current_resnr != 0:
                            if atomtype == 'N':
                                #x_n = float(line[28:].split()[0])
                                x_n = float(line[30:38])
                                #y_n = float(line[28:].split()[1])
                                y_n = float(line[38:46])
                                #z_n = float(line[28:].split()[2])
                                z_n = float(line[46:54])
                                distance = np.sqrt((x_n - x_c) ** 2 + (y_n - y_c) ** 2 + (z_n - z_c) ** 2)
                                if distance > 2.2:
                                    if not found_ter:
                                        newPdb.write('%20s\n' % 'GAP')

                            #elif atomtype != 'O':
                            #    newPdb.write('%20s\n' % 'GAP')
                        current_resnr = int(line[22:26])
                        current_resletter = line[26:27]
                        resnr += 1
                        found_ter = False
                    elif (current_resletter != line[26:27]) and (current_resletter != 'C'):
                        current_resletter = line[26:27]
                        resnr += 1
                #xcord = float(line[28:].split()[0])
                xcord = float(line[30:38])
                #ycord = float(line[28:].split()[1])
                ycord = float(line[38:46])
                #zcord = float(line[28:].split()[2])
                zcord = float(line[46:54])
                if atomtype == 'C':
                    x_c, y_c, z_c = xcord, ycord, zcord
                newline = '%6s%5d  %4s%4s%5d    %8.3f%8.3f%8.3f' % ('ATOM  ', atomnr, atomtype.ljust(4), resname.ljust(4), resnr, xcord, ycord, zcord )
                if resname in residues:
                    if previousAtomtype != atomtype:
                        previousAtomtype = atomtype
                        newPdb.write('%s\n' % newline)
                else:
                    newPdb.write('%s\n' % newline)
            elif line.startswith('TER'):
                newPdb.write('%20s\n' % 'GAP')
                found_ter = True

        except:
            continue

    newPdb.close()
    return newName


def findTerminals(pdb):
    """
    Function to automatically find N-terminals and C-terminals i pdb file for first residue and
    before/after GAPS or if OXT is found in residue. Returns list of residue numbers and N-term/C-term resname:
    """

    residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TYR', 'TRP', 'MET', 'CYS', 'CYX', 'THR',
                'GLU', 'GL-', 'ASP', 'AS-', 'LYS', 'LY+', 'ARG', 'AR+', 'HID', 'HIE', 'HIP', 'SER', 'PRO', 'ASN', 'GLN']

    nt_residues = ['NGLY', 'NALA', 'NVAL', 'NLEU', 'NILE', 'NPHE', 'NTYR', 'NTRP', 'NMET', 'NCYS', 'NCYX', 'NTHR',
                'NGLU', 'NGL-', 'NASP', 'NAS-', 'NLYS', 'NLY+', 'NARG', 'NAR+', 'NHID', 'NHIE', 'NHIP', 'NSER', 'NPRO', 'NASN', 'NGLN']

    ct_residues = ['CGLY', 'CALA', 'CVAL', 'CLEU', 'CILE', 'CPHE', 'CTYR', 'CTRP', 'CMET', 'CCYS', 'CCYX', 'CTHR',
                'CGLU', 'CGL-', 'CASP', 'CAS-', 'CLYS', 'CLY+', 'CARG', 'CAR+', 'CHID', 'CHIE', 'CHIP', 'CSER', 'CPRO', 'CASN', 'CGLN']

    pdbfile = open(pdb, 'r').readlines()

    new_nter = True
    nter_nr = []
    nter_res = []
    cter_nr = []
    cter_res = []

    current_resnr = 'nothing'

    for i in range(len(pdbfile)):
        try:
            if 'ATOM' in pdbfile[i].split()[0] or 'HETATM' in pdbfile[i].split()[0]:
                res = pdbfile[i][17:21].strip()
                resnr = pdbfile[i][21:26]
                try:
                    if pdbfile[i + 1]:
                        pass
                except:
                    if not new_nter:
                        if res in residues:
                            cter_nr.append(resnr)
                            cter_res.append(res)
                            new_nter = True

                if current_resnr != resnr:
                    current_resnr = resnr
                    if new_nter:
                        if res in nt_residues:
                            nter_nr.append(current_resnr)
                            nter_res.append(res)
                            new_nter = False
                        elif res in residues:
                            nter_nr.append(resnr)
                            nter_res.append(res)
                            new_nter = False

                    if res in ct_residues:
                        cter_nr.append(current_resnr)
                        cter_res.append(res)
                        new_nter = True

            elif 'GAP' in pdbfile[i] or 'TER' in pdbfile[i]:
                res = pdbfile[i - 1][17:21].strip()
                resnr = pdbfile[i - 1][21:26]
                if res in residues:
                    cter_nr.append(resnr)
                    cter_res.append(res)
                new_nter = True
        except:
            continue



    return nter_nr, nter_res, cter_nr, cter_res


def convertCterminal(pdb, res_nr):
    """
    Takes a list of residue numbers and a list of new terminal
    residue names and changes it in the curren pdb file
    """

    residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TYR', 'TRP', 'MET', 'CYS', 'CYX', 'THR',
                'GLU', 'GL-', 'ASP', 'AS-', 'LYS', 'LY+', 'ARG', 'AR+', 'HID', 'HIE', 'HIP', 'SER', 'PRO', 'ASN', 'GLN']


    ct_residues = ['CGLY', 'CALA', 'CVAL', 'CLEU', 'CILE', 'CPHE', 'CTYR', 'CTRP', 'CMET', 'CCYS', 'CCYX', 'CTHR',
                'CGLU', 'CGL-', 'CASP', 'CAS-', 'CLYS', 'CLY+', 'CARG', 'CAR+', 'CHID', 'CHIE', 'CHIP', 'CSER', 'CPRO', 'CASN', 'CGLN']

    pdbfile = open(pdb, 'r').readlines()

    residue = []
    oxt_exist = False

    for line in pdbfile:
        try:
            if 'ATOM' in line or 'HETATM' in line:
                current_resnr = line[21:26]
                if current_resnr == res_nr:
                    residue.append(line)
                    if 'OXT' in line:
                        oxt_exist = True
        except:
            continue

    resname = residue[1][17:21].strip()

    create_OT2 = True

    if resname in residues:
        old = residues
        new = ct_residues
    else:
        old = ct_residues
        new = residues
        create_OT2 = False

    newresname = 'CTER'
    for i in range(len(old)):
        if resname == old[i]:
            newresname = new[i]

    if newresname == 'CTER':
        print('Could not toggle C-terminal state ...')
        return


    #Generate coordinates for OT2 and change O to OT1
    if create_OT2:
        if not oxt_exist:
            for i in range(len(residue)):
                if residue[i].split()[2] == 'O':
                    ox = float(residue[i].split()[5])
                    oy = float(residue[i].split()[6])
                    oz = float(residue[i].split()[7])
                    #newline = '%s %4s %s' % (residue[i][0:11],' OT1',residue[i][17:])
                    #del residue[i]
                    #residue.insert(i, newline)
                if residue[i].split()[2] == 'C':
                    cx = float(residue[i].split()[5])
                    cy = float(residue[i].split()[6])
                    cz = float(residue[i].split()[7])
            try:
                co_x = ox - cx
                co_y = oy - cy
                co_z = oz - cz
            except:
                print('something went wrong with terminal generations')
                return

            r1 = np.sqrt(co_x**2 + co_y**2 + co_z**2)
            theta = np.arccos(co_z / r1)
            psi = np.arctan(co_y / co_x)

            theta2 = 2 * theta
            psi2 = 2 * psi

            x = (r1 * np.sin(theta2) * np.cos(psi2)) + cx
            y = (r1 * np.sin(theta2) * np.sin(psi2)) + cy
            z = (r1 * np.cos(theta2)) + cz

            atomnr = int(residue[-1].split()[1]) + 1

            residue.append('%s%5d  OXT %4s%5s    %8.3f%8.3f%8.3f\n' % ('ATOM  ', atomnr, newresname.ljust(4), res_nr,x,y,z ))

    if not create_OT2:
        for i in range(len(residue)):
            if 'OXT' in residue[i]:
                del residue[i]

    #Write new pdbfile:
    newfile = open(pdb, 'w')
    inserted_terminal = False
    atomnr = 0
    for line in pdbfile:
        try:
            current_resnr = line[21:26]
            if current_resnr != res_nr:
                if 'GAP' in line:
                    newfile.write(line)
                else:
                    atomnr += 1
                    newfile.write('%s%5d%s' % (line[0:6], atomnr,line[11:]))
            if current_resnr == res_nr:
                if not inserted_terminal:
                    for newline in residue:
                        atomnr += 1
                        newfile.write('%s%5d%s%4s%s' % (newline[0:6], atomnr,newline[11:17], newresname.ljust(4), newline[21:]))
                    inserted_terminal = True
                else:
                    continue
        except:
            continue

    newfile.close()

    return newresname


def convertNterminal(pdb, res_nr = ' '):
    """
    Takes a  residue number and convrts the N-terminal on/off
    """

    residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TYR', 'TRP', 'MET', 'CYS', 'CYX', 'THR',
                'GLU', 'GL-', 'ASP', 'AS-', 'LYS', 'LY+', 'ARG', 'AR+', 'HID', 'HIE', 'HIP', 'SER', 'PRO', 'ASN', 'GLN']

    nt_residues = ['NGLY', 'NALA', 'NVAL', 'NLEU', 'NILE', 'NPHE', 'NTYR', 'NTRP', 'NMET', 'NCYS', 'NCYX', 'NTHR',
                'NGLU', 'NGL-', 'NASP', 'NAS-', 'NLYS', 'NLY+', 'NARG', 'NAR+', 'NHID', 'NHIE', 'NHIP', 'NSER', 'NPRO', 'NASN', 'NGLN']


    pdbfile = open(pdb, 'r').readlines()
    newfile = open(pdb, 'w')
    residue = []

    resname = 'NTER'
    for line in pdbfile:
        current_resnr = line[21:26]
        if current_resnr == res_nr:
            resname = line[17:21].strip()

    delete_H = False
    if resname in residues:
        old = residues
        new = nt_residues
    else:
        old = nt_residues
        new = residues
        delete_H = True

    newresname = 'NTER'
    for i in range(len(old)):
        if old[i] == resname:
            newresname = new[i]
    if newresname == 'NTER':
        print('Could not toggle N-terminal state ...')
        return

    for line in pdbfile:
        current_resnr = line[21:26]
        if current_resnr == res_nr:
            if delete_H:
                if 'H' not in line[13:14]:
                    newfile.write('%s%4s%s' % (line[0:17], newresname.ljust(4), line[21:]))
            if not delete_H:
                newfile.write('%s%4s%s' % (line[0:17], newresname.ljust(4), line[21:]))
        else:
            newfile.write(line)

    newfile.close()
    return newresname


def getLib(libs = ['Qoplsaa.lib']):
    """
    Reads the currently specified Library files
    and returns list of entries
    """
    lib_entries = []

    for libfiles in libs:
        try:
            libfile = open(libfiles, 'r').readlines()
            for line in libfile:
                if re.search('{.*}', line):
                    lib_entries.append(re.search('{.*}', line).group(0).strip('{*}'))
        except:
            print(('Could not find %s' % libfiles))
            continue

    return lib_entries


def checkPdbLib(res_list, resnr_list, lib_list):
    """
    Checks if residues in res_list exist in lib_list. If they do not,
    residue and residue number is returned.
    """
    missingres = []
    missingresnr = []

    for i in range(len(res_list)):
        if res_list[i] not in lib_list:
            if res_list[i] not in missingres:
                missingres.append(res_list[i])
                missingresnr.append(resnr_list[i])

    return missingres, missingresnr


def toggleAllCharges(pdbFile, charge_dict, atomtype_dict, toggle='off', simrad = 99.9, xc=0, yc=0, zc=0):
    try:
        file = open(pdbFile, 'r').readlines()
    except:
        sys.exit('Could not open %s' % pdbFile)

    simrad = float(simrad) * 0.85
    xc = float(xc)
    yc = float(yc)
    zc = float(zc)

    atomtype = 'CA'
    currentResnr = 0
    j_start = 0
    tmpfile = []

    for i in range(len(file)):
        currentRes = file[i][17:21].strip()
        if file[i][21:26] != currentResnr:
            j_start = i
        currentResnr = file[i][21:26]
        if currentRes in list(charge_dict.keys()):
            if file[i][12:14].strip() != 'H':
                if toggle == 'on':
                    atomtype = atomtype_dict[currentRes]
                    j = j_start
                    while file[j][21:26] == currentResnr:
                        if file[j][12:16].strip() == atomtype:
                            x = float(file[j][26:].split()[0])
                            y = float(file[j][26:].split()[1])
                            z = float(file[j][26:].split()[2])
                            radius = np.sqrt((xc - x)**2 + (yc - y)**2 + (zc - z)**2)
                            if radius <= simrad:
                                newRes = charge_dict[currentRes]
                                tmpfile.append('%s%4s%s' % (file[i][0:17], newRes.ljust(4), file[i][21:]))

                            else:
                                tmpfile.append(file[i])
                        j += 1
                else:
                    newRes = charge_dict[currentRes]
                    tmpfile.append('%s%4s%s' % (file[i][0:17], newRes.ljust(4), file[i][21:]))

            else:
                continue
        else:
            tmpfile.append(file[i])

    pdb_new = open(pdbFile, 'w')

    for line in tmpfile:
        pdb_new.write(line)

    pdb_new.close()


def findSS(pdbfile):
    """
   Functions to find potential S-S bridges. Returns:
   1. Arrays with atomnumber pairs to be used in Qprep5: makebond atom_i atom_j
   2. Array with residue number pairs.
   3. List with distances between S-S below treshold.
   """

    try:
        file = open(pdbfile, 'r').readlines()
    except:
        sys.exit('Could not find %s' % pdbfile)

    #Distance treshold for S-S bridge:
    ssTreshold = 2.20
    #List with atom numbers for SG atoms
    #Used for calculting the distances
    sgList = []
    #Array with x,y,z coordinates for sgList:
    sgCoordinates = []
    #list with CYS residue number
    cysRes = []
    #List with residues where CYS --> CYX
    cyxRes = []
    #Atom number for Qprep make bond command:
    cyxAtoms = []
    #List with distances of S-S below treshold
    cyxDist = []

    #Find all CYS residues and store coordinates:
    for line in file:
        try:
            if 'ATOM' in line.split()[0] or 'HETATM' in line.split()[0]:
                if line[13:17].strip() == 'SG':
                    tmp = []
                    sgList.append(line[6:11])
                    cysRes.append(line[22:26])
                    tmp.append(float(line[27:].split()[0]))
                    tmp.append(float(line[27:].split()[1]))
                    tmp.append(float(line[27:].split()[2]))
                    sgCoordinates.append(tmp)
        except:
            continue

    #Calculate distance between all CYS residues SG atoms
    #Store all below treshold distance:
    for atom1 in range(0, len(sgList) - 1):
        x1, y1, z1 = sgCoordinates[atom1][0:]
        for atom2 in range(atom1 + 1, len(sgList)):
            x2, y2, z2 = sgCoordinates[atom2][0:]
            distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
            if distance < ssTreshold:
                tmp = []
                tmp.append(sgList[atom1].strip())
                tmp.append(sgList[atom2].strip())
                cyxAtoms.append(tmp)
                cyxRes.append(cysRes[atom1])
                cyxRes.append(cysRes[atom2])
                cyxDist.append(distance)
    return cyxAtoms, cyxRes, cyxDist


def getResnrs(pdb, residue='RES'):
    """
    Function that returns list of residue numbers
    for selected residue name
    """
    try:
        pdbFile = open(pdb, 'r').readlines()
    except:
        sys.exit('Could not open %s' % pdb)

    resnrs = []
    current_resnr = 0
    for line in pdbFile:
        try:
            if 'ATOM' in line.split()[0] or 'HETATM' in line.split()[0]:
                if residue == line[17:20]:
                    if current_resnr == 0:
                        current_resnr = line[22:26]
                        resnrs.append(current_resnr)
                    elif current_resnr != line[22:26]:
                        current_resnr = line[22:26]
                        resnrs.append(current_resnr)
        except:
            continue
    return resnrs

def calcDist(pdb,resnr=1, atomtype='C', xc=0.0, yc=0.0, zc=0.0):
    """
    Calculates distance between a residue atom and a given x,y,z.
    Returns distance
    """
    pdbfile = open(pdb,'r').readlines()
    xc = float(xc)
    yc = float(yc)
    zc = float(zc)
    radius = 0.0

    for line in pdbfile:
        if 'ATOM' in line or 'HETATM' in line:
            if int(line[21:26]) == int(resnr):
                if line[12:16].strip() == atomtype:
                    x = float(line[28:].split()[0])
                    y = float(line[28:].split()[1])
                    z = float(line[28:].split()[2])
                    radius = np.sqrt(((xc - x)**2) + ((yc - y)**2) + ((zc - z)**2))

    return radius


def convertOldNew(pdb='structure.pdb', resnrs=[], newres='NEW'):
    """
   Function to change residue name in list of residues (residue numbers).
    :param pdb:
    :param resnrs:
    :param newres:
   """
    exclude_H = ['GLU', 'GL-', 'ASP', 'AS-', 'HID','HIP','HIE','LYS','LY+','ARG','AR+',
                 'NGLU', 'NGL-', 'NASP', 'NAS-', 'NHID','NHIP','NHIE','NLYS','NLY+','NARG','NAR+',
                 'CGLU', 'CGL-', 'CASP', 'CAS-', 'CHID','CHIP','CHIE','CLYS','CLY+','CARG','CAR+']
    try:
        pdbFile = open(pdb, 'r').readlines()
    except:
        sys.exit('Could not open %s' % pdb)

    newres += ' '
    newPdb = open(pdb, 'w')
    if len(resnrs) > 0:
        for line in pdbFile:
            if 'GAP' in line:
                newPdb.write(line)
                continue
            try:
                resnr = line[22:26]
                oldres = line[17:21].strip()
                if len(oldres) == 4:
                    newres_insert = '%1s%3s' % (oldres[0], newres.strip())
                else:
                    newres_insert = newres
                if resnr in resnrs:
                    if oldres in exclude_H:
                        if 'H' not in line[13:14]:
                            newPdb.write('%s%s%s' % (line[0:17], newres_insert, line[21:]))
                    else:
                        newPdb.write('%s%s%s' % (line[0:17], newres_insert, line[21:]))
                else:
                    newPdb.write(line)
            except:
                continue
    else:
        print('No residues speciefied for conversion')


def centerofmass(pdb, use_mass=False):
    """
   Returns center of mass x,y,z coordinates.
   use_mass controls if actual masses are to
   be used or not. If False, all atoms will be treated
   as having the same mass.
   """
    try:
        pdbfile = open(pdb, 'r').readlines()
    except:
        print(('Could not open %s' % pdb.split('/')[-1]))
    masses = []
    x, y, z = [], [], []

    for line in pdbfile:
        try:
            if 'ATOM' == line.split()[0]:
                if not use_mass:
                    masses.append(1.00)
                if use_mass:
                    atomtype = line[13]
                    if atomtype == 'H':
                        mass = 1.008
                    elif atomtype == 'C':
                        mass = 12.011
                    elif atomtype == 'N':
                        mass = 14.007
                    elif atomtype == 'O':
                        mass = 15.999
                    elif atomtype == 'S':
                        mass = 32.06
                    else:
                        mass = 12.011
                    masses.append(mass)
                #x.append(float(line[28:].split()[0]))
                x.append(float(line[30:38]))
                #y.append(float(line[28:].split()[1]))
                y.append(float(line[38:46]))
                #z.append(float(line[28:].split()[2]))
                z.append(float(line[46:54]))
        except:
            continue
    tot_mass = np.sum(masses)
    x_weights = []
    y_weights = []
    z_weights = []

    for atom in range(len(masses)):
        x_weights.append((masses[atom] * x[atom]) / tot_mass)
        y_weights.append((masses[atom] * y[atom]) / tot_mass)
        z_weights.append((masses[atom] * z[atom]) / tot_mass)

    xc = '%.3f' % np.sum(x_weights)
    yc = '%.3f' % np.sum(y_weights)
    zc = '%.3f' % np.sum(z_weights)

    return xc, yc, zc


def findRadius(pdb):
    x, y, z = centerofmass(pdb)
    x = float(x)
    y = float(y)
    z = float(z)
    pdbfile = open(pdb, 'r').readlines()
    max_r = 0.00

    for line in pdbfile:
        try:
            if 'ATOM' == line.split()[0]:
                xi = float(line[28:].split()[0])
                yi = float(line[28:].split()[1])
                zi = float(line[28:].split()[2])
                r = np.sqrt((xi - x) ** 2 + (yi - y) ** 2 + (zi - z) ** 2)
                if r > max_r:
                    max_r = r
        except:
            continue
    return max_r


def pdb_info(pdb):
    """
   Reads Q PDB file and returns 7 lists:
   atom number, atom types, residues, residue numbers,
   x, y, and z coordinates.
   """
    try:
        pdbfile = open(pdb, 'r').readlines()
    except:
        print(('Could not open pdb file %s' % pdb))

    atomnrs = []
    atomtypes = []
    residues = []
    residuenrs = []
    x = []
    y = []
    z = []

    for line in pdbfile:
        try:
            if 'ATOM' in line.split()[0] or 'HETATM' in line.split()[0]:
                atomnrs.append(int(line.split()[1]))
                atomtypes.append(line[12:17].strip())
                residues.append(line[17:21].strip())
                residuenrs.append(int(line[22:26]))
                x = float(line[28:].split()[0])
                y = float(line[28:].split()[1])
                z = float(line[28:].split()[2])
        except:
            continue

    return atomnrs, atomtypes, residues, residuenrs, x, y, z
