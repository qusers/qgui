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
import sys
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
                        libfiles.append(lib)

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




