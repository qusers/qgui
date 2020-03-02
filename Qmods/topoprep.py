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

from tkinter import  Radiobutton, Spinbox, StringVar, Entry, Text, Label, Frame, Button, Scrollbar, Toplevel, \
    Checkbutton, Listbox, DISABLED, NORMAL, END, GROOVE, LEFT, RIGHT, IntVar, PhotoImage

# -*- coding: utf-8 -*-

import os
import shutil
import tkinter.font
import pickle
import prepareTopology as pt
from select_xyz import AtomSelect
from edit_file import FileEdit
import time
import numpy as np
from subprocess import Popen, PIPE
import copy
import signal
import sys
import linecache
import random as rng


class TopologyPrepare(Toplevel):
    """Implements a dialog-box when Prepare -> Topology is chosen from menubar.
    Has got methods to modify pdf file/structure such that it can be used by Q.
    TopologyPrepare has got reference to the MainWindow-class."""

    def __init__(self, app, root, pdbfile, prm, lib, qgui_parent=True):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.pdbfile = pdbfile
        self.check_variable = IntVar()
        self.neutralize = IntVar()
        self.xvar = StringVar()
        self.yvar = StringVar()
        self.zvar = StringVar()
        self.qgui_parent = qgui_parent

        #Pymol session:
        self.sync_pymol = IntVar()
        self.session = None

        #Size of simulation sphere
        self.sphere_radius = StringVar()

        #Pymol zoom buffer when residue is clicked:
        self.pymol_zoom = 10

        if not self.qgui_parent:
            self.app.log = self.app.app.log
            self.app.errorBox = self.app.app.errorBox
            self.app.workdir = self.app.app.workdir

        if not pt.checkPdb(self.pdbfile):
            pdb_id_q = pt.convertPdb(self.pdbfile, self.app.workdir, 0, 0)
            self.app.log('info','%s converted to Q format pdb (%s).' % (self.pdbfile.split('/')[-1], pdb_id_q.split('/')[-1]))
            self.pdbfile = self.app.workdir + '/' + pdb_id_q

        self.app.log('info', 'Topology Prepare session started.')
        self.dialog_box()
        self.set_structure_entry()


        #Auto insert name for topology and topology pdb
        self.topology_entry.insert(0.0,self.pdbfile.split('/')[-1].split('.')[0]+'.top')
        if '_top' not in self.pdbfile.split('/')[-1]:
            self.topo_pdb_entry.insert(0.0,self.pdbfile.split('/')[-1].split('.')[0]+'_top.pdb')
        else:
            self.topo_pdb_entry.insert(0.0, self.pdbfile.split('/')[-1])

        #Generate SS?
        self.makeSS = False

        #XYZ sim center:
        self.set_sim_center()

        #Trace coordinate changes to charge radiuses:
        self.xvar.trace('w', self.simcenter_changed)
        self.yvar.trace('w', self.simcenter_changed)
        self.zvar.trace('w', self.simcenter_changed)

        #Sphere radius
        system_radius = pt.findRadius(self.pdbfile)
        self.sphere_entry.delete(0, END)
        self.sphere_entry.insert(0, '%d' % int(round(system_radius + 5)))

        #If pymol is active, adjust sphere size interactively
        self.sphere_radius.trace('w', self.adjust_simulation_sphere)


        #Set default solvation
        self.check_solvation.set(1)

        #TODO change the HIS module!
        #Check if HIS in file and convert to HID:
        self.his = pt.getResnrs(self.pdbfile, 'HIS')

        if len(self.his) > 0:
            self.app.log('info', 'Found residue HIS in file. Autoconverting to HID/HIE\n')
            #Convert to HID as a first guess. This will be corrected if protons exist later.
            for entry in self.his:
                pt.convertOldNew(self.pdbfile, entry, 'HID')

        #ParameterList (get this from settings)
        self.prm = prm

        #Library list (get this from settings)
        self.lib = lib

        #These are filled from self.checkLib
        #{RES: charge}
        self.res_charge = dict()
        #{RES:[atoms]}
        self.res_atoms = dict()
        #{RES: color}
        self.res_colors = dict()
        #{RES: {nr: {x,y,z, dist} } }
        self.resname_nr_dist = dict()
        #{RES: ATOM}
        self.resname_distatom = dict()

        #Default toggle_res. Will be overwritten if [toggle_residues] in lib file:
        self.toggle_res = {'ARG': 'ARN',
                           'ARN': 'ARG',
                           'LYS': 'LYN',
                           'LYN': 'LYS',
                           'ASP': 'ASH',
                           'ASH': 'ASP',
                           'GLU': 'GLH',
                           'GLH': 'GLU',
                           'HIP': 'HID',
                           'HID': 'HIE',
                           'HIE': 'HIP',
                           }

        #Dictionary describing atoms that changes for new residue in toggle residue:
        self.toggle_res_atoms = {'ARG': '+H12 NH1',
                                 'ARN': '-H12',
                                 'LYS': '+HZ3 NZ',
                                 'LYN': '-HZ3',
                                 'ASP': '-HD2',
                                 'ASH': '+HD2 OD2',
                                 'GLU': '-HE2',
                                 'GLH': '+HE2 OE2',
                                 'HIP': '+HD1 ND1',
                                 'HID': '-HE2',
                                 'HIE': '-HD2 +HE2 NE2'}

        #Check if residues in pdb file exist in library and generate:
        self.checkLib(False)

        self.set_total_charge()

        self.updateList()

    def set_sim_center(self):
        """
        finds the center of current loaded pdb file and
        write xyz into labels.
        """
        self.xc, self.yc, self.zc = pt.centerofmass(self.pdbfile)
        self.center_x_entry.delete(0, END)
        self.center_x_entry.insert(0, self.xc)
        self.center_y_entry.delete(0, END)
        self.center_y_entry.insert(0, self.yc)
        self.center_z_entry.delete(0, END)
        self.center_z_entry.insert(0, self.zc)

    def simcenter_changed(self, *args):
        """
        Updates list with chargable residues and radii whenever
        simulation center is changed
        """
        try:
            self.updateList()
            if self.session:
                self.update_simulation_sphere()
        except:
            return

    def checkLib(self, write_pdb=True):
        """
        #Check if pdb residues exist in current library:
        """
        if write_pdb:
            self.write_pdb()

        status = 'Lib status: OK'

        #Get what atoms to modify upon terminal toggling
        toggle_term_atom = {'n': '-HT3',
                            'N': '+HT3 N',
                            'c': '+HO2 OT2',
                            'C': '-HO2'}

        #Terminals start with:
        term_start = ['N','n','C','c']

        #collect all residues potential terminal residues and fix toggle later:
        term_res = list()

        #Collect temp CRES, cRES, nRES and NRES and delete them when done
        del_list = list()

        #Collect info from lib files
        for lib in self.lib:
            with open(lib, 'r') as lib:
                found_res = False
                found_atoms = False
                found_toggle_res = False

                charges = list()
                res = 'Hello'

                for line in lib:
                    if found_toggle_res:
                        if ' QGUI_TOGGLE ' in line or ' QGUI_N-TERM ' in line or ' QGUI_C-TERM ' in line:
                            try:
                                origial_res = line.split()[2]
                                new_res = line.split()[3]
                                atoms = ' '.join(line.split('radius_to')[0].split()[4:])
                                r_atom = line.split('radius_to')[1].split('#')[0].strip('\n\r ')

                                self.toggle_res[origial_res] = new_res
                                self.toggle_res_atoms[new_res] = atoms
                                if new_res not in self.resname_distatom:
                                    self.resname_distatom[new_res] = r_atom
                                    self.resname_nr_dist[new_res] = dict()
                            except:
                                print('Unexpected line in [toggle_residues]')
                                print(line)

                        if ' QGUI_N-TERM ' in line or ' QGUI_C-TERM ' in line:
                            toggle_term_atom[line.split()[3][0]] = ' '.join(line.split()[4:]).strip('\n')
                            origial_res = line.split()[2]
                            new_res = line.split()[3]
                            for to_remove in [origial_res, new_res]:
                                if to_remove not in del_list:
                                    del_list.append(to_remove)

                    if found_res:
                        if '[bonds]' in line or '*-------' in line or '[charge_groups]' in line:
                            found_toggle_res = False
                            tot_charge = float(np.sum(charges))
                            self.res_charge[res] = round(tot_charge, 1)
                            print(('%4s: %.1f' % (res, round(tot_charge, 1))))
                            if tot_charge > 0.499:
                                self.res_colors[res] = 'blue'
                            elif tot_charge < -0.499:
                                self.res_colors[res] = 'red'
                            elif res == 'CYX':
                                self.res_colors[res] = 'yellow'
                            else:
                                self.res_colors[res] = 'white'

                            found_atoms = False
                            found_res = False

                        if found_atoms:

                            if len(line.split()) > 3:
                                if line.split()[0] != '!':
                                    self.res_atoms[res].append(line.split()[1])
                                    charges.append(float(line.split()[3]))
                        if '[atoms]' in line:
                            found_atoms = True

                    if '{' in line and '}' in line:
                        found_toggle_res = False
                        res = line.split('{')[1].split('}')[0]

                        #Check if residue can be a terminal residue:
                        if res[0] in term_start and len(res) == 4:
                            if res not in term_res:
                                term_res.append(res)

                        self.res_atoms[res] = list()
                        self.res_charge[res] = 0.0
                        charges = list()
                        found_res = True
                    if '[toggle_residues]' in line:
                        print('Found [toggle_residues] in library file. Original will be overwritten.')
                        found_toggle_res = True
                        self.toggle_res = dict()

        #Add N- and C-terminals to toggle_residues
        for term in term_res:
            res = term[1:]
            if res in list(self.res_charge.keys()):
                if term[0].islower():
                    toggle_to = term[0].capitalize() + term[1:]
                else:
                    toggle_to = term[0].lower() + term[1:]

                if toggle_to in list(self.res_charge.keys()):
                    self.toggle_res[term] = toggle_to
                    self.toggle_res_atoms[toggle_to] = toggle_term_atom[toggle_to[0]]
                    self.resname_distatom[term] = copy.deepcopy(self.resname_distatom[term[0]+'RES'])
                    self.resname_distatom[toggle_to] = copy.deepcopy(self.resname_distatom[toggle_to[0]+'RES'])
                else:
                    print(('(Can not toggle %4s to %4s: %4s not found in lib.)' % (term, toggle_to, toggle_to)))

        #Remove temp NRES and CRES definitions for terminals
        for to_remove in del_list:
            del self.toggle_res[to_remove]
            del self.resname_distatom[to_remove]
            del self.resname_nr_dist[to_remove]

        #Get info from pdb file and compare to existing lib entries
        pdb_atoms = list()

        find_CA = False
        c_term_found = False
        create_cterm = False
        create_nterm = False
        res = None

        #count line number
        line_number = 0

        #Comput distance from previous C to current N to check for terminals
        c_xyz = False
        n_xyz = False

        with open(self.pdbfile, 'r') as pdb:
            res_nr = 0
            for line in pdb:
                line_number += 1
                #Check for GAP or TER statements (C-terminals)
                if line.startswith('TER') or 'GAP' in line:
                    c_term_found = True
                    create_nterm = True

                if 'ATOM' in line or 'HETATM' in line:
                    atom_name = line[13:17].strip()

                    #Check if residue is N-terminal:
                    if atom_name == 'N':
                        n_xyz = list()
                        #n_xyz = map(float, line[30:].split()[0:3])
                        n_xyz.append(float(line[30:38]))
                        n_xyz.append(float(line[38:46]))
                        n_xyz.append(float(line[46:54]))
                        if not c_xyz:
                            create_nterm = True

                        elif c_xyz and n_xyz:
                            r = np.sqrt((n_xyz[0] - c_xyz[0])**2 + (n_xyz[1] - c_xyz[1])**2 + (n_xyz[2] - c_xyz[2])**2)
                            if r > 2.2:
                                c_term_found = True
                                create_nterm = True

                    #Collect coordinates for carbonyl carbon
                    elif atom_name == 'C':
                        c_xyz = list()
                        #c_xyz = map(float, line[30:].split()[0:3])
                        c_xyz.append(float(line[30:38]))
                        c_xyz.append(float(line[38:46]))
                        c_xyz.append(float(line[46:54]))

                    #collect atom names for given residue:
                    if res_nr == int(line[21:26]):
                        pdb_atoms.append(atom_name)
                    elif res in list(self.res_atoms.keys()):
                        #Check if the next residue is a aa or not. If not, current res could be C-term:
                        next_res = line[17:21].strip()
                        if next_res not in list(self.res_atoms.keys()):
                            if res in list(self.res_atoms.keys()):
                                c_term_found = True
                        elif 'N' not in self.res_atoms[next_res] and 'CA' not in self.res_atoms[next_res]:
                            if 'C' in self.res_atoms[res] and 'CA' in self.res_atoms[res]:
                                c_term_found = True


                    #New residue number starts:
                    if res_nr != int(line[21:26]):
                        #Check if any heavy atoms in residue are missing in pdb file:
                        self.check_missing_atoms(res, res_nr, pdb_atoms)

                        #Check if previous residue was potential C-terminal
                        if c_term_found:
                            create_cterm = True

                            #Check if C-terminal name is valid or if it is possible to generate one
                            if len(res) > 3 and res[0].lower() == 'c':
                                if res in list(self.res_atoms.keys()):
                                    print(('%4s %3d is assumed a valid C-terminal' % (res, res_nr)))
                                    #No need to generate it, it is already defined
                                    create_cterm = False
                                    #Add the togglable terminal:
                                    if res not in list(self.resname_nr_dist.keys()):
                                        self.resname_nr_dist[res] = dict()
                                    if res_nr not in list(self.resname_nr_dist[res].keys()):
                                        self.resname_nr_dist[res][res_nr] = dict()
                                    if res not in list(self.resname_distatom.keys()):
                                        self.resname_distatom[res] = 'C'
                                    if c_xyz:
                                        self.resname_nr_dist[res][res_nr]['x'] = c_xyz[0]
                                        self.resname_nr_dist[res][res_nr]['y'] = c_xyz[1]
                                        self.resname_nr_dist[res][res_nr]['z'] = c_xyz[2]

                                else:
                                    print(('%4s %3d is not found in lib!' % (res, res_nr)))
                                    self.app.log(' ','WARNING: %4s %3d is not found in lib!\n' % (res, res_nr) )
                                    create_cterm = False

                        if create_cterm:
                            cres = self.check_c_term(res, res_nr, line_number)
                            self.check_missing_atoms(cres, res_nr, pdb_atoms)
                            c_term_found = False
                            create_cterm = False
                            c_xyz = False
                            if not create_nterm:
                                n_xyz = False

                        del pdb_atoms[:]
                        pdb_atoms.append(atom_name)
                        res = line[17:21].strip()
                        res_nr = int(line[21:26])
                        find_CA = False

                        if create_nterm:
                            if res in list(self.res_atoms.keys()):
                                if len(res) > 3 and res[0].lower() == 'n':
                                    print(('%4s %3d is assumed a valid N-terminal' % (res, res_nr)))
                                    #Add the togglable terminal if possible:
                                    if res not in list(self.resname_nr_dist.keys()):
                                        self.resname_nr_dist[res] = dict()
                                    if res_nr not in list(self.resname_nr_dist[res].keys()):
                                        self.resname_nr_dist[res][res_nr] = dict()
                                        if res not in list(self.resname_distatom.keys()):
                                            self.resname_distatom[res] = 'N'
                                        if n_xyz:
                                            self.resname_nr_dist[res][res_nr]['x'] = n_xyz[0]
                                            self.resname_nr_dist[res][res_nr]['y'] = n_xyz[1]
                                            self.resname_nr_dist[res][res_nr]['z'] = n_xyz[2]

                                elif 'N' in self.res_atoms[res]:
                                    if n_xyz:
                                        res = self.check_n_term(res, res_nr, n_xyz,
                                                                 line_number)
                                        #self.check_missing_atoms(nres, res_nr, pdb_atoms)
                            create_nterm = False
                            n_xyz = False

                        if 'CYS' in list(self.toggle_res.keys()):
                            if res == self.toggle_res['CYS']:
                                if not self.makeSS:
                                    self.check_variable.set(1)
                                    self.makeSS = True

                        if not res in list(self.res_atoms.keys()):
                            self.app.log(' ','WARNING: %4s %5d not found in lib entries!\n' % (res, res_nr))
                            if status != 'Lib entries missing':
                                status = 'Lib entries missing'

                        #Check if residue has atom defined for distance to sim center defined
                        #Apply default first atom as default. Change to CA if it exist!
                        if not res in list(self.resname_nr_dist.keys()):
                            if res in list(self.res_charge.keys()):
                                if abs(self.res_charge[res]) > 0.4999:
                                    print(('Added charged residue which is not togglable: %s  ' % res))
                                    self.resname_distatom[res] = atom_name
                                    self.resname_nr_dist[res] = dict()
                                    find_CA = True

                        if res in list(self.resname_nr_dist.keys()):
                            self.resname_nr_dist[res][res_nr] = dict()

                    #Check if atomname matches defined residue in lib file:
                    if res in list(self.res_atoms.keys()):
                        if not atom_name in self.res_atoms[res]:
                            missing_atom = True
                            #Check if residue can be toggled and toggle to correct libname
                            if res in list(self.toggle_res.keys()):
                                if atom_name in self.res_atoms[self.toggle_res[res]]:
                                    toggle_res = self.toggle_res[res]
                                    self.app.log('', 'Changed residue %s to %s (atom %s present)\n' %
                                                     (res, toggle_res, atom_name ))
                                    missing_atom = False
                                    #TODO write code for actual change - Thinking about best way to handle...
                                    res_new = self.toggle_res[res]
                                    if res_new not in list(self.resname_nr_dist.keys()):
                                        self.resname_nr_dist[res_new] = dict()

                                    self.resname_nr_dist[res_new][res_nr] = \
                                        copy.deepcopy(self.resname_nr_dist[res][res_nr])

                                    del self.resname_nr_dist[res][res_nr]
                                    #some dictionary with {atom_nr: True/False} ?
                                    res = res_new

                            if missing_atom:
                                self.app.log(' ','WARNING: Atomname %4s not found in lib for %4s %5d!\n' %
                                                 (atom_name, res, res_nr))
                                if status != 'Lib entries missing':
                                    status = 'Lib entries missing'

                    #Collect x,y,z for residues that can be toggled (compute distances to these):
                    if res in list(self.resname_nr_dist.keys()):
                        #Take coordinates of 1st atom in case distance atom is not found:
                        if 'x' not in list(self.resname_nr_dist[res][res_nr].keys()):
                            #x, y, z = map(float, line[30:].split()[0:3])
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            self.resname_nr_dist[res][res_nr]['x'] = x
                            self.resname_nr_dist[res][res_nr]['y'] = y
                            self.resname_nr_dist[res][res_nr]['z'] = z
                        #Find x,y,z for distatom
                        if find_CA:
                            distatom = 'CA'
                        else:
                            if res not in list(self.resname_distatom.keys()):
                                distatom = 'C'
                            else:
                                distatom = self.resname_distatom[res]

                        #If defined atom to compute distance to is found, add it:
                        if distatom == atom_name:

                            #x, y, z = map(float, line[30:].split()[0:3])
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            self.resname_nr_dist[res][res_nr]['x'] = x
                            self.resname_nr_dist[res][res_nr]['y'] = y
                            self.resname_nr_dist[res][res_nr]['z'] = z

            else:
                if res in list(self.toggle_res.keys()):
                    if res.lower() == self.toggle_res[res].lower():
                        print(('Valid C-terminal residue name found: %s %3d' % (res, res_nr)))
                    else:
                        self.check_c_term(res, res_nr, line_number)

                print('END OF FILE REACHED')

        #If CYX in pdb file, turn on autogenerate S-S bridges:
        if self.makeSS:
            self.app.log(' ', '\nFound residue %s in file. Generating S-S bonds:\n'
                                                  % self.toggle_res['CYS'])
            self.find_ss_bonds()

        self.lib_status.config(state=NORMAL)
        self.lib_status.delete(0.0, END)
        self.lib_status.insert(0.0, status.rjust(19))
        self.lib_status.config(state=DISABLED)

    def check_missing_atoms(self, res, res_nr, pdb_atoms):
        """
        Takes a list of atoms (pdb_atoms) found in pdb file for a residue (res) with residue number (res_nr)
        and checks if there are any atoms defined in the lib entry missing in the pdb file for that residue.
        """
        if res in list(self.res_atoms.keys()):
            for libatom in self.res_atoms[res]:
                if libatom not in pdb_atoms:
                    if libatom[0] != 'H':
                        self.app.log(' ','WARNING: Heavy atom %s missing in %4s %5d \n' % (libatom, res, res_nr))

    def check_n_term(self, res, res_nr, n_xyz, line_number):

        nterm = 'N'+res

        forward = False
        if not nterm in list(self.toggle_res.keys()):
            print(('Generating untogglable standard N-terminal from %s %3d --> %4s\n' % (res, res_nr, nterm)))

        if not nterm in list(self.res_charge.keys()):
            self.app.log(' ', 'WARNING: Could not generate N-terminal from %s %3d\n'  % (res, res_nr))
            self.app.log(' ', 'Residue %s %3d not defined!\n' % (nterm, res_nr))

        else:
            forward = True
            if nterm not in list(self.resname_nr_dist.keys()):
                self.resname_nr_dist[nterm] = dict()

            #Check if information exist about distance atom from original residue:
            if res in list(self.resname_nr_dist.keys()):
                if res_nr in list(self.resname_nr_dist[res].keys()):
                    self.resname_nr_dist[nterm][res_nr] = copy.deepcopy(self.resname_nr_dist[res][res_nr])
                    print(('Created %s %d' % (nterm, res_nr)))
                    del self.resname_nr_dist[res][res_nr]
                    print(('Deleted %s %d' % (res, res_nr)))
                    forward = False
                    res = nterm

        if forward:
            res = nterm
            self.resname_nr_dist[nterm][res_nr] = dict()
            #Insert N as distant atom (default).
            self.resname_nr_dist[nterm][res_nr]['x'] = n_xyz[0]
            self.resname_nr_dist[nterm][res_nr]['y'] = n_xyz[1]
            self.resname_nr_dist[nterm][res_nr]['z'] = n_xyz[2]

            #read pdb file and look for toggle atom for sim sphere distance radius:
            distatom2 = 'N'

            if nterm in list(self.resname_distatom.keys()):
                distatom2 = self.resname_distatom[nterm]

            for i in range(0, 30):
                pdb_line = linecache.getline(self.pdbfile, (line_number + i))

                if 'ATOM' in pdb_line or 'HETATM' in pdb_line:
                    print(pdb_line)
                    if res_nr != int(pdb_line[21:26]):
                        break

                    atom_name2 = pdb_line[13:17].strip()

                    if atom_name2 == distatom2:
                        #x, y, z = map(float, pdb_line[30:].split()[0:3])
                        x = float(pdb_line[30:38])
                        y = float(pdb_line[38:46])
                        z = float(pdb_line[46:54])
                        self.resname_nr_dist[nterm][res_nr]['x'] = x
                        self.resname_nr_dist[nterm][res_nr]['y'] = y
                        self.resname_nr_dist[nterm][res_nr]['z'] = z

                        print(('Found correct distance atom for N-term %s' % nterm))
                        break

        return res

    def check_c_term(self, res, res_nr, line_number):
        cterm = 'C'+res

        rewind = False
        if not cterm in list(self.res_charge.keys()):
            self.app.log(' ', 'WARNING: Could not generate C-terminal from %s %3d\n'  % (res, res_nr))
            self.app.log(' ', 'Residue %s %3d not defined!\n' % (cterm, res_nr))

        #elif not cterm in self.toggle_res.keys():
        #    print 'Could not generate C-terminal from %s %3d\n' % (res, res_nr)

        else:
            rewind = True
            if cterm not in list(self.resname_nr_dist.keys()):
                self.resname_nr_dist[cterm] = dict()

            if cterm not in list(self.resname_distatom.keys()):
                self.resname_distatom[cterm] = 'C'

            if res in list(self.resname_nr_dist.keys()):
                if res_nr in list(self.resname_nr_dist[res].keys()):
                    self.resname_nr_dist[cterm][res_nr] = copy.deepcopy(self.resname_nr_dist[res][res_nr])
                    del self.resname_nr_dist[res][res_nr]

                    rewind = False
                    res = cterm

        if rewind:
            res = cterm
            self.resname_nr_dist[cterm][res_nr] = dict()
            #Rewind pdb file and look for toggle atom for sim sphere distance radius:
            distatom2 = self.resname_distatom[cterm]
            for i in range(1, 30):
                pdb_line = linecache.getline(self.pdbfile, (line_number - i))

                if 'ATOM' in pdb_line or 'HETATM' in pdb_line:
                    if res_nr != int(pdb_line[21:26]):
                        break
                    atom_name2 = pdb_line[13:17].strip()

                    if atom_name2 == distatom2 or atom_name2 == 'C':
                        x, y, z = list(map(float, pdb_line[30:].split()[0:3]))
                        self.resname_nr_dist[cterm][res_nr]['x'] = x
                        self.resname_nr_dist[cterm][res_nr]['y'] = y
                        self.resname_nr_dist[cterm][res_nr]['z'] = z
                    if atom_name2 == distatom2:
                        print(('Found correct distance atom for C-term %s (%d)' % (cterm, res_nr)))
                        break


        return res

    def write_pdb(self):
        """
        Updates the loaded pdb file with the correct modifications made in the topology prepare tool.
        """
        old_pdb = open(self.pdbfile, 'r').readlines()
        print('Writing new pdb file from toggle')

        #Residue numbers to modify:
        #res nr : resname
        res_to_mod = dict()
        for res in list(self.resname_nr_dist.keys()):
            for nr in list(self.resname_nr_dist[res].keys()):
                res_to_mod[nr] = res

        #Write new pdb file to existing name:
        new_pdb = open(self.pdbfile, 'w')

        #Add/delete atoms using:
        #self.toggle_res_atoms[new_res]
        atomnr = 0
        for i in range(len(old_pdb)):
            try:
                res_nr = int(old_pdb[i][21:26])
                orig_res = old_pdb[i][17:21].strip()

                if res_nr in list(res_to_mod.keys()):
                    new_res = res_to_mod[res_nr]
                    if new_res != orig_res:
                        print(('Residue %3d %4s --> %4s' % (res_nr, orig_res, new_res)))
                        atomname = old_pdb[i][12:17].strip()

                        #Get atoms to modify, if some:
                        if new_res in list(self.toggle_res_atoms.keys()):
                            modatoms = self.toggle_res_atoms[new_res]
                        else:
                            modatoms = str()

                        #So far we only need to delete H-atoms. Qprep will add missing hydrogens.
                        del_atoms = list()
                        for atom in modatoms.split():
                            if '-' in atom:
                                del_atoms.append(atom.strip('-'))

                        if atomname in del_atoms:
                            print(('Atom %s was deleted from %s %d' % (atomname, orig_res, res_nr)))
                        else:
                            atomnr += 1
                            print(('%s%5d  %s%4s%s' %
                                            (old_pdb[i][0:6], atomnr, old_pdb[i][13:17], new_res.ljust(4), old_pdb[i][21:])))

                            new_pdb.write('%s%5d  %s%4s%s' %
                                            (old_pdb[i][0:6], atomnr, old_pdb[i][13:17], new_res.ljust(4), old_pdb[i][21:]))

                    else:
                        atomnr += 1
                        new_pdb.write('%s%5d  %s' % (old_pdb[i][0:6], atomnr, old_pdb[i][13:]))

                else:
                    atomnr += 1
                    new_pdb.write('%s%5d  %s' % (old_pdb[i][0:6], atomnr, old_pdb[i][13:]))

            except Exception as e:
                print(e)
                print('Oups! I just tried to magically convert nothing to something')
                print('...this was probably just an empty line, GAP or TER line. Nothing to worry about!\n')
                new_pdb.write(old_pdb[i])

        new_pdb.close()
        self.app.log('info','PDB file updated from topolgy prepare!')

    def writeTopology(self):
        """
        Writes the input files for Qprep. If there are several parameter files defined, these will be merged
        into one file!
        """
        self.write_pdb()

        qprepinp_name = self.pdbfile.split('/')[-1].split('.')[0]+'_Qprep.inp'
        qprepinp = open(self.app.workdir + '/' + qprepinp_name,'w')
        self.topname = self.topology_entry.get(0.0,END).strip()
        self.top_pdbname = self.topo_pdb_entry.get(0.0,END).strip()
        for entry in self.lib:
            qprepinp.write('readlib %s\n' % entry)

        #Can only use 1 prm file. If more than 1 is in settings, merge them!
        if len(self.prm) > 1:
            self.app.log(' ','\n****\nMore than 1 parameter file is set in Settings.\nThese files will'
                                     ' be merged. \n--> Force fields must be of same type.'
                                     '\n--> If atom types with the same name exist,'
                                     'they will be redefined based on which comes last.\n****\n')
            atom_index = []
            bond_index = []
            angle_index = []
            torsion_index = []
            improper_index = []

            merged_name = 'merged.prm'
            merged_file = open('%s/%s' % (self.app.workdir, merged_name), 'w')

            #TODO
            #Add all types to avoid dublicates
            types_added = list()

            #Go through files and find indices:
            prmfiles = []
            for entry in self.prm:
                print((entry,self.prm))
                prmfile = open(entry,'r').readlines()
                for i in range(len(prmfile)):
                    if '[atom_types]' in prmfile[i]:
                        atom_index.append(i)
                    elif '[bonds]' in prmfile[i]:
                        bond_index.append(i)
                    elif '[angles]' in prmfile[i]:
                        angle_index.append(i)
                    elif '[torsions]' in prmfile[i]:
                        torsion_index.append(i)
                    elif '[impropers]' in prmfile[i]:
                        improper_index.append(i)
                    else:
                        continue
                prmfiles.append(prmfile)

            #Merge atom types
            for i in range(0, bond_index[0]):
                merged_file.write(prmfiles[0][i])
            for entry in range(1, len(prmfiles)):
                for j in range(atom_index[entry] + 1, bond_index[entry]):
                    if prmfiles[entry][j].startswith('*'):
                        continue
                    else:
                        merged_file.write(prmfiles[entry][j])
            #Merge bonds
            try:
                for i in range(bond_index[0], angle_index[0]):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(bond_index[entry] + 1 ,angle_index[entry]):
                        if prmfiles[entry][j].startswith('*'):
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            #Merge angles
            try:
                for i in range(angle_index[0], torsion_index[0]):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(angle_index[entry] + 1 ,torsion_index[entry]):
                        if prmfiles[entry][j].startswith('*'):
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            #Merge torsions
            try:
                for i in range(torsion_index[0], improper_index[0]):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(torsion_index[entry] + 1 ,improper_index[entry]):
                        if prmfiles[entry][j].startswith('*'):
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            #Merge impropers
            try:
                for i in range(improper_index[0], len(prmfiles[0])):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(improper_index[entry] + 1 , len(prmfiles[entry])):
                        if prmfiles[entry][j].startswith('*'):
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            merged_file.close()
            qprepinp.write('readprm %s/%s\n' % (self.app.workdir, merged_name))
        else:
            qprepinp.write('readprm %s\n' % self.prm[0])

        qprepinp.write('readpdb %s\n' % self.pdbfile)

        if self.makeSS:
            for cys_i in sorted(self.cys_residues.keys()):
                qprepinp.write('addbond %s:SG %s:SG\n' % (cys_i, self.cys_residues[cys_i]))

        radius = self.sphere_entry.get()

        xc = self.center_x_entry.get()
        yc = self.center_y_entry.get()
        zc = self.center_z_entry.get()

        qprepinp.write('boundary 1 %s %s %s %s\n' % (xc,yc,zc,radius))

        if self.check_solvation.get() != 0:
            if self.check_solvation.get() == 1:
                solvation = 'HOH'
            else:
                solvation = 'SPC'
            qprepinp.write('solvate %s %s %s %s 1 %s\n' % (xc,yc, zc,radius,solvation))

        qprepinp.write('maketop %s topology\n' % self.pdbfile.split('.')[0])
        qprepinp.write('writetop %s' % self.topname)
        qprepinp.write('\n')
        qprepinp.write('writepdb %s' % self.top_pdbname)
        qprepinp.write('\n')
        qprepinp.write('y\n')
        qprepinp.write('q\n')
        qprepinp.close()

        self.app.log('info','Qprep input file written: %s' % qprepinp_name)

    def run_qprep(self):
        q_settings = pickle.load(open(self.app.settings_path + '/Qsettings','rb'))
        qprepinp_name = self.pdbfile.split('/')[-1].split('.')[0]+'_Qprep.inp'
        qprepout_name = self.pdbfile.split('/')[-1].split('.')[0]+'_Qprep.log'

        #Move to workdir
        current_dir = os.getcwd()
        os.chdir(self.app.workdir)

        self.writeTopology()

        Popen(args= "%s < %s/%s > %s/%s" % (q_settings[ 'executables' ][0],self.app.workdir, qprepinp_name, self.app.workdir, qprepout_name),shell=True).wait()

        if os.path.isfile(self.app.workdir + '/' + qprepout_name):
            with open(self.app.workdir + '/' + qprepout_name, 'r') as qpreplog:
                for line in qpreplog:
                    print(line)
                    self.app.log(' ', line)
                    if 'Topology successfully generated' in line:
                            #Qprep run OK: Update entry fields in parent (Qgui main or LIE)
                        if not self.qgui_parent:
                            if self.app.add_title == 'complex':
                                self.app.complex_pdb = self.app.workdir + '/' + self.top_pdbname
                                self.app.complex_top = self.app.workdir + '/' + self.topname
                            elif self.app.add_title == 'ligand':
                                self.app.ligand_pdb = self.app.workdir + '/' + self.top_pdbname
                                self.app.ligand_top = self.app.workdir + '/' + self.topname
                                self.app.update_progress()
                            else:
                                self.app.pdb_id = self.app.workdir + '/' + self.top_pdbname
                                self.app.top_id = self.app.workdir + '/' + self.topname
                                self.app.update_pdb_id_entryfield()
                                self.app.main_window.set_topology_entryfield(self.topname.split('/')[-1])

                        self.app.log('info','Topology successfully generated')
                        break
                    elif 'ERROR:' in line:
                        self.app.log('info','Topology generation failed.')
                        self.app.errorBox('Warning','Topology contains errors. Check log file and parameters.')
                        break
                    elif 'topology is incomplete' in line:
                        self.app.log('info', 'Topology contains errors! Please check log file!')
                        self.app.errorBox('Warning','Topology contains errors. Check log file and parameters.')
                        break
                    elif 'Correct the PDB file' in line:
                        self.app.log('info', 'Topology contains errors! Please check log file!')
                        self.app.errorBox('Warning','Topology contains errors. Check log file and parameters.')
                        break


        #Change back to original path:
        os.chdir(current_dir)

    def set_structure_entry(self):
        """Inserts the name of the current file to structure_entry field. """
        self.structure_entry.config(state = NORMAL)
        self.structure_entry.insert(0.0, self.pdbfile.split('/')[-1].rjust(19))
        self.structure_entry.config(state = DISABLED)

    def set_total_charge(self):

        """Gets the overall charge from countCharges() function from prepareTopology
        file. Inserts the charge into total_c_entry field. """
        self.total_charge = 0

        self.app.log(' ', '\n____________________________________________\n')
        self.app.log(' ', 'Total charge per residue type:')
        self.app.log(' ', '\n____________________________________________\n')
        for res in self.resname_nr_dist:
            if len(self.resname_nr_dist[res]) > 0:
                charge = self.res_charge[res]
                sum_charge = (charge * len(self.resname_nr_dist[res]))
                self.total_charge += sum_charge
                if abs(sum_charge) > 0:
                    self.app.log(' ', '    %4s   %7.2f\n' % (res, sum_charge))

        self.app.log(' ', '____________________________________________\n')
        self.app.log(' ', 'Sum charge %7.2f\n' % self.total_charge)
        self.app.log(' ', '============================================\n')

        self.total_c_entry.config(state=NORMAL)
        self.total_c_entry.delete(0.0, END)
        self.total_c_entry.insert(0.0, self.total_charge)
        self.total_c_entry.config(state=DISABLED)

    def turn_off_charges_button_pressed(self, off=True):
        """This will toggle on/off all charges within 5/6*r where r is the sim. sphere r."""
        charge_on = True

        if off:
            charge_on = False

        #Residues that are not automatically charged with this function call:
        no_toggle = ['HIE', 'HID', 'HIP']

        xc = float(self.center_x_entry.get())
        yc = float(self.center_y_entry.get())
        zc = float(self.center_z_entry.get())

        simrad = float(self.sphere_entry.get()) * 0.85

        #Go through all charged and togglable residues:
        for res in list(self.resname_nr_dist.keys()):
            #Is the residue defined to be togglable?
            if res in list(self.toggle_res.keys()):
                #Check if residue needs to be toggled:
                toggle_res = True

                #check if residue is terminal residue:
                if len(res) == 4 and res not in no_toggle:
                    if res[0].lower() in ['n','c']:
                        res = self.toggle_terminal(res, xc, yc, zc, simrad, charge_on)

                if res in no_toggle:
                    toggle_res = False

                #Turn on all charges?
                elif charge_on:
                    if abs(self.res_charge[res]) >= abs(self.res_charge[self.toggle_res[res]]):
                        toggle_res = False
                #Turn off all charges?
                else:
                    if abs(self.res_charge[res]) <= abs(self.res_charge[self.toggle_res[res]]):
                        toggle_res = False
                print((res, self.toggle_res[res]))
                print((self.res_charge[res], self.res_charge[self.toggle_res[res]]))

                #Toggle all residues of the given residue type:
                if toggle_res:
                    for res_nr in list(self.resname_nr_dist[res].keys()):
                        toggle_res_nr = True

                        #If charges on. Check distance to simulation sphere boundary
                        if charge_on:
                            x = self.resname_nr_dist[res][res_nr]['x']
                            y = self.resname_nr_dist[res][res_nr]['y']
                            z = self.resname_nr_dist[res][res_nr]['z']
                            r = np.sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2)

                            if r > simrad:
                                toggle_res_nr = False

                        if toggle_res_nr:
                            new_res = self.toggle_res[res]
                            if new_res not in list(self.resname_nr_dist.keys()):
                                self.resname_nr_dist[new_res] = dict()
                            self.resname_nr_dist[new_res][res_nr] = copy.deepcopy(self.resname_nr_dist[res][res_nr])
                            del self.resname_nr_dist[res][res_nr]
            else:
                print(('Residue %s is not defined as togglable in lib' % res))

        self.set_total_charge()
        self.updateList()

        self.highlight_charged()

    def toggle_terminal(self, res, xc, yc, zc, simrad, charge_on=False):
        """
        Special function to take care of terminal residues with chargable resiudes (arg, lys, asp, glu)
        """
        new_res = res
        org_res = res[1:]
        terminal = res[0]

        toggle_res = True

        #Is the residue chargable?
        if org_res in list(self.toggle_res.keys()):
            #Turn on charges?
            if charge_on:
                if abs(self.res_charge[org_res]) >= abs(self.res_charge[self.toggle_res[org_res]]):
                    toggle_res = False
            #Turn off charges?
            else:
                if abs(self.res_charge[org_res]) <= abs(self.res_charge[self.toggle_res[org_res]]):
                    toggle_res = False
        else:
            toggle_res = False

        if toggle_res:
            for res_nr in list(self.resname_nr_dist[res].keys()):
                toggle_res_nr = True

                #If charges on. Check distance to simulation sphere boundary
                if charge_on:
                    x = self.resname_nr_dist[res][res_nr]['x']
                    y = self.resname_nr_dist[res][res_nr]['y']
                    z = self.resname_nr_dist[res][res_nr]['z']
                    r = np.sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2)

                    if r > simrad:
                        toggle_res_nr = False

                if toggle_res_nr:
                    new_res = terminal + self.toggle_res[org_res]
                    #Double check that new residue actually is defined in LIB, even though it is put in the header to be togglable:
                    if new_res not in list(self.res_charge.keys()):
                        self.app.log(' ', 'WARNING: Tried to toggle %s --> %s\n-------> %s not found in LIB!'
                                          % (res, new_res, new_res))
                        new_res = res
                    else:
                        if new_res not in list(self.resname_nr_dist.keys()):
                            self.resname_nr_dist[new_res] = dict()
                        self.resname_nr_dist[new_res][res_nr] = copy.deepcopy(self.resname_nr_dist[res][res_nr])
                        del self.resname_nr_dist[res][res_nr]

        return new_res

    def toggle_charge_selected(self):
        """
        Function to toggle a selected residue. The toggle must be defined in the header of lib file.
        """

        try:
            list_index = int(self.listbox.curselection()[0])
        except:
            return

        oldRes, resnr = self.listbox.get(list_index).split()[0:2]

        if oldRes not in list(self.toggle_res.keys()):
            self.app.log(' ', 'Residue %s is not defined as togglable in lib file!\n' % oldRes)
            return

        radius = self.listbox.get(list_index).split('|')[-1]
        #TODO N/C-terminals with chargable residues (4 options)
        newRes = self.toggle_res[oldRes]

        if newRes not in list(self.res_atoms.keys()):
            self.app.log(' ', 'Residue %s is not defined in lib file!\n' % newRes)
            return

        resnr = int(resnr)

        if newRes not in list(self.resname_nr_dist.keys()):
            self.resname_nr_dist[newRes] = dict()

        self.resname_nr_dist[newRes][resnr] = copy.deepcopy(self.resname_nr_dist[oldRes][resnr])
        del self.resname_nr_dist[oldRes][resnr]

        #Update pymol window
        if self.session:
            self.session.stdin.write('set dot_color, %s, i. %s\n' % (self.res_colors[newRes], resnr))

        self.app.log('info','%s %s --> %s %s' % (oldRes,resnr,newRes,resnr))
        self.listbox.insert(list_index, '%4s %4d    |%5.1f' % (newRes, int(resnr), float(radius)))
        self.listbox.delete(list_index + 1 )
        self.listbox.selection_set(list_index)

        self.set_total_charge()

    def find_ss_bonds(self):
        """
        This funciton goes through all CYS residues and check for
        potential S-S bond partners
        """
        all_cys = self.resname_nr_dist['CYS'].copy()
        all_cys.update(self.resname_nr_dist[self.toggle_res['CYS']])

        #S-S bond length is around 2.05 angstrom
        r_max = 2.1

        self.cys_residues = dict()
        for nr1 in sorted(all_cys.keys()):
            for nr2 in sorted(all_cys.keys()):
                if nr1 != nr2:
                    x1 = all_cys[nr1]['x']
                    y1 = all_cys[nr1]['y']
                    z1 = all_cys[nr1]['z']

                    x2 = all_cys[nr2]['x']
                    y2 = all_cys[nr2]['y']
                    z2 = all_cys[nr2]['z']

                    r = float(np.sqrt((x2 -x1)**2 + (y2 -y1)**2 + (z2 -z1)**2))

                    #if SG distance is less than treshold, S-S bond is possible:
                    if r < r_max:
                        self.app.log(' ', 'Possible S-S for CYS %3d - CYS %3d  (r = %.2f A)\n' % (nr1, nr2, r))

                        self.cys_residues[nr1] = nr2

            del all_cys[nr1]

    def create_ss_bonds_checkbutton_pressed(self):
        """Find potential S-S bridges. Returns:
        1. Arrays with atomnumber pairs to be used in Qprep5: makebond atom_i atom_j
        2. Array with residue number pairs.
        3. List with distances between S-S below treshold. """
        state = self.check_variable.get()
        if state == 0:
            self.makeSS = False
        if state == 1:
            self.makeSS = True

        self.find_ss_bonds()

        if self.makeSS:
            print("Create S-S bonds selected:")
            if len(list(self.cys_residues.keys())) > 0:

                self.app.log('info','Generating S-S bonds for:')
                self.app.log(' ','    -------------------------\n')

                for cys1 in sorted(self.cys_residues.keys()):
                    cys2 = self.cys_residues[cys1]

                    self.app.log(' ','    CYX %3d - CYX %3d\n' % (cys1, cys2))

                    if cys2 in list(self.resname_nr_dist['CYS'].keys()):
                        self.resname_nr_dist['CYX'][cys2] = copy.deepcopy(self.resname_nr_dist['CYS'][cys2])
                        del self.resname_nr_dist['CYS'][cys2]

                    if cys1 in list(self.resname_nr_dist['CYS'].keys()):
                        self.resname_nr_dist['CYX'][cys1] = \
                            copy.deepcopy(self.resname_nr_dist['CYS'][cys1])
                        del self.resname_nr_dist['CYS'][cys1]

                self.app.log(' ','    -------------------------\n')

        else:
            if 'CYX' in list(self.resname_nr_dist.keys()):
                for nr in list(self.resname_nr_dist['CYX'].keys()):
                    self.resname_nr_dist['CYS'][nr] = copy.deepcopy(self.resname_nr_dist['CYX'][nr])
                    del self.resname_nr_dist['CYX'][nr]

        self.updateList()

        if self.session:
            self.highlight_charged()

    def get_charge_coordinates(self,charge,pdb,centre,rmin,rmax):
        """
        If checked, adds counter ions to neutralize system
        """
        xc,yc,zc = centre
        xc = np.float(xc)
        yc = np.float(yc)
        zc = np.float(zc)
        r = (np.float(rmin)+np.float(rmax))/2
        n = int(abs(charge))

        x_cart=[]
        y_cart=[]
        z_cart=[]

        coords = set()


        k=0.1 + 1.2*n
        if n == 1:
            start=(-1.+1./(2*n-1))
            increment = 1
        elif n == 2:
            start=(-1.+1./(n-1))
            increment = 0.5
        else:
            start=(-1.+1./(n-1))
            increment = (2.0-2.0/(n-1.))/(n-1.)
        delta = 0.01*r
        j=0
        while j < n:
            s=start+j*increment
            dr = rng.uniform(-1*delta,delta)
            theta = s*k
            phi = np.pi/2.*np.copysign(1,s)*(1.-np.sqrt(1.-abs(s)))
            x = (r)*np.cos(theta)*np.cos(phi)+xc
            y = (r)*np.sin(theta)*np.cos(phi)+yc
            z = (r)*np.sin(phi)+zc
            if j==0:
                coords.add((x,y,z))
                x_cart.append(x)
                y_cart.append(y)
                z_cart.append(z)
                j+=1
            else:
                if (x,y,z) not in coords:
                    coords.add((x,y,z))
                    x_cart.append(x)
                    y_cart.append(y)
                    z_cart.append(z)
                    j+=1
                else:
                    k = 0.1+2.2*(n+j)

        return x_cart,y_cart,z_cart

    def neutralize_system(self):

        # distribute a number of ions to counter a charged protein centred on the protein centre
        # and placed at midpoint between the protein radius and simulation sphere radius
        centre_of_mass = pt.centerofmass(self.pdbfile)
        protein_radius = pt.findRadius(self.pdbfile)
        sim_radius = self.sphere_entry.get()
        self.x_ion,self.y_ion,self.z_ion = self.get_charge_coordinates(self.total_charge,self.pdbfile,centre_of_mass,protein_radius,sim_radius)
        atom_type = "ATOM  "
        if self.total_charge < 0:
            atom_name = "Na  "
            res_name = "Na+ "
        elif self.total_charge > 0:
            atom_name = "Cl  "
            res_name = "Cl- "
        else:
            print("System already neutralized")

        print('Writing new pdb with neutralizing ions')
        atomnr = 0
        resnr = 0

        # write ions to a new pdb file
        file_name_and_extension = self.pdbfile.split('.')
        new_file = file_name_and_extension[0]+'_neut.'+file_name_and_extension[1]
        with open(self.pdbfile,'r') as old_pdb:
            lines = old_pdb.readlines()
            old_pdb.seek(0)
            for line in old_pdb:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM': # find out last current atom and residue number in file
                    atomnr = int(line[6:11])
                    resnr =  int(line[22:26])
            with open(new_file,'w') as new_pdb:
                for item in lines:
                    new_pdb.write(item)
                for i in range(len(self.x_ion)):
                    atomnr+=1
                    resnr+=1
                    x = round(self.x_ion[i],3)
                    y = round(self.y_ion[i],3)
                    z = round(self.z_ion[i],3)
                    new_line = '{:<6}{:>5}  {:<4}{:<4}{:>5}    {:> 8.3f}{:> 8.3f}{:> 8.3f}\n'.format(atom_type,atomnr,atom_name,res_name,resnr,x,y,z)
                    new_pdb.write(new_line)
        old_pdb.close()
        new_pdb.close()

        # move HOH entries to end of file and renumbers the atoms/residues
        HOH_lines = []
        with open(new_file,'r+') as new_pdb:
            all_lines = new_pdb.readlines()
            new_pdb.seek(0)
            for item in all_lines:
                if item[17:20] == 'HOH':
                    HOH_lines.append(item)
                else:
                    new_pdb.write(item)
            for item in HOH_lines:
                new_pdb.write(item)

            new_pdb.seek(0)


        if self.session:
            self.session.stdin.write('delete %s \n' % self.pdbfile)
            self.session.stdin.write('load %s \n' % new_file)

        self.pdbfile = new_file
        self.app.main_window.set_entryfield(new_file.split('/')[-1])
        self.app.update_pdb_id_entryfield()

        #update topo name fields
        self.checkLib()
        self.updateList()
        self.set_total_charge()
        self.topology_entry.delete(1.0, END)
        self.topo_pdb_entry.delete(1.0, END)
        self.topology_entry.insert(0.0,self.pdbfile.split('/')[-1].split('.')[0]+'.top')
        if '_top' not in self.pdbfile.split('/')[-1]:
            self.topo_pdb_entry.insert(0.0,self.pdbfile.split('/')[-1].split('.')[0]+'_top.pdb')
        else:
            self.topo_pdb_entry.insert(0.0, self.pdbfile.split('/')[-1])


    def updateList(self):
        """
        Function to update listbox with charged/neutral residues
        for toggle state selection.
        """

        self.listbox.delete(0, END)

        #Get current x,y,z simulation center:
        xc = float(self.center_x_entry.get())
        yc = float(self.center_y_entry.get())
        zc = float(self.center_z_entry.get())

        #Find N and C terminals: #TODO remove global variable when fixed toggle charge!
        #self.nterm_nr, self.nterm_res, self.cterm_nr, self.cterm_res = pt.findTerminals(self.pdbfile)

        #for residue in sorted(self.toggle_res.keys(), key=lambda s: s.lower()):
        for residue in sorted(list(self.resname_nr_dist.keys()), key=lambda s: s.lower()):
            #if residue in self.resname_nr_dist.keys():
            for nr in sorted(self.resname_nr_dist[residue].keys()):
                x = self.resname_nr_dist[residue][nr]['x']
                y = self.resname_nr_dist[residue][nr]['y']
                z = self.resname_nr_dist[residue][nr]['z']
                r = np.sqrt((xc - x)**2 + (yc - y)**2 + (zc - z)**2)
                self.listbox.insert(END, '%4s %4d    |%5.1f' % (residue, nr, r))

        #TODO write new function to generate terminals. Perhaps a manual option as well?

    def open_atom_select(self):
        """
        Opens a new window to select atom in pdb file as simulation center
        """
        self.select_pdb = AtomSelect(self, self.root, self.pdbfile,
                                     self.center_x_entry, self.center_y_entry, self.center_z_entry)
        self.select_pdb.configure(bg = self.main_color)
        self.select_pdb.title('Select simulation center')
        self.select_pdb.resizable()

    def not_available(self):
        """
        Function not available yet pop-up
        """
        self.app.errorBox('Warning','Sorry, this function is not implemented yet!')
        self.neutralize.set(0)

    def editPdb(self):
        """
        Opens up pdb file for editing
        """
        self.fileEdit = FileEdit(self, self.pdbfile)
        self.fileEdit.config(bg=self.main_color)
        self.fileEdit.title('Edit file')
        self.fileEdit.resizable()

    def toggle_pymol(self):
        """
        Start/close pymol when checkbutton is toggled
        """
        if self.sync_pymol.get() == 1:
            if self.session:
                try:
                    os.killpg(self.session.pid, signal.SIGTERM)
                except:
                    pass
            self.start_pymol()

        elif self.sync_pymol.get() == 0:
            try:
                os.killpg(self.session.pid, signal.SIGTERM)
            except:
                pass

    def start_pymol(self):
        """
        Start syncing of Q-atom selection to pymol
        """
        #QPyMol default settings ('set valence, 0.1' removed)
        self.pymol_settings = ['space cmyk', 'set sphere_scale, 0.4',
                               'set sphere_transparency, 0.7', 'set sphere_color, lightblue',
                               'set sphere_quality, 2', 'set stick_radius, 0.17', 'set stick_quality, 10',
                               'set defer_builds_mode, 3', 'set surface_quality, 1',
                               'set spec_power, 250', 'set spec_reflect, 2',
                               'set cartoon_fancy_helices, 1']

        self.app.pymol_running = True

        tmpfile = open(self.app.workdir+'/.tmpfile','wb')
        if 'darwin' in sys.platform:
            self.session = Popen(["pymol", "-p -x -i", "%s" % self.pdbfile], stdout=tmpfile, stdin=PIPE, universal_newlines=True, bufsize=0, preexec_fn=os.setsid)
        else:
            self.session = Popen(["pymol", "-p", "%s" % self.pdbfile], stdout=tmpfile, stdin=PIPE,  universal_newlines=True, bufsize=0, preexec_fn=os.setsid)

        self.update()
        len_log = 35

        self.session.stdin.flush()

        #Move pymol to workdir
        self.session.stdin.write('cd %s\n' % self.app.workdir)

        #Load default Q-PyMol settings:
        for settings in self.pymol_settings:
            self.session.stdin.write('%s\n' % settings)

        #Color all C-atoms gray:
        self.session.stdin.write('color gray, name C*\n')

        self.session.stdin.write('show cartoon, all\n')
        #Remove internal gui
        self.session.stdin.write('set internal_gui=0\n')

        #Set vdw to 1
        self.session.stdin.write('alter *, vdw=1\n')

        #Show simulation sphere
        self.update_simulation_sphere()

        #Show chargable/toggle residues
        self.highlight_charged()

        while self.session.poll() is None:
            try:
                if self.session:
                    self.update()
                else:
                    break
            except:
                break
            lines = 0
            with open(self.app.workdir + '/.tmpfile', 'r') as pymol_out:
                for line in pymol_out:
                    lines += 1
                    if lines > len_log:
                        len_log = lines
                        if 'You clicked' in line:
                            self.app.log(' ', line)
                            selected = ' '.join(line.split('/')[-2:])
                            if '`' in selected:
                                selected = ' '.join(selected.split('`'))
                            selected = selected.split('\n')[0]
                            print(selected)

                        if 'cmd.id_atom:' in line:
                            pass
                            #self.session.stdin.write('select none\n')
                            #atomnr = line.split()[2].split(')')[0]
                            #qatoms_list = self.qatoms_listbox.get(0, END)
                            #for q in range(len(qatoms_list)):
                            #    if int(atomnr) == int(qatoms_list[q].split()[1]):
                            #        self.qatoms_listbox.select_set(q)
                            #        self.list_q_atoms_event()
                            #        self.qatoms_listbox.yview(q)

            time.sleep(0.2)

        self.sync_pymol.set(0)
        self.sync_check.update()
        self.session = None
        self.app.pymol_running = False
        try:
            os.killpg(self.session.pid, signal.SIGTERM)
        except:
            pass

    def update_simulation_sphere(self, *args):
        if not self.session:
            return

        self.session.stdin.write('delete simsphere\n')

        x = self.center_x_entry.get()
        y = self.center_y_entry.get()
        z = self.center_z_entry.get()

        r = self.sphere_entry.get()

        self.session.stdin.write('pseudoatom simsphere, pos=(%s, %s, %s)\n' % (x, y, z))
        self.session.stdin.write('set nonbonded_transparency, 1, simsphere\n')
        self.session.stdin.write('alter simsphere, vdw=1\n')
        self.session.stdin.write('alter simsphere, resi=99999\n')
        self.session.stdin.write('set sphere_scale, %s, simsphere\n' % r)
        self.session.stdin.write('show spheres, simsphere\n')

    def adjust_simulation_sphere(self, *args):
        if not self.session:
            return

        r = self.sphere_entry.get()

        self.session.stdin.write('set sphere_scale, %s, simsphere\n' % r)

    def update_pymol_ions(self):
        return

    def highlight_charged(self):
        """
        This funciton shows all chargable sidechain as sticks with colored spheres to indicate defined charge.
        """

        if not self.session:
            return

        self.session.stdin.write('hide dots\n')

        color_cmd = {'red': 'set dot_color, red,',
                     'blue': 'set dot_color, blue,',
                     'white': 'set dot_color, white,',
                     'yellow': 'set dot_color, yellow,'}

        show_dots = 'show dots,'

        for i in self.listbox.get(0, END):
            resname = i.split()[0]
            color = self.res_colors[resname]
            resnr = i.split()[1]

            color_cmd[color] += ' i. %s or' % resnr
            show_dots += ' i. %s or' % resnr

        #set dot width
        self.session.stdin.write('set dot_width, 1\n')

        #set dot density
        self.session.stdin.write('set dot_density, 3\n')

        #Show dots
        self.session.stdin.write('%s\n' % show_dots[:-2])

        #Hide dots from backbone
        self.session.stdin.write('hide dots, name C or name N or name O or name H or name CA\n')

        #Color spheres
        for col in list(color_cmd.keys()):
            if len(color_cmd[col].split('i.')) > 1:
                self.session.stdin.write('%s\n' % color_cmd[col][:-2])

    def zoom_out(self):
        """
        Adjust automatic zoom buffer in pymol.
        """
        self.pymol_zoom += 1

        print(('PyMOL zoom buffer = %d' % self.pymol_zoom))

        if not self.session:
            return

        selections = list(map(int, self.listbox.curselection()))

        if len(selections) < 1:
            return

        self.session.stdin.write('zoom i. %s, buffer=%d\n' %
                                 (self.listbox.get(selections[0]).split()[1], self.pymol_zoom))

    def zoom_in(self):
        """
        Adjust automatic zoom buffer in pymol.
        """
        self.pymol_zoom -= 1

        if self.pymol_zoom < 0:
            self.pymol_zoom = 0

        print(('PyMOL zoom buffer = %d' % self.pymol_zoom))

        if not self.session:
            return

        selections = list(map(int, self.listbox.curselection()))

        if len(selections) < 1:
            return

        self.session.stdin.write('zoom i. %s, buffer=%d\n' %
                                 (self.listbox.get(selections[0]).split()[1], self.pymol_zoom))

    def listbox_clicked(self, *args):
        """
        If pymol is active, zoom selection and display atom names
        """
        if not self.session:
            return

        self.session.stdin.write('select none\n')
        selections = list(map(int, self.listbox.curselection()))

        self.session.stdin.write('hide labels\n')

        for selected in selections:
            self.session.stdin.write('label i. %s, name\n' % self.listbox.get(selected).split()[1])
            #self.session.stdin.write('select i. %s\n' % self.listbox.get(selected).split()[1])

        self.session.stdin.write('orient i. %s\n' % self.listbox.get(selections[0]).split()[1])

        self.session.stdin.write('zoom i. %s, buffer=%d\n' %
                                 (self.listbox.get(selections[0]).split()[1], self.pymol_zoom))

    def dialog_box(self):
        """Defines the outlook of Topology Prepare window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('Topology Prepare')
        self.config(background=self.main_color)

        # Define frames
        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=10, padx=(10,0))

        right_frame = Frame(self, bg = self.main_color)
        right_frame.pack(side = RIGHT)


        # Define elements in the left_frame
        sphere_label = Label(left_frame, text = 'Sphere:')
        sphere_label.grid(row = 0, column = 0, sticky='e')
        sphere_label.config(background = self.main_color)

        self.sphere_entry = Spinbox(left_frame, width = 5, highlightthickness = 0, relief = GROOVE, from_=1, to=200,
                                    textvariable=self.sphere_radius)
        self.sphere_entry.grid(row = 0, column = 1)

        aangstrom_label = Label(left_frame, text='%s' % '\xc5')
        aangstrom_label.grid(row=0, column=2, sticky='w')
        aangstrom_label.config(background=self.main_color)

        #Sync pymol
        sync_pymol = Label(left_frame, text='PyMOL:')
        sync_pymol.grid(row=0, column=3)
        sync_pymol.config(background=self.main_color)

        self.sync_check = Checkbutton(left_frame, variable=self.sync_pymol, command=self.toggle_pymol)
        self.sync_check.grid(row = 0, column = 4, sticky = 'w')
        self.sync_check.config(background = self.main_color)

        centre_label = Label(left_frame, text = 'Simulation centre:')
        centre_label.grid(row = 1, column = 0, sticky = 'e')
        centre_label.config(background = self.main_color)

        self.center_x_entry = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE,
                                    textvariable=self.xvar)
        self.center_x_entry.grid(row = 1, column = 1)


        self.center_y_entry = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE,
                                    textvariable=self.yvar)
        self.center_y_entry.grid(row = 1, column = 2)

        self.center_z_entry = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE,
                                    textvariable=self.zvar)
        self.center_z_entry.grid(row = 1, column = 3)

        center_x_label = Label(left_frame, text = 'x')
        center_x_label.grid(row = 2, column = 1, sticky='n')
        center_x_label.config(background = self.main_color)

        center_y_label = Label(left_frame, text = 'y')
        center_y_label.grid(row = 2, column = 2, sticky='n')
        center_y_label.config(background = self.main_color)

        center_z_label = Label(left_frame, text = 'z')
        center_z_label.grid(row = 2, column = 3, sticky='n')
        center_z_label.config(background = self.main_color)

        sim_change_button = Button(left_frame, text = 'Change', command = self.open_atom_select)
        sim_change_button.grid(row = 1, column = 4)
        sim_change_button.config(highlightbackground = self.main_color)

        sim_center_button = Button(left_frame, text = 'Center', command = self.set_sim_center)
        sim_center_button.grid(row = 2, column = 4)
        sim_center_button.config(highlightbackground = self.main_color)

        solvate_label = Label(left_frame, text = 'Solvate:')
        solvate_label.grid(row = 3, column = 0, sticky = 'e')
        solvate_label.config(background = self.main_color)

        tip3p_label = Label(left_frame, text = ' TIP3P')
        tip3p_label.grid(row = 4, column = 0, sticky = 'e')
        tip3p_label.config(background = self.main_color)

        self.check_solvation = IntVar()
        self.tip3p_radiobutton = Radiobutton(left_frame, variable = self.check_solvation, value=1)
        self.tip3p_radiobutton.grid(row = 4, column = 1, sticky = 'w')
        self.tip3p_radiobutton.config(background = self.main_color)

        spc_label = Label(left_frame, text = 'SPC')
        spc_label.grid(row = 5, column = 0, sticky = 'e')
        spc_label.config(background = self.main_color)

        self.scp_radiobutton = Radiobutton(left_frame, variable = self.check_solvation, value=2)
        self.scp_radiobutton.grid(row = 5, column = 1, sticky = 'w')
        self.scp_radiobutton.config(background = self.main_color)

        none_label = Label(left_frame, text = 'None')
        none_label.grid(row = 6, column = 0, sticky = 'e')
        none_label.config(background = self.main_color)

        self.none_radiobutton = Radiobutton(left_frame, variable = self.check_solvation, value = 0)
        self.none_radiobutton.grid(row = 6, column = 1, sticky = 'w')
        self.none_radiobutton.config(background = self.main_color)

        ssbond_label = Label(left_frame, text = 'Create S-S bonds:')
        ssbond_label.grid(row = 7, column = 0, sticky = 'e')
        ssbond_label.config(background = self.main_color)

        ssbond_checkbutton = Checkbutton(left_frame, variable = self.check_variable, command = self.create_ss_bonds_checkbutton_pressed)
        ssbond_checkbutton.grid(row = 7, column = 1, sticky = 'w')
        ssbond_checkbutton.config(background = self.main_color)

        total_c_label = Label(left_frame, text = 'Total charge:')
        total_c_label.grid(row = 8, column = 0, sticky = 'e')
        total_c_label.config(background = self.main_color)

        self.total_c_entry = Text(left_frame, width = 5, height = 1)
        self.total_c_entry.grid(row = 8, column = 1)
        self.total_c_entry.config(state = DISABLED, highlightthickness = 0)

        charge_all = Label(left_frame, text = 'Toggle all charges:', bg=self.main_color)
        charge_all.grid(row=9, column=0, sticky='e')

        charge_on_button = Button(left_frame, text = 'ON ',
                                  command = (lambda: self.turn_off_charges_button_pressed(False)))
        charge_on_button.grid(row = 9, column = 1, sticky = 'e')
        charge_on_button.config(highlightbackground = self.main_color)

        charge_off_button = Button(left_frame, text = 'OFF',
                                   command = (lambda: self.turn_off_charges_button_pressed(True)))
        charge_off_button.grid(row = 9, column = 2, sticky = 'w')
        charge_off_button.config(highlightbackground = self.main_color)


        toggle_sel = Label(left_frame, text = 'Toggle state:', bg=self.main_color)
        toggle_sel.grid(row=10, column=0, sticky='e')

        toggle_state_button = Button(left_frame, text = 'Selected', command = self.toggle_charge_selected)
        toggle_state_button.grid(row = 10, column = 1, columnspan=2)
        toggle_state_button.config(highlightbackground = self.main_color)

        neutralize_label = Label(left_frame, text = 'Neutralize system with NaCl:')
        neutralize_label.grid(row = 11, column = 0, sticky = 'e')
        neutralize_label.config(background = self.main_color)

        neutralize_checkbutton = Checkbutton(left_frame, variable=self.neutralize, command = self.neutralize_system)
        neutralize_checkbutton.grid(row = 11, column = 1, sticky = 'w')
        neutralize_checkbutton.config(background = self.main_color)

        nacl_label = Label(left_frame, text  ='(Na+ Cl-)')
        nacl_label.grid(row = 12, column = 0, sticky = 'e')
        nacl_label.config(background = self.main_color)

        topology_label = Label(left_frame, text = 'Topology name:')
        topology_label.grid(row = 13, column = 0, sticky = 'e',pady=(10,0))
        topology_label.config(background = self.main_color)

        self.topology_entry = Text(left_frame, width = 40, height = 1)
        self.topology_entry.grid(row = 13, column = 1, columnspan = 4, pady=(10,0), sticky = 'w')
        self.topology_entry.config(state = NORMAL, highlightthickness = 0, relief = GROOVE)

        topo_pdb_label = Label(left_frame, text = 'Topology PDB:')
        topo_pdb_label.grid(row = 14, column = 0, sticky = 'e')
        topo_pdb_label.config(background = self.main_color)

        self.topo_pdb_entry = Text(left_frame, width = 40, height = 1)
        self.topo_pdb_entry.grid(row = 14, column = 1, columnspan = 4, sticky='w')
        self.topo_pdb_entry.config(state = NORMAL, highlightthickness = 0, relief = GROOVE)

        edit_pdb = Button(left_frame, text='Edit', width=5, command=self.editPdb)
        edit_pdb.grid(row=4, column=4)
        edit_pdb.config(highlightbackground=self.main_color)

        self.structure_entry = Text(left_frame, width = 20, height = 1)
        self.structure_entry.grid(row = 4, column = 2, columnspan = 2, sticky='e')
        self.structure_entry.config(state = DISABLED, bg = self.main_color, highlightthickness = 0, relief = GROOVE)

        self.lib_status = Text(left_frame, width=20, height=1)
        self.lib_status.grid(row=5, column=2, columnspan=2, sticky='e')
        self.lib_status.config(state = DISABLED, bg = self.main_color, highlightthickness = 0, relief = GROOVE)

        checklib = Button(left_frame, text='Check', command=self.checkLib)
        checklib.grid(row=5, column=4)
        checklib.config(highlightbackground=self.main_color)

        #zoom button
        zoom_image = PhotoImage(file=self.app.qgui_path+ '/Qmods/zoom-in.gif')
        zoom_image.img = zoom_image

        zoom_out_image = PhotoImage(file=self.app.qgui_path+ '/Qmods/zoom-out.gif')
        zoom_out_image.img = zoom_out_image

        zoom_button = Button(left_frame, image=zoom_image, command=self.zoom_in, width=30, height=30)
        zoom_button.grid(row=7, column=3, sticky='e')

        zoom_out_button = Button(left_frame, image=zoom_out_image, command=self.zoom_out, width=30, height=30)
        zoom_out_button.grid(row=7, column=4, sticky='w')

        listbox_scroll = Scrollbar(left_frame)
        listbox_scroll.grid(row = 8, rowspan = 5, column = 7, sticky = 'nsw')
        self.listbox = Listbox(left_frame, yscrollcommand = listbox_scroll.set, highlightthickness = 0, relief = GROOVE)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 8, rowspan = 5, column = 3, columnspan = 4, sticky = 'e')
        self.listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        self.listbox.bind('<<ListboxSelect>>', self.listbox_clicked)


        create_top_button = Button(left_frame, text = 'Write', command = self.writeTopology)
        create_top_button.grid(row = 15, column = 2, pady=(20,10), sticky = 'e')
        create_top_button.config(highlightbackground = self.main_color)

        run_top_button = Button(left_frame, text = ' Run ', command = self.run_qprep)
        run_top_button.grid(row = 15, column = 3,  pady=(20,10))
        run_top_button.config(highlightbackground = self.main_color)

        cancel_button = Button(left_frame, text = 'Close', command = self.destroy)
        cancel_button.grid(row = 15, column = 4, pady=(20,10))
        cancel_button.config(highlightbackground = self.main_color)
