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

from Tkinter import Label, TOP, Button, Listbox, Scrollbar, EXTENDED, Spinbox, Entry, Text, Frame, \
    Toplevel, DISABLED, END, GROOVE, NORMAL, BOTH, OptionMenu, IntVar, StringVar, Checkbutton, HORIZONTAL, LabelFrame

from select_atoms import AtomSelectRange
from edit_file import FileEdit
from edit_evb import EditEvbNotes, ImportParameters, EditParameters, EditBondParameters, EditAngleParameters, \
    EditTorsionParameters, EditImproperParameters
from setup_md import SetupMd
from tkFileDialog import askopenfilename
import tkFont
import copy
import cPickle
import shutil
import random
import os
from subprocess import Popen, PIPE
import time
import signal
import sys
from subprocess import call
import numpy as np


class SetupFEP(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root

        self.topology = self.app.top_id
        self.pdbfile = self.app.pdb_id

        #Mass dictionary
        self.massDict = {'H': 1.008, 'C': 12.01, 'N': 14.01,'O': 16.00,'F': 19.00,'Na': 22.98, 'Mg': 24.305, 'P': 30.97,
                        'S': 32.07, 'Cl': 35.45, 'Ca': 40.078, 'Fe':55.847, 'Zn': 65.38, 'Br':79.90,'I':126.90}

        #Forcefields availiable in ffld_server
        self.forcefields = ('2005', '2001')

        #Control if new inputfiles are to be made and old overwritten or not
        self.overwrite = True
        self.check_overwrite = IntVar()
        self.check_overwrite.set(1)

        #Automatically chop Q atom region on Auto assign, or if checked, select H-add atoms
        self.select_h = IntVar()
        self.select_h.set(0)

        #Use softcoret max potential
        self.use_soft_max_pt = IntVar()
        self.use_soft_max_pt.set(0)

        #Input files written?
        self.files_written = False

        #Number of EVB states to include:
        self.evb_states = IntVar()
        self.evb_states.set(2)
        self.set_evb_states = StringVar()

        #Some functions and buttons are only available in state3 and state4:
        self.state3_enabled = list()
        self.state4_enabled = list()

        #Get libraries and parameter files:
        self.libs = self.app.libs
        self.prms = self.app.prms

        self.fep_written = False

        #Molecular dynamics settings:
        self.md_settings = {'simtime': 0.01,
                            'stepsize': 1.0,
                            'inputfiles': 51,
                            'bath_coupling': 10,
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
                            'trajectory': 100,
                            'trajectory atoms': 'not excluded',
                            'seq_rest': [],
                            'atom_rest': [],
                            'dist_rest': [],
                            'wall_rest': []}


        #Initialize EVB status and progress:
        self.qstatus = {'Q-atoms':'No Q-atoms selected',
                        'Topology pdb': 'NA',
                        u"\N{GREEK SMALL LETTER LAMDA}-steps/run": '51',
                        'Total simulation time (ns)': '%.6f' % (float(self.md_settings['simtime']) * 51.00)}


        #Pymol session:
        self.session = None

        #Initialize global dictionaries
        #Keep track of Q-atom nr to atom nr/names/types etc
        self.q_atom_nr = dict()
        self.q_atom_res = dict()
        self.q_atom_name = dict()
        self.q_notes = dict()

        # {Qi : [type state1, type state2, type state3, type state4]}
        self.q_atomtypes = dict()

        #self.atomtype_prm[atomtype] = [Ri, Ei, Ci, ai, Ri(1-4), Ei(1-4), mass]
        self.atomtype_prm = dict()

        # {Qi : [state1, state2, state3, state4]}
        self.q_charges = dict()

        # {Qi : [[state1 Qjs],[state22 Qjs], [state3 Qjs], [state4 Qjs]]}
        self.q_bonds = dict()

        # {'type1 type2' : [De, alpha, R, Kb]}
        self.bond_prm = dict()

        # {'atom_i':[state1, state2, state3, state4]}
        #Use atom numbers instead of Q nr to avoid errors upon deleting Q-atoms during session ...
        self.change_bonds = dict()
        self.fep_atoms = dict()

        # {'type1 type2 type3: [K, theta]}
        self.angle_prm = dict()

        # {'type1 type2 type3 type4': [[Kt, phase, paths], [Kt, phase, paths], [Kt, phase, paths]]}
        # ==> {'type1 type2 type3 type4': [[min1], [min2], [min3]]
        self.torsion_prm = dict()

        # {q2 : [q1,q3,q4]}    q1             q1
        #                      |              ||
        #                   q3-q2-q4       q3-q2-q4
        # q2 is sp2 and connected to q1, q3 and q4 <-- This definition makes it easy to use with q_bonds and
        # connect it to bonds being broken or formed.
        self.q_impropers = dict()

        #{'type1 type2 type3 type4': [K, phase]}
        self.improper_prm = dict()

        #Soft-pair list (for additional pairs added by user)
        # {Qi : Qj}
        self.softpairs_added = dict()

        # {'atom1 atom2' : [s1, s2, s3, s4]}
        self.excluded_pairs = dict()

        # {'q1 q2' : scale}
        self.qpair_elscale = dict()

        # {Qi: [s1, s2, s3, s4]}
        self.q_softcore = dict()

        # {step nr: [lambda1, lambda2, lambda3, lambda4]}
        self.step_lambdavalues = dict()

        # {Qi : [s1, s2...]}
        self.q_element = dict()

        # {1 : [atom_i, atom_j,...], 2:[..], ...}
        self.monitor_groups = dict()

        # [[group_i, group_j], [group_k, group_l], [..], ...]
        self.monitor_pairs = list()

        #Initialize default step_lambdavalues:
        step = 1
        l1 = 1.00
        l2 = 0.00
        l3 = 0.00
        l4 = 0.00
        self.step_lambdavalues[step] = [l1, l2, l3, l4]
        while l1 > 0.00:
            step += 1
            l1 -= 0.02
            l2 += 0.02
            self.step_lambdavalues[step] = [abs(l1), abs(l2), l3, l4]

        #Temperature(s)
        # self.temp_runs = {T : Runs}
        self.temp_runs = {'300.00' : 1}

        self.setup_evb = StringVar()
        self.sync_pymol = IntVar()
        #When something in pymol is clicked, selection will be removed when Q-atom listbox is clicked:
        self.hide_selected_pymol = False

        self.setup_evb.set('Assign')
        self.setup_evb.trace('w', self.show_var_frame)

        #Choose if all angles are shown or only changing angles
        self.all_angles = IntVar()
        self.all_angles.set(0)
        self.show_changing_angles = False

        #Choose if all torsions are shown or only changing torsions
        self.all_torsions = IntVar()
        self.all_torsions.set(0)
        self.show_changing_torsion = False

        #Choose if all impropers are shown or only changing impropers
        self.all_impropers = IntVar()
        self.all_impropers.set(0)
        self.show_changing_impropers = False

        #Lambda values to trace
        self.l1_start = StringVar()
        self.l1_end = StringVar()
        self.l2_start = StringVar()
        self.l2_end = StringVar()
        self.l3_start = StringVar()
        self.l3_end = StringVar()
        self.l4_start = StringVar()
        self.l4_end = StringVar()

        #Choose what lambda values do vary in lambda step setup:
        self.lock_lambda1 = IntVar()
        self.lock_lambda1.set(0)
        self.lock_lambda2 = IntVar()
        self.lock_lambda2.set(0)
        self.lock_lambda3 = IntVar()
        self.lock_lambda3.set(1)
        self.lock_lambda4 = IntVar()
        self.lock_lambda4.set(1)

        self.lambdastep = StringVar()
        self.lambdastep.set(0.02)

        self.dialog_window()
        if self.app.pdb_id:
            self.pdbfile = self.app.pdb_id
            self.qstatus['Topology pdb'] = 'OK'
            self.sync_check.config(state=NORMAL)
        else:
            self.pdbfile = None
            self.qstatus['Topology pdb'] = 'No topology PDB file loaded'
            self.sync_check.config(state=DISABLED)

        #Insert default lambda-step values
        self.update_lambdasteps()

        #Trace if number of EVB states are changed:
        self.set_evb_states.trace('w', self.evb_states_changed)

        #Trace changes in lambda values
        self.l1_start.trace('w', self.sum_lambda_values)
        self.l2_start.trace('w', self.sum_lambda_values)
        self.l3_start.trace('w', self.sum_lambda_values)
        self.l4_start.trace('w', self.sum_lambda_values)
        self.l1_end.trace('w', self.sum_lambda_values)
        self.l2_end.trace('w', self.sum_lambda_values)
        self.l3_end.trace('w', self.sum_lambda_values)
        self.l4_end.trace('w', self.sum_lambda_values)

        self.update_status()

    def evb_states_changed(self, *args):
        """
        Update windows and layout based on if 2,3 or 4- state EVB is selected
        """
        self.evb_states.set(int(self.set_evb_states.get().split()[0]))

        self.app.log(' ', 'Switched to %d states FEP.\n' % self.evb_states.get())

        #Delete charges for state3 and 4 in change charges
        self.charge3.delete(0, END)
        self.charge4.delete(0, END)

        #Remove pymol states if existing
        if self.session:
            self.session.stdin.write('disable state3\ndisable state4\n')

        for state_funcs in (self.state3_enabled, self.state4_enabled):
            for func in state_funcs:
                func.config(state=DISABLED)

        if self.evb_states.get() == 3:
            for func in self.state3_enabled:
                func.config(state=NORMAL)
            if self.session:
                self.session.stdin.write('enable state3\n')
        if self.evb_states.get() == 4:
            for state_funcs in (self.state3_enabled, self.state4_enabled):
                for func in state_funcs:
                    func.config(state=NORMAL)
            if self.session:
                self.session.stdin.write('enable state3\n')
                self.session.stdin.write('enable state4\n')

        self.update_all()
        self.update_excluded_pairs()
        self.update_lambdasteps()
        self.show_var_frame()

    def sum_lambda_values(self, *args):
        try:
            l1_start = float(self.start1.get())
            l1_end = float(self.end1.get())

            l2_start = float(self.start2.get())
            l2_end = float(self.end2.get())

            l3_start = float(self.start3.get())
            l3_end = float(self.end3.get())

            l4_start = float(self.start4.get())
            l4_end = float(self.end4.get())

            sum_start = (l1_start + l2_start + l3_start + l4_start)
            sum_end = (l1_end + l2_end + l3_end + l4_end)

            self.sum_start.delete(0, END)
            self.sum_start.insert(0, '%.3f' % sum_start)

            self.sum_end.delete(0, END)
            self.sum_end.insert(0, '%.3f' % sum_end)
        except:
            return

    def toggle_show_ang(self):
        """
        Decide if all angles are to be included, or only those involved in changing Q-atoms
        """
        if self.all_angles.get() == 1:
            self.show_changing_angles = True
        else:
            self.show_changing_angles = False

        self.fep_written = False
        self.update_angles()

    def toggle_show_tor(self):
        """
        Decide if all torsions are to be included, or only those involved in changing Q-atoms
        """
        if self.all_torsions.get() == 1:
            self.show_changing_torsion = True
        else:
            self.show_changing_torsion = False

        self.fep_written = False
        self.update_torsions()

    def toggle_show_imp(self):
        """
        Decide if all impropers are to be included, or only those involved in changing Q-atoms
        """
        if self.all_impropers.get() == 1:
            self.show_changing_impropers = True
        else:
            self.show_changing_impropers = False

        self.fep_written = False
        self.update_impropers()

    def add_qatoms(self):
        """
        Select qatoms from pdb file and add to qatoms list
        Fills:  self.q_atom_nr
                self.q_atom_res
                self.q_atom_name
        """
        if not self.pdbfile:
            #Ask user to select pdb file from topology:
            pdbfile = askopenfilename(parent = self, initialdir = self.app.workdir, filetypes=(("pdb", "*.pdb"),("All files","*.*")))

            if pdbfile:
                self.pdbfile = pdbfile
                self.update_status()
                self.sync_check.config(state=NORMAL)
            else:
                print 'No pdb loaded'
                self.sync_check.config(state=DISABLED)
                return

        self.q_index = None

        self.select_qatoms = AtomSelectRange(self, self.root, self.pdbfile, self.q_atom_nr, 'qatoms')


    def del_qatoms(self):
        """
        removes selected q-atsoms from list
        """
        selected = sorted(map(int, self.qatoms_listbox.curselection()))
        selected.reverse()

        stored = [self.q_atom_nr, self.q_atom_res, self.q_atom_name, self.q_notes,
                  self.q_atomtypes, self.q_charges]

        for selection in selected:
            q = int(self.qatoms_listbox.get(selection).split()[0])

            for qlist in stored:
                del qlist[q]

            #Remove Q-atoms from self.q_bonds and q_impropers
            if q in self.q_bonds.keys():
                del self.q_bonds[q]
            for qi in self.q_bonds.keys():
                for state in range(4):
                    if q in self.q_bonds[qi][state]:
                        del self.q_bonds[qi][state][self.q_bonds[qi][state].index(q)]

            if q in self.q_impropers.keys():
                del self.q_impropers[q]

            for qi in self.q_impropers.keys():
                if q in self.q_impropers[qi]:
                    del self.q_impropers[qi]

        #Update Q-atom numbers:
        q_old_new = dict ()
        qnr = 0
        for q_old in sorted(self.q_atom_nr.keys()):
            qnr += 1
            if q_old != qnr:
                q_old_new[q_old] = qnr
                for qlist in stored:
                    qlist[qnr] = qlist[q_old]
                    del qlist[q_old]

        #UPDATE Q-atom number in bonds and impropers
        for q_old in sorted(q_old_new.keys()):
            if q_old in self.q_bonds.keys():
                q_new = q_old_new[q_old]
                self.q_bonds[q_new] = copy.deepcopy(self.q_bonds[q_old])
                del self.q_bonds[q_old]

            for qj in self.q_bonds.keys():
                for state in range(4):
                    print self.q_bonds[qj][state]
                    if q_old in self.q_bonds[qj][state]:
                        print 'befor: %s' % self.q_bonds[qj][state]
                        self.q_bonds[qj][state][self.q_bonds[qj][state].index(q_old)] = q_old_new[q_old]
                        print 'after: %s' % self.q_bonds[qj][state]

            if q_old in self.q_impropers.keys():
                q_new = q_old_new[q_old]
                self.q_impropers[q_new] = copy.deepcopy(self.q_impropers[q_old])
                del self.q_impropers[q_old]
                for i in range(len(self.q_impropers[q_new])):
                    if self.q_impropers[q_new][i] in q_old_new.keys():
                        self.q_impropers[q_new][i] = q_old_new[self.q_impropers[q_new][i]]

        self.update_q_atoms()
        self.update_all()

    def note_qatoms(self):
        """
        Add a note to selected q-atom
        """
        sel_index = map(int, self.qatoms_listbox.curselection())

        if len(sel_index) == 1:
            note = self.qatoms_listbox.get(sel_index[0]).split('!')[-1]
            q = int(self.qatoms_listbox.get(sel_index[0]).split()[0])
            self.edit_note = EditEvbNotes(self, self.root, q, sel_index[0], note)
        else:
            self.app.errorBox('Info','Select exactly one Q-atom entry.')
            return

    def update_q_atoms(self):
        """
        Whenever new Q-atoms are added/deleted, the Q-atoms listbox is updated
        """

        #Check if Q notes exist, if not, make it atomname and res i:
        if len(self.q_notes) == 0:
            for q in self.q_atom_name.keys():
                self.q_notes[q] = '%4s %s' % (self.q_atom_name[q].ljust(4), self.q_atom_res[q])

        #Fill Q atoms listbox
        self.qatoms_listbox.delete(0, END)
        for q in sorted(self.q_atom_nr.keys()):
            self.qatoms_listbox.insert(END, '%2d %6d !%s' % (q, self.q_atom_nr[q], self.q_notes[q]))

    def read_lib(self):
        """
        Reads library files defined in settings and updates:
            self.q_charges
            self.q_bonds
            self.q_impropers
            self.q_atomtypes
        NOTE: Everything is being set up to do up to 4-state EVB
        """
        #For each residue, find all bonds, atomtypes, charges and impropers
        residues = []
        res_nr = []

        # {'resi': [atomnames]}
        improper_names = dict()

        for q in sorted(self.q_atom_res.keys()):
            resi = self.q_atom_res[q]
            if resi not in residues:
                residues.append(resi)
                res_nr.append(int(resi.split()[-1]))

        for resi in residues:
            #The  residue name without the residue number:
            res = resi.split()[0].strip()
            #Collect all Q-atoms and atomnames selected within residue:
            qatoms = []
            atomnames = []
            for q in sorted(self.q_atom_res.keys()):
                if self.q_atom_res[q] == resi:
                    qatoms.append(q)
                    atomnames.append(self.q_atom_name[q])

            for libfile in self.libs:
                found_res = False
                found_atoms = False
                found_bonds = False
                found_connections = False
                found_impropers = False

                with open(libfile, 'r') as lib:
                    for line in lib:
                        if found_res:
                            #### IMPROPERS ####
                            if found_impropers:
                                imp_atoms = True
                                q_imp = [0,0,0,0]
                                if len(line.split()) == 4:
                                    atoms = line.split()[0:]
                                    for atom in range(len(atoms)):
                                        if '+' not in atoms[atom] and '-' not in atoms[atom]:
                                            if atoms[atom] not in atomnames:
                                                imp_atoms = False
                                            else:
                                                q_imp[atom] = qatoms[atomnames.index(atoms[atom])]
                                        if '+' in atoms[atom]:
                                            imp_atoms = False
                                            next_atom = atoms[atom].split('+')[-1]
                                            #Need to check that next_atom exist in next res selection
                                            current_nr = int(resi.split()[-1])
                                            if (current_nr + 1) in res_nr:
                                                next_res = residues[res_nr.index(current_nr + 1)]
                                                #Get all Q-atoms in next residue:
                                                for q_ in self.q_atom_res.keys():
                                                    if self.q_atom_res[q_] == next_res:
                                                        if self.q_atom_name[q_] == next_atom:
                                                            q_imp[atom] = q_
                                                            imp_atoms = True

                                        if '-' in atoms[atom]:
                                            imp_atoms = False
                                            prev_atom = atoms[atom].split('-')[-1]
                                            #Need to check that prev_atom exist in previous res selection
                                            current_nr = int(resi.split()[-1])
                                            if (current_nr - 1) in res_nr:
                                                prev_res = residues[res_nr.index(current_nr - 1)]
                                                #Get all Q-atoms in previous residue:
                                                for q_ in self.q_atom_res.keys():
                                                    if self.q_atom_res[q_] == prev_res:
                                                        if self.q_atom_name[q_] == prev_atom:
                                                            q_imp[atom] = q_
                                                            imp_atoms = True

                                    if imp_atoms:
                                        q1, q2, q3, q4 = q_imp[0:]
                                        self.q_impropers[q2] = [q1, q3, q4]

                            #### CONNECTIONS ####
                            if found_connections:
                                if 'head' in line:
                                    #Connect to previous residues Tail
                                    #Check if Head is in Q-atom selection
                                    atom_head = line.split()[1]
                                    if atom_head in atomnames:
                                        q_head = qatoms[atomnames.index(atom_head)]
                                        #Need to check if previous residue exist:
                                        current_nr = int(resi.split()[-1])
                                        if (current_nr - 1) in res_nr:
                                            #Check if previous residue contains Tail
                                            prev_res = residues[res_nr.index(current_nr - 1)]

                                            #Get all Q-atoms in previous residue:
                                            q_in_prev = []
                                            for q_ in self.q_atom_res.keys():
                                                if self.q_atom_res[q_] == prev_res:
                                                    q_in_prev.append(q_)

                                            #Check if any q_in_prev are + (Tail)
                                            found_tail = False
                                            for q_ in q_in_prev:
                                                #Make sure that Q-atom is involved in bond (it may as well not be..)
                                                if q_ in self.q_bonds.keys():
                                                    for state in range(4):
                                                        if '+' in self.q_bonds[q_][state]:
                                                            q_tail = q_
                                                            found_tail = True

                                                #Append connecting Q-atoms to self.q_bonds
                                                if found_tail:
                                                    if q_head not in self.q_bonds.keys():
                                                        self.q_bonds[q_head] = [[q_tail], [q_tail], [q_tail], [q_tail]]
                                                    else:
                                                        for state in range(4):
                                                            if q_tail not in self.q_bonds[q_head][state]:
                                                                self.q_bonds[q_head][state].append(q_tail)

                                                    if q_tail not in self.q_bonds.keys():
                                                        self.q_bonds[q_tail] = [[q_head], [q_head], [q_head], [q_head]]
                                                    else:
                                                        for state in range(4):
                                                            if q_head not in self.q_bonds[q_tail][state]:
                                                                self.q_bonds[q_tail][state].append(q_head)


                                if 'tail' in line:
                                    atom_tail = line.split()[1]
                                    if atom_tail in atomnames:
                                        q_tail = qatoms[atomnames.index(atom_tail)]
                                        if q_tail not in self.q_bonds.keys():
                                            self.q_bonds[q_tail] = [['+'], ['+'], ['+'], ['+']]
                                        else:
                                            for state in range(4):
                                                if '+' not in self.q_bonds[q_tail][state]:
                                                    self.q_bonds[q_tail][state].append('+')
                            #### BONDS ####
                            if found_bonds:
                                #if not len(line.split()[0]) > 4:
                                if not '[' in line and not ']' in line:
                                    try:
                                        atom1, atom2 = line.split()[0:2]
                                    except:
                                        continue
                                    #Check if atom1 and atom2 are selected Q-atoms.
                                    if atom1 in atomnames and atom2 in atomnames:
                                        #Bond exist in Q-atom selection, get Q-atom nr and append bond:
                                        qi = qatoms[atomnames.index(atom1)]
                                        qj = qatoms[atomnames.index(atom2)]
                                        if qi not in self.q_bonds.keys():
                                            self.q_bonds[qi] = [[qj], [qj], [qj], [qj]]
                                        else:
                                            for state in range(4):
                                                if qj not in self.q_bonds[qi][state]:
                                                    self.q_bonds[qi][state].append(qj)
                                        if qj not in self.q_bonds.keys():
                                            self.q_bonds[qj] = [[qi], [qi], [qi], [qi]]
                                        else:
                                            for state in range(4):
                                                if qi not in self.q_bonds[qj][state]:
                                                    self.q_bonds[qj][state].append(qi)

                            if found_atoms:
                                if not '[bonds]' in line:
                                    try:
                                        if line.split()[0].isdigit():
                                            atomname = line.split()[1]
                                            if atomname in atomnames:
                                                #Atomname in libfile is selected as Q-atom.
                                                #Find Q-atom nr and assign Q-atom type and charge:
                                                qi = qatoms[atomnames.index(atomname)]
                                                type1 = line.split()[2]
                                                charge = float(line.split()[3])
                                                self.q_atomtypes[qi] = [type1, type1, type1, type1]
                                                self.q_charges[qi] = [charge, charge, charge, charge]
                                    except:
                                        continue



                            if '[atoms]' in line:
                                found_atoms = True
                            if '[bonds]' in line:
                                found_atoms = False
                                found_bonds = True
                            if '[connections]' in line:
                                found_bonds = False
                                found_connections = True
                            if '[impropers]' in line:
                                found_impropers = True
                                found_atoms = False
                                found_bonds = False
                                found_connections = False
                            if '[charge_groups]' in line:
                                found_impropers = False
                                found_atoms = False
                                found_bonds = False
                                found_connections = False
                            if '*-----' in line:
                                found_res = False
                                break

                        if '{' in line and res in line:
                            found_res = True


        #Go through self.q_bonds and remove all '+' (Tail symbols)
        for q in self.q_bonds.keys():
            for state in range(4):
                if '+' in self.q_bonds[q][state]:
                    del self.q_bonds[q][state][self.q_bonds[q][state].index('+')]

        #Update all lists:
        self.update_all()

    def get_atomtype_parameters(self):
        """
        Goes through parameters files and collects atomtype parameters for atomtypes in self.q_atomtypes
        ==> updates self.atomtype_prm (does not overwrite existing parameters!)
        """
        if len(self.q_atomtypes) == 0:
            self.app.log(' ', 'No atomtypes found. Can not update.')
            return

        #Make a list with atomtypes (no redundancies)
        atomtypes_to_find = []
        for q in self.q_atomtypes.keys():
            for state in range(self.evb_states.get()):
                atomtype = self.q_atomtypes[q][state]
                if atomtype not in atomtypes_to_find:
                    atomtypes_to_find.append(atomtype)

        not_in_line = ['*','-','!', '']
        for parameterfile in self.prms:
            found_atomtypes = False
            with open(parameterfile, 'r') as prm:
                for line in prm:
                    if found_atomtypes:
                        try:
                            if line.split()[0].strip()[0] not in not_in_line:
                                atomtype = line.split()[0].strip()
                                if atomtype in atomtypes_to_find:
                                    #Check that atomtype is not already in list (may be added/edited by user):
                                    if atomtype not in self.atomtype_prm.keys():
                                        ri = float(line.split()[1].strip())
                                        ei = float(line.split()[3].strip())
                                        ci = 0.0
                                        if atomtype[0] == 'H':
                                            ci = 0.0
                                        ai = 0.0
                                        ri1_4 = float(line.split()[4].strip())
                                        ei_1_4 = float(line.split()[5].strip())
                                        mass = float(line.split()[6].strip())
                                        self.atomtype_prm[atomtype] = [ri, ei, ci, ai, ri1_4, ei_1_4, mass]
                        except:
                            continue
                    if '[atom_types]' in line:
                        found_atomtypes = True
                    if '[bonds]' in line:
                        break
        print self.atomtype_prm

    def update_atom_parameters(self):
        """
        Updates atomtype listboxes after reading parameter file(s) or editing parameters.
        ==> self.atomtypes_listbox
        ==> self.changetypes_listbox
        """
        #Update atomtype listboxes for all n (n = 2,3 or 4) states:
        prm_nr = dict()
        nr = 0
        self.changetypes_listbox.delete(0, END)
        for q in sorted(self.q_atomtypes.keys()):
            state_types = [' ', ' ', ' ', ' ']
            for state in range(self.evb_states.get()):
                atomtype = self.q_atomtypes[q][state]
                state_types[state] = atomtype
                if atomtype not in prm_nr.keys():
                    nr += 1
                    prm_nr[atomtype] = nr

            t1, t2, t3, t4 = state_types[0:]
            self.changetypes_listbox.insert(END, '%2d  %4s %4s %4s %4s' % (q, t1.ljust(4), t2.ljust(4),
                                                                               t3.ljust(4), t4.ljust(4)))

        print prm_nr

        #Insert atomtypes to self.atomtypes_listbox
        self.atomtypes_listbox.delete(0, END)
        for prm in prm_nr.keys():
            if prm in self.atomtype_prm.keys():
                ri, ei, ci, ai, ri1_4, ei1_4, mass  = self.atomtype_prm[prm][0:]
                self.atomtypes_listbox.insert(END, '%4s %7.2f %5.2f %5.2f %4.2f %7.2f %5.2f %5.2f' %
                                               (prm.ljust(4), ri, ei, ci, ai, ri1_4, ei1_4, mass))
            else:
                self.atomtypes_listbox.insert(END,'%4s %7s %5s %5s %4s %7s %5s %5s' %
                                               (prm.ljust(4), '??', '??', '??', '??', '??', '??', '??'))

    def update_charges(self):
        """
        Takes all Q-atom charges from self.q_charges and updates each state in listbox
        ==> self.charge_listbox
        """
        self.charge_listbox.delete(0, END)
        charge_sum = [0,0,0,0]
        for q in sorted(self.q_charges.keys()):
            charges = [' ', ' ', ' ', ' ']
            types = [' ', ' ', ' ', ' ']
            for state in range(self.evb_states.get()):
                charges[state] = self.q_charges[q][state]
                types[state] = self.q_atomtypes[q][state]
                try:
                    charge_sum[state] += self.q_charges[q][state]
                except:
                    continue
            q1, q2, q3, q4 = charges[0:]
            t1, t2, t3, t4 = types[0:]
            self.charge_listbox.insert(END, '%3d %6s %6s %6s %6s !%4s %4s %4s %4s' %
                                            (q, str(q1), str(q2), str(q3), str(q4), t1.ljust(4), t2.ljust(4),
                                             t3.ljust(4), t4.ljust(4)))

        #Insert sum of charges for each state
        charge_sum = map(lambda x: round(x, 3), charge_sum)

        sum1, sum2, sum3, sum4 = charge_sum[0:]
        if self.evb_states.get() < 4:
            sum4 = ' '
        if self.evb_states.get() < 3:
            sum3 = ' '
        if len(self.q_atom_nr) > 0:
            self.charge_listbox.insert(END, '--------------------------------')
            self.charge_listbox.insert(END, 'SUM %6s %6s %6s %6s' % (sum1, sum2, sum3, sum4))

        #CHECK SUM CHARGES
        charge_sum = [sum1, sum2, sum3, sum4]
        check_charge = []
        for charge in charge_sum:
            if charge != ' ':
                charge = round(charge, 3)
                if charge not in check_charge:
                    check_charge.append(charge)
        if len(check_charge) > 1:
            self.qstatus['Charge sum'] = 'WARNING! Not identical between states!'
        else:
            self.qstatus['Charge sum'] = 'OK'

        self.update_status()

    def update_softcore(self):
        """
        Update listbox from self.q_softcore dict
        """
        self.softcore_listbox.delete(0, END)

        for q in sorted(self.q_softcore.keys()):
            soft = [' ', ' ', ' ', ' ']
            for state in range(self.evb_states.get()):
                soft[state] = self.q_softcore[q][state]
            s1, s2, s3, s4 = soft[0:]
            self.softcore_listbox.insert(END, '%3d %6s %6s %6s %6s ' %
                                            (q, str(s1), str(s2), str(s3), str(s4)))

    def get_bond_parameters(self):
        """
        Goes through parameters files and collects bond parameters
        ==> Update self.bond_prm
        """
        #Find all unique bonds:
        bonds = []

        for q1 in self.q_bonds.keys():
            for state in range(self.evb_states.get()):
                for q2 in self.q_bonds[q1][state]:
                    t1 = self.q_atomtypes[q1][state]
                    t2 = self.q_atomtypes[q2][state]
                    bond = '%s %s' % (t1, t2)
                    bond_rev = '%s %s' % (t2, t1)
                    if bond not in bonds and bond_rev not in bonds:
                        #If parameters already exist, do not look for them:
                        if bond not in self.bond_prm.keys() and bond_rev not in self.bond_prm.keys():
                            bonds.append(bond)
        print bonds
        if len(bonds) == 0:
            return

        #Go through parameter files and try to find bond parameters:
        for parameterfile in self.prms:
            found_bonds = False
            with open(parameterfile, 'r') as prm:
                for line in prm:
                    if '[angles]' in line:
                        break
                    if found_bonds:
                        try:
                            type1 = line.split()[0]
                            type2 = line.split()[1]
                            prm_bond = '%s %s' % (type1, type2)
                            prm_bond_rev = '%s %s' % (type2, type1)

                            if prm_bond in bonds or prm_bond_rev in bonds:
                                kb = float(line.split()[2])
                                rb = float(line.split()[3])
                                de = (kb / 8.0)
                                alpha = 2.0
                                self.bond_prm[prm_bond] = [de, alpha, rb, kb]
                        except:
                            continue
                    if '[bonds]' in line:
                        found_bonds = True
        print self.bond_prm

    def find_non_dummy(self,states, du='du'):
        for state in states:
            if state != du:
                du = state
                break

        return du

    def update_q_bonds(self):
        """
        Collects all bonds for all states and updates listboxes.
        ==> self.bondtypes_listbox
        ==> self.changebond_listbox
        """
        self.unique_bonds = copy.deepcopy(self.q_bonds)
        bond_nr_prm = dict()
        prm_nr = 0

        self.bondtypes_listbox.delete(0, END)
        self.changebond_listbox.delete(0, END)
        for q1 in sorted(self.unique_bonds.keys()):
            state_types = [[],[],[],[]]
            state_qbonds = [[],[],[],[]]
            for state in range(self.evb_states.get()):
                for q2 in self.unique_bonds[q1][state]:
                    #Get atomtypes:
                    t1 = self.q_atomtypes[q1][state]
                    t2 = self.q_atomtypes[q2][state]
                    if 'du' in t1:
                        t1 = self.find_non_dummy(self.q_atomtypes[q1], t1)
                    if 'du' in t2:
                        t2 = self.find_non_dummy(self.q_atomtypes[q2], t2)
                    bond = '%s %s' % (t1, t2)
                    bond_rev = '%s %s' % (t2, t1)
                    state_types[state].append(bond)
                    state_qbonds[state].append('%d %d' % (q1, q2))
                    #Append bondtype with number to bond_prm_nr
                    if bond not in bond_nr_prm.values() and bond_rev not in bond_nr_prm.values():
                        prm_nr += 1
                        bond_nr_prm[prm_nr] = bond

                    #Delete q1 from q2 keys to avoid redundancies
                    try:
                        del self.unique_bonds[q2][state][self.unique_bonds[q2][state].index(q1)]
                    except:
                        continue

            #Insert bond types for each state ==> self.changebond_listbox
            s = 0
            max_ = self.evb_states.get()
            while s < max_:
                if s == 0:
                    insert_bonds = True
                else:
                    insert_bonds = False

                for i in range(len(state_qbonds[s])):
                    type_nr = [' ', ' ', ' ', ' ']
                    qpair = state_qbonds[s][i]
                    qtype = state_types[s][i]

                    qtype_rev = '%s %s' % (qtype.split()[1], qtype.split()[0])
                    #Find prm for bond type in bond_nr_prm:
                    for nr in bond_nr_prm:
                        if bond_nr_prm[nr] == qtype:
                            type_nr[s] = nr
                        elif bond_nr_prm[nr] == qtype_rev:
                            type_nr[s] = nr

                    #if state > 1, check if bond does not exist in previous state
                    if s > 0:
                        insert_bonds = True
                        for prev_s in range(0, s):
                            if qpair not in state_qbonds[prev_s]:
                                type_nr[prev_s] = 0
                            if qpair in state_qbonds[prev_s]:
                                insert_bonds = False

                    #Check if bond exist in next state(s)
                    next_s = s + 1
                    while next_s < max_:
                        if qpair in state_qbonds[next_s]:
                            type_next_s = state_types[next_s][state_qbonds[next_s].index(qpair)]
                            type_next_s_rev = '%s %s' % (type_next_s.split()[1], type_next_s.split()[0])
                            #Find prm for bond type in bond_nr_prm:
                            for nr_next in bond_nr_prm:
                                if bond_nr_prm[nr_next] == type_next_s:
                                    type_nr[next_s] = nr_next
                                elif bond_nr_prm[nr_next] == type_next_s_rev:
                                    type_nr[next_s] = nr_next

                        else:
                            type_nr[next_s] = 0
                        next_s += 1

                    if insert_bonds:
                        bt1, bt2, bt3, bt4 = type_nr[0:]
                        q1_, q2_ = qpair.split()
                        a1_, a2_ = self.q_atom_nr[int(q1_)], self.q_atom_nr[int(q2_)]
                        self.changebond_listbox.insert(END, '%6s %6s  %2s  %2s  %2s  %2s' %
                                                            (a1_, a2_, bt1, bt2, bt3, bt4))

                s += 1

        #Insert relevant bond parameter in ==> self.bondtypes_listbox
        for prm in sorted(bond_nr_prm.keys()):
            bond = bond_nr_prm[prm]
            bond_rev = '%s %s' % (bond.split()[1], bond.split()[0])
            if bond in self.bond_prm.keys():
                de, alpha, rb = self.bond_prm[bond][0:3]
                self.bondtypes_listbox.insert(END, '%2d %6.2f %5.2f %5.2f !%s' % (prm, de, alpha, rb, bond))
            elif bond_rev in self.bond_prm.keys():
                de, alpha, rb = self.bond_prm[bond_rev][0:3]
                self.bondtypes_listbox.insert(END, '%2d %6.2f %5.2f %5.2f !%s' % (prm, de, alpha, rb, bond))
            else:
                de, alpha, rb = '??', '??','??'
                self.bondtypes_listbox.insert(END, '%2d %6s %5s %5s !%s' % (prm, de, alpha, rb, bond))

    def update_fep_atoms(self):
        """
        Takes self.change_bonds and updates current states to listbox
        ==> self.fepatoms_listbox
        """
        self.fepatoms_listbox.delete(0, END)
        if len(self.fep_atoms) == 0:
            return

        evb_states = self.evb_states.get()
        for atom in self.fep_atoms.keys():

            s1, s2, s3, s4 = self.fep_atoms[atom]
            atom = int(atom)

            #Find Q-atom number corresponding to atom1 and atom2
            q1 = None
            for q in self.q_atom_nr.keys():
                if self.q_atom_nr[q] == atom:
                    q1 = q

            if evb_states < 3:
                s3, s4 = ' ', ' '
            if evb_states < 4:
                s4 = ' '

            if not q1:
                del self.fep_atoms[atom]
            else:
                self.fepatoms_listbox.insert(END, '%2d  %3s %3s  %3s %3s' %
                                                (q1, str(s1), str(s2), str(s3), str(s4)))

    def update_angles(self):
        """
        Uses self.q_bonds to generate all angle types.
        If angle types are not present in self.bond_prm, a search through parameter files will be performed.
        Updates angle listboxes
        ==> self.unique_bonds
        ==> self.angle_prm
        ==> self.angletypes_listbox
        ==> self.changeangle_listbox
        """
        self.changeangle_listbox.delete(0, END)
        self.angletypes_listbox.delete(0, END)
        change_q_bonds = dict()
        angle_nr_prm = dict()
        prm_nr = 0

        #List with Q-atoms to include for angle terms
        q_atoms = self.q_atom_nr.keys()

        if self.show_changing_angles:
            insert_angle = False
            q_atoms[:] = []
            #Collect q atoms that are involved in bond forming/breaking
            if len(self.change_bonds) > 0:
                for atompair in self.change_bonds.keys():
                    atom1, atom2 = map(int, atompair.split())

                    #Find corresponding Q-atom:
                    for q in self.q_atom_nr.keys():
                        if int(self.q_atom_nr[q]) == atom1 or int(self.q_atom_nr[q]) == atom2:
                            if q not in q_atoms:
                                q_atoms.append(q)

            #collect Q atoms that change atomtype
            for q in self.q_atomtypes.keys():
                atomtypes = []
                for s in range(self.evb_states.get()):
                    t = self.q_atomtypes[q][s]
                    if t not in atomtypes:
                        #Do not include dummy atoms!
                        if not t.startswith('du'):
                            atomtypes.append(t)
                if len(atomtypes) > 1:
                    q_atoms.append(q)


        elif not self.show_changing_angles:
            return

        #Collect all angle types:
        q_angles = [[], [], [], []]
        for q1 in sorted(self.q_bonds.keys()):
            insert_angle = False
            state_types = [[],[],[],[]]
            state_qangles = [[],[],[],[]]
            for state in range(self.evb_states.get()):
                for q2 in self.q_bonds[q1][state]:
                    if q2 != q1:
                        for q3 in self.q_bonds[q2][state]:
                            if q3 != q2 and q3 != q1:
                                if [q3, q2, q1] not in q_angles[state] or [q1, q2, q3] not in q_angles[state]:
                                    q_angles[state].append([q3, q2, q1])
                                    q_angles[state].append([q1, q2, q3])
                                    if (q1 in q_atoms) or (q2 in q_atoms) or (q3 in q_atoms):
                                        insert_angle = True
                                        t1 = self.q_atomtypes[q1][state]
                                        t2 = self.q_atomtypes[q2][state]
                                        t3 = self.q_atomtypes[q3][state]

                                        angle = '%s %s %s' % (t1, t2, t3)
                                        angle_rev = '%s %s %s' % (t3, t2, t1)

                                        state_types[state].append(angle)
                                        state_qangles[state].append('%d %d %d' % (q1, q2, q3))

                                        #Append angletype with number to bond_prm_nr
                                        if angle not in angle_nr_prm.values() and \
                                                        angle_rev not in angle_nr_prm.values():
                                            prm_nr += 1
                                            angle_nr_prm[prm_nr] = angle
            if insert_angle:
                #Insert angle types for each state ==> self.changeangle_listbox
                s = 0
                max_ = self.evb_states.get()
                while s < max_:
                    if s == 0:
                        insert_angle = True
                    else:
                        insert_angle = False

                    for i in range(len(state_qangles[s])):
                        type_nr = [' ', ' ', ' ', ' ']
                        qpair = state_qangles[s][i]
                        qtype = state_types[s][i]
                        qtype_rev = '%s %s %s' % (qtype.split()[2], qtype.split()[1], qtype.split()[0])

                        #Find prm for angle type in angle_nr_prm:
                        for nr in angle_nr_prm:
                            if angle_nr_prm[nr] == qtype:
                                type_nr[s] = nr
                            elif angle_nr_prm[nr] == qtype_rev:
                                type_nr[s] = nr

                        #if state > 1, check if bond does not exist in previous state
                        if s > 0:
                            insert_angle = True
                            for prev_s in range(0, s):
                                if qpair not in state_qangles[prev_s]:
                                    type_nr[prev_s] = 0
                                if qpair in state_qangles[prev_s]:
                                    insert_angle = False

                        #Check if bond exist in next state(s)
                        next_s = s + 1
                        while next_s < max_:
                            if qpair in state_qangles[next_s]:
                                type_next_s = state_types[next_s][state_qangles[next_s].index(qpair)]
                                type_next_s_rev = '%s %s %s' % (type_next_s.split()[2], type_next_s.split()[1],
                                                               type_next_s.split()[0])
                                #Find prm for bond type in bond_nr_prm:
                                for nr_next in angle_nr_prm:
                                    if angle_nr_prm[nr_next] == type_next_s:
                                        type_nr[next_s] = nr_next
                                    elif angle_nr_prm[nr_next] == type_next_s_rev:
                                        type_nr[next_s] = nr_next

                            else:
                                type_nr[next_s] = 0

                            next_s += 1

                        if insert_angle:
                            ang1, ang2, ang3, ang4 = type_nr[0:]
                            q1_, q2_, q3_ = map(int, qpair.split())
                            atom1, atom2, atom3 = map(int, [self.q_atom_nr[q1_], self.q_atom_nr[q2_], self.q_atom_nr[q3_]])
                            self.changeangle_listbox.insert(END, '%6d %6d %6d %2s  %2s %2s  %2s' %
                                                                 (atom1, atom2, atom3, ang1, ang2, ang3, ang4))

                    s += 1

        #Check if any parameters are missing, and collect missing from parameter file(s) if possible:
        prms = angle_nr_prm.values()
        missing_prms = []
        for prm in prms:
            prm_rev = '%s %s %s' % (prm.split()[2], prm.split()[1], prm.split()[0] )
            if prm not in self.angle_prm.keys() and prm_rev not in self.angle_prm.keys():
                missing_prms.append(prm)

        if len(missing_prms) > 0:
            for parameterfile in self.prms:
                found_angles = False
                with open(parameterfile, 'r') as prm:
                    for line in prm:
                        if found_angles:
                            if '[' in line and ']' in line:
                                found_angles = False
                            try:
                                type1 = line.split()[0]
                                type2 = line.split()[1]
                                type3 = line.split()[2]
                                angle = '%s %s %s' % (type1, type2, type3)
                                angle_rev = '%s %s %s' % (type3, type2, type1)
                                found_prm = False
                                if angle in missing_prms:
                                    found_prm = True
                                    del missing_prms[missing_prms.index(angle)]
                                elif angle_rev in missing_prms:
                                    found_prm = True
                                    del missing_prms[missing_prms.index(angle_rev)]
                                if found_prm:
                                    #If parameters already exist, do not overwrite (could be edited by user)
                                    k = float(line.split()[3])
                                    theta = float(line.split()[4])
                                    self.angle_prm[angle] = [k, theta]
                                    self.angle_prm[angle_rev] = [k, theta]

                            except:
                                continue
                        if '[angles]' in line:
                            found_angles = True

        #If parameters still are missing, the do not exist in parameters, append ?? to values
        if len(missing_prms) > 0:
            for angle in missing_prms:
                angle_rev = '%s %s %s' % (angle.split()[2], angle.split()[1], angle.split()[0])
                del missing_prms[missing_prms.index(angle)]
                if angle not in self.angle_prm.keys():
                    self.angle_prm[angle] = ['??', '??']
                if angle_rev not in self.angle_prm.keys():
                    self.angle_prm[angle_rev] = ['??', '??']

        #Insert angle types to listbox
        for prm in sorted(angle_nr_prm.keys()):
            angle = angle_nr_prm[prm]
            angle_rev = '%s %s %s' % (angle.split()[2], angle.split()[1], angle.split()[0])
            if angle in self.angle_prm.keys():
                k, theta = self.angle_prm[angle][0:]
            elif angle_rev in self.angle_prm:
                k, theta = self.angle_prm[angle_rev][0:]
            else:
                k, theta = '??', '??'
            self.angletypes_listbox.insert(END, '%2s %6s %6s !%s' % (prm, k, theta, angle))

    def update_torsions(self):
        """
        Uses self.q_bonds to generate all torsion types.
        If angle torsion types are not present in self.torsion_prm, a search through parameter files will be performed.
        Updates torsion listboxes
        ==> self.q_bonds
        ==> self.torsion_prm
        ==> self.torsiontypes_listbox
        ==> self.changetorsion_listbox
        """
        self.changetorsion_listbox.delete(0, END)
        self.torsiontypes_listbox.delete(0, END)
        torsion_nr_prm = dict()
        prm_nr = 1

        #List with Q-atoms to include for angle terms
        q_atoms = self.q_atom_nr.keys()

        if self.show_changing_torsion:
            insert_torsion = False
            q_atoms[:] = []
            #Collect q atoms that are involved in bond forming/breaking
            if len(self.change_bonds) > 0:
                for atompair in self.change_bonds:
                    atom1, atom2 = map(int, atompair.split())

                    #Find corresponding Q-atom:
                    for q in self.q_atom_nr.keys():
                        if int(self.q_atom_nr[q]) == atom1 or int(self.q_atom_nr[q]) == atom2:
                            if q not in q_atoms:
                                q_atoms.append(q)

            #collect Q atoms that change atomtype
            for q in self.q_atomtypes.keys():
                atomtypes = []
                for s in range(self.evb_states.get()):
                    t = self.q_atomtypes[q][s]
                    if t not in atomtypes:
                        #Do not include dummy atoms!
                        if not t.startswith('du'):
                            atomtypes.append(t)
                if len(atomtypes) > 1:
                    q_atoms.append(q)
        elif not self.show_changing_torsion:
            return

        #Collect all torsion types:
        q_torsions = [[], [], [], []]
        for q1 in sorted(self.q_bonds.keys()):
            insert_torsion = False
            state_types = [[],[],[],[]]
            state_qtorsions = [[],[],[],[]]
            for state in range(self.evb_states.get()):
                for q2 in self.q_bonds[q1][state]:
                    if q2 != q1:
                        for q3 in self.q_bonds[q2][state]:
                            if q3 != q2 and q3 != q1:
                                for q4 in self.q_bonds[q3][state]:
                                    if q4 != q3 and q4 != q2 and q4 != q1:
                                        if [q4, q3, q2, q1] not in q_torsions[state] or \
                                                [q1, q2, q3, q4] not in q_torsions[state]:
                                            q_torsions[state].append([q4, q3, q2, q1])
                                            q_torsions[state].append([q1, q2, q3, q4])
                                            if (q1 in q_atoms) or (q2 in q_atoms) or (q3 in q_atoms) or (q4 in q_atoms):
                                                insert_torsion = True
                                                t1 = self.q_atomtypes[q1][state]
                                                t2 = self.q_atomtypes[q2][state]
                                                t3 = self.q_atomtypes[q3][state]
                                                t4 = self.q_atomtypes[q4][state]

                                                torsion = '%s %s %s %s' % (t1, t2, t3, t4)
                                                torsion_rev = '%s %s %s %s' % (t4, t3, t2, t1)

                                                state_types[state].append(torsion)
                                                state_qtorsions[state].append('%d %d %d %d' % (q1, q2, q3, q4))

                                                #Append angletype with number to bond_prm_nr
                                                if torsion not in torsion_nr_prm.values() and \
                                                        torsion_rev not in torsion_nr_prm.values():
                                                    torsion_nr_prm[prm_nr] = torsion
                                                    prm_nr += 3
            if insert_torsion:
                #Insert torsion types for each state ==> self.changtorsion_listbox
                s = 0
                max_ = self.evb_states.get()
                while s < max_:
                    if s == 0:
                        insert_torsion = True
                    else:
                        insert_torsion = False

                    for i in range(len(state_qtorsions[s])):
                        type_nr = [' ', ' ', ' ', ' ']
                        qpair = state_qtorsions[s][i]
                        qtype = state_types[s][i]
                        qtype_rev = '%s %s %s %s' % \
                                    (qtype.split()[3], qtype.split()[2], qtype.split()[1], qtype.split()[0])

                        #Find prm for angle type in angle_nr_prm:
                        for nr in torsion_nr_prm:
                            if torsion_nr_prm[nr] == qtype:
                                type_nr[s] = nr
                            elif torsion_nr_prm[nr] == qtype_rev:
                                type_nr[s] = nr

                        #if state > 1, check if bond does not exist in previous state
                        if s > 0:
                            insert_torsion = True
                            for prev_s in range(0, s):
                                if qpair not in state_qtorsions[prev_s]:
                                    type_nr[prev_s] = 0
                                if qpair in state_qtorsions[prev_s]:
                                    insert_torsion = False

                        #Check if bond exist in next state(s)
                        next_s = s + 1
                        while next_s < max_:
                            if qpair in state_qtorsions[next_s]:
                                type_next_s = state_types[next_s][state_qtorsions[next_s].index(qpair)]
                                type_next_s_rev = '%s %s %s %s' % (type_next_s.split()[3], type_next_s.split()[2],
                                                                   type_next_s.split()[1], type_next_s.split()[0])
                                #Find prm for bond type in bond_nr_prm:
                                for nr_next in torsion_nr_prm:
                                    if torsion_nr_prm[nr_next] == type_next_s:
                                        type_nr[next_s] = nr_next
                                    elif torsion_nr_prm[nr_next] == type_next_s_rev:
                                        type_nr[next_s] = nr_next

                            else:
                                type_nr[next_s] = 0

                            next_s += 1

                        if insert_torsion:
                            tor1, tor2, tor3, tor4 = type_nr[0:]
                            q1_, q2_, q3_, q4_ = map(int, qpair.split())
                            atom1, atom2, atom3, atom4 = map(int, [self.q_atom_nr[q1_], self.q_atom_nr[q2_],
                                                                   self.q_atom_nr[q3_],  self.q_atom_nr[q4_]])
                            for z in range(3):
                                self.changetorsion_listbox.insert(END, '%6d %6d %6d %6d %3s  %3s  %3s  %3s' %
                                                                       (atom1, atom2, atom3, atom4,
                                                                        str(tor1), str(tor2), str(tor3), str(tor4)))
                                if tor1 != 0:
                                    tor1 += 1
                                if tor2 != 0:
                                    tor2 += 1
                                try:
                                    if tor3 != 0:
                                        tor3 += 1
                                except:
                                    continue
                                try:
                                    if tor4 != 0:
                                        tor4 += 1
                                except:
                                    continue

                    s += 1

        #Check if any parameters are missing, and collect missing from parameter file(s) if possible:
        prms = torsion_nr_prm.values()
        missing_prms = []
        for prm in prms:
            prm_rev = '%s %s %s %s' % (prm.split()[3], prm.split()[2], prm.split()[1], prm.split()[0])
            if prm not in self.torsion_prm.keys() and prm_rev not in self.torsion_prm.keys():
                missing_prms.append(prm)

        if len(missing_prms) > 0:
            for parameterfile in self.prms:
                found_torsions = False
                with open(parameterfile, 'r') as prm:
                    for line in prm:
                        if found_torsions:
                            if '[' in line and ']' in line:
                                found_torsions = False
                            try:
                                type1 = line.split()[0]
                                type2 = line.split()[1]
                                type3 = line.split()[2]
                                type4 = line.split()[3]

                                torsion = '%s %s %s %s' % (type1, type2, type3, type4)
                                torsion_rev = '%s %s %s %s' % (type4, type3, type2, type1)
                                found_prm = False

                                if torsion in missing_prms:
                                    found_prm = True
                                elif torsion_rev in missing_prms:
                                    found_prm = True

                                if found_prm:
                                    kt = float(line.split()[4])
                                    minima = int(float(line.split()[5]))
                                    phase = float(line.split()[6])
                                    paths = int(float(line.split()[7]))

                                    if torsion not in self.torsion_prm.keys():
                                        self.torsion_prm[torsion] = [[0,0.0,1],[0, 180.0,1],[0, 0.0,1]]
                                    if torsion_rev not in self.torsion_prm.keys():
                                        self.torsion_prm[torsion_rev] = [[0,0.0,1],[0, 180.0,1],[0, 0.0,1]]
                                    #Make sure to note overwrite existing parameters (check that kt > 0)
                                    if self.torsion_prm[torsion][abs(minima) - 1][0] == 0:
                                        self.torsion_prm[torsion][abs(minima) - 1] = [kt, phase, paths]
                                        self.torsion_prm[torsion_rev][abs(minima) - 1] = [kt, phase, paths]


                                    if minima > 0:
                                        del missing_prms[missing_prms.index(torsion)]

                            except:
                                continue
                        if '[torsions]' in line:
                            found_torsions = True

        #If parameters still are missing, the do not exist in parameters, append ?? to values
        if len(missing_prms) > 0:
            for torsion in missing_prms:
                torsion_rev = '%s %s %s %s' % (torsion.split()[2], torsion.split()[2], torsion.split()[1],
                                               torsion.split()[0])
                del missing_prms[missing_prms.index(torsion)]
                if torsion not in self.torsion_prm.keys() and torsion_rev not in self.torsion_prm.keys():
                    self.torsion_prm[torsion] = [['??', 0.00, 1], ['??', 180.00, 1], ['??', 0.00, 1]]
                    self.torsion_prm[torsion_rev] = [['??', 0.00, 1], ['??', 180.00, 1], ['??', 0.00, 1]]


        #Insert torsion types to listbox
        for prm in sorted(torsion_nr_prm.keys()):
            torsion = torsion_nr_prm[prm]
            torsion_rev = '%s %s %s %s' % (torsion.split()[3], torsion.split()[2], torsion.split()[1],
                                           torsion.split()[0])
            for j in range(3):
                list_nr = prm + j
                insert_tor = torsion
                minima = str(-(j + 1))
                if int(minima) == -3:
                    minima = str(3)

                if torsion in self.torsion_prm.keys():
                    k, phase, path = map(str, self.torsion_prm[torsion][j][0:])
                elif torsion_rev in self.torsion_prm.keys():
                    k, phase, path = map(str, self.torsion_prm[torsion_rev][j][0:])
                    insert_tor = torsion_rev
                else:
                    k, phase, path = '??', '0.00', '1'
                    if j == 1:
                        phase = '180.00'

                self.torsiontypes_listbox.insert(END, '%3d %7s %6s %6s %6s !%s' %
                                                      (list_nr, k, minima, phase, path, insert_tor))

    def update_impropers(self):
        """
        Uses information from self.q_impropers (from the library files) to generate all improper types.
        self.q_bonds is checked against each state to verify if bond exist.

        ==> self.q_impropers  {q2 : [q1,q3,q4]} --> search for q1 q2 q3 q4
        ==> self.q_bonds
        ==> self.improper_prm
        ==> self.impropertypes_listbox
        ==> self.changeimproper_listbox
        """
        self.changeimproper_listbox.delete(0, END)
        self.impropertypes_listbox.delete(0, END)

        #List with Q-atoms to include for angle terms
        q_atoms = self.q_atom_nr.keys()

        if self.show_changing_impropers:
            q_atoms[:] = []
            #Collect q atoms that are involved in bond forming/breaking
            if len(self.change_bonds) > 0:
                for atompair in self.change_bonds.keys():
                    atom1, atom2 = map(int, atompair.split())

                    #Find corresponding Q-atom:
                    for q in self.q_atom_nr.keys():
                        if int(self.q_atom_nr[q]) == atom1 or int(self.q_atom_nr[q]) == atom2:
                            if q not in q_atoms:
                                q_atoms.append(q)

            #collect Q atoms that change atomtype
            for q in self.q_atomtypes.keys():
                atomtypes = []
                for s in range(self.evb_states.get()):
                    t = self.q_atomtypes[q][s]
                    if t not in atomtypes:
                        #Do not include dummy atoms!
                        if not t.startswith('du'):
                            atomtypes.append(t)
                if len(atomtypes) > 1:
                    q_atoms.append(q)

        elif not self.show_changing_impropers:
            return

        #Write change impropers lisbox
        prm_nr = 0
        nr_type = dict()
        for q2 in self.q_impropers.keys():
            q1, q3, q4 = self.q_impropers[q2][0:]
            if q1 in q_atoms or q2 in q_atoms or q3 in q_atoms or q4 in q_atoms:
                states = [' ',' ',' ',' ']
                for state in range(self.evb_states.get()):
                    bonds_q2 = self.q_bonds[q2][state]
                    if q1 not in bonds_q2 or q3 not in bonds_q2 or q4 not in bonds_q2:
                        states[state] = 0
                    #If true improper, q2 is bonded to exactly 3 atoms!
                    elif len(self.q_bonds[q2][state]) != 3:
                        states[state] = 0
                    #q3 and q4 can not both be sp3!
                    elif len(self.q_bonds[q3][state]) > 3 and len(self.q_bonds[q4][state]) > 3:
                        states[state] = 0
                    else:
                        t1 = self.q_atomtypes[q1][state]
                        t2 = self.q_atomtypes[q2][state]
                        t3 = self.q_atomtypes[q3][state]
                        t4 = self.q_atomtypes[q4][state]
                        imp_type = '%s %s %s %s' % (t1, t2, t3, t4)
                        imp_type_rev = '%s %s %s %s' % (t4, t3, t2, t1)
                        if imp_type not in nr_type.values() and imp_type_rev not in nr_type.values():
                            prm_nr += 1
                            nr_type[prm_nr] = imp_type
                            states[state] = str(prm_nr)
                        elif imp_type in nr_type.values() or imp_type_rev in nr_type.values():
                            for nr in nr_type.keys():
                                if nr_type[nr] == imp_type or nr_type[nr] == imp_type_rev:
                                    states[state] = str(nr)

                try:
                    s1, s2, s3, s4 = states[0:]
                    a1 = self.q_atom_nr[q1]
                    a2 = self.q_atom_nr[q2]
                    a3 = self.q_atom_nr[q3]
                    a4 = self.q_atom_nr[q4]

                    self.changeimproper_listbox.insert(END, '%6d %6d %6d %6d  %3s %3s %3s %3s' %
                                                   (a1, a2, a3, a4, s1, s2, s3, s4))
                except:
                    continue
        #Check if any parameters are missing, and collect missing from parameter file(s) if possible:
        prms = nr_type.values()
        missing_prms = []
        for prm in prms:
            prm_rev = '%s %s %s %s' % (prm.split()[3], prm.split()[2], prm.split()[1], prm.split()[0])
            if prm not in self.improper_prm.keys() and prm_rev not in self.improper_prm.keys():
                missing_prms.append(prm)

        if len(missing_prms) > 0:
            for parameterfile in self.prms:
                found_impropers = False
                with open(parameterfile, 'r') as prm:
                    for line in prm:
                        if found_impropers:
                            if '[' in line and ']' in line:
                                found_impropers = False
                            try:
                                type1 = line.split()[0]
                                type2 = line.split()[1]
                                type3 = line.split()[2]
                                type4 = line.split()[3]

                                improper = '%s %s %s %s' % (type1, type2, type3, type4)
                                improper_rev = '%s %s %s %s' % (type4, type3, type2, type1)
                                found_prm = False

                                if improper in missing_prms:
                                    found_prm = True
                                elif improper_rev in missing_prms:
                                    found_prm = True
                                    improper = improper_rev

                                if found_prm:
                                    kt = float(line.split()[4])
                                    phase = float(line.split()[5])

                                    self.improper_prm[improper] = [kt, phase]
                                    del missing_prms[missing_prms.index(improper)]
                            except:
                                continue
                        if '[impropers]' in line:
                            found_impropers = True

        if len(missing_prms) > 0:
            for imp in missing_prms:
                self.improper_prm[imp] = ['??', '??']

        #Insert imroper types to listbox
        for nr in sorted(nr_type.keys()):
            type_ = nr_type[nr]
            type_rev = '%s %s %s %s' % (type_.split()[3], type_.split()[2], type_.split()[1], type_.split()[0])
            if type_ in self.improper_prm.keys():
                kt, phase = map(str, self.improper_prm[type_][0:])
            elif type_rev in self.improper_prm.keys():
                kt, phase = map(str, self.improper_prm[type_rev][0:])

            self.impropertypes_listbox.insert(END, '%3d %6s %6s !%s' % (nr, kt, phase, type_))

        print self.improper_prm

    def update_soft_pairs(self):
        """
        All bonds broken or formed are automatically added..
        """
        self.softpairs_listbox.delete(0, END)
        #Go through self.change_bonds and set softpairs for bonds broken/formed
        for atompair in self.change_bonds.keys():
            q1 = None
            q2 = None

            states = []
            atom1, atom2 = map(int, atompair.split()[0:])

            #Check that bondings are changed between states
            for state in range(self.evb_states.get()):
                s = self.change_bonds[atompair][state]
                if s not in states:
                    states.append(s)
            if len(states) > 1:
                #Bond is broken/formed. Find Q nr for atom1 and atom2
                for q in self.q_atom_nr.keys():
                    if int(self.q_atom_nr[q]) == atom1:
                        q1 = q
                    if int(self.q_atom_nr[q]) == atom2:
                        q2 = q

                if q1 and q2:
                    note = ''
                    for state in range(self.evb_states.get()):
                        qt1 = self.q_atomtypes[q1][state]
                        qt2 = self.q_atomtypes[q2][state]
                        if int(self.change_bonds[atompair][state]) == 0:
                            note += '%s + %s' % (qt1, qt2)
                        else:
                            note += '%s-%s' % (qt1, qt2)
                        if state != self.evb_states.get() - 1:
                                note += ' --> '
                    self.softpairs_listbox.insert(END, '%3d %3d  !%s' % (q1, q2, note))

        #Insert additional soft-pairs (if any)
        if len(self.softpairs_added) > 0:
            for q1 in self.softpairs_added.keys():
                q2 = self.softpairs_added[q1]
                bond_change = False
                bonded = False

                #find correct atompair:
                atom1 = self.q_atom_nr[q1]
                atom2 = self.q_atom_nr[q2]
                atompair = '%d %d' % (atom1, atom2)
                atompair_rev = '%d %d' % (atom2, atom1)

                if atompair in self.change_bonds:
                    bond_change = True
                elif atompair_rev in self.change_bonds:
                    atompair = atompair_rev
                    bond_change = True
                else:
                    if q2 in self.q_bonds[q1]:
                        bonded = True
                note = ''
                for state in range(self.evb_states.get()):
                    qt1 = self.q_atomtypes[q1][state]
                    qt2 = self.q_atomtypes[q2][state]
                    if bond_change:
                        if int(self.change_bonds[atompair][state]) == 0:
                            note += '%s + %s' % (qt1, qt2)
                        else:
                            note += '%s-%s' % (qt1, qt2)
                    if not bond_change:
                        if bonded:
                            note += '%s-%s' % (qt1, qt2)
                        else:
                            note += '%s + %s' % (qt1, qt2)

                    if state != self.evb_states.get() - 1:
                            note += ' --> '
                self.softpairs_listbox.insert(END, '%3d %3d  !%s' % (q1, q2, note))

    def update_couplings(self):
        """
        Angles, torsions and impropers dependent on bonds formed/broken are automatically added
        """
        if len(self.change_bonds) == 0:
            return

        self.angle_couplings_listbox.delete(0, END)
        self.torsion_couplings_listbox.delete(0, END)
        self.improper_couplings_listbox.delete(0, END)

        #Get bonds broken/formed and enumerate
        bonds_to_couple = []
        nr_bonds = []
        nr = 0
        bond_list = self.changebond_listbox.get(0, END)
        for line in bond_list:
            nr += 1
            a1, a2 = map(int, line.split()[0:2])
            if 0 in map(int, line.split()[2:]):
                q1, q2 = None, None
                for q in self.q_atom_nr.keys():
                    if int(self.q_atom_nr[q]) == a1:
                        q1 = q
                    if int(self.q_atom_nr[q]) == a2:
                        q2 = q
                if q1 and q2:
                    bonds_to_couple.append([q1,q2])
                    nr_bonds.append(nr)

        #Enumerate angles and insert to listbox
        nr = 0
        angle_list = self.changeangle_listbox.get(0, END)
        for line in angle_list:
            nr += 1
            atom1, atom2, atom3 = map(int, line.split()[0:3])
            q1, q2, q3 = 0, 0, 0
            for q in self.q_atom_nr.keys():
                if self.q_atom_nr[q] == atom1:
                    q1 = q
                if self.q_atom_nr[q] == atom2:
                    q2 = q
                if self.q_atom_nr[q] == atom3:
                    q3 = q

            for i in range(len(bonds_to_couple)):
                bond = bonds_to_couple[i]
                insert_coupling = False

                if q1 in bond:
                    if q2 in bond or q3 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                elif q2 in bond:
                    if q1 in bond or q3 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                if insert_coupling:
                    self.angle_couplings_listbox.insert(END, '%3d %3d' % (nr, bnd_nr))

        #Get torsion couplings
        nr = 0
        torsion_list = self.changetorsion_listbox.get(0, END)
        for line in torsion_list:
            nr += 1
            atom1, atom2, atom3, atom4 = map(int, line.split()[0:4])
            q1, q2, q3, q4 = 0, 0, 0, 0

            for q in self.q_atom_nr.keys():
                if self.q_atom_nr[q] == atom1:
                    q1 = q
                if self.q_atom_nr[q] == atom2:
                    q2 = q
                if self.q_atom_nr[q] == atom3:
                    q3 = q
                if self.q_atom_nr[q] == atom4:
                    q4 = q

            for i in range(len(bonds_to_couple)):
                bond = bonds_to_couple[i]
                insert_coupling = False

                if q1 in bond:
                    if q2 in bond or q3 in bond or q4 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                elif q2 in bond:
                    if q1 in bond or q3 in bond or q4 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                elif q3 in bond:
                    if q1 in bond or q2 in bond or q4 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                if insert_coupling:
                    self.torsion_couplings_listbox.insert(END, '%3d %3d' % (nr, bnd_nr))


        #Get improper couplings
        nr = 0
        improper_list = self.changeimproper_listbox.get(0, END)
        for line in improper_list:
            nr += 1
            atom1, atom2, atom3, atom4 = map(int, line.split()[0:4])
            q1, q2, q3, q4 = 0, 0, 0, 0

            for q in self.q_atom_nr.keys():
                if self.q_atom_nr[q] == atom1:
                    q1 = q
                if self.q_atom_nr[q] == atom2:
                    q2 = q
                if self.q_atom_nr[q] == atom3:
                    q3 = q
                if self.q_atom_nr[q] == atom4:
                    q4 = q

            for i in range(len(bonds_to_couple)):
                bond = bonds_to_couple[i]
                insert_coupling = False

                if q1 in bond:
                    if q2 in bond or q3 in bond or q4 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                elif q2 in bond:
                    if q1 in bond or q3 in bond or q4 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                elif q3 in bond:
                    if q1 in bond or q2 in bond or q4 in bond:
                        bnd_nr = nr_bonds[i]
                        insert_coupling = True
                if insert_coupling:
                    self.improper_couplings_listbox.insert(END, '%3d %3d' % (nr, bnd_nr))


    def update_excluded_pairs(self):
        """
        Updates excluded pairs listbox
        """
        self.excludedpairs_listbox.delete(0, END)

        for atompair in self.excluded_pairs.keys():
            states = ' '
            for state in range(self.evb_states.get()):
                states += '%3d ' % self.excluded_pairs[atompair][state]
            self.excludedpairs_listbox.insert(END, '%s %s' % (atompair, states))

    def update_lambdasteps(self):
        self.lambdasteps_listbox.delete(0, END)

        for i in sorted(self.step_lambdavalues.keys()):
            l1, l2, l3, l4 = map(float, self.step_lambdavalues[i][0:])
            l1, l2, l3, l4 = '%03.3f' % l1, '%03.3f' % l2, '%03.3f' % l3, '%03.3f' % l4
            if self.evb_states.get() < 3:
                l3 = ' '
                l4 = ' '
            if self.evb_states.get() < 4:
                l4 = ' '
            self.lambdasteps_listbox.insert(END, '%4d %s %s %s %s' % (i, l1, l2, l3, l4))

        try:
            lambda_steps = max(self.step_lambdavalues.keys())
        except:
            lambda_steps = 0

        self.qstatus[u"\N{GREEK small LETTER LAMDA}-steps/run"] = '%s' % str(lambda_steps)
        status_list = self.qstatus_listbox.get(0, END)

        for i in range(len(status_list)):
            if u"\N{GREEK small LETTER LAMDA}-steps/run" in status_list[i]:
                self.qstatus_listbox.delete(i)
                self.qstatus_listbox.insert(i, u"\N{GREEK small LETTER LAMDA}-steps/run: %s" % str(lambda_steps))

        #Update temperatures to get correct total simulation time:
        self.update_temperatures()

    def update_temperatures(self):
        self.temp_listbox.delete(0, END)
        runs = 0
        for temp in sorted(self.temp_runs.keys()):
            self.temp_listbox.insert(END, '%6s  %3d' % (temp, self.temp_runs[temp]))
            runs += self.temp_runs[temp]

        #Calculate total simulation time and update status
        try:
            lambda_steps = max(self.step_lambdavalues.keys())
        except:
            lambda_steps = 0

        simtime_step = self.md_settings['simtime']
        self.md_settings['inputfiles'] = lambda_steps * runs
        tot_simtime = round(float(lambda_steps) * float(runs) * float(simtime_step), 6)
        self.qstatus['Total simulation time (ns)'] = '%s' % str(tot_simtime)

        #Update Status listbox:
        status_list = self.qstatus_listbox.get(0, END)

        for i in range(len(status_list)):
            if 'Total simulation time (ns)' in status_list[i]:
                self.qstatus_listbox.delete(i)
                self.qstatus_listbox.insert(i, 'Total simulation time (ns): %.6f ' % tot_simtime)

    def update_all(self):
        """
        Updates all lists without reading libraries again. Takes care of Q-atom changes, types, bonds, charges etc.
        """

        self.get_atomtype_parameters()
        self.update_atom_parameters()

        self.update_charges()

        self.get_bond_parameters()
        self.update_q_bonds()

        self.update_fep_atoms()

        self.update_angles()

        self.update_torsions()

        self.update_impropers()

        self.update_soft_pairs()

        self.update_temperatures()

        self.update_status()

        self.update_change_atoms_listbox()

        if self.session:
            self.highlight_q_atoms()

    def update_status(self, print_all=False):
        """
        Update status list on EVB setup progress
        """
        self.fep_written = False
        self.files_written = False
        #Go through listboxes and check if parameters are missing
        missing_check = {'Charges' : self.charge_listbox, 'Atom types' : self.atomtypes_listbox,
                         'Bond parameters' : self.bondtypes_listbox, 'Angle parameters' : self.angletypes_listbox,
                         'Torsion parameters' : self.torsiontypes_listbox,
                         'Improper parameters' : self.impropertypes_listbox}

        for term in missing_check.keys():
            listbox_ = missing_check[term].get(0, END)
            listbox = []
            for i in listbox_:
                listbox.append(i)
            count = 0
            miss = False
            for line in listbox:
                if '??' in line:
                    count += 1
                    miss = True
            if miss:
                prm = 'parameters'
                if count == 1:
                    prm = 'parameter'
                self.qstatus[term] = 'WARNING! Missing %d %s!' % (count, prm)
            else:
                self.qstatus[term] = 'OK'

        #Check if Q-atoms and bond change is defined:
        exist_check = {'Q-atoms': self.qatoms_listbox, 'Form/Break bonds': self.fepatoms_listbox}
        for term in exist_check.keys():
            listbox_ = exist_check[term].get(0, END)
            listbox = []
            for i in listbox_:
                listbox.append(i)

            if len(listbox) == 0:
                self.qstatus[term] = 'WARNING! None selected.'
            else:
                self.qstatus[term] = 'OK'

        #Update statues listbox
        self.qstatus_listbox.delete(0, END)
        do_not_print = []
        if not print_all:
            do_not_print = ['OK','ok','na','NA','']

        for i in self.qstatus.keys():
            if self.qstatus[i] not in do_not_print:
                self.qstatus_listbox.insert(END, '%s: %s' % (i, self.qstatus[i]))

    def import_fepfile(self):
        """
        Import existing FEP file
        """
        pass
        self.update_status()

    def auto_setup(self):
        """
        Use Impact to get charges, bonds, angles, torsions and impropers for Q-atoms in state 1 and 2.
        """
        self.app.errorBox('Info', 'Sorry, not implemented yet.')
        self.update_status()

    def import_atom_prm(self):
        """
        opens up a window with list of all atom-parameters in parameter file.
        """
        self.import_atomtypes = ImportParameters(self, self.root, self.atomtypes_listbox)

    def import_bond_prm(self):
        self.import_bondprm = ImportParameters(self, self.root, self.bondtypes_listbox, '[bonds]')

    def import_angle_prm(self):
        self.import_angleprm = ImportParameters(self, self.root, self.angletypes_listbox, '[angles]')

    def import_torsion_prm(self):
        self.import_torsionprm = ImportParameters(self, self.root, self.angletypes_listbox, '[torsions]')

    def import_improper_prm(self):
        self.import_improperprm = ImportParameters(self, self.root, self.angletypes_listbox, '[impropers]')

    def set_atomtype_state(self, state):
        """
        Takes selected atomtype from list and appends it to Q-atom state 1,2,3 or 4
        """
        type_selected = map(int, self.atomtypes_listbox.curselection())
        if len(type_selected) != 1:
            self.app.errorBox('Error', 'Select exactly 1 atomtype to assign to Q-atoms.')
            return

        newtype = self.atomtypes_listbox.get(type_selected[0]).split()[0].strip()

        q_atoms_selections = map(int, self.changetypes_listbox.curselection())
        for selected in q_atoms_selections:
            qi = int(self.changetypes_listbox.get(selected).split()[0])
            self.q_atomtypes[qi][state] = newtype

        #If new atom parameter is not in atomtype:parameter list, append it:
        if newtype not in self.atomtype_prm.keys():
            ri, ei, ci, ai, ri1_4, ei_1_4, mass = map(float, self.atomtypes_listbox.get(type_selected[0]).split()[1:])
            self.atomtype_prm[newtype] = [ri, ei, ci, ai, ri1_4, ei_1_4, mass]

        #Go through parameter files and try to find bond parameters for new atomtyp
        self.get_bond_parameters()

        self.update_all()

        #Highlight selections again after updating lists:
        self.atomtypes_listbox.selection_set(type_selected[0])
        for selected in q_atoms_selections:
            self.changetypes_listbox.selection_set(selected)

    def edit_atom_prm(self):
        """
        Edit selected atom type
        """
        selections = map(int, self.atomtypes_listbox.curselection())
        if len(selections) != 1:
            return

        self.add_prm = EditParameters(self, self.root, self.atomtypes_listbox, selections[0], True)

    def add_atom_prm(self):
        """
        Manually add new atom parameter
        """
        self.add_prm = EditParameters(self, self.root, self.atomtypes_listbox, END, False)

    def delete_atom_prm(self):
        """
        Deletes selected atom parameter from list
        """
        selections = map(int, self.atomtypes_listbox.curselection())
        if len(selections) != 1:
            return

        del_type = self.atomtypes_listbox.get(selections[0]).split()[0]
        del self.atomtype_prm[del_type]

        #Check if any Q atoms have atomtype, and change:
        for atomtype in self.q_atomtypes.keys():
            for state in range(self.evb_states.get()):
                if del_type in self.q_atomtypes[atomtype][state]:
                    #Check if atomtype can be adapted from other states for Q-atom
                    found_replacement = False
                    for others in range(self.evb_states.get()):
                        if del_type not in self.q_atomtypes[atomtype][others]:
                            found_replacement = True
                            self.q_atomtypes[atomtype][state] = self.q_atomtypes[atomtype][others]
                    if not found_replacement:
                        self.q_atomtypes[atomtype][state] = '??'

        self.update_all()

    def edit_bond_prm(self):
        """
        Edit existing bond parameter
        """
        selections = map(int, self.bondtypes_listbox.curselection())
        if len(selections) != 1:
            return

        self.editbond = EditBondParameters(self, self.root, self.bondtypes_listbox, selections[0], True)

    def add_bond_prm(self):
        """
        Manually add new bond parameter
        """
        self.editbond = EditBondParameters(self, self.root, self.bondtypes_listbox, END, False)

    def add_fep_atom(self):
        """
        Add selected Q-atom pair to bond change list
        """
        sel_index = map(int, self.qatoms_listbox.curselection())

        if len(sel_index) < 1:
            return
        
        for selected in sel_index:
            qatom = int(self.qatoms_listbox.get(selected).split()[0])
            atom = self.q_atom_nr[qatom]
            if atom not in self.fep_atoms.keys():
                self.fep_atoms[atom] = [1, 1, 1, 1]

        self.qatoms_listbox.select_clear(0, END)
        self.update_pymol_bonds()
        self.update_all()

    def set_fep_state(self, state):
        """
        Turns bond for selected state on/off (1/0) and updated self.change_bonds
        ==> self.fepatoms_listbox
        """

        selections = map(int, self.fepatoms_listbox.curselection())
        if len(selections) == 0:
            return

        dummyDict = {1.01:'duH', 12.01:'duC', 14.01: 'duN', 16.00: 'duO', 19.00: 'duF', 24.305: 'duMg', 30.97: 'duP',
                     32.07: 'duS', 35.45: 'duCl', 40.08: 'duCa', 55.85: 'duFe', 65.38: 'duZn', 79.90:'duBr',
                     126.90: 'duI'}

        for selected in selections:
            atomchange = self.fepatoms_listbox.get(selected)
            q = int(atomchange.split()[0])
            atom = self.q_atom_nr[q]
            states = map(int, atomchange.split()[1:])

            #Change selected state:
            states[state] = abs(states[state] - 1)
            fep_atom = int(states[state])
            if fep_atom == 0:
                mass = round(self.atomtype_prm[self.q_atomtypes[q][state]][-1], 2)
                print mass
                dummyType = dummyDict[mass]
                if not dummyType in self.atomtype_prm.keys():
                    self.atomtype_prm[dummyType] = [0, 0, 0, 0, 0, 0, mass]
                self.q_atomtypes[q][state] = dummyType
                self.q_charges[q][state] = 0
            elif fep_atom == 1:
                for i in range(4):
                    if i != state:
                        if self.q_atomtypes[q][i] != self.q_atomtypes[q][state]:
                            self.q_atomtypes[q][state] = self.q_atomtypes[q][i]
                        if self.q_charges[q][i] != self.q_charges[q][state]:
                            self.q_charges[q][state] = self.q_charges[q][i]

            for s in range(len(states)):
                self.fep_atoms[atom][s] = states[s]


            self.update_pymol_bonds()

        self.update_all()


        #Highlight selections again after updating lists:
        for selected in selections:
            self.fepatoms_listbox.selection_set(selected)

    def del_fep_atom(self):
        """
        removes bond change from list
        ==> self.change_bonds
        """
        sel_index = map(int, self.fepatoms_listbox.curselection())

        if len(sel_index) == 0:
            return

        for selected in sel_index:
            bondchange = self.fepatoms_listbox.get(selected)
            q = int(bondchange.split()[0])
            atom = self.q_atom_nr[q]
            del self.fep_atoms[atom]
            self.fepatoms_listbox.delete(selected)

            #Needs special attention for pymol session
            if self.session:
                self.session.stdin.write('show sticks, id %s\n' % atom)

        #self.update_q_atoms()
        self.update_all()

    def edit_angle_prm(self):
        """
        Edit existing angle parameter
        """
        selections = map(int, self.angletypes_listbox.curselection())
        if len(selections) != 1:
            return
        self.editangle = EditAngleParameters(self, self.root, self.angletypes_listbox, selections[0], True)

    def add_angle_prm(self):
        """
        Manually add new angle parameter
        """
        self.editangle = EditAngleParameters(self, self.root, self.angletypes_listbox, END, False)

    def edit_torsion_prm(self):
        """
        edit existing torsion parameter
        """
        selections = map(int, self.torsiontypes_listbox.curselection())
        if len(selections) != 1:
            return
        self.edittorsion = EditTorsionParameters(self, self.root, self.torsiontypes_listbox, selections[0], True)

    def add_torsion_prm(self):
        """
        Manually add new torsion parameter
        """
        self.addtorsion = EditTorsionParameters(self, self.root, self.torsiontypes_listbox, END, False)

    def edit_improper_prm(self):
        """
        edit existing improper parameter
        """
        selections = map(int, self.impropertypes_listbox.curselection())
        if len(selections) != 1:
            return
        self.editimproper = EditImproperParameters(self, self.root, self.impropertypes_listbox, selections[0], True)

    def add_improper_prm(self):
        """
        Manually add new improper parameter
        """
        self.addimproper = EditImproperParameters(self, self.root, self.impropertypes_listbox, END, False)

    def add_improper(self):
        """
        Add selection of 4 Q-atoms to improper list
        """
        try:
            selections = map(int, self.qatoms_listbox.curselection())
        except:
            return

        if len(selections) != 4:
            self.app.errorBox('Warning', 'Select exactly 4 atoms to define new improper!')
            return
        q_atoms = []
        for selected in selections:
            q_atoms.append(int(self.qatoms_listbox.get(selected).split()[0]))

        #Find out how q atoms are bonded to setup correct {q2: [q1, q3,q4]}
        # q2 is bonded to all and q1 is only bonded to q2. q3 and q4 can be bonded to others or not.
        q2 = None
        q1 = None
        s = None
        bonds = dict()
        for state in range(self.evb_states.get()):
            #Check if improper is possible in state, else clear dict and continue! (len=3)
            check_len = []
            for q in q_atoms:
                bonds[q] = []
                for qj in self.q_bonds[q][state]:
                    if qj in q_atoms:
                        if qj not in bonds[q]:
                            bonds[q].append(qj)
                check_len.append(len(bonds[q]))
                if len(bonds[q]) == 3:
                    q2 = q
            if 3 in check_len:
                s = state
                break
        if not s:
            self.app.errorBox('Error', 'Could not add improper. Maybe it is not supposed to be improper?')
            return
        del q_atoms[q_atoms.index(q2)]

        #Find Q atom that is not bonded to q2 in state n --> assign this as q1:
        for state in range(self.evb_states.get()):
            if state != s:
                for qj in q_atoms:
                    if qj != q2:
                        if qj not in self.q_bonds[q2][state]:
                            q1 = qj
                            break

        print bonds
        print 'q2 is %d in state %d' % (q2, s+1)

        if not q1:
            self.app.errorBox('Error', 'Could not add new improper. '
                                       'It seems selection is not unique and therefore already existing.')
            return
        print 'q1 is %d' % q1
        del q_atoms[q_atoms.index(q1)]

        #Find q3 and q4 based on atomnumbers:
        q3, q4 = q_atoms[0], q_atoms[1]

        if int(self.q_atom_nr[q3]) > int(self.q_atom_nr[q4]):
            q3,q4 = q4, q3

        self.q_impropers[q2] = [q1, q3, q4]
        self.update_impropers()

    def del_improper(self):
        """
        delete selected improper
        """
        selections = map(int, self.changeimproper_listbox.curselection())

        for selected in selections:
            line = self.changeimproper_listbox.get(selected)
            atom2 = line.split()[1]
            #Find Q atom number for atom2:
            for q in self.q_atom_nr.keys():
                if int(self.q_atom_nr[q]) == int(atom2):
                    if q in self.q_impropers.keys():
                        del self.q_impropers[q]

        self.update_impropers()
        self.update_status()

    def add_softpair(self):
        """
        add additional softpair (that is not autogenerated)
        """
        selections = map(int, self.qatoms_listbox.curselection())

        if len(selections) != 2:
            self.app.errorBox('Warning','Select exactly two Q-atoms for soft pairs!')
            return

        for i in range(len(selections)):
            if i == 0:
                q1 = int(self.qatoms_listbox.get(selections[i]).split()[0])
            else:
                q2 = int(self.qatoms_listbox.get(selections[i]).split()[0])

        if q1 not in self.softpairs_added.keys():
            self.softpairs_added[q1] = q2
        elif q2 not in self.softpairs_added.keys():
            if self.softpairs_added[q1] != q2:
                self.softpairs_added[q2] = q1
        else:
            print 'Could not add softpair'
        self.update_soft_pairs()

    def del_softpair(self):
        """
        Delete manually added softpair (it is not possible to delete required soft-pairs)
        """
        selections = map(int, self.softpairs_listbox.curselection())

        if len(selections) == 0:
            return

        for selected in selections:
            q1, q2 = map(int, self.softpairs_listbox.get(selected).split()[0:2])
            if q1 in self.softpairs_added.keys():
                if self.softpairs_added[q1] == q2:
                    del self.softpairs_added[q1]
            if q2 in self.softpairs_added.keys():
                if self.softpairs_added[q2] == q1:
                    del self.softpairs_added[q2]
        self.update_soft_pairs()

    def add_excludedpair(self):
        selections = map(int, self.qatoms_listbox.curselection())

        if len(selections) != 2:
            self.app.errorBox('Warning','Select exactly two atoms for excluded pairs!')
            return

        for i in range(len(selections)):
            if i == 0:
                atom1 = int(self.qatoms_listbox.get(selections[i]).split()[1])
            else:
                atom2 = int(self.qatoms_listbox.get(selections[i]).split()[1])

        atompair = '%6d %6d'  % (atom1, atom2)
        atompair_rev = '%6d %6d' % (atom2, atom1)

        if atompair not in self.excluded_pairs.keys() and atompair_rev not in self.excluded_pairs.keys():
            self.excluded_pairs[atompair] = [1,1,1,1]
            self.update_excluded_pairs()

    def del_excludedpair(self):
        selections = map(int, self.excludedpairs_listbox.curselection())

        if len(selections) == 0:
            return

        for selected in selections:
            atom1, atom2 = map(int, self.excludedpairs_listbox.get(selected).split()[0:2])
            atompair = '%6d %6d' % (atom1, atom2)
            atompair_rev = '%6d %6d' % (atom2, atom1)
            if atompair in self.excluded_pairs.keys():
                del self.excluded_pairs[atompair]
            elif atompair_rev in self.excluded_pairs.keys():
                del self.excluded_pairs[atompair_rev]

        self.update_excluded_pairs()

    def excluded_state(self, state):
        """
        Toggle excluded pair ON/OFF in state
        """
        selections = map(int, self.excludedpairs_listbox.curselection())
        if len(selections) == 0:
            return

        indexes = []
        for selected in selections:
            indexes.append(selected)
            atom1, atom2 = map(int, self.excludedpairs_listbox.get(selected).split()[0:2])
            excluded_state = int(self.excludedpairs_listbox.get(selected).split()[state + 2])
            excluded_state = abs(excluded_state - 1)
            atompair = '%6d %6d' % (atom1, atom2)
            atompair_rev = '%6d %6d' % (atom2, atom1)
            if atompair_rev in self.excluded_pairs.keys():
                atompair = atompair_rev
            self.excluded_pairs[atompair][state] = excluded_state

        self.update_excluded_pairs()
        for i in indexes:
            self.excludedpairs_listbox.select_set(indexes[i])

    def add_nonq_excludedpair(self):
        """
        Select atoms from pdb file for excluded pair
        """
        pass

    def add_temperature(self):
        try:
            temp = '%.2f' % float(self.temperature.get())
            runs = int(self.runs.get())
        except:
            print 'Invalid values for temperature and runs'
            return

        if temp not in self.temp_runs.keys():
            self.temp_runs[temp] = runs
        else:
            self.temp_runs[temp] += runs

        self.update_temperatures()

    def del_temperature(self):
        selections = map(int, self.temp_listbox.curselection())

        if len(selections) == 0:
            return

        for selected in selections:
            temp = self.temp_listbox.get(selected).split()[0]
            try:
                del self.temp_runs[temp]
            except:
                continue

        self.update_temperatures()

    def edit_charge(self):
        """
        Takes new charges from entry fields, and replaces original charges
        ==> self.q_charges
        ==> self.charge_listbox
        """
        selections = map(int, self.charge_listbox.curselection())
        if len(selections) != 1:
            return
        try:
            charge1 = float(self.charge1.get())
            charge2 = float(self.charge2.get())
            charge3 = charge1
            charge4 = charge1
            if self.evb_states.get() > 2:
                charge3 = float(self.charge3.get())
            if self.evb_states.get() == 4:
                charge4 = float(self.charge4.get())
            charges = [charge1, charge2, charge3, charge4]

        except:
            print 'Can not change selection to specified charge!'
            return

        try:
            qi = int(self.charge_listbox.get(selections[0]).split()[0])
        except:
            return

        for state in range(len(charges)):
            self.q_charges[qi][state] = charges[state]

        self.update_charges()
        self.charge_listbox.selection_set(selections[0])

    def add_softcore(self):
        """
        Adds selected q atoms from Q-atoms listbox to softcore list
        """
        selections = map(int, self.qatoms_listbox.curselection())

        if len(selections) < 1:
            return

        for selected in selections:
            try:
                q = int(self.qatoms_listbox.get(selected).split()[0])
                self.q_softcore[q] = [0, 0, 0, 0]
            except:
                continue

        self.update_softcore()

    def del_softcore(self):
        """
        Deletes selected softcores in softore listbox
        """
        selections = map(int, self.softcore_listbox.curselection())

        if len(selections) < 1:
            return

        for selected in reversed(selections):
            q = int(self.softcore_listbox.get(selected).split()[0])
            del self.q_softcore[q]

        self.update_softcore()

    def edit_softcore(self):
        """
        Edit softcore value in softcore listbox
        """
        selections = map(int, self.softcore_listbox.curselection())

        if len(selections) < 1:
            return

        soft = [0, 0, 0, 0]

        try:
            s1 = float(self.softcore1.get())
            s2 = float(self.softcore2.get())
            s3, s4 = s1, s1
            if self.evb_states.get() > 2:
                s3 = float(self.softcore3.get())
            if self.evb_states.get() == 4:
                s4 = float(self.softcore4.get())

            soft = [s1, s2, s3, s4]

        except:
            print 'Please specify valid number'
            return

        for selected in selections:
            try:
                q = int(self.softcore_listbox.get(selections[selected]).split()[0])
            except:
                print 'Could not get Q atom number ...'

            for state in range(4):
                self.q_softcore[q][state] = soft[state]

        self.update_softcore()
        self.softcore_listbox.selection_set(selections[0])


    def set_elscale(self):
        """
        Takes new scaling factor from self.elscale field and replaces original value
        """
        try:
            selections = map(int, self.elscale_listbox.curselection())
        except:
            return

        if len(selections) != 1:
            return

        q1, q2 = map(int, self.elscale_listbox.get(selections[0]).split()[0:2])
        scale = self.elscale.get()

        qpair = '%3d %3d' % (q1, q2)
        self.qpair_elscale[qpair] = scale
        self.update_elscale()

    def update_elscale(self):
        self.elscale_listbox.delete(0, END)

        for i in self.qpair_elscale.keys():
            scale = str(self.qpair_elscale[i])
            self.elscale_listbox.insert(END, '%7s   %s' % (i, scale))

    def add_elscale(self):
        try:
            selections = map(int, self.qatoms_listbox.curselection())
            if len(selections) != 2:
                self.app.errorBox('Warning', 'Selected exactly 2 atoms for el scale!')
                return
        except:
            return

        for i in range(len(selections)):
            if i == 0:
                q1 = int(self.qatoms_listbox.get(selections[i]).split()[0])
            else:
                q2 = int(self.qatoms_listbox.get(selections[i]).split()[0])

        qpair = '%3d %3d' % (q1, q2)
        if qpair not in self.qpair_elscale.keys():
            self.qpair_elscale[qpair] = 1.00

        self.update_elscale()

    def del_elscale(self):
        try:
            selections = map(int, self.elscale_listbox.curselection())
        except:
            return

        for selected in selections:
            q1, q2 = map(int, self.elscale_listbox.get(selected).split()[0:2])
            qpair = '%3d %3d' % (q1, q2)
            if qpair in self.qpair_elscale.keys():
                del self.qpair_elscale[qpair]
            else:
                qpair = '%3d %3d' % (q2, q1)
                if qpair in self.qpair_elscale.keys():
                    del self.qpair_elscale[qpair]
                else:
                    self.app.errorBox('Error', 'An error has occured in el scale lists. Please report problem!')

        self.update_elscale()


    def sigmoid_point(self, beta, gamma, x):
        """
        Returns lambda value for point x when
        lambda(x) = 1 / (beta + exp(-10x + 5)) + gamma
        """
        if round(x, 3) == 1:
            lamda = 1./beta + gamma
        elif round(x, 3) == 0:
            lamda = gamma
        else:
            lamda = 1. / (beta + np.exp(-10.0 * x + 5.)) + gamma

        return lamda

    def add_lambda_steps(self):
        """
        Add sequence of lambda steps to state 1-4 based on selection
        """
        #Check that sum lambda = 1 or break!
        sum_start = float(self.sum_start.get())
        sum_end = float(self.sum_end.get())

        #Check that two states are locked
        lock_lambda = [self.lock_lambda1.get(), self.lock_lambda2.get(),
                       self.lock_lambda3.get(),self.lock_lambda4.get()]
        if sum(lock_lambda) != 2:
            self.app.errorBox('Warning', '2 lambda values must be locked!')
            return

        if sum_start != float(1) or sum_end != float(1):
            self.app.errorBox('Error', 'Invalid lambda sum. Must be equal to 1!')
            return

        #Check if lambda steps exist and get step number:
        if len(self.step_lambdavalues) > 0:
            step = max(self.step_lambdavalues.keys())
        else:
            step = 0

        l1_start = float(self.start1.get())
        l1_end = float(self.end1.get())
        self.start1.delete(0, END)
        self.start1.insert(0, '%.2f' % l1_end)

        l2_start = float(self.start2.get())
        l2_end = float(self.end2.get())
        self.start2.delete(0, END)
        self.start2.insert(0, '%.2f' % l2_end)

        l3_start = float(self.start3.get())
        l3_end = float(self.end3.get())
        self.start3.delete(0, END)
        self.start3.insert(0, '%.2f' % l3_end)

        l4_start = float(self.start4.get())
        l4_end = float(self.end4.get())
        self.start4.delete(0, END)
        self.start4.insert(0, '%.2f' % l4_end)

        #Check that values are changing:
        if (l1_start == l1_end) and (l2_start == l2_end):
            if (l3_start == l3_end) and (l4_start == l4_end):
                return

        start_list = [l1_start, l2_start, l3_start, l4_start]
        end_list = [l1_end, l2_end, l3_end, l4_end]

        #Find out lambda values to vary:
        starts = []
        ends = []
        indexes = []
        for i in range(len(lock_lambda)):
            if lock_lambda[i] == 0:
                starts.append(start_list[i])
                ends.append(end_list[i])
                indexes.append(i)

        increase_ind = indexes[1]
        decrease_ind = indexes[0]

        increase_start = start_list[increase_ind]
        increase_end = end_list[increase_ind]

        decrease_start = start_list[decrease_ind]
        decrease_end  = end_list[decrease_ind]

        if abs(round(decrease_start - decrease_end, 3)) != abs(round(increase_start - increase_end, 3)):
            self.app.errorBox('Warning', 'Illegal lambda sums for varying lambda pair.')

        if increase_start > increase_end:
            decrease_start, decrease_end, increase_start, increase_end, decrease_ind, increase_ind= \
                increase_start, increase_end, decrease_start, decrease_end, increase_ind, decrease_ind

        #get step_size
        step_size  = float(self.lambdastep.get())
        steps = (increase_end - increase_start) / step_size + 1
        x_step = 1. / (steps - 1)

        beta = 1.00 / (increase_end - increase_start)
        gamma = increase_start

        print 'BETA GAMMA'
        print (beta, gamma)

        #nitialize
        x = 0
        step += 1
        self.step_lambdavalues[step] = [l1_start, l2_start, l3_start, l4_start]
        self.step_lambdavalues[step][increase_ind] = increase_start
        self.step_lambdavalues[step][decrease_ind] = decrease_start

        #l1 decreases until end, while l2 increases
        loop_count = 0
        while True:
            loop_count += 1
            x += x_step
            if self.lambda_sampling.get() == 'Linear':
                decrease_start = round(abs(decrease_start - step_size), 3)
                increase_start = round(increase_start + step_size, 3)
            else:
                decrease_start = round(self.sigmoid_point(beta, gamma, 1. - x), 3)
                increase_start = round(self.sigmoid_point(beta, gamma, x), 3)


            l_values = [l1_start, l2_start, l3_start, l4_start]
            l_values[increase_ind] = increase_start
            l_values[decrease_ind] = decrease_start

            if l_values not in self.step_lambdavalues.values():
                step += 1
                self.step_lambdavalues[step] = l_values


            print (decrease_start, increase_start)

            if float(abs(increase_start)) >= float(increase_end):
                break

            if loop_count > 2999:
                print 'WOOOPS!! You just reached 3000 lambda steps!!'
                print 'Something probably went wrong!? Or did you really want this?'
                print 'Either way, please report this problem!'
                break

        self.update_lambdasteps()

    def del_lambda_steps(self):
        """
        Deletes selected lambda steps from list, and renumbers list:
        """
        try:
            selections = map(int, self.lambdasteps_listbox.curselection())
        except:
            return

        for selected in selections:
            step = int(self.lambdasteps_listbox.get(selected).split()[0])
            del self.step_lambdavalues[step]

        self.lambdasteps_listbox.delete(0, END)
        tmp = dict()
        #Renumber lambdavalues
        if len(self.step_lambdavalues) > 0:
            step = 0
            for old in sorted(self.step_lambdavalues.keys()):
                step += 1
                tmp[step] = self.step_lambdavalues[old]

        self.step_lambdavalues = tmp
        self.update_lambdasteps()

    def reverse_lambda_list(self):
        nr = 0
        rev_lambda = dict()

        for old_nr in reversed(sorted(self.step_lambdavalues.keys())):
            nr += 1
            rev_lambda[nr] = self.step_lambdavalues[old_nr]

        self.step_lambdavalues = rev_lambda

        self.update_lambdasteps()

    def update_change_atoms_listbox(self):
        self.atoms_listbox.delete(0, END)
        for q in sorted(self.q_element.keys()):
            a1, a2 = self.q_element[q][0:2]
            a3, a4 = ' ', ' '
            if self.evb_states.get() > 2:
                a3 = self.q_element[q][2]
            if self.evb_states.get() == 4:
                a4 = self.q_element[q][3]
            self.atoms_listbox.insert(END, '%3d  %3s  %3s  %3s  %3s' % (q, a1, a2, a3, a4))

    def add_atom_change(self):
        """
        Add selected atoms in Qatoms listbox to the change atom list
        """
        

        selections = map(int, self.qatoms_listbox.curselection())

        if len(selections) < 1:
            return

        for selected in selections:
            q = int(self.qatoms_listbox.get(selected).split()[0])
            if q not in self.q_element.keys():
                types = []
                #Get atomtype(s) for Q atom:
                for state in range(4):
                    types.append(self.q_atomtypes[q][state])

                #Try to figure out what element Qi is based on its weight
                atoms = []
                for at in types:
                    dmass = 1000.00
                    atom = '??'
                    mass = float(self.atomtype_prm[at][-1])
                    for elem in self.massDict.keys():
                        if abs(self.massDict[elem] - mass) < dmass:
                            dmass = abs(self.massDict[elem] - mass)
                            atom = elem
                    atoms.append(atom)
                self.q_element[q] = atoms

        self.update_change_atoms_listbox()

    def del_atom_change(self):
        """
        Deletes selection from the change atom listbox
        """
        pass

    def change_atom_state(self, state=0):
        """
        Takes selected atom in element list, and replaces the selected state atom
        """

        elem_ind = map(int, self.elem_listbox.curselection())

        if len(elem_ind) != 1:
            print 'Select exactly one element!'

        elem = self.elem_listbox.get(elem_ind[0]).strip()

        selection = map(int, self.atoms_listbox.curselection())

        if len(selection) < 1:
            print 'No Q-atoms selected!'
            return

        for selected in selection:
            q = int(self.atoms_listbox.get(selected).split()[0])
            self.q_element[q][state] = elem
            old_type = self.q_atomtypes[q][state]
            new_type = '.%s.' % elem
            mass = self.massDict[elem]
            if not new_type in self.atomtype_prm.keys():
                self.atomtype_prm[new_type] = copy.deepcopy(self.atomtype_prm[old_type])
                self.atomtype_prm[new_type][-1] = mass
            self.q_atomtypes[q][state] = new_type

        self.update_all()
        self.update_pymol_bonds()
            

    def list_q_atoms_event(self, *args):
        """
        Highlight Q-atoms in pymol when Q-atoms are selected
        """
        if not self.session:
            return

        self.session.stdin.write('hide spheres\n')
        atoms = []
        try:
            selections = map(int, self.qatoms_listbox.curselection())
        except:
            return

        if len(selections) == 0:
            self.highlight_q_atoms()
            return

        for selected in selections:
            atoms.append(int(self.qatoms_listbox.get(selected).split()[1]))

        pymol_sel_string = ''

        for atom in range(len(atoms)):
            #Make sure list has not ended:
            if (atom + 1) != len(atoms):
                pymol_sel_string += 'id %d or ' % int(atoms[atom])
            else:
                #End of atom list reached
                pymol_sel_string += 'id %d' % int(atoms[atom])

        self.session.stdin.write('show spheres, %s\n' % pymol_sel_string)

    def list_fepatoms_event(self, *args):
        try:
            selections = map(int, self.fepatoms_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        for selected in selections:

            q = int(self.fepatoms_listbox.get(selected).split()[0])

            self.qatoms_listbox.selection_set(q - 1)
            self.qatoms_listbox.yview(q - 1)
        self.list_q_atoms_event()

    def list_charge_event(self, *args):
        self.charge1.delete(0, END)
        self.charge2.delete(0, END)
        self.charge3.delete(0, END)
        self.charge4.delete(0, END)

        try:
            selections = map(int, self.charge_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        indexes = []
        for selected in selections:
            try:
                q1 = int(self.charge_listbox.get(selected).split()[0])
                _ind = q1-1
                self.qatoms_listbox.selection_set(_ind)
                indexes.append(_ind)

            except:
                return
        self.qatoms_listbox.yview(min(indexes))
        if len(selections) == 1:
            try:
                #8 (4) 7 (3) 6 (2)
                charges = self.charge_listbox.get(selections[0]).split()[1:]
                self.charge1.insert(0, charges[0])
                self.charge2.insert(0, charges[1])
                if len(charges) > 5:
                    self.charge3.insert(0, charges[2])
                if len(charges) > 7:
                    self.charge4.insert(0, charges[3])

            except:
                pass
        self.list_q_atoms_event()

    def list_softcore_event(self, *args):
        self.softcore1.delete(0, END)
        self.softcore2.delete(0, END)
        self.softcore3.delete(0, END)
        self.softcore4.delete(0, END)

        try:
            selections = map(int, self.softcore_listbox.curselection())
        except:
            return

        if len(selections) == 0:
            return

        self.qatoms_listbox.selection_clear(0, END)
        indexes = []
        for selected in selections:
            try:
                q1 = int(self.charge_listbox.get(selected).split()[0])
                _ind = q1-1
                self.qatoms_listbox.selection_set(_ind)
                indexes.append(_ind)

            except:
                return
        self.qatoms_listbox.yview(min(indexes))
        if len(selections) == 1:
            try:
                #8 (4) 7 (3) 6 (2)
                charges = self.softcore_listbox.get(selections[0]).split()[1:]
                self.softcore1.insert(0, charges[0])
                self.softcore2.insert(0, charges[1])
                if len(charges) > 3:
                    self.softcore3.insert(0, charges[2])
                if len(charges) > 4:
                    self.softcore4.insert(0, charges[3])

            except:
                pass
        self.list_q_atoms_event()

    def list_changetypes_event(self, *args):
        try:
            selections = map(int, self.changetypes_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        atomtypes = []
        indexes = []
        for selected in selections:
            try:
                q1 = int(self.changetypes_listbox.get(selected).split()[0])
                _ind = q1 - 1
                indexes.append(_ind)
                self.qatoms_listbox.selection_set(_ind)
                types = self.changetypes_listbox.get(selected).split()[1:]

                for atomtype in types:
                    if atomtype not in atomtypes:
                        atomtypes.append(atomtype)
            except:
                continue
        self.qatoms_listbox.yview(min(indexes))

        if len(atomtypes) > 0:
            self.atomtypes_listbox.selection_clear(0, END)
            list_types = self.atomtypes_listbox.get(0,END)
            ind_list = []
            for i in range(len(list_types)):
                if list_types[i].split()[0] in atomtypes:
                    ind_list.append(i)

            for index in ind_list:
                self.atomtypes_listbox.selection_set(index)
        self.list_q_atoms_event()

    def list_changebond_event(self, *args):
        try:
            selections = map(int, self.changebond_listbox.curselection())
        except:
            return

        qatoms = self.qatoms_listbox.get(0, END)

        self.qatoms_listbox.selection_clear(0, END)
        self.bondtypes_listbox.selection_clear(0, END)
        bonds = self.bondtypes_listbox.get(0, END)
        indexes = []
        for selected in selections:
            atom1, atom2 = map(int, self.changebond_listbox.get(selected).split()[0:2])
            states = map(str, self.changebond_listbox.get(selected).split()[2:])

            for ind in range(len(qatoms)):
                if int(qatoms[ind].split()[1]) == atom1 or int(qatoms[ind].split()[1]) == atom2:
                    self.qatoms_listbox.selection_set(ind)
                    indexes.append(ind)

            for bond_sele in states:
                try:
                    bond = int(bond_sele)
                    for i in range(len(bonds)):
                        if bond == int(bonds[i].split()[0]):
                            self.bondtypes_listbox.selection_set(i)
                except:
                    continue
        self.qatoms_listbox.yview(min(indexes))
        self.list_q_atoms_event()

    def list_changeangle_event(self, *args):
        try:
            selections = map(int, self.changeangle_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        self.angletypes_listbox.selection_clear(0, END)
        angles = self.angletypes_listbox.get(0, END)
        indexes = []
        for selected in selections:
            a1, a2, a3 = map(int, self.changeangle_listbox.get(selected).split()[0:3])
            states = map(str, self.changeangle_listbox.get(selected).split()[3:])

            for i in self.q_atom_nr.keys():
                if int(self.q_atom_nr[i]) == a1:
                    q1 = i
                if int(self.q_atom_nr[i]) == a2:
                    q2 = i
                if int(self.q_atom_nr[i]) == a3:
                    q3 = i
            for q in (q1 - 1, q2 - 1, q3 - 1):
                self.qatoms_listbox.selection_set(q)
                indexes.append(q)

            for prm in range(len(angles)):
                if angles[prm].split()[0] in states:
                    self.angletypes_listbox.selection_set(prm)

        self.qatoms_listbox.yview(min(indexes))
        self.list_q_atoms_event()

    def list_changetorsion_event(self, *args):
        try:
            selections = map(int, self.changetorsion_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        self.torsiontypes_listbox.selection_clear(0, END)
        torsions = self.torsiontypes_listbox.get(0, END)
        indexes = []
        for selected in selections:
            a1, a2, a3, a4 = map(int, self.changetorsion_listbox.get(selected).split()[0:4])
            states = map(str, self.changetorsion_listbox.get(selected).split()[4:])

            for i in self.q_atom_nr.keys():
                if int(self.q_atom_nr[i]) == a1:
                    q1 = i
                if int(self.q_atom_nr[i]) == a2:
                    q2 = i
                if int(self.q_atom_nr[i]) == a3:
                    q3 = i
                if int(self.q_atom_nr[i]) == a4:
                    q4 = i

            for q in (q1 - 1, q2 - 1, q3 - 1, q4 - 1):
                self.qatoms_listbox.selection_set(q)
                indexes.append(q)

            for prm in range(len(torsions)):
                if torsions[prm].split()[0] in states:
                    self.torsiontypes_listbox.selection_set(prm)

        self.qatoms_listbox.yview(min(indexes))
        self.list_q_atoms_event()

    def list_changeimproper_event(self, *args):
        try:
            selections = map(int, self.changeimproper_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        self.impropertypes_listbox.selection_clear(0, END)
        torsions = self.impropertypes_listbox.get(0, END)
        indexes = []
        for selected in selections:
            a1, a2, a3, a4 = map(int, self.changeimproper_listbox.get(selected).split()[0:4])
            states = map(str, self.changeimproper_listbox.get(selected).split()[4:])

            for i in self.q_atom_nr.keys():
                if int(self.q_atom_nr[i]) == a1:
                    q1 = i
                if int(self.q_atom_nr[i]) == a2:
                    q2 = i
                if int(self.q_atom_nr[i]) == a3:
                    q3 = i
                if int(self.q_atom_nr[i]) == a4:
                    q4 = i

            for q in (q1 - 1, q2 - 1, q3 - 1, q4 - 1):
                self.qatoms_listbox.selection_set(q)
                indexes.append(q)

            for prm in range(len(torsions)):
                if torsions[prm].split()[0] in states:
                    self.impropertypes_listbox.selection_set(prm)

        self.qatoms_listbox.yview(min(indexes))
        self.list_q_atoms_event()

    def list_softpairs_event(self, *args):
        try:
            selections = map(int, self.softpairs_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        indexes = []

        for selected in selections:
            q1, q2 = map(int, self.softpairs_listbox.get(selected).split()[0:2])
            indexes.append(q1 - 1)
            indexes.append(q2 - 1)
            self.qatoms_listbox.selection_set(q1 - 1)
            self.qatoms_listbox.selection_set(q2 - 1)

        self.qatoms_listbox.yview(min(indexes))
        self.list_q_atoms_event()

    def list_excludedpairs_event(self, *args):
        try:
            selections = map(int, self.excludedpairs_listbox.curselection())
        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        indexes = []

        for selected in selections:
            a1, a2 = map(int, self.excludedpairs_listbox.get(selected).split()[0:2])

            for i in self.q_atom_nr.keys():
                if int(self.q_atom_nr[i]) == a1 or int(self.q_atom_nr[i]) == a2:
                    indexes.append(i - 1)
                    self.qatoms_listbox.selection_set(i-1)

        self.qatoms_listbox.yview(min(indexes))
        self.list_q_atoms_event()

    def list_elscale_event(self, *args):

        try:
            selections = map(int, self.elscale_listbox.curselection())

        except:
            return

        self.qatoms_listbox.selection_clear(0, END)
        indexes = []

        for selected in selections:
            q1, q2 = map(int, self.elscale_listbox.get(selected).split()[0:2])
            scale = float(self.elscale_listbox.get(selected).split()[2])
            indexes.append(q1 - 1)
            indexes.append(q2 - 1)
            self.qatoms_listbox.selection_set(q1 - 1)
            self.qatoms_listbox.selection_set(q2 - 1)
            self.elscale.delete(0, END)
            self.elscale.insert(0, scale)

        self.qatoms_listbox.yview(min(indexes))
        self.list_q_atoms_event()

    def list_atomchange_event(self, *args):
        """
        Do something when change atoms are selected... hmm
        """
        selections = map(int, self.atoms_listbox.curselection())

        if len(selections) < 1:
            return
        self.qatoms_listbox.selection_clear(0, END)
        self.elem_listbox.selection_clear(0, END)

        elemlist = list(self.elem_listbox.get(0, END))

        for selected in selections:
            q = int(self.atoms_listbox.get(selected).split()[0])
            self.qatoms_listbox.selection_set(q - 1)
            elements = self.atoms_listbox.get(selected).split()[1:]
            for i in range(len(elemlist)):
                if elemlist[i].strip() in elements:
                    self.elem_listbox.selection_set(i)

        self.list_q_atoms_event()

    def list_monitorgroups_event(self, *args):
        self.group_atoms_listbox.delete(0, END)
        selections = map(int, self.groups_listbox.curselection())

        if len(selections) != 1:
            return

        for i in selections:
            group = int(self.groups_listbox.get(i).split()[1])
            for atom in self.monitor_groups[group]:
                self.group_atoms_listbox.insert(END, '%6s' % atom)


    def toggle_pymol(self, *args):
        """
        Start/close pymol when checkbutton is toggled
        """
        if self.sync_pymol.get() == 1:
            self.auto_setup.config(state=NORMAL)
            self.optimize_button.config(state=NORMAL)
            self.autocharge.config(state=NORMAL)
            self.force_field.config(state=NORMAL)
            self.selH_check.config(state=NORMAL)
            if self.session:
                try:
                    os.killpg(self.session.pid, signal.SIGTERM)
                except:
                    pass
            self.start_pymol()

        elif self.sync_pymol.get() == 0:
            self.auto_setup.config(state=DISABLED)
            self.autocharge.config(state=DISABLED)
            self.optimize_button.config(state=DISABLED)
            self.force_field.config(state=DISABLED)
            self.selH_check.config(state=DISABLED)
            try:
                os.killpg(self.session.pid, signal.SIGTERM)
            except:
                pass

    def geom_opt(self):
        """
        Geometry optimization of selected atoms and last clicked state
        """
        if not self.session:
            return

        try:
            selections = map(int, self.qatoms_listbox.curselection())
        except:
            return

        pml_atom = []
        for selected in selections:
            atom = self.qatoms_listbox.get(selected).split()[1]
            pml_atom.append('id %s' % atom)

        pml_string = ' or '.join(pml_atom)

        #Get last selected state from logfile:
        state = 'state2'
        with open(self.app.workdir+'/.tmpfile', 'r') as log:
            for line in log:
                if 'You clicked' in line:
                    state = line.split('/')[1]

        print 'Fixing geometry for %s: %s' % (state, pml_string)
        self.session.stdin.write('clean %s and (%s)\n' % (state, pml_string))

    def update_pymol_bonds(self):
        """
        Updates atoms visible and not in pymol in each state
        """
        if not self.session:
            return

        hide_show = {0: 'hide', 1:'show'}

        atom_color = {'H': 'white',
                      'C': 'chartreuse',
                      'N': 'blue',
                      'O': 'red',
                      'P': 'orange',
                      'S': 'yellow',
                      'F': 'green',
                      'Cl': 'green',
                      'Br': 'firebrick',
                      'Na': 'violet',
                      'Mg': 'forest',
                      'K': 'violet',
                      'Ca': 'forest',
                      'Fe': 'orange',
                      'Zn': 'pink'}

        for atom in self.fep_atoms.keys():
            for state in range(self.evb_states.get()):
                cmd = hide_show[self.fep_atoms[atom][state]]
                self.session.stdin.write('%s sticks, state%d and id %s\n' % (cmd, state + 1, atom))

        #Go through change atoms list and make changes:
        for q in self.q_element.keys():
            atom_id = self.q_atom_nr[q]
            for state in range(self.evb_states.get()):
                elem = self.q_element[q][state].strip()
                color = atom_color[elem]
                self.session.stdin.write('alter state%d and id %s, elem="%s"\n' % (state + 1, atom_id, elem))
                self.session.stdin.write('color %s, state%d and id %s\n' % (color, state + 1, atom_id))


    def start_pymol(self):
        """
        Start syncing of Q-atom selection to pymol
        """
        #QPyMol default settings ('set valence, 0.1' removed)
        self.pymol_settings = ['space cmyk', 'set sphere_scale, 0.4',
                               'set sphere_transparency, 0.5', 'set sphere_color, lightblue',
                               'set sphere_quality, 2', 'set stick_radius, 0.17', 'set stick_quality, 10',
                               'set defer_builds_mode, 3', 'set surface_quality, 1',
                               'set spec_power, 250', 'set spec_reflect, 2',
                               'set cartoon_fancy_helices, 1']

        self.app.pymol_running = True

        tmpfile = open(self.app.workdir+'/.tmpfile','wb')
        if 'darwin' in sys.platform:
            self.session = Popen(["MacPyMol", "-p -x -i", "%s" % self.pdbfile], stdout=tmpfile, stdin=PIPE, preexec_fn=os.setsid)
        else:
            self.session = Popen(["pymol", "-p", "%s" % self.pdbfile], stdout=tmpfile, stdin=PIPE, preexec_fn=os.setsid)

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

        #Set atom selection mode
        self.session.stdin.write('set mouse_selection_mode, 0\n')

        #For some reason the internal pymol gui needs to be loaded to initialize grid mode..
        self.session.stdin.write('set internal_gui=1\n')

        #Set grid mode on in pymol (display changes in states)
        self.session.stdin.write('set grid_mode, 1\n')

        #Copy pdb file to state1-4
        pdbname = self.pdbfile.split('/')[-1].split('.')[0]
        for i in range(4):
            self.session.stdin.write('copy state%d, %s\n' % ((i + 1), pdbname))
            if (i + 1) > self.evb_states.get():
                self.session.stdin.write('disable state%d\n' % (i + 1))

        #Remove original pdb name from pymol:
        self.session.stdin.write('delete %s\n' % pdbname)

        #Remove internal gui
        self.session.stdin.write('set internal_gui=0\n')

        #Highlight Q-atoms
        self.highlight_q_atoms()

        #Update fep atoms
        self.update_pymol_bonds()

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
                            print selected
                            #Get atomnumer of clicked atom:
                            self.session.stdin.write('id_atom sele\n')
                        if 'cmd.id_atom:' in line:
                            self.session.stdin.write('select none\n')
                            atomnr = line.split()[2].split(')')[0]
                            qatoms_list = self.qatoms_listbox.get(0, END)
                            for q in range(len(qatoms_list)):
                                if int(atomnr) == int(qatoms_list[q].split()[1]):
                                    self.qatoms_listbox.select_set(q)
                                    self.list_q_atoms_event()
                                    self.qatoms_listbox.yview(q)

            time.sleep(0.2)

        self.sync_pymol.set(0)
        self.sync_check.update()
        self.session = None
        self.app.pymol_running = False
        try:
            os.killpg(self.session.pid, signal.SIGTERM)
        except:
            pass

    def highlight_q_atoms(self):
        """
        Highlights Q-atoms in pymol and zoom
        """
        if not self.session:
            return
        self.session.stdin.write('hide sticks, all\n')
        self.session.stdin.write('hide spheres, all\n')
        self.session.stdin.write('color gray, name C*\n')
        if len(self.q_atom_nr) == 0:
            return

        atoms= []
        for qnr in sorted(self.q_atom_nr.keys()):
            atoms.append(int(self.q_atom_nr[qnr]))

        pymol_sel_string = ''

        for atom in range(len(atoms)):
            #Make sure list has not ended:
            if (atom + 1) != len(atoms):
                pymol_sel_string += 'id %d or ' % int(atoms[atom])
            else:
                #End of atom list reached
                pymol_sel_string += 'id %d' % int(atoms[atom])

        self.session.stdin.write('color chartreuse, (%s) and name C* \n' % pymol_sel_string)
        self.session.stdin.write('show sticks, %s \n' % pymol_sel_string)
        self.session.stdin.write('hide lines, %s \n' % pymol_sel_string)
        self.session.stdin.write('orient %s\n' % pymol_sel_string)
        self.session.stdin.write('zoom %s, 5\n' % pymol_sel_string)

        self.update_pymol_bonds()

    def write_fep(self):
        self.fep_written = True
        fepname = self.pdbfile.split('/')[-1].split('.')[0] + '.fep'

        if not os.path.isdir('%s/inputfiles' % self.app.workdir):
            os.mkdir('%s/inputfiles' % self.app.workdir)

        fepfile = open('%s/inputfiles/%s' % (self.app.workdir, fepname), 'w')

        #Write header with info about topology (for load FEP)
        fepfile.write('!%s\n' % self.pdbfile.split('/')[-1])
        try:
            fepfile.write('!%s\n' % self.topology.split('/')[-1])
        except:
            pass

        # [FEP]
        fepfile.write('[FEP]\nstates %d\n' % self.evb_states.get())
        if self.use_soft_max_pt.get() == 1:
            fepfile.write('softcore_use_max_potential on\n')

        # [atoms]
        fepfile.write('\n[atoms]\n')
        qatoms = self.qatoms_listbox.get(0, END)
        for qatom in qatoms:
            fepfile.write('%s\n' % qatom)
        fepfile.write('\n')

        # [change_charges]
        fepfile.write('[change_charges]\n')
        charges = self.charge_listbox.get(0, END)
        for charge in charges:
            if '-------' not in charge and 'SUM' not in charge:
                fepfile.write('%s\n' % charge)
        fepfile.write('\n')

        # [atom_types]
        fepfile.write('[atom_types]\n')
        atomtypes = self.atomtypes_listbox.get(0, END)
        for atomtype in atomtypes:
            fepfile.write('%s\n' % atomtype)
        fepfile.write('\n')

        # [change_atoms]
        fepfile.write('[change_atoms]\n')
        changetypes = self.changetypes_listbox.get(0, END)
        for changetype in changetypes:
            fepfile.write('%s\n' % changetype)
        fepfile.write('\n')

        # [soft_pairs]
        softpairs = self.softpairs_listbox.get(0, END)
        if len(softpairs) > 0:
            fepfile.write('[soft_pairs]\n')

            for softpair in softpairs:
                fepfile.write('%s\n' % softpair)
            fepfile.write('\n')

        # [excluded_pairs]
        expairs = self.excludedpairs_listbox.get(0, END)
        if len(expairs) > 0:
            fepfile.write('[excluded_pairs]\n')
            for expair in expairs:
                fepfile.write('%s \n' % expair)
            fepfile.write('\n')

        # [el_scale]
        if len(self.elscale_listbox.get(0, END)) > 0:
            fepfile.write('[el_scale]\n')
            for line in self.elscale_listbox.get(0, END):
                fepfile.write(line)

        # [bond_types]
        fepfile.write('[bond_types]\n')
        bondtypes = self.bondtypes_listbox.get(0, END)
        for bondtype in bondtypes:
            fepfile.write('%s\n' % bondtype)
        fepfile.write('\n')

        # [change_bonds]
        fepfile.write('[change_bonds]\n')
        changebonds = self.changebond_listbox.get(0, END)
        for changebond in changebonds:
            fepfile.write('%s\n' % changebond)
        fepfile.write('\n')

        # [angle_types]
        angletypes = self.angletypes_listbox.get(0, END)
        if len(angletypes) > 0:
            fepfile.write('[angle_types]\n')
            for angletype in angletypes:
                fepfile.write('%s\n' % angletype)
            fepfile.write('\n')

        # [change_angles]
        changeangles = self.changeangle_listbox.get(0, END)
        if len(changeangles) > 0:
            fepfile.write('[change_angles]\n')
            for changeangle in changeangles:
                fepfile.write('%s\n' % changeangle)
            fepfile.write('\n')

        # [torsion_types]
        torsiontypes = self.torsiontypes_listbox.get(0, END)
        if len(torsiontypes) > 0:
            fepfile.write('[torsion_types]\n')
            for torsiontype in torsiontypes:
                fepfile.write('%s\n' % torsiontype)
            fepfile.write('\n')

        # [change_torsions]
        changetorsions = self.changetorsion_listbox.get(0, END)
        if len(changetorsions) > 0:
            fepfile.write('[change_torsions]\n')
            for changetorsion in changetorsions:
                fepfile.write('%s\n' % changetorsion)
            fepfile.write('\n')

        # [improper_types]
        impropertypes = self.impropertypes_listbox.get(0, END)
        if len(impropertypes) > 0:
            fepfile.write('[improper_types]\n')
            for impropertype in impropertypes:
                fepfile.write('%s\n' % impropertype)
            fepfile.write('\n')

        # [change_impropers]
        changeimpropers = self.changeimproper_listbox.get(0, END)
        if len(changeimpropers) > 0:
            fepfile.write('[change_impropers]\n')
            for changeimproper in changeimpropers:
                fepfile.write('%s\n' % changeimproper)
            fepfile.write('\n')

        # [angle_couplings]
        anglecouplings = self.angle_couplings_listbox.get(0, END)
        if len(anglecouplings) > 0:
            fepfile.write('[angle_couplings]\n')
            for anglecoupling in anglecouplings:
                fepfile.write('%s\n' % anglecoupling)
            fepfile.write('\n')

        # [torsion_couplings]
        torsioncouplings = self.torsion_couplings_listbox.get(0, END)
        if len(torsioncouplings) > 0:
            fepfile.write('[torsion_couplings]\n')
            for torsioncoupling in torsioncouplings:
                fepfile.write('%s\n' % torsioncoupling)
            fepfile.write('\n')

        # [improper_couplings]
        impropercouplings = self.improper_couplings_listbox.get(0, END)
        if len(impropercouplings) > 0:
            fepfile.write('[improper_couplings]\n')
            for impropercoupling in impropercouplings:
                fepfile.write('%s\n' % impropercoupling)
            fepfile.write('\n')

        # [softcore]
        softcores = self.softcore_listbox.get(0, END)
        if len(softcores) > 0:
            fepfile.write('[softcore]\n')
            for softcore in softcores:
                fepfile.write('%s\n' % softcore)
            fepfile.write('\n')

        #[monitor_groups]
        if len(self.monitor_groups) > 0:
            fepfile.write('[monitor_groups]\n')
            for group in sorted(self.monitor_groups.keys()):
                atoms = ' '.join(self.monitor_groups[group])
                fepfile.write('%s\n' % atoms)
            fepfile.write('\n')

        #[monitor_group_pairs]
        if len(self.monitor_pairs) > 0:
            fepfile.write('[monitor_group_pairs]\n')
            for pair in self.monitor_pairs:
                fepfile.write('%s\n' % ' '.join(pair))
            fepfile.write('\n')


        fepfile.close()
        self.app.log('info','FEP file (inputfiles/%s) written' % fepname.split('/')[-1])

    def auto_evb(self, only_charges=False):
        """
        Generates structures from pymol (must be enabled) and use ffld_server to get parameters
        """
        if self.sync_pymol.get() == 0:
            return

        #Remove tmp pdb and prm files if existing:
        files_to_delete = ['qevb_org1.pdb','qevb_org2.pdb', 'qevb_org3.pdb', 'qevb_org4.pdb',
                           'qevb_tmp1.pdb', 'qevb_tmp2.pdb', 'qevb_tmp3.pdb', 'qevb_tmp4.pdb',
                           'prm_tmp1.out', 'prm_tmp2.out', 'prm_tmp3.out', 'prm_tmp4.out']
        for tmpfile in files_to_delete:
            if os.path.isfile(self.app.workdir + '/%s' % tmpfile):
                os.remove(self.app.workdir + '/%s' % tmpfile)

        self.session.stdin.write('set internal_gui=1\n')

        #Auto expand and add H around Q-atom region, or manualy select atoms to be chopped.
        h_add_all = True
        if self.select_h.get() == 1:
            h_add_all = False


        pdb = open(self.pdbfile, 'r').readlines()

        #Get q-atoms that are changin bonds
        q_changing = []
        for pair in self.change_bonds.keys():
            atom1, atom2 = map(int, pair.split()[0:2])
            for q in self.q_atom_nr.keys():
                if (self.q_atom_nr[q] == atom1) or (self.q_atom_nr[q] == atom2):
                    q_changing.append(q)

        h_add_list = []
        residues = []
        if h_add_all:
            #Get Q-atoms and exapand around selection if possible (include entire residue):
            resnames = []
            res_atoms = []
            tmp_atoms = []
            for q in sorted(self.q_atom_res.keys()):
                if self.q_atom_res[q].split()[0] not in resnames:
                    if len(tmp_atoms) > 0:
                        res_atoms.append(tmp_atoms)
                        tmp_atoms = []
                    resnames.append(self.q_atom_res[q].split()[0])
                    residues.append('i. %d' % int(self.q_atom_res[q].split()[-1]))
                tmp_atoms.append(self.q_atom_nr[q])
            res_atoms.append(tmp_atoms)

            #Go through libraries and find connecting atoms for residues:
            connect_atoms = []
            for resname in resnames:
                tmp = []
                for lib in self.libs:
                    with open(lib, 'r') as libfile:
                        found_res = False
                        found_connections = False

                        for line in libfile:
                            if found_res:
                                if found_connections:
                                    if '[' in line or '*-------' in line:
                                        found_connections = False
                                        found_res = False
                                    elif 'head' in line or 'tail' in line:
                                        if line.split()[1] not in tmp:
                                            tmp.append(line.split()[1])

                                if '[connections]' in line:
                                    found_connections = True
                                if '{' in line and '}' in line:
                                    break
                            try:
                                if line.split()[0] == '{%s}' % resname:
                                    found_res = True
                            except:
                                continue
                connect_atoms.append(tmp)

            #Generate pml selection string for h_add on terminals that are chopped:
            for residue in range(len(residues)):
                if len(connect_atoms[residue]) > 0:
                    for atom in connect_atoms[residue]:
                        h_add_list.append('%s and name %s' % (residues[residue], atom))

            #Generate pymol selection strings


        elif not h_add_all:
            for q in sorted(self.q_atom_nr.keys()):
                residues.append('id %d' % self.q_atom_nr[q])
            try:
                selections = map(int, self.qatoms_listbox.curselection())
                for selected in selections:
                    h_add_list.append('id %s' % self.qatoms_listbox.get(selected).split()[1])
            except:
                pass

        pml_select = ' or '.join(residues)
        h_add_string = ' or '.join(h_add_list)
        print h_add_string

        #Make templates for ffld_server in pymol:
        for state in range(1, self.evb_states.get() + 1):
            self.session.stdin.write('create state%d_auto, state%d and (%s) \n' % (state, state, pml_select))
            for atomnr in sorted(self.fep_atoms.keys()):
                if int(self.fep_atoms[atomnr][state - 1]) == 0:
                    self.session.stdin.write('remove state%d_auto and id %d\n' % (state, atomnr))
            if not h_add_all:
                self.session.stdin.write('fix_chemistry state%d_auto and (%s)\n' % (state, h_add_string))
            if len(h_add_string) > 0:
                self.session.stdin.write('h_add state%d_auto and (%s)\n' % (state, h_add_string))
            if not h_add_all and len(h_add_string) > 0:
                self.session.stdin.write('clean state%d_auto\n' % state)

            self.session.stdin.write('save qevb_org%d.pdb, state%d_auto\n' % (state, state))
            self.session.stdin.write('disable state%d\n' % state)

        time.sleep(1)
        self.session.stdin.write('set grid_mode, 1\n')
        self.session.stdin.write('set internal_gui=0 \n')

        #Create state list with dictionaries for later tracing of coordinates back to original Q-atoms:
        #[{'x y z' : Qi (s1)}, {'x y z' : Qi (s2}, ...]

        org_pdb = []
        for line in pdb:
            try:
                if int(line.split()[1]) in self.q_atom_nr.values():
                    org_pdb.append(line)
            except:
                continue

        cord_q_state = list()
        for state in range(1, self.evb_states.get() + 1):
            cord_q = dict()

            check_count = 0
            while True:
                if os.path.isfile('%s/qevb_org%d.pdb' % (self.app.workdir, state)):
                    print 'Found file'
                    break
                else:
                    check_count += 1
                    time.sleep(0.1)
                if check_count >= 50:
                    print 'Could not find file'
                    break

            q_state = open('%s/qevb_org%d.pdb' % (self.app.workdir, state), 'r').readlines()
            for line in q_state:
                atomtype = line[13:17]
                res = line[17:26]
                cord = line[26:55].strip()
                for orgline in org_pdb:
                    orgtype = orgline[13:17]
                    orgres = orgline[17:26]
                    if orgtype == atomtype and orgres == res:
                        atomnr = int(orgline.split()[1])
                        for q in self.q_atom_nr.keys():
                            if self.q_atom_nr[q] == atomnr:
                                cord_q[cord] = q

            cord_q_state.append(cord_q)

            #Rename atoms and insert occupancy and atom and write pdb files for ffld_server
            response = self.rename_atoms('%s/qevb_org%d.pdb' % (self.app.workdir, state),
                                          '%s/qevb_tmp%d.pdb' % (self.app.workdir,state), 'QS%d' % state)
            if not response:
                self.app.errorBox('Error', 'Could not translate all atomtypes. Pleas verify pdb file.')
                return
            print 'ffld_server templet qevb_tmp%d.pdb generated' % state

            #Run ffld_server and generate parameter for state
            ff = 'OPLS'+self.force_field.get()
            ffld_failed = self.run_ffld_server('qevb_tmp%d.pdb' % state, ff,
                                               '%s/prm_tmp%d.out' % (self.app.workdir, state))
            if ffld_failed:
                self.app.errorBox('Error', 'Parameter assignment with ffld_server failed!')
                return

        #Bring back original pymol structures again and get parameters
        for state in range(1, self.evb_states.get() + 1):
            self.session.stdin.write('enable state%d\n' % state)
            self.session.stdin.write('delete state%d_auto\n' % state)
            self.get_auto_prm(cord_q_state[state - 1], '%s/qevb_tmp%d.pdb' % (self.app.workdir, state),
                               '%s/prm_tmp%d.out' % (self.app.workdir, state), state - 1, only_charges)

        self.update_all()
        self.app.log('info','Parameter assignment for %d FEP states completed.' % self.evb_states.get())
        self.app.log(' ', 'Always control and verify parameters. Use with care!\n')

    def check_imp(self,q1,q2,q3,q4):
        """
        Goes through q atoms in improper, and puts them in correct order for Q and Qgui prm dictionaries.
        """
        qlist = [q1,q2,q3,q4]
        qbonds = [[],[],[],[]]
        for state in range(self.evb_states.get()):
            for i in range(len(qlist)):
                qi = qlist[i]
                for qj in self.q_bonds[qi][state]:
                    if qj in qlist:
                        if qj not in qbonds[i]:
                            qbonds[i].append(qj)

        #Find q2 and q1
        found_q2 = False
        done = False
        while not done:
            if not found_q2:
                for i in range(len(qlist)):
                    qi = qlist[i]
                    tmp_q = []
                    #Get the other 3 Q atoms
                    for j in range(len(qlist)):
                        if qlist[j] != qi:
                            tmp_q.append(qlist[j])
                    qk, ql, qm = tmp_q[0:]
                    #Find q2 (connected to q1, q3 and q4)
                    if qk in qbonds[i] and ql in qbonds[i]:
                        if qm in qbonds[i]:
                            q2 = qi
                            found_q2 = True

            #Find q1 (only connected to q2)
            if found_q2:
                del qbonds[qlist.index(q2)]
                del qlist[qlist.index(q2)]

                for i in range(len(qlist)):
                    qi = qlist[i]
                    #Get the two remaining q atoms
                    tmp_q = []
                    for j in range(len(qlist)):
                        if qlist[j] != qi:
                            tmp_q.append(qlist[j])
                    if len(tmp_q) < 2:
                        done = True
                        break
                    ql, qm = tmp_q[0:]
                    if q2 in qbonds[i]:
                        if not ql in qbonds[i] and not qm in qbonds[i]:
                            q1 = qi
                            done = True
                            break
                    if done:
                        break
        del qbonds[qlist.index(q1)]
        del qlist[qlist.index(q1)]

        print qlist
        q3, q4 = qlist[0:]

        #Check if improper exist and rearange accordingly:
        if q2 in self.q_impropers.keys():
            if q1 == self.q_impropers[q2][0]:
                #Improper exists in dict. Check that q3 and q4 are in correct order:
                if q3 == self.q_impropers[q2][-1]:
                    q3, q4 = q4, q3

        return q1, q2, q3, q4

    def get_auto_prm(self,cord_q, qpdb, qprm, state, only_charges=True):

        #Read in lines from qpdb:
        cord = list()
        with open(qpdb, 'r') as rqpdb:
            for line in rqpdb:
                if 'ATOM' in line or 'HETATM' in line:
                    cord.append(line[26:55].strip())

        #Get newnames, atomtypes, charges and vdw prm.
        found_vdw = False
        found_stretch = False
        found_bending = False
        found_torsion = False
        found_improper = False


        count = 0

        newname_q = dict()
        newname_atomtype = dict()

        type_nr = 0
        with open(qprm, 'r') as rprm:
            for line in rprm:
                if found_improper:
                    if len(line.split()) < 5:
                        found_improper = False
                    else:
                        atom1, atom2, atom3, atom4 = line.split()[0:4]
                        k = 0.5 * float(line.split()[4])
                        if atom1 in newname_q.keys() and atom2 in newname_q.keys():
                            if atom3 in newname_q.keys() and atom4 in newname_q.keys():
                                q1 = newname_q[atom1]
                                q2 = newname_q[atom2]
                                q3 = newname_q[atom3]
                                q4 = newname_q[atom4]

                                tmp_dict = {q1: atom1, q2: atom2,
                                            q3: atom3, q4: atom4}

                                q1, q2, q3, q4 = self.check_imp(q1,q2,q3,q4)
                                self.q_impropers[q2] = [q1,q3,q4]

                                atom1 = tmp_dict[q1]
                                atom2 = tmp_dict[q2]
                                atom3 = tmp_dict[q3]
                                atom4 = tmp_dict[q4]

                                quart = '%s %s %s %s' % (newname_atomtype[atom1], newname_atomtype[atom2],
                                                        newname_atomtype[atom3], newname_atomtype[atom4])
                                quart_rev = '%s %s %s %s' % (newname_atomtype[atom4], newname_atomtype[atom3],
                                                              newname_atomtype[atom2], newname_atomtype[atom1])
                                self.improper_prm[quart] = [k, 180.00]
                                self.improper_prm[quart_rev] = [k, 180.0]

                if found_torsion:
                    if len(line.split()) < 5:
                        found_torsion = False
                    else:
                        atom1, atom2, atom3, atom4 = line.split()[0:4]
                        k1, k2, k3 = map(float, line.split()[4:7])
                        k1 *= 0.5
                        k2 *= 0.5
                        k3 *= 0.5
                        if atom1 in newname_q.keys() and atom2 in newname_q.keys():
                            if atom3 in newname_q.keys() and atom4 in newname_q.keys():
                                quart = '%s %s %s %s' % (newname_atomtype[atom1], newname_atomtype[atom2],
                                                         newname_atomtype[atom3], newname_atomtype[atom4])
                                quart_rev = '%s %s %s %s' % (newname_atomtype[atom1], newname_atomtype[atom2],
                                                             newname_atomtype[atom4], newname_atomtype[atom3])
                                self.torsion_prm[quart] = [[k1, 0.0, 1], [k2, 180.0, 1], [k3, 0.0, 1]]
                                self.torsion_prm[quart_rev] = [[k1, 0.0, 1], [k2, 180.0, 1], [k3, 0.0, 1]]


                if found_bending:
                    if len(line.split()) < 5:
                        found_bending = False
                    else:
                        atom1, atom2, atom3 = line.split()[0:3]
                        if atom1 in newname_q.keys() and atom2 in newname_q.keys():
                            if atom3 in newname_q.keys():
                                trip = '%s %s %s' % (newname_atomtype[atom1], newname_atomtype[atom2],
                                                     newname_atomtype[atom3])
                                trip_rev = '%s %s %s' % (trip.split()[2], trip.split()[1], trip.split()[0])
                                k = 2.0 * float(line.split()[3])
                                r = float(line.split()[4])
                                self.angle_prm[trip] = [k, r]
                                self.angle_prm[trip_rev] = [k, r]

                if found_stretch:
                    if len(line.split()) < 5:
                        found_stretch = False
                    else:
                        atom1, atom2 = line.split()[0:2]
                        if atom1 in newname_q.keys() and atom2 in newname_q.keys():
                            #{'type1 type2' : [De, alpha, R, Kb]}
                            pair = '%s %s' % (newname_atomtype[atom1], newname_atomtype[atom2])
                            pair_rev = '%s %s' % (pair.split()[1], pair.split()[0])

                            kb = 2 * float(line.split()[2])
                            r = float(line.split()[3])
                            alpha = 2.00
                            De = kb/8.00

                            self.bond_prm[pair] = [De, alpha, r, kb]
                            self.bond_prm[pair_rev] = [De, alpha, r, kb]

                if found_vdw:
                    if len(line.split()) < 6:
                        found_vdw = False
                    elif '--------------' in line:
                        found_vdw = False
                    else:
                        newname = line.split()[0]
                        atomtype = 'Q'+line.split()[3]

                        charge = float(line.split()[4])
                        sigma = float(line.split()[5])
                        epsilon = float(line.split()[6])
                        ri, avdw2, ri1_4, ei, ei_1_4 = self.calcNBpar(sigma, epsilon)
                        ai = 1.60
                        ci = 70.0
                        if newname[0] == 'H':
                            ci = 7.0
                        mass = self.massDict[newname[0]]
                        if cord[type_nr] in cord_q.keys():
                            q = cord_q[cord[type_nr]]
                            if not only_charges:
                                newname_atomtype[newname] = atomtype
                                newname_q[newname] = q
                                self.atomtype_prm[atomtype] = [ri, ei, ci, ai, ri1_4, ei_1_4, mass]
                                self.q_atomtypes[q][state] = atomtype
                            if q in self.q_charges.keys():
                                self.q_charges[q][state] = charge

                        type_nr += 1


                if '--------------' in line:
                    count += 1
                if count == 2:
                    found_vdw = True
                if not only_charges:
                    if 'Stretch' in line:
                        found_vdw = False
                        found_stretch = True
                    if 'Bending' in line:
                        found_stretch = False
                        found_bending = True
                    if 'proper Torsion' in line:
                        found_bending = False
                        found_torsion = True
                    if 'improper Torsion' in line:
                        found_torsion = False
                        found_improper = True

    def calcNBpar(self, sigma = 0, epsilon = 0):
        """
        >> npbar(sigma, epsilon)

        Returns: avdw1,avdw2,avdw3,bvdw1,bvdw2
        """
        eps = float(epsilon)
        sig = float(sigma)

        avdw1 = np.sqrt(4*eps*sig**(12))
        avdw2 = avdw1
        avdw3 = np.sqrt(0.5) * avdw1

        bvdw1 = np.sqrt(4 * eps * sig**(6) )
        bvdw23 = np.sqrt(0.5) * bvdw1

        return (avdw1, avdw2, avdw3, bvdw1, bvdw23)

    def rename_atoms(self, original_file, new_file, resname = None, resnr = None):
        """
        Function to rename atomtypes when merging multiple residues to singular residues
        """
        atom_count = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0}
        pdb_created = True

        inputfile = open(original_file,'r').readlines()
        outputfile = open(new_file, 'w')

        if not resnr:
            resnr = inputfile[1][21:26]
        if not resname:
            resname = inputfile[1][17:21]

        for line in inputfile:
            try:
                j = 13
                found_atom = False
                while not found_atom:
                    atom = line[j]
                    if atom not in atom_count.keys():
                        if j > 16:
                            pdb_created = False
                            break
                        j += 1
                    else:
                        found_atom = True

                atom_count[atom] += 1
                outputfile.write('%s%1s%3s%s%s%s' % (line[0:13], atom, str(atom_count[atom]).ljust(3), resname, resnr,
                                                 line[26:]))
            except:
                continue

        return pdb_created

    def run_ffld_server(self, ipdb=None, ff='OPLS2005', tmpname='.tmpfile'):
        """
        Uses $SCHRODINGER/utilities/ffld_server -ipdb -print_parameters -version 11 or version 14
        """
        if not ipdb:
            self.app.errorBox('Error', 'No input pdb given for ffld_server!')
            return

        ffld_failed = False

        ffld_path = self.app.q_settings[6]

        print ff

        if ff == 'OPLS2001':
            version = '11'
        else:
            version = '14'

        tmpfile = open(tmpname, 'wb')

        call('%s/utilities/ffld_server -ipdb %s/%s -print_parameters -version %s' %
             (ffld_path, self.app.workdir, ipdb, version), shell=True, stdout=tmpfile, stderr=tmpfile)

        self.update()

        done = False

        self.app.log('info','Generating parameters ...')
        job_id = open(tmpname, 'r').readlines()
        for line in job_id:
            self.app.main_window.update_txt(line)
            if 'OPLSAA FORCE FIELD TYPE ASSIGNED' in line:
                done = True
            if 'exception' in line:
                done = True
                ffld_failed = True
            if 'FATAL' in line:
                ffld_failed = True
                done = True

        if not done:
            #Wait for ffld_server to finish if not done:
            self.app.main_window.txt.config(state=NORMAL)
            self.app.main_window.txt.insert(END, '\nWaiting for parameter assignment to finish\n')
            logline = len(self.app.main_window.txt.get(0.0, END).split('\n')) - 2
            logline = str(logline) + '.0'

        dotcount = 1
        count = 0
        while not done:
            job_id = open(self.app.workdir + '/.tmpfile','r').readlines()
            for line in job_id:
                if 'OPLSAA FORCE FIELD TYPE ASSIGNED' in line:
                    done = True
            #Update monitor so that user can see that something is happening:
            self.app.main_window.txt.delete(logline, END)
            self.app.main_window.txt.insert(END, '\nWaiting for parameter assignment to finish' + dotcount * '.' + '\n')
            self.update()
            time.sleep(0.2)
            dotcount += 1
            count += 1
            if dotcount >= 40:
                dotcount = 1
            if count > 100:
                ffld_failed = True
                break

        return ffld_failed

    def open_file(self):
        """
        opens the FEP file for editing
        """

        if not self.fep_written:
            self.write_fep()

        fepname = self.pdbfile.split('/')[-1].split('.')[0] + '.fep'
        fepname = 'inputfiles/%s' % fepname
        if not os.path.isfile(fepname):
            fepname = '%s/%s' % (self.app.workdir, fepname)
            if not os.path.isfile(fepname):
                print 'Could not find %s' % fepname
                return

        self.app.log(' ', 'Editing FEP file: \n'
                          'Saved changes will be overwritten if settings are toggled later in Setup EVB')

        self.fileEdit = FileEdit(self, fepname)
        self.fileEdit.config(bg=self.main_color)
        self.fileEdit.title('Edit FEP file')
        self.fileEdit.resizable()

    def config_md(self):
        setup_md_ = SetupMd(self,self.root, self.pdbfile, self.topology, False, fep=True)
        setup_md_.configure(background = self.main_color)
        setup_md_.title('Configure MD settings for FEP')

    def write_evb(self, feedback=True):
        """
        Write button clicked. All directories with inputfiles will be written with submit script, but not submitted.
        """
        inputfiles = '%s/inputfiles' % self.app.workdir
        if self.check_overwrite.get() == 0:
            self.overwrite = False
            if not os.path.exists(inputfiles):
                self.write_inputfiles()
            else:
                self.app.log('info','FEP inputfile templates exist and will not be overwritten.')
        else:
            self.write_inputfiles()
        self.copy_inputfiles(feedback)

    def run_evb(self):
        """
        Writes inputfiles and submits EVB run(s)
        """
        inputfiles = '%s/inputfiles' % self.app.workdir
        #Input files have not been written in current session:
        if not self.files_written:
            #Check if files should be overwritten
            if self.check_overwrite.get() == 0:
                self.overwrite = False
                if not os.path.exists(inputfiles):
                    self.app.errorBox('Warning', 'Could not find inputfiles templates. Generating new templates.')
                    self.write_evb(feedback=False)
                else:
                    self.app.log(' ','Using existing inputfiles, only changing temperature and runs')
            else:
                self.overwrite = True
                self.write_evb(feedback=False)
        else:
            os.chdir(self.app.workdir)
            tmpfile = open('.tmpfile', 'w')
            #os.system('bash runLIE.sh')
            call('bash EVB_run.sh', shell=True, stdout=tmpfile, stderr=tmpfile)
            job_id = open(self.app.workdir + '/.tmpfile','r').readlines()
            self.app.log('info','Submitting FEP jobs ...')
            for line in job_id:
                self.app.main_window.update_txt(line)

            self.app.errorBox('Info','Jobs submitted!')

    def write_inputfiles(self):
        """
        Makes inputfiles directory with FEP file, equlibration inputs and MD inputs.
        These files are templates for the subdirectory hierarchy (temperatures --> runs)
        """
        if not self.fep_written:
            self.write_fep()

        fepname = self.pdbfile.split('/')[-1].split('.')[0] + '.fep'

        q_settings = cPickle.load(open(self.app.settings_path + '/.Qsettings','rb'))


        lambda_list = list(self.lambdasteps_listbox.get(0, END))

        inputfiles = '%s/inputfiles' % self.app.workdir


        #Make inputfiles directory and move to directory:
        if not os.path.exists(inputfiles):
            os.makedirs(inputfiles)
        try:
            shutil.copy(self.pdbfile, inputfiles)
        except:
            pass
        try:
            shutil.copy(self.topology, inputfiles)
        except:
            pass

        os.chdir(inputfiles)

        #Open submission script and write head:
        self.submit = 'run.sh'
        self.submitcommand = q_settings[4][1]
        submitfile = open(self.submit,'w')
        if int(q_settings[4][0]) == 1:
            if os.path.isfile(self.app.settings_path + '/.qsubmit'):
                submissionscipt = open(self.app.settings_path + '/.qsubmit','r').readlines()
            elif os.path.isfile(self.app.workdir + '/' + '.qsubmit'):
                submissionscipt = open(self.app.workdir + '/' + '.qsubmit','r').readlines()
            else:
                submissionscipt = ['#!/bin/bash\n#Qdyn I/O\n']
                print 'submission script not found! Please edit this in settings'
            for line in submissionscipt:
                if '#Qdyn I/O' in line:
                    break
                else:
                    submitfile.write(line)

        #Check if Qdyn is MPI run or not:
        qdyn = q_settings[5][1]
        if qdyn[-1] == 'p':
            qdyn = 'mpirun %s' % qdyn

        #Find atom nrs for solute and solvent:
        pdbfile = open(self.pdbfile,'r').readlines()
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

        output_int = self.md_settings['ene_summary']
        trj_int = self.md_settings['trajectory']
        ene_int = self.md_settings['ene_file']
        non_bond_int = self.md_settings['nonbond_list']
        radial_force = self.md_settings['radial_force']
        polarisation = self.md_settings['pol_force']
        shell_force = self.md_settings['shell_force']
        shell_radius = self.md_settings['shell_rad']
        md_steps = int(round((float(self.md_settings['simtime']) * 1000000.00)/ float(self.md_settings['stepsize'])))
        md_stepsize = self.md_settings['stepsize']
        md_temp = 'T_VAR'
        md_bath = self.md_settings['bath_coupling']
        shake_solvent = on_off[self.md_settings['shake_solvent']]
        shake_solute = on_off[self.md_settings['shake_solute']]
        shake_hydrogens = on_off[self.md_settings['shake_hydrogens']]
        lrf = on_off[self.md_settings['lrf']]
        lrf_cutoff = self.md_settings['lrf_cut']
        use_pol = on_off[self.md_settings['polarisation']]

        #Write default equilibration procedure
        #Leave base_name as md (names are long enough with the lambda values included)
        base_name = 'eq'
        random_seed = 'SEED_VAR'
        count = 0
        for i in range(len(q_settings[3])):
            count += 1
            inputfile =  base_name + '%d.inp' % count
            logfile =  base_name + '%d.log' % count
            eq_file = open(inputfile, 'w')

            if q_settings[3][i][0] == 'End':
                temp = md_temp
            else:
                temp = q_settings[3][i][0]

            if int(q_settings[4][0]) == 1:
                inputfile = base_name + '%d.inp' % count
                logfile = base_name + '%d.log' % count

            submitfile.write('%s %s > %s\n' % (qdyn, inputfile, logfile))

            eq_file.write('[MD]\n')
            eq_file.write('%25s %s\n' % ('steps'.ljust(25), q_settings[3][i][5]))
            eq_file.write('%25s %s\n' % ('stepsize'.ljust(25), q_settings[3][i][4]))
            eq_file.write('%25s %s\n' % ('temperature'.ljust(25), temp))
            eq_file.write('%25s %s\n' % ('bath_coupling'.ljust(25), q_settings[3][i][1]))
            if count == 1:
                eq_file.write('%25s %s\n' % ('random_seed'.ljust(25), random_seed))
                eq_file.write('%25s %s\n' % ('initial_temperature'.ljust(25), q_settings[3][i][0]))
                eq_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), 'on'))
            if count > 1:
                eq_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), shake_solvent))

            eq_file.write('%25s %s\n' % ('shake_hydrogens'.ljust(25), shake_hydrogens))
            eq_file.write('%25s %s\n' % ('shake_solute'.ljust(25), shake_solute))
            eq_file.write('%25s %s\n' % ('lrf'.ljust(25), lrf))

            eq_file.write('\n[cut-offs]\n')
            eq_file.write('%25s %s\n' % ('solute_solvent'.ljust(25), self.md_settings['solute_solvent_cut']))
            eq_file.write('%25s %s\n' % ('solute_solute'.ljust(25), self.md_settings['solute_solute_cut']))
            eq_file.write('%25s %s\n' % ('solvent_solvent'.ljust(25), self.md_settings['solvent_solvent_cut']))
            eq_file.write('%25s %s\n' % ('q_atom'.ljust(25), self.md_settings['q_atoms_cut']))
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

            eq_file.write('\n[files]\n')
            if int(q_settings[4][0]) == 1:
                topology = self.topology.split('/')[-1]
                fepname = fepname.split('/')[-1]
            else:
                topology = self.topology
            eq_file.write('%25s %s\n' % ('topology'.ljust(25), topology))
            eq_file.write('%25s %s%d.dcd\n' % ('trajectory'.ljust(25), base_name, count))
            if count != 1:
                eq_file.write('%25s %s\n' % ('restart'.ljust(25), self.restart_file))
            eq_file.write('%25s %s%d.re\n' % ('final'.ljust(25), base_name, count))
            eq_file.write('%25s %s\n' % ('fep'.ljust(25), fepname))
            self.restart_file = logfile.split('.')[0]+'.re'

            eq_file.write('\n[trajectory_atoms]\n')
            eq_file.write('%s\n' % self.md_settings['trajectory atoms'])

            eq_file.write('\n[lambdas]\n')
            eq_file.write('%s\n' % ' '.join(lambda_list[0].split()[1:]))

            eq_file.write('\n[sequence_restraints]\n')
            force = float(q_settings[3][i][3])
            if q_settings[3][i][2] != 'None':
                if q_settings[3][i][2] == 'All':
                    atomlist = all_atoms
                    #eq_file.write(' not excluded  %4.1f 0  0\n' % force)
                elif q_settings[3][i][2] == 'Solute':
                    atomlist = solute
                elif q_settings[3][i][2] == 'Solvent':
                    atomlist = solvent
                try:
                    atom_i = atomlist[0]
                    atom_j = atomlist[-1]
                    eq_file.write('%6s %6s %4.1f 0  0\n' % (atom_i, atom_j, force))
                except:
                    continue

            if len(self.md_settings['seq_rest']) > 0:
                for restraint in self.md_settings['seq_rest']:
                    eq_file.write('%s\n' % restraint)

            if len(self.md_settings['dist_rest']) > 0:
                eq_file.write('\n[distance_restraints]\n')
                for restraint in self.md_settings['dist_rest']:
                    eq_file.write('%s\n' % restraint)

            if len(self.md_settings['atom_rest']) > 0:
                eq_file.write('\n[atom_restraints]\n')
                for restraint in self.md_settings['atom_rest']:
                    eq_file.write('%s\n' % restraint)

            if len(self.md_settings['wall_rest']) > 0:
                eq_file.write('\n[wall_restraints]\n')
                for restraint in self.md_settings['wall_rest']:
                    eq_file.write('%s\n' % restraint)

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
            md_file = open(inputfile, 'w')

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
            md_file.write('%25s %s\n' % ('solute_solvent'.ljust(25), self.md_settings['solute_solvent_cut']))
            md_file.write('%25s %s\n' % ('solute_solute'.ljust(25), self.md_settings['solute_solute_cut']))
            md_file.write('%25s %s\n' % ('solvent_solvent'.ljust(25), self.md_settings['solvent_solvent_cut']))
            md_file.write('%25s %s\n' % ('q_atom'.ljust(25), self.md_settings['q_atoms_cut']))
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
            if int(q_settings[4][0]) == 1:
                topology = self.topology.split('/')[-1]
                fepname = fepname.split('/')[-1]
            else:
                topology = self.topology
            md_file.write('%25s %s\n' % ('topology'.ljust(25), topology))
            md_file.write('%25s %s.dcd\n' % ('trajectory'.ljust(25), base_name))
            md_file.write('%25s %s\n' % ('restart'.ljust(25), self.restart_file))
            md_file.write('%25s %s.en\n' % ('energy'.ljust(25), base_name))
            md_file.write('%25s %s.re\n' % ('final'.ljust(25), base_name))
            md_file.write('%25s %s\n' % ('fep'.ljust(25), fepname))
            self.restart_file = logfile.split('.')[0]+'.re'

            md_file.write('\n[trajectory_atoms]\n')
            md_file.write('%s\n' % self.md_settings['trajectory atoms'])

            md_file.write('\n[lambdas]\n')
            md_file.write('%s\n' % lambda_value)

            if len(self.md_settings['seq_rest']) > 0:
                md_file.write('\n[sequence_restraints]\n')
                for restraint in self.md_settings['seq_rest']:
                    md_file.write('%s\n' % restraint)

            if len(self.md_settings['dist_rest']) > 0:
                md_file.write('\n[distance_restraints]\n')
                for restraint in self.md_settings['dist_rest']:
                    md_file.write('%s\n' % restraint)

            if len(self.md_settings['atom_rest']) > 0:
                md_file.write('\n[atom_restraints]\n')
                for restraint in self.md_settings['atom_rest']:
                    md_file.write('%s\n' % restraint)

            if len(self.md_settings['wall_rest']) > 0:
                md_file.write('\n[wall_restraints]\n')
                for restraint in self.md_settings['wall_rest']:
                    md_file.write('%s\n' % restraint)
            md_file.close()

        #If use submission script, check for end statements (comes after #Qdyn I/O):
        if int(q_settings[4][0]) == 1:
            write_end = False
            submissionscipt = open(self.app.settings_path + '/.qsubmit','r').readlines()
            for k in range(len(submissionscipt)):
                if '#Qdyn I/O' in submissionscipt[k]:
                    end_statements_start = k + 1
                    write_end = True
            if write_end:
                for line in range(end_statements_start, len(submissionscipt)):
                    submitfile.write(submissionscipt[line])
        submitfile.close()

        os.chdir(self.app.workdir)
        self.app.log('info', 'FEP templates written to "/inputfiles"')

    def copy_inputfiles(self, feedback=False):
        """
        Copies md files from inputfiles and changes random seed and temperature
        If feedback=False, jobs will be submitted on the fly.
        """
        self.app.log('info', 'Generating FEP inputfiles from templates ...')
        self.update()

        submitfile = open(self.app.workdir+'/FEP_run.sh','w')
        submitfile.write('#!/bin/bash\n')

        #If no feedback, jobs are to be submitted
        if not feedback:
            os.chdir(self.app.workdir)


        if not os.path.exists('inputfiles'):
            self.app.errorBox('Error','Directory "inputfiles" not found. Please make sure that you are in correct'
                                      ' workdir')
            return

        #Get list of template files to copy
        evb_files = os.listdir('inputfiles')
        for i in range(len(evb_files)):
            if evb_files[i].startswith('.'):
                del evb_files[i]

        temp_runs = list(self.temp_listbox.get(0, END))
        for i in temp_runs:
            T = i.split()[0]
            runs = i.split()[1]
            if not os.path.isdir(T):
                os.mkdir(T)

            start = 1
            end = int(runs)

            #If not overwrite, find latest subrun in T dir:
            if not self.overwrite:
                subdirs = os.listdir(T)
                last_run = 0
                for sub in subdirs:
                    if sub.isdigit():
                        if int(sub) > last_run:
                            last_run = int(sub)

                start = last_run + 1
                end += last_run
            submitfile.write('for i in {%d..%d}\ndo\ncd %s/$i\n%s run.sh\n' % (start, end, T,self.app.q_settings[4][1]))
            submitfile.write('cd ../../\ndone')
            for run in range(start, end + 1):
                print 'T=%s   run=%d' % (T, run)
                if not os.path.exists('%s/%d' % (T, run)):
                    os.makedirs('%s/%d' % (T, run))
                for inputfile in evb_files:
                    shutil.copy('inputfiles/' + inputfile, '%s/%d' % (T, run))

                #Change T_VAR and SEED_VAR in files:
                for inputfile in evb_files:
                    file_ = open('%s/%d/%s' % (T, run, inputfile), 'r').readlines()
                    newfile = open('%s/%d/%s' % (T, run, inputfile), 'w')
                    for line in file_:
                        if 'SEED_VAR' not in line and 'T_VAR' not in line:
                            newfile.write(line)
                        elif 'SEED_VAR' in line:
                            random_seed = random.randrange(1000, 9999)
                            newfile.write('%25s %s\n' % ('random_seed'.ljust(25), random_seed))
                        elif 'T_VAR' in line:
                            newfile.write('%25s %s\n' % ('temperature'.ljust(25), T))

                    newfile.close()

                #If feedback=False, jobs are submitted as they are created
                if not feedback:
                    os.chdir('%s/%d' % (T, run))
                    tmpfile = open('.tmpfile', 'w')
                    call('bash run.sh', shell=True, stdout=tmpfile, stderr=tmpfile)
                    job_id = open(self.app.workdir + '/.tmpfile','r').readlines()

                    for line in job_id:
                        self.app.main_window.update_txt(line)
                    tmpfile.close()
                    os.chdir(self.app.workdir)

        submitfile.close()
        if feedback:
            self.app.errorBox('Info','FEP inputfiles written. Use FEP_run.sh to submit jobs from command line.')
        self.app.log('info', 'FEP inputfiles written.')

        self.files_written = True

    def add_monitor_group(self):
        if len(self.monitor_groups) == 0:
            self.monitor_groups[1] = list()
        else:
            self.monitor_groups[max(self.monitor_groups.keys()) + 1] = list()

        self.update_monitor_listboxes()

    def del_monitor_group(self):
        selections = map(int, self.groups_listbox.curselection())

        if len(selections) < 1:
            return

        for selected in selections:
            group = int(self.groups_listbox.get(selected).split()[1])
            del self.monitor_groups[group]
            for i in range(len(self.monitor_pairs)):
                if str(group) in self.monitor_pairs[i]:
                    del self.monitor_pairs[i]

        self.update_monitor_listboxes()

    def add_monitor_atoms(self):
        sel_group = map(int, self.groups_listbox.curselection())

        if len(sel_group) != 1:
            print 'Select exactly one group to append monitor atoms to!'
            return

        group = int(self.groups_listbox.get(sel_group[0]).split()[1])

        sel_atoms = map(int, self.qatoms_listbox.curselection())

        if len(sel_atoms) < 1:
            print 'No atoms selected!'
            return

        for i in sel_atoms:
            atom = self.qatoms_listbox.get(i).split()[1]
            if atom not in self.monitor_groups[group]:
                self.monitor_groups[group].append(atom)

        self.update_monitor_listboxes()

    def del_monitor_atoms(self):
        sel_group = map(int, self.groups_listbox.curselection())

        if len(sel_group) != 1:
            print 'Select exactly one group to remove monitor atoms from!'
            return

        group = int(self.groups_listbox.get(sel_group[0]).split()[1])

        sel_atoms = map(int, self.group_atoms_listbox.curselection())

        for i in sel_atoms:
            atom = self.group_atoms_listbox.get(i).strip()
            del self.monitor_groups[group][self.monitor_groups[group].index(atom)]

        self.update_monitor_listboxes()

    def add_monitor_pair(self):
        sel_groups = map(int, self.groups_listbox.curselection())

        if len(sel_groups) != 2:
            print 'Select exactly to groups to monitor!'
            return

        groups = list()

        for i in sel_groups:
            group = self.groups_listbox.get(i).split()[1]
            groups.append(group)

        if groups not in self.monitor_pairs:
            self.monitor_pairs.append(groups)
        self.update_monitor_listboxes()

    def del_monitor_pair(self):
        sel_pairs = map(int, self.monitor_pairs_listbox.curselection())

        if len(sel_pairs) < 1:
            return

        for i in sel_pairs:
            atom1, atom2 = self.monitor_pairs_listbox.get(i).split()[0:2]
            for j in range(len(self.monitor_pairs)):
                pair = self.monitor_pairs[j]
                if atom1 in pair and atom2 in pair:
                    del self.monitor_pairs[j]

        self.update_monitor_listboxes()

    def update_monitor_listboxes(self):
        self.groups_listbox.delete(0, END)
        self.group_atoms_listbox.delete(0, END)
        self.monitor_pairs_listbox.delete(0, END)

        for group in sorted(self.monitor_groups.keys()):
            self.groups_listbox.insert(END, 'Group %d' % group)
            for atom in self.monitor_groups[group]:
                self.group_atoms_listbox.insert(END, '%6s' % atom)

        for pair in self.monitor_pairs:
            self.monitor_pairs_listbox.insert(END,'%2s  %2s' % (pair[0], pair[1]))

    def show_var_frame(self, *args):
        frames = {'Charges': self.charge_frame,
                  'Atomtypes': self.atomtype_frame,
                  'Atoms': self.atoms_frame,
                  'Bonds': self.bond_frame,
                  'Angles': self.angle_frame,
                  'Torsions': self.torsion_frame,
                  'Impropers': self.improper_frame,
                  'Monitor groups': self.monitorgroups_frame,
                  #'Couplings': self.coupling_frame,
                  #'Soft pairs': self.softpair_frame,
                  #'El-scale': self.elscale_frame,
                  #'Excluded pairs': self.excludedpairs_frame,
                  'Softcore': self.softcore_frame,
                   u"\N{GREEK SMALL LETTER LAMDA}-steps": self.lambdasteps_frame,
                   'Temperature(s)': self.temperature_frame}

        for i in frames.keys():
            frames[i].grid_forget()
        try:
            frames[self.setup_evb.get()].grid(row=2, column=0, columnspan=2)
        except:
            pass

    def dialog_window(self):

        self.title('FEP setup')

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Qatoms frame
        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        #Change bonds frame
        frame2 = Frame(mainframe, bg=self.main_color)
        frame2.grid(row=0, column=1, columnspan=2, padx=(10,0), sticky='n')

        #Status frame
        frame3 = Frame(mainframe, bg=self.main_color)
        frame3.grid(row=1, column=0, columnspan=2)

        #Variable frames
        self.charge_frame = Frame(mainframe, bg=self.main_color)
        self.atomtype_frame = Frame(mainframe, bg=self.main_color)
        self.atoms_frame = Frame(mainframe, bg=self.main_color)
        self.bond_frame = Frame(mainframe, bg=self.main_color)
        self.angle_frame = Frame(mainframe, bg=self.main_color)
        self.torsion_frame = Frame(mainframe, bg=self.main_color)
        self.improper_frame = Frame(mainframe, bg=self.main_color)
        self.coupling_frame = Frame(mainframe, bg=self.main_color)
        self.softpair_frame = Frame(mainframe, bg=self.main_color)
        self.elscale_frame = Frame(mainframe, bg=self.main_color)
        self.excludedpairs_frame = Frame(mainframe, bg=self.main_color)
        self.softcore_frame = Frame(mainframe, bg=self.main_color)
        self.lambdasteps_frame = Frame(mainframe, bg=self.main_color)
        self.temperature_frame = Frame(mainframe, bg=self.main_color)
        self.monitorgroups_frame = Frame(mainframe, bg=self.main_color)

        #Save/close frame
        frame4 = Frame(mainframe, bg=self.main_color)
        frame4.grid(row=3, column=0, columnspan=2, pady=(20, 0))


        #Q-ATOMS
        qatoms_label = LabelFrame(frame1, text='Q-atoms', bg=self.main_color)
        qatoms_label.pack(fill=BOTH, expand=1)

        q_label = Label(frame1, text='Qi   Atom    !Note', bg=self.main_color)
        q_label.grid(in_=qatoms_label, row=0, column=0, sticky='w')

        qatoms_yscroll = Scrollbar(frame1)
        qatoms_yscroll.grid(in_=qatoms_label, row = 1, rowspan=10, column = 1, sticky = 'nsw', padx=(0,10))
        qatoms_xscroll = Scrollbar(frame1, orient=HORIZONTAL)
        qatoms_xscroll.grid(in_=qatoms_label, row=11, column=0, sticky='we')

        self.qatoms_listbox = Listbox(frame1, yscrollcommand = qatoms_yscroll.set,
                                      xscrollcommand=qatoms_xscroll.set,
                                      width=20, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        qatoms_yscroll.config(command=self.qatoms_listbox.yview)
        qatoms_xscroll.config(command=self.qatoms_listbox.xview)
        self.qatoms_listbox.grid(in_=qatoms_label, row=1, rowspan=10, column = 0, sticky='e')
        self.qatoms_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.qatoms_listbox.bind('<<ListboxSelect>>', self.list_q_atoms_event)


        add_qatoms = Button(frame1, text='+', highlightbackground=self.main_color, command = self.add_qatoms)
        add_qatoms.grid(in_=qatoms_label, row=3, column=2)

        del_qatoms = Button(frame1, text='-', highlightbackground=self.main_color, command = self.del_qatoms)
        del_qatoms.grid(in_=qatoms_label, row=3, column=3)

        note_qatoms = Button(frame1, text='Note', highlightbackground=self.main_color, command=self.note_qatoms)
        note_qatoms.grid(in_=qatoms_label, row=4, column=2, columnspan=2)

        self.optimize_button = Button(frame1, text='Fix geom.', highlightbackground=self.main_color, command=self.geom_opt)
        self.optimize_button.grid(in_=qatoms_label, row=5, column=2, columnspan=2)
        self.optimize_button.config(state=DISABLED)

        #CHANGE BONDS
        change_bonds = LabelFrame(frame2, text='Remove atoms', bg=self.main_color)
        change_bonds.grid(sticky='ns')

        bond_label = Label(frame2, text=u" Qi     \N{GREEK SMALL LETTER PHI}1   \N{GREEK SMALL LETTER PHI}2"
                                        u"   \N{GREEK SMALL LETTER PHI}3   \N{GREEK SMALL LETTER PHI}4",
                           bg=self.main_color)
        bond_label.grid(in_=change_bonds, row=0, column=0, sticky='w')

        bondchange_yscroll = Scrollbar(frame2)
        bondchange_yscroll.grid(in_=change_bonds, row = 1, rowspan=5, column = 1, sticky = 'nsw', padx=(0,10))
        self.fepatoms_listbox = Listbox(frame2, yscrollcommand = qatoms_yscroll.set,
                                      width=25, height=5, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        bondchange_yscroll.config(command=self.fepatoms_listbox.yview)
        self.fepatoms_listbox.grid(in_=change_bonds, row=1, rowspan=5, column = 0, sticky = 'e')
        self.fepatoms_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.fepatoms_listbox.bind('<<ListboxSelect>>', self.list_fepatoms_event)

        add_bondchange = Button(frame2, text='+', highlightbackground=self.main_color, command=self.add_fep_atom)
        add_bondchange.grid(in_=change_bonds, row=0, column=2)

        del_bondchange = Button(frame2, text='-', highlightbackground=self.main_color, command=self.del_fep_atom)
        del_bondchange.grid(in_=change_bonds, row=0, column=3)

        alter_bond1 = Button(frame2, text=u"\N{GREEK SMALL LETTER PHI}1", highlightbackground=self.main_color,
                              command = lambda: self.set_fep_state(0), width=2)
        alter_bond1.grid(in_=change_bonds,row=1, column=2)

        alter_bond2 = Button(frame2, text=u"\N{GREEK SMALL LETTER PHI}2", highlightbackground=self.main_color,
                              command = lambda: self.set_fep_state(1), width=2)
        alter_bond2.grid(in_=change_bonds, row=1, column=3)

        self.alter_bond3 = Button(frame2, text=u"\N{GREEK SMALL LETTER PHI}3", highlightbackground=self.main_color,
                              command = lambda: self.set_fep_state(2), width=2)
        self.alter_bond3.grid(in_=change_bonds, row=2, column=2)
        self.alter_bond3.config(state=DISABLED)
        self.state3_enabled.append(self.alter_bond3)

        self.alter_bond4 = Button(frame2, text=u"\N{GREEK SMALL LETTER PHI}4", highlightbackground=self.main_color,
                              command = lambda: self.set_fep_state(3), width=2)
        self.alter_bond4.grid(in_=change_bonds, row=2, column=3)
        self.alter_bond4.config(state=DISABLED)
        self.state4_enabled.append(self.alter_bond4)


        #Setup
        setup_label = LabelFrame(frame2, text='Setup', bg=self.main_color)
        setup_label.grid()

        set_evb_states = OptionMenu(frame2, self.set_evb_states, '2 states FEP', '3 states FEP', '4 states FEP')
        set_evb_states.config(bg=self.main_color, highlightbackground=self.main_color, width=15)
        set_evb_states.grid(in_=setup_label, row=3, column=0)
        self.set_evb_states.set('2 states FEP')

        sync_label = Label(frame2, text='Sync pyMol:', bg=self.main_color)
        sync_label.grid(in_=setup_label, row=4, column=1, sticky='w')

        self.sync_check = Checkbutton(frame2, bg=self.main_color, variable=self.sync_pymol, command=self.toggle_pymol)
        self.sync_check.grid(in_=setup_label, row=4, column=1, sticky='e')

        selH = Label(frame2, text='selH:', bg=self.main_color)
        selH.grid(in_=setup_label, row=4, column=2, sticky='w')

        self.selH_check = Checkbutton(frame2, bg=self.main_color, variable=self.select_h)
        self.selH_check.grid(in_=setup_label, row=4, column=2, sticky='e')
        self.selH_check.config(state=DISABLED)

        self.setup_menu = OptionMenu(frame2, self.setup_evb,
                                   'Charges',
                                   'Atoms',
                                   'Atomtypes',
                                   'Bonds',
                                   'Angles',
                                   'Torsions',
                                    'Impropers',
                                   # 'Couplings',
                                   #'Soft pairs',
                                   #'El-scale',
                                   #'Excluded pairs',
                                   'Monitor groups',
                                   'Softcore',
                                   u"\N{GREEK SMALL LETTER LAMDA}-steps",
                                   'Temperature(s)')
        self.setup_menu.config(highlightbackground=self.main_color, bg=self.main_color, width=15)
        self.setup_menu.grid(in_=setup_label, row=4, column=0)

        self.auto_setup = Button(frame2, text='Auto assign', highlightbackground=self.main_color,
                                 command=lambda: self.auto_evb(only_charges=False))
        self.auto_setup.grid(in_=setup_label, row=3, column=1)
        self.auto_setup.config(state=DISABLED)

        self.force_field = Spinbox(frame2, width=5, highlightthickness=0, relief=GROOVE, values=self.forcefields)
        self.force_field.grid(in_=setup_label, row=3, column=2)
        #self.force_field.config(font=tkFont.Font(family="Courier", size=10))
        self.force_field.config(state=DISABLED)

        #Status frame
        status_label = LabelFrame(frame3, text='Status', bg=self.main_color)
        status_label.grid()

        qstatus_yscroll = Scrollbar(frame3)
        qstatus_yscroll.grid(in_=status_label, row = 1, rowspan=3, column = 1, sticky = 'nsw', padx=(0,10))
        self.qstatus_listbox = Listbox(frame3, yscrollcommand = qstatus_yscroll.set,
                                      width=77, height=3, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED)
        qstatus_yscroll.config(command=self.qstatus_listbox.yview)
        self.qstatus_listbox.grid(in_=status_label, row=1, rowspan=3, column = 0, sticky = 'e')
        self.qstatus_listbox.config(font=tkFont.Font(family="Courier", size=12))

        #Variable frames from dropdown menu (.grid_remove)
        ######### Charges ##########
        charge_label = Label(self.charge_frame, text=u" Qi       \N{GREEK SMALL LETTER PHI}1    "
                                                     u"    \N{GREEK SMALL LETTER PHI}2"
                                                     u"       \N{GREEK SMALL LETTER PHI}3"
                                                     u"       \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        charge_label.grid(row=0, column=0, sticky='w')

        charge_yscroll = Scrollbar(self.charge_frame)
        charge_yscroll.grid(row = 1, rowspan=10, column = 1, sticky = 'nsw', padx=(0,10))
        self.charge_listbox = Listbox(self.charge_frame, yscrollcommand = charge_yscroll.set,
                                      width=55, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        charge_yscroll.config(command=self.charge_listbox.yview)
        self.charge_listbox.grid(row=1, rowspan=10, column = 0, sticky = 'e')
        self.charge_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.charge_listbox.bind('<<ListboxSelect>>', self.list_charge_event)

        self.autocharge = Button(self.charge_frame, text='Auto Charge', highlightbackground=self.main_color,
                            command=lambda: self.auto_evb(only_charges=True))
        self.autocharge.grid(row=0, column=3, columnspan=2)
        self.autocharge.config(state=DISABLED)

        s1 = Label(self.charge_frame, text=u"\N{GREEK SMALL LETTER PHI}1", bg=self.main_color)
        s1.grid(row=3, column=3)

        s2 = Label(self.charge_frame, text=u"\N{GREEK SMALL LETTER PHI}2", bg=self.main_color)
        s2.grid(row=3, column=4)

        self.charge1 = Spinbox(self.charge_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-3.00, to=3.00, increment=0.01)
        self.charge1.grid(row=4, column=3)
        self.charge1.delete(0, END)
        self.charge1.insert(0, '')

        self.charge2 = Spinbox(self.charge_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-3.00, to=3.00, increment=0.01)
        self.charge2.grid(row=4, column=4)
        self.charge2.delete(0, END)
        self.charge2.insert(0, '')

        self.charge3 = Spinbox(self.charge_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-3.00, to=3.00, increment=0.01)
        self.charge3.grid(row=6, column=3)
        self.charge3.delete(0, END)
        self.charge3.insert(0, '')
        self.charge3.config(state=DISABLED)
        self.state3_enabled.append(self.charge3)

        self.charge4 = Spinbox(self.charge_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-3.00, to=3.00, increment=0.01)
        self.charge4.grid(row=6, column=4)
        self.charge4.delete(0, END)
        self.charge4.insert(0, '')
        self.charge4.config(state=DISABLED)
        self.state4_enabled.append(self.charge4)

        s3 = Label(self.charge_frame, text=u"\N{GREEK SMALL LETTER PHI}3", bg=self.main_color)
        s3.grid(row=5, column=3)

        s4 = Label(self.charge_frame, text=u"\N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        s4.grid(row=5, column=4)

        edit_charges = Button(self.charge_frame, text='Change charge', highlightbackground=self.main_color,
                              command=self.edit_charge)
        edit_charges.grid(row=7, column=3, columnspan=2)

        ######## Atomtypes #########
        atomtypes_label = Label(self.atomtype_frame, text="Name      Ri       Ei        Ci      ai    Ri(1-4)  "
                                                          "Ei(1-4)  Mass", bg=self.main_color)
        atomtypes_label.grid(row=0, column=0, columnspan=6, sticky='w')

        atomtypes_yscroll = Scrollbar(self.atomtype_frame)
        atomtypes_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.atomtypes_listbox = Listbox(self.atomtype_frame, yscrollcommand = atomtypes_yscroll.set,
                                      width=50, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        atomtypes_yscroll.config(command=self.atomtypes_listbox.yview)
        self.atomtypes_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.atomtypes_listbox.config(font=tkFont.Font(family="Courier", size=12))
        #self.atomtypes_listbox.bind('<<ListboxSelect>>', self.list_charge_event)

        import_atom_prm = Button(self.atomtype_frame, text='Import', highlightbackground=self.main_color,
                            command=self.import_atom_prm)
        import_atom_prm.grid(row=11, column=1)

        edit_atom_prm = Button(self.atomtype_frame, text='Edit', highlightbackground=self.main_color,
                               command=self.edit_atom_prm)
        edit_atom_prm.grid(row=11, column=2)

        add_atom_prm = Button(self.atomtype_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_atom_prm)
        add_atom_prm.grid(row=11, column=3)

        del_atom_prm = Button(self.atomtype_frame, text='-', highlightbackground=self.main_color,
                              command=self.delete_atom_prm)
        del_atom_prm.grid(row=11, column=4)

        changetypes_label = Label(self.atomtype_frame, text=u"Qi  \N{GREEK SMALL LETTER PHI}1"
                                                            u"    \N{GREEK SMALL LETTER PHI}2"
                                                            u"    \N{GREEK SMALL LETTER PHI}3"
                                                            u"    \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        changetypes_label.grid(row=0, column=7, columnspan=4, sticky='w')

        changetypes_yscroll = Scrollbar(self.atomtype_frame)
        changetypes_yscroll.grid(row = 1, rowspan=10, column = 11, sticky = 'nsw', padx=(0,10))
        self.changetypes_listbox = Listbox(self.atomtype_frame, yscrollcommand = changetypes_yscroll.set,
                                      width=27, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        changetypes_yscroll.config(command=self.changetypes_listbox.yview)
        self.changetypes_listbox.grid(row=1, rowspan=10, column = 7, columnspan=4, sticky = 'e')
        self.changetypes_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.changetypes_listbox.bind('<<ListboxSelect>>', self.list_changetypes_event)

        alter_state1 = Button(self.atomtype_frame, text=u"\N{GREEK SMALL LETTER PHI}1", highlightbackground=self.main_color,
                              command = lambda: self.set_atomtype_state(0))
        alter_state1.grid(row=11, column=7)

        alter_state2 = Button(self.atomtype_frame, text=u"\N{GREEK SMALL LETTER PHI}2", highlightbackground=self.main_color,
                              command = lambda: self.set_atomtype_state(1))
        alter_state2.grid(row=11, column=8)

        self.alter_state3 = Button(self.atomtype_frame, text=u"\N{GREEK SMALL LETTER PHI}3", highlightbackground=self.main_color,
                              command = lambda: self.set_atomtype_state(2))
        self.alter_state3.grid(row=11, column=9)
        self.alter_state3.config(state=DISABLED)
        self.state3_enabled.append(self.alter_state3)

        self.alter_state4 = Button(self.atomtype_frame, text=u"\N{GREEK SMALL LETTER PHI}4", highlightbackground=self.main_color,
                              command = lambda: self.set_atomtype_state(3))
        self.alter_state4.grid(row=11, column=10)
        self.alter_state4.config(state=DISABLED)
        self.state4_enabled.append(self.alter_state4)


        ####### Bond types ########
        bondtypes_label = Label(self.bond_frame, text=u"#      De       \N{GREEK SMALL LETTER ALPHA}        R0      "
                                , bg=self.main_color)
        bondtypes_label.grid(row=0, column=0, columnspan=6, sticky='w')

        bondtypes_yscroll = Scrollbar(self.bond_frame)
        bondtypes_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.bondtypes_listbox = Listbox(self.bond_frame, yscrollcommand = bondtypes_yscroll.set,
                                      width=35, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        bondtypes_yscroll.config(command=self.bondtypes_listbox.yview)
        self.bondtypes_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.bondtypes_listbox.config(font=tkFont.Font(family="Courier", size=12))

        import_bond_prm = Button(self.bond_frame, text='Import', highlightbackground=self.main_color,
                            command=self.import_bond_prm)
        import_bond_prm.grid(row=11, column=1)

        edit_bond_prm = Button(self.bond_frame, text='Edit', highlightbackground=self.main_color,
                               command=self.edit_bond_prm)
        edit_bond_prm.grid(row=11, column=2)

        add_bond_prm = Button(self.bond_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_bond_prm)
        add_bond_prm.grid(row=11, column=3)

        changetypes_label = Label(self.bond_frame, text=u"     i            j       \N{GREEK SMALL LETTER PHI}1"
                                                            u"  \N{GREEK SMALL LETTER PHI}2"
                                                            u"  \N{GREEK SMALL LETTER PHI}3"
                                                            u"  \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        changetypes_label.grid(row=0, column=7, sticky='w')

        changebond_yscroll = Scrollbar(self.bond_frame)
        changebond_yscroll.grid(row = 1, rowspan=10, column = 8, sticky = 'nsw', padx=(0,10))
        self.changebond_listbox = Listbox(self.bond_frame, yscrollcommand = changebond_yscroll.set,
                                      width=35, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        changebond_yscroll.config(command=self.changebond_listbox.yview)
        self.changebond_listbox.grid(row=1, rowspan=10, column = 7, sticky = 'e')
        self.changebond_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.changebond_listbox.bind('<<ListboxSelect>>', self.list_changebond_event)

        #update_bonds = Button(self.bond_frame, text='Update', highlightbackground=self.main_color,
        #                      command = self.update_status)
        #update_bonds.grid(row=11, column=7)

        #### ANGLE TYPES ####
        angle_label = Label(self.angle_frame, text=u"#      K          \N{GREEK CAPITAL LETTER THETA}  ",
                            bg=self.main_color)
        angle_label.grid(row=0, column=0, columnspan=6, sticky='w')

        angletypes_yscroll = Scrollbar(self.angle_frame)
        angletypes_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.angletypes_listbox = Listbox(self.angle_frame, yscrollcommand = angletypes_yscroll.set,
                                      width=30, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        angletypes_yscroll.config(command=self.angletypes_listbox.yview)
        self.angletypes_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.angletypes_listbox.config(font=tkFont.Font(family="Courier", size=12))

        import_angle_prm = Button(self.angle_frame, text='Import', highlightbackground=self.main_color,
                            command=self.import_angle_prm)
        import_angle_prm.grid(row=11, column=1)

        edit_angle_prm = Button(self.angle_frame, text='Edit', highlightbackground=self.main_color,
                               command=self.edit_angle_prm)
        edit_angle_prm.grid(row=11, column=2)

        add_angle_prm = Button(self.angle_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_angle_prm)
        add_angle_prm.grid(row=11, column=3)

        changetypes_label = Label(self.angle_frame, text=u"   i             j           k     "
                                                         u"\N{GREEK SMALL LETTER PHI}1"
                                                            u"   \N{GREEK SMALL LETTER PHI}2"
                                                            u"   \N{GREEK SMALL LETTER PHI}3"
                                                            u"   \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        changetypes_label.grid(row=0, column=7, sticky='w')

        changeangle_yscroll = Scrollbar(self.angle_frame)
        changeangle_yscroll.grid(row = 1, rowspan=10, column = 8, sticky = 'nsw', padx=(0,10))
        self.changeangle_listbox = Listbox(self.angle_frame, yscrollcommand = changeangle_yscroll.set,
                                      width=40, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        changeangle_yscroll.config(command=self.changeangle_listbox.yview)
        self.changeangle_listbox.grid(row=1, rowspan=10, column = 7, sticky = 'e')
        self.changeangle_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.changeangle_listbox.bind('<<ListboxSelect>>', self.list_changeangle_event)

        include_ang = Label(self.angle_frame, text='Include changing angles:', bg=self.main_color)
        include_ang.grid(row=11, column=7)

        all_angles_check = Checkbutton(self.angle_frame, bg=self.main_color, variable=self.all_angles, command=self.toggle_show_ang)
        all_angles_check.grid(row=11, column=7, sticky='e')

        #### TORSION TYPES ####
        torsion_label = Label(self.torsion_frame, text=u"#      K          min.    phase     paths  ",
                            bg=self.main_color)
        torsion_label.grid(row=0, column=0, columnspan=6, sticky='w')

        torsiontypes_yscroll = Scrollbar(self.torsion_frame)
        torsiontypes_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.torsiontypes_listbox = Listbox(self.torsion_frame, yscrollcommand = torsiontypes_yscroll.set,
                                      width=45, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        torsiontypes_yscroll.config(command=self.torsiontypes_listbox.yview)
        self.torsiontypes_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.torsiontypes_listbox.config(font=tkFont.Font(family="Courier", size=12))

        import_torsion_prm = Button(self.torsion_frame, text='Import', highlightbackground=self.main_color,
                            command=self.import_torsion_prm)
        import_torsion_prm.grid(row=11, column=1)

        edit_torsion_prm = Button(self.torsion_frame, text='Edit', highlightbackground=self.main_color,
                               command=self.edit_torsion_prm)
        edit_torsion_prm.grid(row=11, column=2)

        add_torsion_prm = Button(self.torsion_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_torsion_prm)
        add_torsion_prm.grid(row=11, column=3)

        changetypes_label = Label(self.torsion_frame, text=u"    i           j           k           l       "
                                                         u"\N{GREEK SMALL LETTER PHI}1"
                                                            u"    \N{GREEK SMALL LETTER PHI}2"
                                                            u"    \N{GREEK SMALL LETTER PHI}3"
                                                            u"    \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        changetypes_label.grid(row=0, column=7, sticky='w')

        changetorsion_yscroll = Scrollbar(self.torsion_frame)
        changetorsion_yscroll.grid(row = 1, rowspan=10, column = 8, sticky = 'nsw', padx=(0,10))
        self.changetorsion_listbox = Listbox(self.torsion_frame, yscrollcommand = changetorsion_yscroll.set,
                                      width=48, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        changetorsion_yscroll.config(command=self.changetorsion_listbox.yview)
        self.changetorsion_listbox.grid(row=1, rowspan=10, column = 7, sticky = 'e')
        self.changetorsion_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.changetorsion_listbox.bind('<<ListboxSelect>>', self.list_changetorsion_event)

        include_ang = Label(self.torsion_frame, text='Include changing torsions:', bg=self.main_color)
        include_ang.grid(row=11, column=7)

        all_torsion_check = Checkbutton(self.torsion_frame, bg=self.main_color, variable=self.all_torsions,
                                        command=self.toggle_show_tor)
        all_torsion_check.grid(row=11, column=7, sticky='e')

        #### IMPROPER TYPES ####
        improper_label = Label(self.improper_frame, text=u"#      K          phase  ",
                            bg=self.main_color)
        improper_label.grid(row=0, column=0, columnspan=6, sticky='w')

        impropertypes_yscroll = Scrollbar(self.improper_frame)
        impropertypes_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.impropertypes_listbox = Listbox(self.improper_frame, yscrollcommand = impropertypes_yscroll.set,
                                      width=30, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        impropertypes_yscroll.config(command=self.impropertypes_listbox.yview)
        self.impropertypes_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.impropertypes_listbox.config(font=tkFont.Font(family="Courier", size=12))

        import_improper_prm = Button(self.improper_frame, text='Import', highlightbackground=self.main_color,
                            command=self.import_improper_prm)
        import_improper_prm.grid(row=11, column=1)

        edit_improper_prm = Button(self.improper_frame, text='Edit', highlightbackground=self.main_color,
                               command=self.edit_improper_prm)
        edit_improper_prm.grid(row=11, column=2)

        add_improper_prm = Button(self.improper_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_improper_prm)
        add_improper_prm.grid(row=11, column=3)

        changetypes_label = Label(self.improper_frame, text=u"    i           j           k           l       "
                                                         u"\N{GREEK SMALL LETTER PHI}1"
                                                            u"    \N{GREEK SMALL LETTER PHI}2"
                                                            u"    \N{GREEK SMALL LETTER PHI}3"
                                                            u"    \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        changetypes_label.grid(row=0, column=7, sticky='w')

        changeimproper_yscroll = Scrollbar(self.improper_frame)
        changeimproper_yscroll.grid(row = 1, rowspan=10, column = 8, sticky = 'nsw', padx=(0,10))
        self.changeimproper_listbox = Listbox(self.improper_frame, yscrollcommand = changeimproper_yscroll.set,
                                      width=48, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        changeimproper_yscroll.config(command=self.changeimproper_listbox.yview)
        self.changeimproper_listbox.grid(row=1, rowspan=10, column = 7, sticky = 'e')
        self.changeimproper_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.changeimproper_listbox.bind('<<ListboxSelect>>', self.list_changeimproper_event)

        add_imp = Button(self.improper_frame, text = 'Add', highlightbackground=self.main_color, command=self.add_improper)
        add_imp.grid(row=11, column=7, sticky='w', padx=(100,0))

        del_imp = Button(self.improper_frame, text = 'Delete', highlightbackground=self.main_color, command=self.del_improper)
        del_imp.grid(row=11, column=7, sticky='e', padx=(0,100))

        include_ang = Label(self.improper_frame, text='Include changing impropers:', bg=self.main_color)
        include_ang.grid(row=12, column=7)

        all_improper_check = Checkbutton(self.improper_frame, bg=self.main_color, variable=self.all_impropers,
                                         command=self.toggle_show_imp)
        all_improper_check.grid(row=12, column=7, sticky='e')


        #### SOFT PAIRS ####
        soft_label = Label(self.softpair_frame, text=u"Qi      Qj ",
                            bg=self.main_color)
        soft_label.grid(row=0, column=0, columnspan=6, sticky='w')

        softpairs_yscroll = Scrollbar(self.softpair_frame)
        softpairs_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.softpairs_listbox = Listbox(self.softpair_frame, yscrollcommand = softpairs_yscroll.set,
                                      width=45, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        softpairs_yscroll.config(command=self.softpairs_listbox.yview)
        self.softpairs_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.softpairs_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.softpairs_listbox.bind('<<ListboxSelect>>', self.list_softpairs_event)

        add_softpair = Button(self.softpair_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_softpair)
        add_softpair.grid(row=4, column=7)

        del_softpair = Button(self.softpair_frame, text='-', highlightbackground=self.main_color,
                              command=self.del_softpair)
        del_softpair.grid(row=4, column=8)

        #### COUPLINGS ####
        angle_coupling = Label(self.coupling_frame, text='Angles',
                            bg=self.main_color)
        angle_coupling.grid(row=0, column=0, columnspan=1, sticky='w')

        angcop_yscroll = Scrollbar(self.coupling_frame)
        angcop_yscroll.grid(row = 1, rowspan=10, column = 1, sticky = 'nsw', padx=(0,10))
        self.angle_couplings_listbox = Listbox(self.coupling_frame, yscrollcommand = angcop_yscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        angcop_yscroll.config(command=self.angle_couplings_listbox.yview)
        self.angle_couplings_listbox.grid(row=1, rowspan=10, column = 0, columnspan=1, sticky = 'e')
        self.angle_couplings_listbox.config(font=tkFont.Font(family="Courier", size=12))

        torsion_coupling = Label(self.coupling_frame, text='Torsions',
                            bg=self.main_color)
        torsion_coupling.grid(row=0, column=2, columnspan=1, sticky='w')

        torcop_yscroll = Scrollbar(self.coupling_frame)
        torcop_yscroll.grid(row = 1, rowspan=10, column = 3, sticky = 'nsw', padx=(0,10))
        self.torsion_couplings_listbox = Listbox(self.coupling_frame, yscrollcommand = torcop_yscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        torcop_yscroll.config(command=self.torsion_couplings_listbox.yview)
        self.torsion_couplings_listbox.grid(row=1, rowspan=10, column = 2, columnspan=1, sticky = 'e')
        self.torsion_couplings_listbox.config(font=tkFont.Font(family="Courier", size=12))

        improper_coupling = Label(self.coupling_frame, text='Impropers',
                            bg=self.main_color)
        improper_coupling.grid(row=0, column=4, columnspan=1, sticky='w')

        impcop_yscroll = Scrollbar(self.coupling_frame)
        impcop_yscroll.grid(row = 1, rowspan=10, column = 5, sticky = 'nsw', padx=(0,10))
        self.improper_couplings_listbox = Listbox(self.coupling_frame, yscrollcommand = impcop_yscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        impcop_yscroll.config(command=self.improper_couplings_listbox.yview)
        self.improper_couplings_listbox.grid(row=1, rowspan=10, column = 4, columnspan=1, sticky = 'e')
        self.improper_couplings_listbox.config(font=tkFont.Font(family="Courier", size=12))


        #### EXCLUDED PAIRS ####
        excluded_label = Label(self.excludedpairs_frame, text=u"     i          j       \N{GREEK SMALL LETTER PHI}1"
                                                            u"   \N{GREEK SMALL LETTER PHI}2"
                                                            u"   \N{GREEK SMALL LETTER PHI}3"
                                                            u"   \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        excluded_label.grid(row=0, column=0, columnspan=6, sticky='w')

        excludedpairs_yscroll = Scrollbar(self.excludedpairs_frame)
        excludedpairs_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.excludedpairs_listbox = Listbox(self.excludedpairs_frame, yscrollcommand = excludedpairs_yscroll.set,
                                      width=35, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        excludedpairs_yscroll.config(command=self.excludedpairs_listbox.yview)
        self.excludedpairs_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.excludedpairs_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.excludedpairs_listbox.bind('<<ListboxSelect>>', self.list_excludedpairs_event)

        add_excludedpair = Button(self.excludedpairs_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_excludedpair)
        add_excludedpair.grid(row=4, column=7)

        #add_non_q_pair = Button(self.excludedpairs_frame, text='Non Q-atoms', highlightbackground=self.main_color,
        #                        command=self.add_nonq_excludedpair)
        #add_non_q_pair.grid(row=5, column=7, columnspan=2)

        del_excludedpair = Button(self.excludedpairs_frame, text='-', highlightbackground=self.main_color,
                              command=self.del_excludedpair)
        del_excludedpair.grid(row=4, column=8)

        state1 = Button(self.excludedpairs_frame, text=u"\N{GREEK SMALL LETTER PHI}1", highlightbackground=self.main_color,
                        command=lambda: self.excluded_state(0))
        state1.grid(row=11, column=2)

        state2 = Button(self.excludedpairs_frame, text=u"\N{GREEK SMALL LETTER PHI}2", highlightbackground=self.main_color,
                        command=lambda: self.excluded_state(1))
        state2.grid(row=11, column=3)

        self.state3 = Button(self.excludedpairs_frame, text=u"\N{GREEK SMALL LETTER PHI}3",
                             highlightbackground=self.main_color, command=lambda: self.excluded_state(2))
        self.state3.grid(row=11, column=4)
        self.state3.config(state=DISABLED)
        self.state3_enabled.append(self.state3)

        self.state4 = Button(self.excludedpairs_frame, text=u"\N{GREEK SMALL LETTER PHI}4",
                             highlightbackground=self.main_color, command=lambda: self.excluded_state(3))
        self.state4.grid(row=11, column=5)
        self.state4.config(state=DISABLED)
        self.state4_enabled.append(self.state4)

        #### EL SCALE ####
        elscale_label = Label(self.elscale_frame, text=u"Qi      Qj  Scaling",
                            bg=self.main_color)
        elscale_label.grid(row=0, column=0, columnspan=6, sticky='w')

        elscale_yscroll = Scrollbar(self.softpair_frame)
        elscale_yscroll.grid(row = 1, rowspan=10, column = 6, sticky = 'nsw', padx=(0,10))
        self.elscale_listbox = Listbox(self.elscale_frame, yscrollcommand = softpairs_yscroll.set,
                                      width=20, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        elscale_yscroll.config(command=self.elscale_listbox.yview)
        self.elscale_listbox.grid(row=1, rowspan=10, column = 0, columnspan=6, sticky = 'e')
        self.elscale_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.elscale_listbox.bind('<<ListboxSelect>>', self.list_elscale_event)

        add_elscale = Button(self.elscale_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_elscale)
        add_elscale.grid(row=4, column=7)

        del_elscale = Button(self.elscale_frame, text='-', highlightbackground=self.main_color,
                              command=self.del_elscale)
        del_elscale.grid(row=4, column=8)

        self.elscale = Spinbox(self.elscale_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.1)
        self.elscale.grid(row=5, column=7, columnspan=2)

        set_elscale = Button(self.elscale_frame, text='Set scale', highlightbackground=self.main_color,
                             command=self.set_elscale)
        set_elscale.grid(row=6, column=7, columnspan=2)


        #### LAMBDA STEPS ####
        lock_label = Label(self.lambdasteps_frame, text='Lock', bg=self.main_color)
        lock_label.grid(row=0, column=0, columnspan=2, sticky='e')

        start_label = Label(self.lambdasteps_frame, text='Start', bg=self.main_color)
        start_label.grid(row=0, column=2)

        end_label = Label(self.lambdasteps_frame, text='END', bg=self.main_color)
        end_label.grid(row=0, column=3)

        l1_label = Label(self.lambdasteps_frame, text=u'\N{GREEK SMALL LETTER LAMDA}1', bg=self.main_color)
        l1_label.grid(row=1, column=0)

        l2_label = Label(self.lambdasteps_frame, text=u'\N{GREEK SMALL LETTER LAMDA}2', bg=self.main_color)
        l2_label.grid(row=2, column=0)

        l3_label = Label(self.lambdasteps_frame, text=u'\N{GREEK SMALL LETTER LAMDA}3', bg=self.main_color)
        l3_label.grid(row=3, column=0)

        l4_label = Label(self.lambdasteps_frame, text=u'\N{GREEK SMALL LETTER LAMDA}4', bg=self.main_color)
        l4_label.grid(row=4, column=0)

        sum_label = Label(self.lambdasteps_frame, text='SUM', bg=self.main_color)
        sum_label.grid(row=5, column=0, columnspan=2)

        lock1 = Checkbutton(self.lambdasteps_frame, bg=self.main_color, variable=self.lock_lambda1)
        lock1.grid(row=1, column=1)

        lock2 = Checkbutton(self.lambdasteps_frame, bg=self.main_color, variable=self.lock_lambda2)
        lock2.grid(row=2, column=1)

        self.lock3 = Checkbutton(self.lambdasteps_frame, bg=self.main_color, variable=self.lock_lambda3)
        self.lock3.grid(row=3, column=1)
        self.state3_enabled.append(self.lock3)
        self.lock3.config(state=DISABLED)

        self.lock4 = Checkbutton(self.lambdasteps_frame, bg=self.main_color, variable=self.lock_lambda4)
        self.lock4.grid(row=4, column=1)
        self.state4_enabled.append(self.lock4)
        self.lock4.config(state=DISABLED)

        self.start1 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l1_start)
        self.start1.grid(row= 1, column=2)
        self.start1.delete(0, END)
        self.start1.insert(0, '0.00')

        self.start2 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l2_start)
        self.start2.grid(row= 2, column=2)
        self.start2.delete(0, END)
        self.start2.insert(0, '1.00')

        self.start3 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l3_start)
        self.start3.grid(row=3, column=2)
        self.state3_enabled.append(self.start3)
        self.start3.config(state=DISABLED)

        self.start4 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l4_start)
        self.start4.grid(row=4, column=2)
        self.state4_enabled.append(self.start4)
        self.start4.config(state=DISABLED)

        self.sum_start = Entry(self.lambdasteps_frame, width=6, highlightthickness=0, relief=GROOVE)
        self.sum_start.grid(row=5, column=2)
        self.sum_start.insert(0, '1.000')

        self.end1 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l1_end)
        self.end1.grid(row= 1, column=3)
        self.end1.delete(0, END)
        self.end1.insert(0, '0.00')

        self.end2 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l2_end)
        self.end2.grid(row= 2, column=3)
        self.end2.delete(0, END)
        self.end2.insert(0, '1.00')

        self.end3 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l3_end)
        self.end3.grid(row=3, column=3)
        self.state3_enabled.append(self.end3)
        self.end3.config(state=DISABLED)

        self.end4 = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.00, increment=0.01, textvariable=self.l4_end)
        self.end4.grid(row=4, column=3)
        self.state4_enabled.append(self.end4)
        self.end4.config(state=DISABLED)

        self.sum_end = Entry(self.lambdasteps_frame, width=6, highlightthickness=0, relief=GROOVE)
        self.sum_end.grid(row=5, column=3)
        self.sum_end.insert(0, '1.000')



        self.lambda_sampling = StringVar()
        self.lambda_sampling.set('Linear')
        lambda_sampling = OptionMenu(self.lambdasteps_frame, self.lambda_sampling, 'Linear', 'Sigmoidal')
        lambda_sampling.configure(bg=self.main_color, width=10)
        lambda_sampling.grid(row=1, column=4, columnspan=2)

        step_label = Label(self.lambdasteps_frame, text=u'\N{GREEK SMALL LETTER LAMDA}-Step', bg=self.main_color)
        step_label.grid(row=2, column=4, sticky='s', padx=10)

        lambdastep = Spinbox(self.lambdasteps_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=0.01, to=0.40, increment=0.01, textvariable=self.lambdastep)
        lambdastep.grid(row=3,rowspan=3, column=4, sticky='n', padx=10)

        add_button = Button(self.lambdasteps_frame, text=u"\u21D2",
                            highlightbackground=self.main_color, command=self.add_lambda_steps)
        add_button.grid(row=2, rowspan=2, column=5, padx=(0,10))
        add_button.config(font=tkFont.Font(family="Courier", size=28))



        lambda_label = Label(self.lambdasteps_frame, text=u"  #      \N{GREEK SMALL LETTER LAMDA}1    "
                                                          u"  \N{GREEK SMALL LETTER LAMDA}2    "
                                                          u"   \N{GREEK SMALL LETTER LAMDA}3    "
                                                          u"   \N{GREEK SMALL LETTER LAMDA}4",
                                                          bg=self.main_color)
        lambda_label.grid(row=0, column=6, columnspan=6, sticky='w')

        lambdasteps_yscroll = Scrollbar(self.lambdasteps_frame)
        lambdasteps_yscroll.grid(row = 1, rowspan=10, column = 12, sticky = 'nsw', padx=(0,10))
        self.lambdasteps_listbox = Listbox(self.lambdasteps_frame, yscrollcommand = lambdasteps_yscroll.set,
                                      width=30, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        lambdasteps_yscroll.config(command=self.lambdasteps_listbox.yview)
        self.lambdasteps_listbox.grid(row=1, rowspan=10, column = 6, columnspan=6, sticky = 'e')
        self.lambdasteps_listbox.config(font=tkFont.Font(family="Courier", size=12))

        del_lamda_steps = Button(self.lambdasteps_frame, text='Delete', highlightbackground=self.main_color,
                                 command=self.del_lambda_steps)
        del_lamda_steps.grid(row=11, column=7)

        reverse_lamda = Button(self.lambdasteps_frame, text='Reverse', highlightbackground=self.main_color,
                               command=self.reverse_lambda_list)
        reverse_lamda.grid(row=11, column=8)

        #### TEMPERATURE(S) ####

        temperature_label = Label(self.temperature_frame, text='Temperature', bg=self.main_color)
        temperature_label.grid(row=4, column=2)

        runs_label = Label(self.temperature_frame, text='Runs', bg=self.main_color)
        runs_label.grid(row=4, column=3)

        self.temperature = Spinbox(self.temperature_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=273, to=373, increment=1)
        self.temperature.grid(row= 5, column=2)
        self.temperature.delete(0, END)
        self.temperature.insert(0, 300)

        self.runs = Spinbox(self.temperature_frame, width=5, highlightthickness=0, relief=GROOVE,
                                   from_=1, to=500, increment=1)
        self.runs.grid(row= 5, column=3)
        self.runs.delete(0, END)
        self.runs.insert(0, 1)

        add_button = Button(self.temperature_frame, text=u"\u21D2",
                            highlightbackground=self.main_color, command=self.add_temperature)
        add_button.grid(row=0, rowspan=10, column=5, padx=(10,10))
        add_button.config(font=tkFont.Font(family="Courier", size=28))



        temp_label = Label(self.temperature_frame, text='  Temp.   Runs',
                                                          bg=self.main_color)
        temp_label.grid(row=0, column=6, columnspan=6, sticky='w')

        temp_yscroll = Scrollbar(self.temperature_frame)
        temp_yscroll.grid(row = 1, rowspan=10, column = 12, sticky = 'nsw', padx=(0,10))
        self.temp_listbox = Listbox(self.temperature_frame, yscrollcommand = temp_yscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        temp_yscroll.config(command=self.temp_listbox.yview)
        self.temp_listbox.grid(row=1, rowspan=10, column = 6, columnspan=6, sticky = 'e')
        self.temp_listbox.config(font=tkFont.Font(family="Courier", size=12))

        del_temp = Button(self.temperature_frame, text='Delete', highlightbackground=self.main_color,
                                 command=self.del_temperature)
        del_temp.grid(row=11, column=6, columnspan=6)

        #SOFTCORE
        ######### Charges ##########
        softcore_label = Label(self.softcore_frame, text=u" Qi       \N{GREEK SMALL LETTER PHI}1    "
                                                     u"    \N{GREEK SMALL LETTER PHI}2"
                                                     u"       \N{GREEK SMALL LETTER PHI}3"
                                                     u"       \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        softcore_label.grid(row=0, column=0, sticky='w')

        softcore_yscroll = Scrollbar(self.softcore_frame)
        softcore_yscroll.grid(row = 1, rowspan=10, column = 1, sticky = 'nsw', padx=(0,10))
        self.softcore_listbox = Listbox(self.softcore_frame, yscrollcommand = softcore_yscroll.set,
                                      width=35, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        softcore_yscroll.config(command=self.softcore_listbox.yview)
        self.softcore_listbox.grid(row=1, rowspan=10, column = 0, sticky = 'e')
        self.softcore_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.softcore_listbox.bind('<<ListboxSelect>>', self.list_softcore_event)

        add_softcore = Button(self.softcore_frame, text = '+', highlightbackground=self.main_color,
                              command=self.add_softcore)
        add_softcore.grid(row=0, column=3, sticky='e')

        del_softcore = Button(self.softcore_frame, text= '-', highlightbackground=self.main_color,
                              command=self.del_softcore)
        del_softcore.grid(row=0, column=4, sticky='w')

        s1 = Label(self.softcore_frame, text=u"\N{GREEK SMALL LETTER PHI}1", bg=self.main_color)
        s1.grid(row=3, column=3)

        s2 = Label(self.softcore_frame, text=u"\N{GREEK SMALL LETTER PHI}2", bg=self.main_color)
        s2.grid(row=3, column=4)

        s3 = Label(self.softcore_frame, text=u"\N{GREEK SMALL LETTER PHI}3", bg=self.main_color)
        s3.grid(row=5, column=3)

        s4 = Label(self.softcore_frame, text=u"\N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        s4.grid(row=5, column=4)

        self.softcore1 = Spinbox(self.softcore_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=0, to=100, increment=1)
        self.softcore1.grid(row=4, column=3)
        self.softcore1.delete(0, END)
        self.softcore1.insert(0, '')

        self.softcore2 = Spinbox(self.softcore_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-3.00, to=3.00, increment=0.01)
        self.softcore2.grid(row=4, column=4)
        self.softcore2.delete(0, END)
        self.softcore2.insert(0, '')

        self.softcore3 = Spinbox(self.softcore_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-3.00, to=3.00, increment=0.01)
        self.softcore3.grid(row=6, column=3)
        self.softcore3.delete(0, END)
        self.softcore3.insert(0, '')
        self.softcore3.config(state=DISABLED)
        self.state3_enabled.append(self.softcore3)

        self.softcore4 = Spinbox(self.softcore_frame, width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-3.00, to=3.00, increment=0.01)
        self.softcore4.grid(row=6, column=4)
        self.softcore4.delete(0, END)
        self.softcore4.insert(0, '')
        self.softcore4.config(state=DISABLED)
        self.state4_enabled.append(self.softcore4)

        edit_softcore = Button(self.softcore_frame, text='Set softcore', highlightbackground=self.main_color,
                              command=self.edit_softcore)
        edit_softcore.grid(row=7, column=3, columnspan=2)

        softcore = Label(self.softcore_frame, text='Use softcore max potential: ', bg=self.main_color)
        softcore.grid(row=11, column=0, columnspan=2)

        use_soft_max_pot = Checkbutton(self.softcore_frame, bg=self.main_color, variable= self.use_soft_max_pt)
        use_soft_max_pot.grid(row=11,column=1, columnspan=2, sticky='w')

        #ATOMS frame
        atoms_label = Label(self.atoms_frame, text=u"   Qi     \N{GREEK SMALL LETTER PHI}1"
                                                            u"    \N{GREEK SMALL LETTER PHI}2"
                                                            u"    \N{GREEK SMALL LETTER PHI}3"
                                                            u"    \N{GREEK SMALL LETTER PHI}4", bg=self.main_color)
        atoms_label.grid(row=0, column=2, columnspan=6, sticky='w')

        atoms_yscroll = Scrollbar(self.atoms_frame)
        atoms_yscroll.grid(row = 1, rowspan=10, column = 8, sticky = 'nsw', padx=(0,10))
        self.atoms_listbox = Listbox(self.atoms_frame, yscrollcommand = atoms_yscroll.set,
                                      width=30, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        atoms_yscroll.config(command=self.atoms_listbox.yview)
        self.atoms_listbox.grid(row=1, rowspan=10, column = 2, columnspan=6, sticky = 'e')
        self.atoms_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.atoms_listbox.bind('<<ListboxSelect>>', self.list_atomchange_event)


        element_label = Label(self.atoms_frame, text='Atom', bg=self.main_color)
        element_label.grid(row=0, column=0, columnspan=2)

        elem_yscroll = Scrollbar(self.atoms_frame)
        elem_yscroll.grid(row = 1, rowspan=10, column = 1, sticky = 'nsw', padx=(0,10))
        self.elem_listbox = Listbox(self.atoms_frame, yscrollcommand = elem_yscroll.set,
                                      width=5, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        elem_yscroll.config(command=self.elem_listbox.yview)
        self.elem_listbox.grid(row=1, rowspan=10, column = 0, sticky = 'e')
        self.elem_listbox.config(font=tkFont.Font(family="Courier", size=12))

        #Insert default elements here
        elements = ['H', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'Na', 'Mg', 'K', 'Ca', 'Fe', 'Zn']
        for i in elements:
            self.elem_listbox.insert(END, '%3s' % i)


        add_atomchange = Button(self.atoms_frame, text='+', highlightbackground=self.main_color,
                              command=self.add_atom_change)
        add_atomchange.grid(row=4, column=9)

        del_atomchange = Button(self.atoms_frame, text='-', highlightbackground=self.main_color,
                              command=self.del_atom_change)
        del_atomchange.grid(row=4, column=10)

        atom_state1 = Button(self.atoms_frame, text=u"\N{GREEK SMALL LETTER PHI}1", highlightbackground=self.main_color,
                        command=lambda: self.change_atom_state(0))
        atom_state1.grid(row=11, column=4)

        atom_state2 = Button(self.atoms_frame, text=u"\N{GREEK SMALL LETTER PHI}2", highlightbackground=self.main_color,
                        command=lambda: self.change_atom_state(1))
        atom_state2.grid(row=11, column=5)

        self.atom_state3 = Button(self.atoms_frame, text=u"\N{GREEK SMALL LETTER PHI}3",
                             highlightbackground=self.main_color, command=lambda: self.change_atom_state(2))
        self.atom_state3.grid(row=11, column=6)
        self.atom_state3.config(state=DISABLED)
        self.state3_enabled.append(self.atom_state3)

        self.atom_state4 = Button(self.atoms_frame, text=u"\N{GREEK SMALL LETTER PHI}4",
                             highlightbackground=self.main_color, command=lambda: self.change_atom_state(3))
        self.atom_state4.grid(row=11, column=7)
        self.atom_state4.config(state=DISABLED)
        self.state4_enabled.append(self.atom_state4)

        #MONITOR GROUPS FRAME

        groups_label = Label(self.monitorgroups_frame, text='Groups',
                            bg=self.main_color)
        groups_label.grid(row=0, column=0, columnspan=1, sticky='w')

        groups_yscroll = Scrollbar(self.monitorgroups_frame)
        groups_yscroll.grid(row = 1, rowspan=10, column = 1, sticky = 'nsw', padx=(0,10))
        self.groups_listbox = Listbox(self.monitorgroups_frame, yscrollcommand = groups_yscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        groups_yscroll.config(command=self.groups_listbox.yview)
        self.groups_listbox.grid(row=1, rowspan=10, column = 0, columnspan=1, sticky = 'e')
        self.groups_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.groups_listbox.bind('<<ListboxSelect>>', self.list_monitorgroups_event)


        group_atoms = Label(self.monitorgroups_frame, text='Atoms',
                            bg=self.main_color)
        group_atoms.grid(row=0, column=2, columnspan=1, sticky='w')

        add_group = Button(self.monitorgroups_frame, text='+', highlightbackground=self.main_color,
                           command=self.add_monitor_group)
        add_group.grid(row=11, column=0, sticky='w')

        del_group = Button(self.monitorgroups_frame, text='-', highlightbackground=self.main_color,
                           command=self.del_monitor_group)
        del_group.grid(row=11, column=0, sticky='e')


        group_atoms_yscroll = Scrollbar(self.monitorgroups_frame)
        group_atoms_yscroll.grid(row = 1, rowspan=10, column = 3, sticky = 'nsw', padx=(0,10))
        self.group_atoms_listbox = Listbox(self.monitorgroups_frame, yscrollcommand = group_atoms_yscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        group_atoms_yscroll.config(command=self.group_atoms_listbox.yview)
        self.group_atoms_listbox.grid(row=1, rowspan=10, column = 2, columnspan=1, sticky = 'e')
        self.group_atoms_listbox.config(font=tkFont.Font(family="Courier", size=12))

        add_monitor_atoms = Button(self.monitorgroups_frame, text='+', highlightbackground=self.main_color,
                           command=self.add_monitor_atoms)
        add_monitor_atoms.grid(row=11, column=2, sticky='w')

        del_monitor_atoms = Button(self.monitorgroups_frame, text='-', highlightbackground=self.main_color,
                           command=self.del_monitor_atoms)
        del_monitor_atoms.grid(row=11, column=2, sticky='e')


        monitor_pairs = Label(self.monitorgroups_frame, text='Monitor pairs',
                            bg=self.main_color)
        monitor_pairs.grid(row=0, column=4, columnspan=1, sticky='w')

        monitor_pairs_yscroll = Scrollbar(self.monitorgroups_frame)
        monitor_pairs_yscroll.grid(row = 1, rowspan=10, column = 5, sticky = 'nsw', padx=(0,10))
        self.monitor_pairs_listbox = Listbox(self.monitorgroups_frame, yscrollcommand = monitor_pairs_yscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        monitor_pairs_yscroll.config(command=self.monitor_pairs_listbox.yview)
        self.monitor_pairs_listbox.grid(row=1, rowspan=10, column = 4, columnspan=1, sticky = 'e')
        self.monitor_pairs_listbox.config(font=tkFont.Font(family="Courier", size=12))

        add_monitor_pair = Button(self.monitorgroups_frame, text='+', highlightbackground=self.main_color,
                           command=self.add_monitor_pair)
        add_monitor_pair.grid(row=11, column=4, sticky='w')

        del_monitor_pair = Button(self.monitorgroups_frame, text='-', highlightbackground=self.main_color,
                           command=self.del_monitor_pair)
        del_monitor_pair.grid(row=11, column=4, sticky='e')


        #Save/close frame3
        md = Button(frame4, text='Configure MD', highlightbackground=self.main_color, command=self.config_md)
        md.grid(row=0, column=0, columnspan=4)


        overwrite_label = Label(frame4, text='Overwrite existing files', bg=self.main_color)
        overwrite_label.grid(row=2, column=0, columnspan=3)

        overwrite = Checkbutton(frame4, bg=self.main_color, variable=self.check_overwrite)
        overwrite.grid(row=2,column=3, sticky='w')

        submit_button = Button(frame4, text='Run', highlightbackground=self.main_color, command=None)
        submit_button.grid(row=1, column=0)

        save_button = Button(frame4, text='Write', highlightbackground=self.main_color,
                             command=lambda: self.write_evb(True))
        save_button.grid(row=1, column=1)

        viewfep = Button(frame4, text = 'Edit FEP', highlightbackground=self.main_color, command=self.open_file)
        viewfep.grid(row=1, column=2)

        close_button= Button(frame4, text='Close', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=1, column=3)


