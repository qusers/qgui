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

from Tkinter import Label, Button, Frame, Toplevel, DISABLED, NORMAL, Scrollbar, GROOVE, Listbox, EXTENDED, END, \
    OptionMenu, StringVar, Spinbox, SINGLE, LabelFrame, MULTIPLE
from subprocess import Popen, PIPE
import tkFont
import time
import os
import signal
import sys
import numpy as np
import shutil
from copy import deepcopy
import build as bld
from tkFileDialog import askopenfilenames


class ViewPyMol(Toplevel):
    def __init__(self, app, root, start_cmd=None, edit_mode=False):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root

        self.pymol_mode = StringVar()

        #Variable to trace
        self.new_bond = StringVar()
        self.new_angle = StringVar()
        self.new_torsion = StringVar()
        self.existing_bond = StringVar()
        self.existing_angle = StringVar()
        self.existing_torsion = StringVar()

        if self.app.pdb_id:
            self.pdbfile = self.app.pdb_id
        else:
            self.pdbfile = '%s/.new.pdb' % self.app.workdir
            #Need to add empty line to avoid pymol complaining about it
            f = open(self.pdbfile, 'w')
            f.write('\n')
            f.close()
            edit_mode = True

        self.app.log('info', 'Launching pyMol with %s' % self.pdbfile.split('/')[-1])

        self.start_cmd = start_cmd
        self.session = None

        #Keep track if pseudo atom is created or not
        self.pseudo_atom = False

        self.pymol_color = StringVar()
        self.pymol_atoms_color = StringVar()
        self.select_mode = StringVar()
        self.pymol_presets = StringVar()
        self.create_bond = StringVar()
        self.bond_atom_name = StringVar()
        self.atoms_fragments = StringVar()
        self.p_atoms_selection = StringVar()


        self.pymol_color.set('Gray')
        self.pymol_atoms_color.set('C')
        self.select_mode.set('Residues')
        self.pymol_presets.set('Presets')
        self.create_bond.set('H')
        self.atoms_fragments.set('ATOMS')
        self.p_atoms_selection.set('P-atoms')

        #Keep track of click sequence in edit mode. Max 4 atoms {number: atomnumber}
        self.selected_atoms = list()

        #Keep track of temporary pdb files for undo function
        self.tmp_pdb = list()

        #Max number of undo (this in principle means how many temporary pdb files to store):
        self.max_undo = 20

        #Store tmp pdb files:
        self.tmp_pdb_dir = '%s/.tmp_pdb_history' % self.app.workdir

        #Get fragments/residues fragments.dat
        self.fragments_dict = self.get_fragments()

        #Get atom info from file:
        atom_file = self.app.qgui_path + '/Qmods/atoms.dat'
        self.atom_dict = self.get_atoms(atom_file)

        self.dialog_window()
        #Trace changes made to parameters in New atoms section
        self.new_bond.trace('w', self.update_pseudo_atom)
        self.new_angle.trace('w', self.update_pseudo_atom)
        self.new_torsion.trace('w', self.update_pseudo_atom)

        #Trace changes made to existing atoms
        self.existing_bond.trace('w', self.existing_bond_changed)
        self.existing_angle.trace('w', self.existing_angle_changed)
        self.existing_torsion.trace('w', self.existing_torsion_changed)

        #Variable determining if existing atoms should be modified or not
        self.modify_atoms = False

        if edit_mode:
            self.pymol_mode.set('Edit')
        else:
            self.pymol_mode.set('Display')

        self.show_var_frame()

        #Edit mode or display mode:
        self.pymol_mode.trace('w', self.show_var_frame)
        #Select atoms, residues...
        self.select_mode.trace('w', self.select_mode_changed)
        self.p_atoms_selection.trace('w', self.p_group_selection_changed)
        #Build by atoms or pre-defined fragments
        self.atoms_fragments.trace('w', self.build_mode_changed)

        #QPyMol default settings:
        self.pymol_settings = ['space cmyk', 'set valence, 0.1', 'set sphere_scale, 0.8',
                               'set sphere_quality, 2', 'set stick_radius, 0.17', 'set stick_quality, 10',
                               'set defer_builds_mode, 3', 'set surface_quality, 1',
                               'set spec_power, 250', 'set spec_reflect, 2',
                               'set cartoon_fancy_helices, 1']

        #QpyMol can be launched with additional start commands:
        #self.start_cmd = ['show cartoon', 'hide lines']

        self.add_atoms_to_list(self.pdbfile)

        #Fill atoms in listbox
        self.fill_build_list("ATOMS")

        if edit_mode:
            self.buildlist.selection_set(0)

        self.start_pymol()

        #Kill pymol if window is closed... though this does not quite work does it? TODO
        root.protocol('WM_DELETE_WINDOW', self.close_pymol)

    def get_p1_group(self):
        """
        Takes all atoms selected in listbox + p_atom 1 (self.selected_atoms)
        and returns a dictionary {atomnumber:{'xyz'[x,y,z]}}
        :return:
        """
        atomnumbers = list()

        #Original selected_atoms list is strings:

        selected_atoms = map(int, self.selected_atoms)

        atomnumbers.append(selected_atoms[0])
        selections = self.listbox.curselection()

        for seleced in selections:
            nr = int(self.listbox.get(seleced).split()[0])
            if nr not in atomnumbers and nr not in selected_atoms:
                atomnumbers.append(nr)
        p1_dict = dict()

        #Open current pdb file and get xyz for atomnumbers
        #Get original/previous pdb file:
        if len(self.tmp_pdb) < 1:
            old_pdb = self.pdbfile
        else:
            old_pdb = self.tmp_pdb[-1]

        with open(old_pdb, 'r') as pdb:
            for line in pdb:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atomnr = int(line.split()[1])
                    if atomnr in atomnumbers:
                        p1_dict[atomnr] = dict()
                        xyz = map(float, line[30:].split()[0:3])
                        p1_dict[atomnr]['xyz'] = xyz
                    if len(p1_dict.keys()) == len(atomnumbers):
                        break

        return p1_dict

    def existing_bond_changed(self, *args):
        """
        Modify existing bond
        :param args:
        :return:
        """
        if len(self.selected_atoms) < 2:
            return

        if not self.modify_atoms:
            return

        bond, angle, torsion = None, None, None

        try:
            bond = float(self.existing_bond.get())
        except:
            #Return silently, this typically happens while user is actually typing in a number
            return

        p1_group = self.get_p1_group()
        p_atoms = self.get_p_atoms()
        #Get xyz for pseudo atom dict from build.py
        new_group = bld.BuildByGroup('b', p_atoms[0], p_atoms[1], p_atoms[2], p_atoms[3], bond, angle, torsion,
                                     None, None, p1_group)

        self.write_tmp_pdb(group=new_group, insert_after_atomnr=None, resnr_add=0, replace=True)

    def existing_angle_changed(self, *args):
        """
        Modify existing angle
        :param args:
        :return:
        """
        if len(self.selected_atoms) < 3:
            return
        if not self.modify_atoms:
            return

        bond, angle, torsion = None, None, None

        try:
            angle = float(self.existing_angle.get())
        except:
            #Return silently, this typically happens while user is actually typing in a number
            return

        p1_group = self.get_p1_group()
        p_atoms = self.get_p_atoms()
        #Get xyz for pseudo atom dict from build.py
        new_group = bld.BuildByGroup('c', p_atoms[0], p_atoms[1], p_atoms[2], p_atoms[3], bond, angle, torsion,
                                     None, None, p1_group)

        self.write_tmp_pdb(group=new_group, insert_after_atomnr=None, resnr_add=0, replace=True)

    def existing_torsion_changed(self, *args):
        """
        Modify existing torsion
        :param args:
        :return:
        """
        if len(self.selected_atoms) < 4:
            return

        if not self.modify_atoms:
            return

        bond, angle, torsion = None, None, None

        try:
            torsion = float(self.existing_torsion.get())
        except:
            #Return silently, this typically happens while user is actually typing in a number
            return

        p1_group = self.get_p1_group()
        p_atoms = self.get_p_atoms()
        #Get xyz for pseudo atom dict from build.py
        new_group = bld.BuildByGroup('d', p_atoms[0], p_atoms[1], p_atoms[2], p_atoms[3], bond, angle, torsion,
                                     None, None, p1_group)

        self.write_tmp_pdb(group=new_group,insert_after_atomnr=None, resnr_add=0, replace=True)


    def get_p_atoms(self):
        """
        Get list of atom dictionaries for p1-p4 and return
        :return: p_atoms = [dict()]
        """

        p_atoms = [None, None, None, None]
        if len(self.selected_atoms) < 1:
            return p_atoms

        pdb_lines = self.get_pdb_atoms(self.selected_atoms)

        for i in range(len(self.selected_atoms)):
            atomnr = self.selected_atoms[i]
            xyz = map(float, pdb_lines[i][30:].split()[0:3])
            name = pdb_lines[i][13:17].strip()
            residue = pdb_lines[i][17:21].strip()
            residue_nr = pdb_lines[i][22:26].strip()
            p_atoms[i] = {'atomnr': int(atomnr),
                          'xyz': xyz,
                          'name': name,
                          'residue': residue,
                          'residuenr': int(residue_nr)}

        return p_atoms

    def add_build(self):
        """
        Triggered when Add button is pressed. Decide if to build by atoms or fragments and go to corresponding
        function.
        """
        bond, angle, torsion = None, None, None

        if self.pseudo_atom:
            self.session.stdin.write('delete new_atom\n')
            self.pseudo_atom = False


        if len(self.selected_atoms) > 0:
            #Get bond length
            try:
                bond = float(self.bond_new.get())
            except:
                self.app.errorBox('Error', 'Invalid bond length value')
                return
        if len(self.selected_atoms) > 1:
            #Get angle
            try:
                angle = float(self.angle_new.get())
            except:
                self.app.errorBox('Error', 'Invalid angle value')
                return
        if len(self.selected_atoms) > 2:
            #Get torsion
            try:
                torsion = float(self.torsion_new.get())
            except:
                self.app.errorBox('Error', 'Invalid dihedral value')
                return


        if len(self.buildlist.curselection()) > 0:
            #Get selected fragment/atom
            fragment = self.buildlist.get(self.buildlist.curselection()[0]).split()[-1]
            #Get list of p atoms (selected atoms ordered 1 - 4)
            p_atoms = self.get_p_atoms()
            if self.atoms_fragments.get() == 'ATOMS':
                self.build_atom(fragment, p_atoms, bond, angle, torsion)
            else:
                self.build_fragment(fragment, p_atoms, bond, angle, torsion)
        else:
            return

    def undo_build(self):
        """
        Unde build moves to the previous pdb structure
        :return:
        """
        if len(self.tmp_pdb) < 1:
            return

        #Delete last pdb
        os.remove(self.tmp_pdb[-1])
        del self.tmp_pdb[-1]

        if len(self.tmp_pdb) < 1:
            new = self.pdbfile
        else:
            new = self.tmp_pdb[-1]

        self.update_pymol_structure(new)
        self.add_atoms_to_list(new)
        self.regulate_edit_controls()

    def delete_atoms(self):
        """
        Deletes selected pdb atoms and create new tmp pdb file.
        """
        if len(self.selected_atoms) < 1:
            return

        old_pdb, new_pdb = self.get_old_new_pdb_name()

        #Add new pdb file to history
        self.add_pdb_history(new_pdb)

        new = open(new_pdb, 'w')
        i = 0
        with open(old_pdb, 'r') as old:
            for line in old:
                if 'ATOM' in line or 'HETATM' in line:
                    atomnr = line.split()[1]
                    if atomnr not in self.selected_atoms:
                        i += 1
                        new.write('%s%5d%s' % (line[0:6], i, line[11:]))

        #If nothing has been written to file, add empty line to avoid pymol complaining!
        if i < 1:
            new.write('\n')

        new.close()

        #make selected atoms zero:
        del self.selected_atoms[:]
        self.listbox.selection_clear(0, END)

        #update pymol with new tmp pdb file
        self.update_pymol_structure(new_pdb)
        self.add_atoms_to_list(new_pdb)

    def get_old_new_pdb_name(self):
        """
        Decide what the current PDB file is and name of new pdb file
        :return: old_pdb new_pdb
        """
        #Get original/previous pdb file:
        if len(self.tmp_pdb) < 1:
            old_pdb = self.pdbfile
        else:
            old_pdb = self.tmp_pdb[-1]

        new_pdb = '%s/tmp%03d.pdb' % (self.tmp_pdb_dir, len(self.tmp_pdb))

        return old_pdb, new_pdb

    def add_pdb_history(self, pdb):
        """
        Adds pdb file to history (for undo function) and controls that history does not exceed the max undo
        :param pdb:
        """
        #Add pdb to history
        self.tmp_pdb.append(pdb)

        #Check that PDB files in history is not exceeding the "max undo"
        if len(self.tmp_pdb) > self.max_undo:
            remove_from_history = self.tmp_pdb.pop(0)
            os.remove(remove_from_history)

    def write_tmp_pdb(self, group, insert_after_atomnr=None, resnr_add=0, replace=False):
        """
        writes a atom to pdb file, update list and pymol
        :param group: {atom_nr: {'atomnnr', 'xyz', 'name', 'residue', 'residuenr'}}
        :return:
        """
        #Get the current pdb file and name for new pdb file
        old_pdb, tmp_name = self.get_old_new_pdb_name()

        #Add new pdb file to history
        self.add_pdb_history(tmp_name)

        atomnr = 0
        inserted_fragment = False
        insert_fragment = False

        if not insert_after_atomnr and not replace:
            insert_fragment = True

        #Write the new pdb file:
        tmp_pdb = open(tmp_name, 'w')
        with open(old_pdb, 'r') as old:
            for line in old:
                if 'ATOM' in line or 'HETATM' in line:
                    atomnr = int(line.split()[1])
                    if inserted_fragment:
                        atomnr += 1
                        resnr = int(line[22:26]) + resnr_add
                        tmp_pdb.write('%s%5d%s%5d%s' % (line[0:6], atomnr, line[11:21], resnr, line[26:]))
                    elif replace:
                        if atomnr in group.keys():
                            tmp_pdb.write('%s' % line[0:30])
                            tmp_pdb.write('%8.3f%8.3f%8.3f\n' % (group[atomnr]['xyz'][0], group[atomnr]['xyz'][1],
                                                                 group[atomnr]['xyz'][2]))
                        else:
                            tmp_pdb.write(line)
                    else:
                        tmp_pdb.write(line)
                        if atomnr == insert_after_atomnr:
                            #insert the new fragment
                            insert_fragment = True

                if insert_fragment:
                    #insert the new fragment
                    for i in sorted(group.keys()):
                        atomnr += 1
                        resnr = group[i]['residuenr']
                        res = group[i]['residue']
                        x, y, z = group[i]['xyz'][0:]
                        name = group[i]['name']
                        tmp_pdb.write('ATOM  %5d  %4s%4s %4d    %8.3f%8.3f%8.3f\n' %
                                      (atomnr, name.ljust(4), res.ljust(4), resnr, x, y, z))
                    inserted_fragment = True
                    insert_fragment = False

        tmp_pdb.close()

        #make selected atoms zero:
        if not replace:
            del self.selected_atoms[:]
            self.listbox.selection_clear(0, END)

        #update pymol with new tmp pdb file
        self.update_pymol_structure(tmp_name)
        self.add_atoms_to_list(tmp_name)

        if replace:
            atoms_select = group.keys()
            #for i in self.selected_atoms:
            #    atoms_select.append(i)
            self.set_selection(atoms_select)
            self.pymol_highlight_edit(self.selected_atoms)
            self.pymol_command('select')

    def update_pymol_structure(self, pdb):
        """
        updates currens structure in pymol to the given pdb
        :param pdb:
        :return: Nada de nada
        """
        self.session.stdin.write('load %s, mol, state=1\n' % pdb)
        self.pymol_highlight_edit(self.selected_atoms)

    def get_new_atomname(self, name, p1_atom):
        """
        Decides proper atom name for new atom (f. ex name C would be C1)
        :param name: atom name of new atom
        :param p1_atom:
        :return: new_name
        """
        res_nr = p1_atom['residuenr']
        nr = 1

        atom_len = len(name)

        #Get correct pdb file to read from:
        pdb_file = self.pdbfile
        if len(self.tmp_pdb) > 0:
            pdb_file = self.tmp_pdb[-1]

        found_res = False
        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                if 'ATOM' in line or 'HETATM' in line:
                    residuenr = int(line[22:26])

                    if residuenr == res_nr:
                        found_res = True

                    if found_res:
                        if residuenr != res_nr:
                            break

                        atom = line[13:13+atom_len]
                        if atom == name:
                            nr += 1

        return '%s%d' % (name, nr)

    def change_atom(self):
        """
        Function to change the selected atom in the PDB to atom type selected in listbox
        :return:
        """
        if len(self.buildlist.curselection()) < 1:
            return

        atom = self.buildlist.get(self.buildlist.curselection()[0]).split()[-1]

        print('Changing atom to %s' % atom)
        print self.selected_atoms

        p_atoms = self.get_p_atoms()
        atomnames = self.get_atom_names(int(p_atoms[0]['residuenr']))

        if atom in atomnames.keys():
            p_atoms[0]['name'] = '%s%d' % (atom, (atomnames[atom] + 1))
        else:
            p_atoms[0]['name'] = '%s%d' % (atom, 1)


        #Write pdb with new atom:
        old_pdb, new_pdb = self.get_old_new_pdb_name()

        #Add new pdb file to history
        self.add_pdb_history(new_pdb)

        new = open(new_pdb, 'w')

        with open(old_pdb, 'r') as old:
            for line in old:
                if 'ATOM' in line or 'HETATM' in line:
                    atomnr = line.split()[1]
                    if atomnr == self.selected_atoms[0]:
                        resnr = p_atoms[0]['residuenr']
                        res = p_atoms[0]['residue']
                        x, y, z = p_atoms[0]['xyz'][0:]
                        name = p_atoms[0]['name']
                        new.write('ATOM  %5d  %4s%4s %4d    %8.3f%8.3f%8.3f\n' %
                                      (int(atomnr), name.ljust(4), res.ljust(4), resnr, x, y, z))
                    else:
                        new.write(line)
                else:
                    new.write(line)

        new.close()

        #make selected atoms zero:
        del self.selected_atoms[:]
        self.listbox.selection_clear(0, END)

        #update pymol with new tmp pdb file
        self.update_pymol_structure(new_pdb)
        self.add_atoms_to_list(new_pdb)

    def get_residue_nr(self):
        """
        Opens current PDB file and finds out how many residues are in file.
        :return: resnr: new residue nr (+1 of max in pdb)
        """
        pdb_file = self.pdbfile
        if len(self.tmp_pdb) > 0:
            pdb_file = self.tmp_pdb[-1]

        resnr = 0

        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                print line
                if 'ATOM' in line or 'HETATM' in line:
                    resnr = int(line[22:26])

        resnr += 1

        return resnr

    def build_atom(self, atomname=None, p_atoms=None, bond=None, angle=None, torsion=None):
        """
        Will add atom to atom 1 with bond length X-1 and
        with optional angle X-1-2 or torsion X-1-2-3.
        Changes are temporary saved to tmp1.pdb - tmp20.pdb for undo function to work.

        Calls build.py BuildByAtom(flag,p1,p2=None,p3=None,p4=None,d=None,thetaA=None,thetaT=None):
        p1 = {'atomnnr', 'xyz', 'name'}}
        """

        new_atom = dict()
        if p_atoms[0]:
            #Adapt properties from p1 atoms:
            new_atom[1] = deepcopy(p_atoms[0])
            atomnames = self.get_atom_names(int(p_atoms[0]['residuenr']))
        else:
            #No atoms selected. Get residue nr
            resnr = self.get_residue_nr()

            new_atom[1] = {'name':atomname+'1', 'atomnr': 1, 'xyz': [0, 0, 0], 'residue': 'UNK', 'residuenr': resnr}

            self.write_tmp_pdb(new_atom)
            return

        #Get xyx dict from build
        q = bld.BuildByAtom('a', p_atoms[0], p_atoms[1], p_atoms[2], p_atoms[3], bond, angle, torsion)

        new_atom[1]['xyz'] = q['xyz']
        #atomtype = ''.join([i for i in deepcopy(q['name']) if not i.isdigit()])
        if atomname in atomnames.keys():
            new_atom[1]['name'] = '%s%d' % (atomname, (atomnames[atomname] + 1))

        #Remove potential None from p_atoms to get last p_atom
        p_atoms = filter(lambda x: x!=None, p_atoms)

        self.write_tmp_pdb(new_atom, p_atoms[-1]['atomnr'])


    def build_fragment(self, fragment=None, p_atoms=None, bond=None, angle=None, torsion=None):
        """
        :param fragment:
        :param p_atoms:
        :param bond:
        :param angle:
        :param torsion:
        :return:
        """
        fragmentlib = '%s/Qmods/fragments.dat' % self.app.qgui_path

        #TODO
        #Insert fragment after P1 residue. Renumber residue number and atom numbers.
        new_atoms = dict()

        #Create a dicitonary with atomnames {C:3, H:6,...}
        atomnames = dict()

        #Get dummy atom position
        if not p_atoms[0]:
            # No atoms selected. Get residue nr
            resnr = self.get_residue_nr()
            p_atoms[0] = {'name':'DUM', 'atomnr': 0, 'xyz': [0, 0, 0], 'residue': 'UNK', 'residuenr': resnr}
            bond=0.1
        else:
            atomnames = self.get_atom_names(p_atoms[0]['residuenr'])

        #Append to existing residue number or add new residue number?
        #TODO

        #{symbol:atom  xyz:[x,y,z]  name:C1}
        frag = bld.BuildByGroup('a', p_atoms[0], p_atoms[1], p_atoms[2], p_atoms[3], bond, angle, torsion, fragment,
                                fragmentlib)

        atomnr = deepcopy(p_atoms[0]['atomnr'])
        for atom in sorted(frag.keys()):
            new_atoms[atom] = deepcopy(p_atoms[0])
            atomnr += 1
            atomtype = ''.join([i for i in deepcopy(frag[atom]['name']) if not i.isdigit()])

            if atomtype in atomnames.keys():
                atomnames[atomtype] += 1
                atomname = '%s%d' % (atomtype, atomnames[atomtype])
            else:
                atomname = deepcopy(frag[atom]['name'])

            new_atoms[atom]['name'] = atomname
            new_atoms[atom]['xyz'] = deepcopy(frag[atom]['xyz'])
            new_atoms[atom]['atomnr'] = atomnr
            print frag[atom]['xyz']
            print new_atoms[atom]['xyz']

        print new_atoms

        self.write_tmp_pdb(new_atoms)


    def get_atom_names(self, res_nr):
        """
        Reads current pdb file and extracts highest number for all atom names in residue number.
        :param res_number: residue number to get atom names from
        :return: atomnames: dictionary {C:3, H:6,....}
        """
        #Get correct pdb file to read from:
        pdb_file = self.pdbfile
        if len(self.tmp_pdb) > 0:
            pdb_file = self.tmp_pdb[-1]

        atomnames = dict()
        found_res = False
        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                if 'ATOM' in line or 'HETATM' in line:
                    if int(line[22:26]) == int(res_nr):
                        found_res = True

                    if found_res:
                        if int(line[22:26]) != int(res_nr):
                            break

                        atom = ''.join([i for i in line[11:17].strip() if not i.isdigit()])
                        if not atom in atomnames.keys():
                            atomnames[atom] = 1
                        else:
                            atomnames[atom] += 1

        return atomnames

    def adjust_group(self):
        #TODO
        pass


    def get_fragments(self):
        """
        :return: dictionary with entries in lib files defined in settings
        {resname: {atomnames, types, bonds, head, tail}}
        #TODO this is only temporary... will read from Qmods/fragments.dat
        """
        fragments = dict()

        fragment_files = ['%s/Qmods/fragments.dat' % self.app.qgui_path]

        for f in fragment_files:
            with open(f,'r') as frags:
                for line in frags:
                    if line.startswith('#'):
                        _type = line.split()[0].strip('#')
                        if _type not in fragments:
                            fragments[_type] = list()
                        fragments[_type].append(line.split()[1])

        return fragments

    def get_atoms(self, atom_file):
        """
        This function reads a text file (atom_list) with atomic information and returns a atom dictionary
        :param atom_list: text file with atom numbers in 1st column, mass in 2nd, name in 3rd and symbol in 4th.
        :return: dictionary {atomnumber: {mass, name, symbol}}
        """
        atoms_dict = dict()
        self.fragments_dict['ATOMS'] = list()

        with open(atom_file) as atoms:
            for line in atoms:
                if len(line.split()) > 6 and not line.startswith('#'):
                    atomnumber = int(line.split()[0])
                    mass = float(line.split()[1])
                    name = line.split()[2]
                    symbol = line.split()[3]
                    cov_r = float(line.split()[6])

                    self.fragments_dict['ATOMS'].append(symbol)

                    atoms_dict[symbol] = dict()
                    atoms_dict[symbol]['mass'] = mass
                    atoms_dict[symbol]['name'] = name
                    atoms_dict[symbol]['atomnumber'] = atomnumber
                    atoms_dict[symbol]['covalent r'] = cov_r

        return atoms_dict

    def fill_build_list(self, _type):
        """
        Fill self.buildlist with atoms or fragments.
        :param: dictionary with integers as keys {name: {atomnumber, mass, name, symbol}} or {number:{symbol, atoms}}
        :return: nothing
        """
        self.buildlist.delete(0, END)

        for i in sorted(self.fragments_dict[_type]):
            self.buildlist.insert(END, '%5s' % i)

        self.buildlist.selection_set(0)

    def select_mode_changed(self, *args):
        """
        Change the pymol selectmode
        """
        pymol_value = 1
        if self.select_mode.get() == 'Atoms':
            pymol_value = 0

        self.session.stdin.write('set mouse_selection_mode, %d\n' % pymol_value)

    def p_group_selection_changed(self, *args):
        """
        Change selection mode in pymol and update edit functionalities
        :param args:
        :return:
        """
        modes = {'Atoms': 0, 'Residues': 1, 'Molecules': 5, 'P-atoms': 0}

        if self.session:
            if self.p_atoms_selection.get() != 'P-atoms':
                self.modify_atoms = True
            else:
                self.modify_atoms = False

            self.session.stdin.write('set mouse_selection_mode, %d\n' % modes[self.p_atoms_selection.get()])

    def start_pymol(self):
        self.app.pymol_running = True
        self.reopen_button.config(state=DISABLED)
        for entry in self.widgetList:
            entry.config(state=NORMAL)

        pymol_log = self.app.workdir+'/.tmpfile'
        tmpfile = open(pymol_log,'wb')

        if 'darwin' in sys.platform:
            try:
                self.session = Popen(["pymol", "-p -x -i"], stdout=tmpfile, stdin=PIPE,
                                     preexec_fn=os.setsid)
            except:
                print 'pymol not found'

        else:
            try:
                self.session = Popen(["pymol", "-p"], stdout=tmpfile, stdin=PIPE, preexec_fn=os.setsid)
            except:
                print 'No pymol version found'

        self.update()
        len_log = 35

        self.session.stdin.flush()

        #Update select mode (Atoms, residues)
        self.select_mode_changed()

        self.session.stdin.write('load %s, mol\n' % self.pdbfile)

        #Load default Q-PyMol settings:
        for settings in self.pymol_settings:
            self.session.stdin.write('%s\n' % settings)

        #Check if other special start-up settings are defined:
        if self.start_cmd:
            for command in self.start_cmd:
                self.session.stdin.write('%s\n' % command)
            self.start_cmd = None

        #Check if file has changed before bottering to read it.
        log_size = os.stat(pymol_log).st_size
        while self.session.poll() is None:
            try:
                if self.session:
                    self.update()
                else:
                    break
            except:
                break

            if os.stat(pymol_log).st_size != log_size:

                with open(pymol_log, 'r') as pymol_out:
                    lines = 0
                    for line in pymol_out:
                        lines += 1
                        if lines > len_log:
                            len_log = lines
                            if 'You clicked' in line:
                                self.app.log(' ', line)
                                #Get atom numbers for selected atoms:
                                atomnumbers = list()
                                self.session.stdin.write('identify sele\n')
                                time.sleep(0.2)
                                for rline in reversed(open(self.app.workdir + '/.tmpfile').readlines()):
                                    if 'cmd.identify' in rline:
                                        print rline
                                        atom = rline.split('id ')[-1]
                                        atom = atom.strip(')\n')
                                        atomnumbers.append(atom)
                                    elif 'identify sele' in rline:
                                        break

                                if self.pymol_mode.get() == 'Edit':
                                    if self.p_atoms_selection.get() == 'P-atoms':
                                        for existing in self.selected_atoms:
                                            if existing in atomnumbers:
                                                del atomnumbers[atomnumbers.index(existing)]
                                            else:
                                                atomnumbers.append(existing)

                                self.set_selection(atomnumbers)

                #Change log_size to current value
                log_size = os.stat(pymol_log).st_size

            time.sleep(0.2)

        self.reopen_button.config(state=NORMAL)

        for entry in self.widgetList:
            entry.config(state=DISABLED)

        self.app.pymol_running = False

    def add_atoms_to_list(self, pdbfile):
        """
        Add atoms from pdbfile to selection list.
        """
        self.listbox.delete(0, END)

        with open(pdbfile) as pdb:
            for line in pdb:
                if 'ATOM' in line or 'HETATM' in line:
                    try:
                        atomnr = int(line.split()[1])
                        atomname = line[12:17].strip()
                        resname = line[17:21].strip()
                        resnr = int(line[21:26])

                        self.listbox.insert(END, '%6d %4s %4s %5d' % (atomnr, atomname, resname, resnr))
                    except:
                        continue

    def close_pymol(self):
        try:
            os.killpg(self.session.pid, signal.SIGTERM)
            self.reopen_button.config(state=NORMAL)
            self.app.pymol_running = False
        except:
            self.app.pymol_running = False
            self.destroy()

        #Remove pdb temp directory
        if os.path.exists(self.tmp_pdb_dir):
            shutil.rmtree(self.tmp_pdb_dir)

    def set_selection(self, atoms='all'):
        """
        Highlights selection from pymol in atomlist.
        """
        if atoms == 'all':
            self.listbox.selection_set(0, END)
            return

        tmp_list = self.listbox.get(0, END)
        self.listbox.selection_clear(0, END)

        zoom_to = None
        for index_ in range(len(tmp_list)):
            if tmp_list[index_].split()[0] in atoms or int(tmp_list[index_].split()[0]) in atoms:
                self.listbox.selection_set(index_)
                if not zoom_to:
                    zoom_to = index_

        self.listbox.yview(zoom_to)

        if self.pymol_mode.get() == 'Edit' and self.p_atoms_selection.get() == 'P-atoms':
            self.atomlist_event()

    def clear_list(self):
        """
        Clears all selection from list
        """
        if self.pymol_mode.get() == 'Edit':
            if self.session:
                self.session.stdin.write('select none\n')
                self.session.stdin.write('hide labels\n')
                for atom in self.selected_atoms:
                    self.session.stdin.write('hide spheres, id %s\n' % atom)
                if self.pseudo_atom:
                    self.session.stdin.write('delete new_atom\n')
                    self.pseudo_atom = False

        del self.selected_atoms[:]
        self.listbox.select_clear(0, END)
        self.regulate_edit_controls()

    def invert_selection(self):
        """
        Unselects current items, and selects all other
        """
        selection = map(int, self.listbox.curselection())
        self.listbox.select_clear(0, END)

        self.listbox.select_set(0, END)
        for index_ in selection:
            self.listbox.select_clear(index_)

    def get_selected_atoms(self):
        """
        Gets selected atoms from atom list with their atom numbers.
        Returns list ['id 1-123','id 321-421'] f.ex
        """
        #Get selected atoms from list:
        selections = map(int, self.listbox.curselection())

        atoms = []
        for selection in selections:
            atoms.append(self.listbox.get(selection).split()[0])

        pymol_sel_list = []
        start = False
        for atom in range(len(atoms)):
            if not start:
                first_atom = int(atoms[atom])
                start = True
            if start:
                #Make sure list has not ended:
                if (atom + 1) != len(atoms):
                    if int(atoms[atom + 1]) - int(atoms[atom]) > 1:
                        last_atom = int(atoms[atom])
                        pymol_sel_list.append('id %d-%d' % (first_atom, last_atom))
                        start = False
                else:
                    #End of atom list reached
                    last_atom = int(atoms[atom])
                    pymol_sel_list.append('id %d-%d' % (first_atom, last_atom))

        return pymol_sel_list

    def pymol_command(self, cmd='', sele='all'):
        """
        General pymol commands
        """
        #Get selected atoms from list:
        selections = map(int, self.listbox.curselection())

        if len(selections) > 0:
            sel_list = self.get_selected_atoms()

            sele = sel_list[0]

            for i in range(1, len(sel_list)):
                sele += ' or %s' % sel_list[i]

            self.session.stdin.write('%s %s\n' % (cmd, sele))
        else:
            if sele != '':
                self.session.stdin.write('%s %s\n' % (cmd, sele))
            else:
                self.session.stdin.write('%s\n' % cmd)
            return

    def pymol_highlight_edit(self, atoms):

        if not self.session:
            return

        #Clear selection
        self.session.stdin.write('select none\nhide labels\n')
        self.session.stdin.write('hide spheres, all\n')

        nr = 0
        for i in range(len(atoms)):
            nr += 1
            self.session.stdin.write('show spheres, id %s\n' % atoms[i])
            self.session.stdin.write('set sphere_scale, 0.3, id %s\n' % atoms[i])
            self.session.stdin.write('set sphere_transparency, 0.3, id %s\n' % atoms[i])

            color='cyan'
            if i == 0:
                color='yellow'
            self.session.stdin.write('set sphere_color, %s, id %s\n' % (color, atoms[i]))

            self.session.stdin.write('label id %s, %d\n' % (atoms[i], nr))

    def pymol_special(self, cmd=''):
        """
        Special types pymol command
        """
        selections = map(int, self.listbox.curselection())
        if len(selections) > 0:
            sele = self.get_selected_atoms()
        else:
            sele = ['']

        #BALL & STICKS
        if cmd.split()[1] == 'b&s':
            for idrange in range(len(sele)):
                if cmd.split()[0] == 'show':
                    scale = '0.3'
                else:
                    scale = '0.8'
                self.session.stdin.write('set sphere_scale, %s, %s\n' % (scale, sele[idrange]))
                self.session.stdin.write('%s spheres, %s\n' % (cmd.split()[0], sele[idrange]))
                self.session.stdin.write('%s sticks, %s\n' % (cmd.split()[0], sele[idrange]))
        #CUTEMOL
        elif cmd.split()[1] == 'QuteMol':
            cmd_list = ['set_color "oxygen", [1.0,0.4,0.4],',
                        'set_color "nitrogen", [0.5,0.5,1.0],',
                        'show spheres,',
                        'util.cbaw, ', 'set sphere_scale, 1,']

            atom_independent = ['unset depth_cue,',
                                'bg_color white,',
                        'set light_count, 10,',
                        'set spec_count, 1,', 'set shininess, 10,',
                        'set specular, 0.25,', 'set ambient, 0,',
                        'set direct, 0,', 'set reflect, 1.5,',
                        'set ray_shadow_decay_factor, 0.1,',
                        'set ray_shadow_decay_range, 2,',
                        'util.cbaw']
            self.session.stdin.write('reinitialize settings\n')


            if cmd.split()[0] == 'show':
                for pycmd in atom_independent:
                    self.session.stdin.write('%s\n' %pycmd)
                for idrange in range(len(sele)):
                    for command in cmd_list:
                        self.session.stdin.write('%s %s\n' % (command, sele[idrange]))
            else:
                self.pymol_special('hide b&s')

        #PlastRay
        elif cmd.split()[1] == 'PlastRay':
            self.session.stdin.write('reinitialize settings\n')
            for settings in self.pymol_settings:
                    self.session.stdin.write('%s\n' % settings)
            if cmd.split()[0] == 'show':

                    atom_independent = ['unset depth_cue,', 'set light_count, 10,', 'set shininess, 10,',
                                        'set specular, 0.25,', 'set ambient, 0,',
                                        'set direct, 0,', 'set reflect, 1.5,',
                                        'set ray_shadow_decay_factor, 0.1,',
                                        'set ray_shadow_decay_range, 2,']
                    for pycmd in atom_independent:
                        self.session.stdin.write('%s\n' %pycmd)

        #Default Q-pymol settings:
        elif cmd.split()[1] == 'Q-pymol':
            self.session.stdin.write('reinitialize settings\n')
            if cmd.split()[0] == 'show':
                for settings in self.pymol_settings:
                    self.session.stdin.write('%s\n' % settings)
            else:
                self.session.stdin.write('util.cbag\n')

        #HIGHLIGHT CHARGED RESIDUES
        elif cmd.split()[1] == 'Charged':
            negative_id = []
            negative = ['AS-', 'GL-']
            positive_id = []
            positive = ['AR+', 'LY+']
            cation_id = []
            cations = ['Ca2']
            anion_id = []
            anions = ['Cl-']


            self.session.stdin.write('reinitialize settings\n')
            if cmd.split()[0] == 'show':
                for settings in self.pymol_settings:
                    self.session.stdin.write('%s\n' % settings)
                self.session.stdin.write('hide all\nshow lines\ncolor gray, name C*\n')
                list_items = self.listbox.get(0, END)
                for element in range(len(list_items)):
                    resname = list_items[element].split()[2]
                    atom_id = list_items[element].split()[3]
                    if resname in negative:
                        self.listbox.selection_set(element)
                        self.listbox.yview(element)
                        if atom_id not in negative_id:
                            negative_id.append(atom_id)
                    if len(resname) == 4 and resname[0] == 'C':
                        self.listbox.selection_set(element)
                        self.listbox.yview(element)
                        if atom_id not in negative_id:
                            negative_id.append(atom_id)
                    if resname in anions:
                        self.listbox.selection_set(element)
                        self.listbox.yview(element)
                        if atom_id not in anion_id:
                            anion_id.append(atom_id)
                    if resname in positive:
                        self.listbox.selection_set(element)
                        self.listbox.yview(element)
                        if atom_id not in positive_id:
                            positive_id.append(atom_id)
                    if len(resname) == 4 and resname[0] == 'N':
                        self.listbox.selection_set(element)
                        self.listbox.yview(element)
                        if atom_id not in positive_id:
                            positive_id.append(atom_id)
                    if resname in cations:
                        self.listbox.selection_set(element)
                        self.listbox.yview(element)
                        if atom_id not in cation_id:
                            cation_id.append(atom_id)

                all_nrs = [negative_id, anion_id, positive_id, cation_id]

                if len(negative_id) > 0:
                    red_sticks = 'i. %s ' % [negative_id[0]]
                    for i in range(1, len(negative_id)):
                        red_sticks += 'or i. %s ' % negative_id[i]
                    self.session.stdin.write('show sticks, %s\n' % red_sticks)
                    self.session.stdin.write('color red, %s\n' % red_sticks)
                if len(anion_id) > 0:
                    red_spheres = 'i. %s ' % anion_id[0]
                    for i in range(1, len(anion_id)):
                        red_spheres += 'or i. %s ' % anion_id[i]
                    self.session.stdin.write('show spheres, %s\n' % red_spheres)
                    self.session.stdin.write('color red, %s\n' % red_spheres)
                if len(positive_id) > 0:
                    blue_sticks = 'i. %s ' % positive_id[0]
                    for i in range(1, len(positive_id)):
                        blue_sticks += 'or i. %s ' % positive_id[i]
                    self.session.stdin.write('show sticks, %s\n' % blue_sticks)
                    self.session.stdin.write('color blue, %s\n' % blue_sticks)
                if len(cation_id) > 0:
                    blue_spheres = 'i. %s ' % cation_id[0]
                    for i in range(1, len(cation_id)):
                        blue_spheres += 'i. %s ' % cation_id[i]
                    self.session.stdin.write('show spheres, %s\n' % blue_spheres)
                    self.session.stdin.write('color blue, %s\n' % blue_spheres)

            else:
                self.session.stdin.write('util.cbag\n')
                self.session.stdin.write('hide sticks\n')


        #'set cartoon_side_chain_helper, 1'

    def pymol_color_atoms_selected(self, color, atomtype):
        """
        Colors selected atoms and atomtypes with color
        """

        if atomtype == 'All':
            atomtype = ''

        if atomtype == 'BG':
            self.session.stdin.write('bg_color %s\n' % color)
            return

        selections = map(int, self.listbox.curselection())

        if len(selections) > 0:
            sel_list = self.get_selected_atoms()

            sele = sel_list[0] + ' and name %s*' % atomtype

            for i in range(1, len(sel_list)):
                sele += ' or (%s and name %s*)' % (sel_list[i], atomtype)

            self.session.stdin.write('color %s,  %s\n' % (color, sele))
        else:
            self.session.stdin.write('color %s, name %s*\n' % (color, atomtype))
            return

    def change_bond_order(self, order):
        """
        Change the bond order between selected pair, or generate bond
        between pair.
        """
        selections = map(int, self.listbox.curselection())

        if len(selections) == 2:
            atom1 = self.listbox.get(selections[0]).split()[0]
            atom2 = self.listbox.get(selections[1]).split()[0]
            self.session.stdin.write('bond id %s, id %s, 1\n' % (atom1, atom2))
            self.session.stdin.write('unbond id %s, id %s\n' % (atom1, atom2))
            if order != 0:
                self.session.stdin.write('bond id %s, id %s, %d\n' % (atom1, atom2, order))
        else:
            self.app.errorBox('Warning', 'Select exactly two atoms to change bond order!')
            return

    def fix_geometry(self):
        """
        Fix the geometry for selected atoms:
        """
        selections = map(int, self.listbox.curselection())
        if len(selections) > 998:
            self.app.errorBox('Warning', 'Too many atoms selected. PyMol can at best do 999 atoms.')
            return
        if len(selections) < 1:
            return

        sel_list = self.get_selected_atoms()
        sele = sel_list[0]
        for i in range(1, len(sel_list)):
            sele += ' or %s' % sel_list[i]

        self.session.stdin.write('clean %s\n' % sele)

    def load_dcd(self):
        dcd_files = list()
        trjs = askopenfilenames(parent = self, initialdir = self.app.workdir,
                                  filetypes=(("Trajectory", "*.dcd"),("All files","*.*")))

        if trjs:
            for dcd in trjs:
                dcd_files.append(dcd)

        else:
            return

        for dcd in dcd_files:
            self.session.stdin.write('load_traj %s, %s\n' % (dcd, 'mol' )) #self.pdbfile.split('/')[-1].split('.')[0]))

        self.session.stdin.write('set internal_gui, 1\n')

    def regulate_edit_controls(self):
        """
        Function to control what buttons and boxes to be active in edit mode
        """
        self.modify_atoms = False
        #Clear all spinboxes:
        spinboxes = [self.bond_exist, self.angle_exist, self.torsion_exist,
                     self.bond_new, self.angle_new, self.torsion_new]

        for spinbox in spinboxes:
            spinbox.delete(0,END)
            spinbox.config(state=DISABLED)

        if self.p_atoms_selection.get() != 'P-atoms':
            self.modify_atoms = True

    def build_mode_changed(self, *args):
        """
        Update fragments list with correct value when changed
        """

        self.fill_build_list(self.atoms_fragments.get())

        
    def show_var_frame(self, *args):

        frames = {'Display': self.display_frame,
                  'Edit': self.edit_frame}

        for i in frames.keys():
            frames[i].grid_forget()
        try:
            frames[self.pymol_mode.get()].grid(row=1, column=1, pady=10, padx=(10, 10), sticky='ns')
        except:
            pass

        self.listbox.selection_clear(0, END)
        del self.selected_atoms[:]

        #Changing to Edit mode
        if self.pymol_mode.get() == 'Edit':
            #Create pdb temp directory:
            if not os.path.exists(self.tmp_pdb_dir):
                os.makedirs(self.tmp_pdb_dir)

            self.select_mode.set('Atoms')
            self.p_atoms_selection.set('P-atoms')
            self.listbox.config(selectmode=MULTIPLE)
            self.regulate_edit_controls()

            if len(self.buildlist.curselection()) == 0:
                self.buildlist.selection_set(0)

        #Changing to Display mode
        else:
            self.listbox.config(selectmode=EXTENDED)

        #Common operations if pymol is active:
        if self.session:
            if self.pseudo_atom:
                self.session.stdin.write('delete new_atom\n')
                self.pseudo_atom = False

            self.session.stdin.write('select none\n')
            self.session.stdin.write('hide spheres\nhide labels\n')

            #Update select mode (Atoms, residues)
            self.select_mode_changed()

        self.update()

    def get_pdb_atoms(self, atoms):
        """
        Reads pdb file and returns lines for atoms
        """
        pdb = list()

        #Make lists same length (also to put entry at correct place in new list)
        for i in range(len(atoms)):
            pdb.append(None)

        pdb_file = self.pdbfile
        if len(self.tmp_pdb) > 0:
            pdb_file = self.tmp_pdb[-1]

        with open(pdb_file) as pdbfile:
            for line in pdbfile:
                if 'ATOM' in line or 'HETATM' in line:
                    atom = line.split()[1]
                    if atom in atoms:
                        pdb[atoms.index(atom)] = line
                #No need to read the rest of the file if all atoms are found:
                if not None in pdb:
                    break

        return pdb

    def meassure(self, atoms):
        """
        Use the meassure function
        """
        #Get xyz for atoms 1-2-3-4
        pdb = self.get_pdb_atoms(atoms)

        xyz = list()
        #xyz = map(float, line[30:].split()[0:3])
        for i in pdb:
            xyz.append(map(float, i[30:].split()[0:3]))

        xyz = np.array(xyz)

        #Get bond
        r = bld.measure(xyz[0], xyz[1])
        print 'BOND = %f' % r

        self.update_spinbox(self.bond_exist, r)

        if len(xyz) > 2:
            angle = bld.measure(xyz[0], xyz[1], xyz[2])
            print 'ANGLE = %f' % angle
            self.update_spinbox(self.angle_exist, angle)

        if len(xyz) > 3:
            torsion = bld.measure(xyz[0], xyz[1], xyz[2], xyz[3])
            print 'TORSION = %f' % torsion
            self.update_spinbox(self.torsion_exist, torsion)

        if len(xyz) < 4:
            self.suggest_buld(atoms, xyz)


    def guess_atom_type(self, pdbline):
        """
        Takes a line from a pdb file, reads the pdb atomname and returns the atom type (H, C, N... etc.)
        """
        atomtype = False

        atom = ''.join(i for i in pdbline[13:17].strip() if not i.isdigit())

        #If atom symbol consist of two letters, convert first to capital and second to lower.
        if len(atom) > 1:
            atom = atom[0].upper()+atom[1].lower()

        #Check if atom exist:
        if atom not in self.atom_dict.keys():
            if atom[0] in self.atom_dict.keys():
                atomtype = atom[0]

        elif atom in self.atom_dict.keys():
            atomtype = atom

        return atomtype

    def suggest_bond_length(self):
        """
        Based on atoms connected, suggest a proper bond lenth
        self.atoms_dict[atomsymbol]['covalent r'] = cov_r
        """

        pdb_lines = self.get_pdb_atoms(self.selected_atoms)
        p1_atom = self.guess_atom_type(pdb_lines[0])

        p1_r = self.atom_dict[p1_atom]['covalent r']

        #Get new atom bonding to p1
        if self.atoms_fragments.get() == 'ATOMS':
            new_atom = self.buildlist.get(self.buildlist.curselection()).strip()

        else:
            #TODO read from fragment list dummy atom type and get atom
            new_atom = 'C'

        new_atom_r = self.atom_dict[new_atom]['covalent r']

        bond = p1_r + new_atom_r

        return bond

    def suggest_buld(self, atoms, xyz):
        #TODO based on connectivety of atoms 1(-2-3) suggest values for bond/angle/torsion
        #need bond/angle/torsion list (dict) for atoms
        #need to think about how to handle fragments... not all can be grown, so they must be placed fex.

        bond, angle, torsion = None, None, None

        #Get selected
        if len(atoms) > 0:
            #Suggest bond-length (X-1-2)
            bond = self.suggest_bond_length()
            self.update_spinbox(self.bond_new, bond)
        else:
            return

        if len(atoms) > 1:
            #Suggest angle (X-1-2)
            angle = 109.50
            self.update_spinbox(self.angle_new, angle)

        if len(atoms) > 2:
            torsion = 180.00
            self.update_spinbox(self.torsion_new, torsion)

        self.show_pseudo_atom(bond, angle, torsion)

    def show_pseudo_atom(self, bond, angle, torsion):
        """
        Display position of new atom as pseudo atom in pymol before actually building it (adjustments and stuff)
        """
        p_atoms = self.get_p_atoms()

        #Get xyz for pseudo atom dict from build.py
        pseudo = bld.BuildByAtom('a', p_atoms[0], p_atoms[1], p_atoms[2], p_atoms[3], bond, angle, torsion)

        xyz = pseudo['xyz']

        if self.pseudo_atom:
            self.alter_atom_xyz('new_atom', xyz)
        else:
            self.session.stdin.write('pseudoatom new_atom, pos=%s\n' % (tuple(xyz),))
            self.pseudo_atom = True

    def update_pseudo_atom(self, *args):
        """
        :param args:
        :return:
        """
        if not self.pseudo_atom:
            return

        bond, angle, torsion = None, None, None

        if len(self.selected_atoms) > 0:
            try:
                bond = float(self.bond_new.get())
            except:
                bond = 1.5
        else:
            return
        if len(self.selected_atoms) > 1:
            try:
                angle = float(self.angle_new.get())
            except:
                angle = 109.5
        if len(self.selected_atoms) > 2:
            try:
                torsion = float(self.torsion_new.get())
            except:
                torsion= 180.0

        p_atoms = self.get_p_atoms()
        #Get xyz for pseudo atom dict from build.py
        pseudo = bld.BuildByAtom('a', p_atoms[0], p_atoms[1], p_atoms[2], p_atoms[3], bond, angle, torsion)

        xyz = pseudo['xyz']
        self.alter_atom_xyz('new_atom', xyz)


    def alter_atom_xyz(self, selection, xyz):
        """
        moves pymol selection to new coordinates xyz:
        """
        self.session.stdin.write('alter_state 1, %s, (x,y,z)=%s\n' % (selection, tuple(xyz),))

    def update_spinbox(self, spinbox, value):
        """
        inserts value in spinbox
        """
        self.modify_atoms = False
        spinbox.config(state=NORMAL)
        spinbox.insert(0, '%.2f' % value)

        if self.p_atoms_selection.get() != 'P-atoms':
            self.modify_atoms = True

    def fragmenlist_event(self, *args):
        """
        Update suggestions for bond/angle/torsion when atom or fragment is changed
        """
        if len(self.selected_atoms) > 0:
            self.regulate_edit_controls()
            self.suggest_buld(self.selected_atoms, [])

    def atomlist_event(self, *args):
        """
        TODO
        """

        #If in display mode, do not do anything here.
        if self.pymol_mode.get() == 'Display' or self.p_atoms_selection.get() != 'P-atoms':
            self.pymol_command('select')
            return

        try:
            selections = map(int, self.listbox.curselection())
        except:
            return

        selected_atoms = list()

        #If more than 4 is selected, remove first existing in self.selected_atoms
        remove_atom = None
        if len(selections) > 4:
            remove_atom = self.selected_atoms[0]
            del self.selected_atoms[0]

        for selected in selections:
            atomnumber = self.listbox.get(selected).split()[0]
            if atomnumber == remove_atom:
                self.listbox.selection_clear(selected, selected)
            elif atomnumber not in self.selected_atoms:
                self.selected_atoms.append(atomnumber)
            selected_atoms.append(atomnumber)

        #Check if user unselected something, and remove that atom
        if len(selected_atoms) < 4:
            for i in self.selected_atoms:
                if i not in selected_atoms:
                    del self.selected_atoms[self.selected_atoms.index(i)]

        #Update existing and new bond/angle/torsion values:
        self.regulate_edit_controls()

        #Update correct selection to pymol
        self.pymol_highlight_edit(self.selected_atoms)

        if len(self.selected_atoms) > 1:
            self.meassure(self.selected_atoms)
        else:
            self.suggest_buld(self.selected_atoms, [])

        #delete pseudo atom if no atoms are selected
        if len(self.selected_atoms) == 0:
            if self.pseudo_atom:
                self.session.stdin.write('delete new_atom\n')
                self.pseudo_atom = False

    #Main window
    def dialog_window(self):
        self.title('Qgui PyMol')
        self.config(background=self.main_color)

        self.widgetList = []

        # Main frame
        frame = Frame(self, bg = self.main_color)
        frame.grid(row=0, column=0, columnspan=2, pady=10, padx=10)

        # Frame with residue/atoms list (static)
        frame2 = Frame(self, bg=self.main_color)
        frame2.grid(row=1,column=0, pady=10, padx=(10,0))

        # Frame with pymol commands (variable)
        self.display_frame = Frame(self, bg=self.main_color)

        # Frame with edit and build commands (variable)
        self.edit_frame = Frame(self, bg=self.main_color)

        # Frame 1
        label_pymol = Label(frame, text='PyMol session (%s)' % self.pdbfile.split('/')[-1],
                            bg=self.main_color)
        label_pymol.grid(row=0, column=0, columnspan=3)

        self.close_button = Button(frame, text='Quit', highlightbackground=self.main_color, command=self.close_pymol)
        self.close_button.grid(row=1,column=0)

        self.reopen_button = Button(frame, text='View again', highlightbackground=self.main_color, command=self.start_pymol)
        self.reopen_button.grid(row=1, column=1)
        
        self.pymol_mode_menu = OptionMenu(frame, self.pymol_mode,
                                   'Display',
                                   'Edit')
        self.pymol_mode_menu.config(highlightbackground=self.main_color, bg=self.main_color, width=15)
        self.pymol_mode_menu.grid(row=2, column=0, columnspan=3)

        self.load_trj = Button(frame, text ='Load trajectory', highlightbackground=self.main_color, command=self.load_dcd)
        self.load_trj.grid(row=1, column=2)

        self.widgetList.append(self.load_trj)

        # frame2 (atom list)
        listbox_scroll = Scrollbar(frame2)
        listbox_scroll.grid(row = 0, rowspan = 10, column = 11, sticky = 'nsw')
        self.listbox = Listbox(frame2, yscrollcommand = listbox_scroll.set, width = 23, height=20,
                               highlightthickness=0, relief=GROOVE, selectmode=EXTENDED, exportselection=False)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 0, rowspan = 10, column = 0, columnspan = 10, sticky = 'w')
        self.listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.listbox.bind('<<ListboxSelect>>', self.atomlist_event)

        sel_all_button = Button(frame2, text='Select all', highlightbackground=self.main_color,
                                command=lambda: self.set_selection('all'))
        sel_all_button.grid(row=11, column=0, columnspan=5)

        clear_button = Button(frame2, text='Clear', highlightbackground=self.main_color, command=self.clear_list)
        clear_button.grid(row=11, column=6, columnspan=5)

        invert_button = Button(frame2, text='Invert selection', highlightbackground=self.main_color,
                               command=self.invert_selection)
        invert_button.grid(row=12, column=0, columnspan=11)

        #Get types from fragments.dat
        bld_types = list()
        for _type in sorted(self.fragments_dict.keys()):
            bld_types.append(_type)

        # edit_frame
        atoms_fragments = OptionMenu(self.edit_frame, self.atoms_fragments, *bld_types)
        atoms_fragments.config(bg=self.main_color, highlightbackground=self.main_color, width=12)
        atoms_fragments.grid(row=0, column=0, columnspan=2)


        p_group_selector = OptionMenu(self.edit_frame, self.p_atoms_selection,
                                      *('P-atoms', 'Atoms', 'Residues', 'Molecules'))
        p_group_selector.config(bg=self.main_color, highlightbackground=self.main_color, width=12)
        p_group_selector.grid(row=0, column=2)

        buildlist_scroll = Scrollbar(self.edit_frame)
        buildlist_scroll.grid(row=1, rowspan=5, column=1, sticky='nsw')
        self.buildlist = Listbox(self.edit_frame, yscrollcommand=buildlist_scroll.set, width=13, height=10,
                                 highlightthickness=0, relief=GROOVE, selectmode=SINGLE, exportselection=False)
        buildlist_scroll.config(command=self.buildlist.yview)
        self.buildlist.grid(row=1, rowspan=5, column=0, sticky='e')
        self.buildlist.config(font=tkFont.Font(family="Courier", size=12))
        self.buildlist.bind('<<ListboxSelect>>', self.fragmenlist_event)

        #information about existing atoms:
        exist_frame_label = LabelFrame(self.edit_frame, text='Existing', bg=self.main_color)
        exist_frame_label.grid(sticky='ns', row=1, rowspan=5, column=2)

        bond_label = Label(self.edit_frame, text='Bond', bg=self.main_color)
        bond_label.grid(in_=exist_frame_label, column=2, row=1)

        angle_label = Label(self.edit_frame, text='Angle', bg=self.main_color)
        angle_label.grid(in_=exist_frame_label, column=2, row=2)

        torsion_label = Label(self.edit_frame, text='Torsion', bg=self.main_color)
        torsion_label.grid(in_=exist_frame_label, column=2, row=3)

        self.bond_exist = Spinbox(self.edit_frame, width=5, highlightthickness=0, relief=GROOVE,
                            from_=0.00, to=100.00, increment=0.1, format="%.2f", textvariable=self.existing_bond)
        self.bond_exist.grid(in_=exist_frame_label, column=3, row=1)

        self.angle_exist = Spinbox(self.edit_frame, width=5, highlightthickness=0, relief=GROOVE,
                            from_=-360.00, to=360.00, increment=1.0, format="%.2f", textvariable=self.existing_angle)
        self.angle_exist.grid(in_=exist_frame_label, column=3, row=2)

        self.torsion_exist = Spinbox(self.edit_frame, width=5, highlightthickness=0, relief=GROOVE,
                            from_=-360.00, to=360.00, increment=1.0, format="%.2f", textvariable=self.existing_torsion)
        self.torsion_exist.grid(in_=exist_frame_label, column=3, row=3)

        #Information for new atom to be added:
        add_frame_label = LabelFrame(self.edit_frame, text='New', bg=self.main_color)
        add_frame_label.grid(sticky='ns', row=1, rowspan=5, column=3)

        self.bond_new = Spinbox(self.edit_frame, width=5, highlightthickness=0, relief=GROOVE,
                            from_=0.00, to=100.00, increment=0.1, format="%.2f", textvariable=self.new_bond)
        self.bond_new.grid(in_=add_frame_label, column=4, row=1)

        self.angle_new = Spinbox(self.edit_frame, width=5, highlightthickness=0, relief=GROOVE,
                            from_=-360.00, to=360.00, increment=1.0, format="%.2f", textvariable=self.new_angle)
        self.angle_new.grid(in_=add_frame_label, column=4, row=2)

        self.torsion_new = Spinbox(self.edit_frame, width=5, highlightthickness=0, relief=GROOVE,
                            from_=-360.00, to=360.00, increment=1.0, format="%.2f", textvariable=self.new_torsion)
        self.torsion_new.grid(in_=add_frame_label, column=4, row=3)

        self.addbutton = Button(self.edit_frame, text='Add', highlightbackground=self.main_color,
                                command=self.add_build)
        self.addbutton.grid(in_=add_frame_label, row=4, column=4, columnspan=1)

        self.deletebutton = Button(self.edit_frame, text='Delete', highlightbackground=self.main_color,
                                   command=self.delete_atoms)
        self.deletebutton.grid(in_=exist_frame_label, row=4, column=2, columnspan=1)

        self.undobutton = Button(self.edit_frame, text='Undo', highlightbackground=self.main_color,
                                 command=self.undo_build)
        self.undobutton.grid(in_=exist_frame_label, row=4, column=3, columnspan=1)

        self.change_atom_button = Button(self.edit_frame, text='Change atom', highlightbackground=self.main_color,
                                  command=self.change_atom)
        self.change_atom_button.grid(in_=exist_frame_label, row=5, column=2, columnspan=2)



        # display_frame (pymol commands)
        selectmode = Label(self.display_frame, text='Select', bg=self.main_color)
        selectmode.grid(row=0, column=0, sticky='e')

        self.select_box = OptionMenu(self.display_frame, self.select_mode, 'Atoms', 'Residues')
        self.select_box.config(bg=self.main_color, highlightbackground=self.main_color, width=15)
        self.select_box.grid(row=0, column=1, columnspan=2)
        self.widgetList.append(self.select_box)

        self.orient_button = Button(self.display_frame, text='Orient', highlightbackground=self.main_color,
                                    command=lambda: self.pymol_command('orient'))
        self.orient_button.grid(row=1, column=1, sticky='e')
        self.widgetList.append(self.orient_button)

        self.zoom_button = Button(self.display_frame, text='Zoom', highlightbackground=self.main_color,
                                  command=lambda: self.pymol_command('zoom'))
        self.zoom_button.grid(row=1, column = 2, sticky='e')
        self.widgetList.append(self.zoom_button)

        self.colorbutton = Button(self.display_frame, text='Color', highlightbackground=self.main_color,
                                  command=lambda: self.pymol_color_atoms_selected(
                                      self.pymol_color.get(), self.pymol_atoms_color.get()))
        self.colorbutton.grid(row=2, column=0, sticky='e')
        self.widgetList.append(self.colorbutton)

        self.color_atoms = OptionMenu(self.display_frame, self.pymol_atoms_color,
                                   'C', 'N', 'O', 'H', 'S', 'P', 'All', 'BG')
        self.color_atoms.config(highlightbackground=self.main_color, bg=self.main_color, width=5)
        self.color_atoms.grid(row=2, column=2, sticky='w')
        self.widgetList.append(self.color_atoms)

        self.colorbox = OptionMenu(self.display_frame, self.pymol_color,
                                   'Brown', 'Black', 'Blue', 'Cyan', 'Orange', 'Gray', 'Green', 'Pink', 'Purple',
                                   'Red', 'Sand', 'Silver', 'White', 'Yellow')
        self.colorbox.config(highlightbackground=self.main_color, bg=self.main_color, width=8)
        self.colorbox.grid(row=2, column=1, sticky='e')
        self.widgetList.append(self.colorbox)

        #Ball and stick:
        show_bs = Label(self.display_frame, text='B & S', bg=self.main_color)
        show_bs.grid(row=3, column=0, sticky='w')
        self.show_bs = Button(self.display_frame, text= 'On', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_special('show b&s'))
        self.show_bs.grid(row=3, column = 1, sticky='e')
        self.widgetList.append(self.show_bs)

        self.hide_bs = Button(self.display_frame, text= 'Off', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_special('hide b&s'))
        self.hide_bs.grid(row=3, column = 2, sticky='w')
        self.widgetList.append(self.hide_bs)


        showcartoon = Label(self.display_frame, text='Cartoon ', bg=self.main_color)
        showcartoon.grid(row=4, column=0, sticky='w')
        self.show_cartoon = Button(self.display_frame, text= 'On', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_command('show cartoon,'))
        self.show_cartoon.grid(row=4, column = 1, sticky='e')
        self.widgetList.append(self.show_cartoon)

        self.hide_cartoon = Button(self.display_frame, text= 'Off', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_command('hide cartoon,'))
        self.hide_cartoon.grid(row=4, column = 2, sticky='w')
        self.widgetList.append(self.hide_cartoon)

        showlines = Label(self.display_frame, text='Lines', bg=self.main_color)
        showlines.grid(row=5, column=0, sticky='w')

        self.show_lines = Button(self.display_frame, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show lines,'))
        self.show_lines.grid(row=5, column=1, sticky='e')

        self.widgetList.append(self.show_lines)
        self.hide_lines = Button(self.display_frame, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide lines,'))
        self.hide_lines.grid(row=5, column=2, sticky='w')
        self.widgetList.append(self.hide_lines)

        showsphere = Label(self.display_frame, text='Spheres', bg=self.main_color)
        showsphere.grid(row=6, column=0, sticky='w')

        self.show_sphere = Button(self.display_frame, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show spheres,'))
        self.show_sphere.grid(row=6, column=1, sticky='e')
        self.widgetList.append(self.show_sphere)

        self.hide_sphere = Button(self.display_frame, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide spheres,'))
        self.hide_sphere.grid(row=6, column=2, sticky='w')
        self.widgetList.append(self.hide_sphere)

        showsticks = Label(self.display_frame, text='Sticks', bg=self.main_color)
        showsticks.grid(row=7, column=0, sticky='w')

        self.show_sticks = Button(self.display_frame, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show sticks,'))
        self.show_sticks.grid(row=7, column=1, sticky='e')
        self.widgetList.append(self.show_sticks)

        self.hide_sticks = Button(self.display_frame, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide sticks,'))
        self.hide_sticks.grid(row=7, column=2, sticky='w')
        self.widgetList.append(self.hide_sticks)

        showlines = Label(self.display_frame, text='Surface', bg=self.main_color)
        showlines.grid(row=8, column=0, sticky='w')

        self.show_surface = Button(self.display_frame, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show surface,'))
        self.show_surface.grid(row=8, column=1, sticky='e')
        self.widgetList.append(self.show_surface)

        self.hide_surface = Button(self.display_frame, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide surface,'))
        self.hide_surface.grid(row=8, column=2, sticky='w')
        self.widgetList.append(self.hide_surface)

        self.special_presets = OptionMenu(self.display_frame, self.pymol_presets,
                                   'QuteMol', 'Charged', 'PlastRay', 'Q-pymol')
        self.special_presets.config(highlightbackground=self.main_color, bg=self.main_color, width=10)
        self.special_presets.grid(row=9, column=0, sticky='w')
        self.widgetList.append(self.special_presets)

        self.show_preset = Button(self.display_frame, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_special('show %s' % self.pymol_presets.get()))
        self.show_preset.grid(row=9, column=1, sticky='e')
        self.widgetList.append(self.show_preset)

        self.hide_preset = Button(self.display_frame, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_special('hide %s' % self.pymol_presets.get()))
        self.hide_preset.grid(row=9, column=2, sticky='w')
        self.widgetList.append(self.hide_preset)



