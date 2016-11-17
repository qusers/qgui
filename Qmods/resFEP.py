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
from tkFileDialog import askopenfilename
import qgui_functions as qf
from select_return import SelectReturn
from subprocess import Popen, PIPE
import tkFont
import time
import os
import signal
import sys
import numpy as np
import shutil
from copy import deepcopy


class ResFEP(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root

        #Store topologies 1..N:
        self.topology_paths = dict()

        #Store total FEP files for all mutations in each topology
        # [topology][nr] = fepdict <-- see read_fepfile in qgui_functions.py for info on fepdict
        self.topology_fep = dict()

        #Store individual FEP files for every residual mutation in each topology
        # [toplogy][res_nr] = fepdict
        self.topology_res_fep = dict()

        #Remember the order of how Residues are added
        self.topology_res_order = dict()

        #Atoms added manually to residue mutations (f.ex counterions) {topology: res nr: atomname: atomnr}
        self.added_atoms = dict()

        self.selected_topology = StringVar()
        self.selected_topology.set('Topology 1')
        self.topology_paths['Topology 1'] = '*.top'

        self.mutate_to = StringVar()
        self.mutate_from = StringVar()

        self.mutate_from.set('start')
        self.mutate_to.set('end')

        #remember residue numbers to mutate in topologies [topology][resnr] = [WT,MUT]
        self.topology_mutation = dict()
        self.topology_mutation['Topology 1'] = dict()

        #Check if a topology is loaded in main window
        if self.app.top_id:
            self.topology_paths['Topology 1'] = self.app.top_id
            self.insert_topology_name()


        #Trace stuff
        self.selected_topology.trace('w', self.topology_changed)

        self.dialog_window()

        #Collect predifined FEP protocols {res_wt: {res_mut: path}}
        self.feps = self.get_fep_protocols()

    def insert_topology_name(self):
        """
        Write topology name on screen
        """
        topology = self.topology_paths[self.selected_topology.get()]
        topname = topology.split('/')[-1]

        if not self.selected_topology.get() in self.topology_fep.keys():
            self.topology_fep[self.selected_topology.get()] = dict()

        self.topology_label.config(text=topname)

    def topology_changed(self, *args):
        """
        Update values in window when topology is toggeld
        """
        self.insert_topology_name()
        self.refresh_residue_list()
        self.refresh_feplist()

    def refresh_feplist(self):
        """
        Show FEP files for selected topology in FEP file(s) listbox
        """
        self.feplist.delete(0, END)

        for fep in sorted(self.topology_fep[self.selected_topology.get()].keys()):
            self.feplist.insert(END, 'FEP%d' % fep)


    def load_topology(self):
        """
        :return:
        """
        filename = askopenfilename(parent=self, initialdir=self.app.workdir,
                                   filetypes=(("TOP", "*.top"), ("All files", '*.*')))
        if filename != '':
            self.topology_paths[self.selected_topology.get()] = filename
            self.topology_fep[self.selected_topology.get()] = dict()
            self.topology_res_fep[self.selected_topology.get()] = dict()

            self.insert_topology_name()

    def add_topology(self):
        """
        Add a new topology
        :return:
        """
        top_nr = 'Topology %d' % (int(self.selected_topology.get().split()[1]) + 1)

        self.topology_paths[top_nr] = '*.top'
        self.topology_mutation[top_nr] = dict()
        self.topology_res_fep[top_nr] = dict()
        self.topology_res_order[top_nr] = dict()


        self.selected_topology.set(top_nr)
        self.update_topology_menu()

    def del_topology(self):
        """
        Remove selected topology
        :return:
        """
        top = self.selected_topology.get()

        del self.topology_paths[top]
        del self.topology_res_fep[top]
        del self.topology_mutation[top]
        del self.topology_res_order[top]
        nr = int(top.split()[1]) - 1

        if nr < 1:
            self.add_topology()
        else:
            self.update_topology_menu()

    def update_topology_menu(self):
        menu = self.select_topology['menu']
        menu.delete(0, END)

        for top in sorted(self.topology_paths.keys()):
            menu.add_command(label=top, command=lambda v=top: self.selected_topology.set(v))

    def dict_to_list(self, adict):
        """
        Converts the keys and values of a dictionary to a list with ['value  key']
        :param adict:
        :return: list
        """
        new_list = list()
        for i in sorted(adict.keys()):
            v = adict[i]
            i = str(i)
            new_list.append('%8s  %8s' % (v.ljust(8), i))

        return new_list

    def add_residue(self):
        """
        Opens a new window to select residue(s) to mutate

        """
        if self.topology_paths[self.selected_topology.get()] == '*.top':
            print('No topology loaded for %s' % self.selected_topology.get())
            return

        #Get topology residues
        residues = qf.read_topology(self.topology_paths[self.selected_topology.get()])[1]

        #Make dict to list to be read by selector window:
        residues = self.dict_to_list(residues)

        self.select_res = SelectReturn(self, self.root, elements=residues, select_title='Select residue to mutate',
                                       Entry=self.reslist)
        self.select_res.configure(bg=self.main_color)
        self.select_res.title('Select residue to mutate')
        self.select_res.resizable()

    def del_residue(self):
        """
        Delete selected residue from mutation list
        :return:
        """
        selection = self.reslist.curselection()
        if len(selection) < 1:
            return

        top = self.selected_topology.get()
        for i in selection:
            res_nr = int(self.reslist.get(i).split()[1])
            del self.topology_mutation[top][res_nr]
            del self.topology_res_fep[top][res_nr]
            for nr in sorted(self.topology_res_order[top].keys()):
                if self.topology_res_order[top][nr] == res_nr:
                    del self.topology_res_order[top][nr]

            self.reslist.delete(i)

        #Need to create res_fep again, because nr of FEP files might have changed.

        #Make FEP for topology in total:
        self.topology_fep[top] = self.make_topology_fep(self.topology_res_fep[top], self.topology_res_order)
        self.refresh_feplist()

        self.refresh_residue_list()
        self.refresh_feplist()

    def mutate_residue(self, residue):
        """
        :param residue: RES 123
        Checks if RES exist in FEP library and loads FEP files. Abort if not existing!
        """
        top = self.selected_topology.get()
        resname = residue.split()[0]

        res_nr = int(residue.split()[-1])

        #Keeping mutate from as a list() in case we later decide to expand some unforseen options here...
        mutate_from = list()
        mutate_to_options = list()
        mutate_to = None

        #Check if resname exist in FEP protocols
        for fepfrom in self.feps.keys():
            if fepfrom == resname:
                mutate_from.append(fepfrom)
                for fepto in self.feps[fepfrom].keys():
                    mutate_to_options.append(fepto)
                    if not '+' in fepto:
                        mutate_to = fepto

        if len(mutate_from) < 1:
            print('Found no FEP protocol for residue: %s' % resname)
            return

        if not mutate_to:
            mutate_to = mutate_to_options[0]

        if not top in self.topology_res_order.keys():
            self.topology_res_order[top] = dict()

        self.topology_res_order[top][len(self.topology_res_order[top].keys())+1] = res_nr

        self.topology_mutation[self.selected_topology.get()][res_nr] = [resname, mutate_to]

        #Get FEP files for default mutation
        self.get_fep_files(resname, mutate_to, res_nr)

        self.mutate_to.set(mutate_to)
        self.mutate_from.set(resname)

        self.update_mutate_to_options(mutate_to_options)
        self.update_mutate_from_options(mutate_from)

        self.refresh_residue_list()

    def check_fep_top(self, fepfile, res_nr):
        """
        :param fepfile:
        :param topolgy:
        :param res_nr:
        :return:
        """
        pdb_atoms_order = list()
        pdb_atomnumbers_order = list()

        #Get residue atom names and atom number from pdb file:
        pdb_res = qf.get_pdb_resnr(qf.create_pdb_from_topology(self.topology_paths[self.selected_topology.get()]),
                                   res_nr)

        for line in pdb_res:
            if 'ATOM' in line:
                atom_nr = int(line.split()[1])
                atom_name = line[12:17].strip()
                pdb_atoms_order.append(atom_name)
                pdb_atomnumbers_order.append(atom_nr)

        #atom nr and atom names from FEP file
        q_atomnr, q_atomname = qf.get_fep_atoms(fepfile)

        #Remove atoms from pdb that is not present in FEP file:
        for i in range(len(pdb_atoms_order) - 1, -1, -1):
            if pdb_atoms_order[i] not in q_atomname.values():
                pdb_atoms_order.pop(i)
                pdb_atomnumbers_order.pop(i)

        missing_atoms = dict()
        wrong_order = False

        #Check if additional atoms (f.ex counter ions) are missing.
        top = self.selected_topology.get()

        if len(pdb_atoms_order) != len(q_atomname.keys()):
            residues = qf.create_pdb_from_topology(self.topology_paths[self.selected_topology.get()],
                                               self.app.q_settings['library'])
            for q in q_atomname.keys():
                atomtype = q_atomname[q]

                if atomtype not in pdb_atoms_order:
                    #Let us add the missing atom to this FEP residue
                    #TODO we probably do not need self.added_atoms ... delete it later if so!
                    if not top in self.added_atoms.keys():
                        self.added_atoms[top] = dict()
                    if not res_nr in self.added_atoms[top].keys():
                        self.added_atoms[top][res_nr] = dict()
                    self.added_atoms[top][res_nr][atomtype] = False

                    #atomnr = self.select_missing_atom(atomtype)
                    print('Q atom %d %s not found in residue %d ' % (q, atomtype, res_nr))
                    self.app.log(' ', 'Q atom %d %s not found in residue %d \n' % (q, atomtype, res_nr))
                    self.app.log(' ', 'Select atom %s from topology.' % atomtype)
                    atomnr = SelectReturn(self, self.root, elements=residues,
                                       select_title='atom %s for FEP' % atomtype, Entry=self.reslist).show()
                    if not atomnr:
                        print()

                    self.added_atoms[top][res_nr][atomtype] = atomnr
                    pdb_atomnumbers_order.append(self.added_atoms[top][res_nr][atomtype])
                    pdb_atoms_order.append(atomtype)


        #{atomnumber in FEP file: offset value making it the same as TOP atom nr}
        atomoffset = dict()

        #Time to make every Q-atom nr correspond to the actual atom nr in *.top! (also check if something is wrong)
        i = 0

        for q in sorted(q_atomnr.keys()):
            atomtype = q_atomname[q]

            if atomtype == pdb_atoms_order[i]:
                atomoffset[q_atomnr[q]] = int(pdb_atomnumbers_order[i]) - q_atomnr[q]
            else:
                if not wrong_order:
                    #Get the atom order in FEP file
                    fep_atoms_order = list()
                    for qi in sorted(q_atomname.keys()):
                        fep_atoms_order.append(q_atomname[qi])
                    print 'FEP atoms order VS topology atoms order:'
                    print(fep_atoms_order)
                    print(pdb_atoms_order)
                    print('This will most likely be handled by Qgui, but verify final FEP files!')

                wrong_order = True
                if not atomtype in pdb_atoms_order:
                    missing_atoms[q] = atomtype
                else:
                    #Find atom number for atomtype
                    for j in range(len(pdb_atoms_order)):
                        if pdb_atoms_order[j] == atomtype:
                            atomoffset[q_atomnr[q]] = int(pdb_atomnumbers_order[j]) - q_atomnr[q]
                            break
            i += 1

        if len(missing_atoms.keys()) > 0:
            self.app.errorBox('warning', 'Could not find all atoms defined in FEP in topology!')
            for q in sorted(missing_atoms.keys()):
                print('Q-atom %d with atom name %s not found in topolgy.' % (q, missing_atoms[q]))
                self.app.log('', 'Q-atom %d with atom name %s not found in topolgy.\n' % (q, missing_atoms[q]))

            return None

        if wrong_order:
            self.app.errorBox('warning', 'Atom order in FEP file and topology does not match!')
            #Todo make some translateion dictionary insted of using offset to handle this

        return atomoffset


    def get_fep_files(self, rest_wt, res_mut, res_nr):
        """
        Gets FEP file for mutation res_wt-->res_mut for residue nr res_nr
        :param rest_wt:
        :param res_mut:
        :return: nothing, but updates self.topology_res_fep dictionary
        """
        feps = dict()
        top = self.selected_topology.get()

        checked_pdb = False
        atomoffset = None

        for fepfile in os.listdir(self.feps[rest_wt][res_mut]):
            if fepfile.endswith('.fep'):
                nr = int(filter(str.isdigit, fepfile))
                fepfile = '%s/%s' % (self.feps[rest_wt][res_mut], fepfile)
                #Check atom numebers from pdb
                if not checked_pdb:
                    atomoffset = self.check_fep_top(fepfile, res_nr)
                    if not atomoffset:
                        return
                    checked_pdb = True

                feps[nr] = qf.read_fep(fepfile, qoffset=0, atomoffset=atomoffset)

        self.topology_res_fep[top][res_nr] = feps

        #Make FEP for topology in total:
        self.topology_fep[top] = self.make_topology_fep(self.topology_res_fep[top], self.topology_res_order)
        self.refresh_feplist()

    def make_topology_fep(self, res_fep, res_order):
        """
        :param: res_fep {res nr: fep nr: nr: section: stuff}
        :param: res_order {nr: res nr}
        :return: fep_tot
        """

        #Order to merge FEP files:
        merge_sequence = list()
        for nr in sorted(res_order[self.selected_topology.get()].keys()):
            merge_sequence.append(nr)

        if len(merge_sequence) < 1:
            return

        #Initialize FEP with first residue
        res = res_order[self.selected_topology.get()][merge_sequence.pop(0)]
        fep_tot = deepcopy(res_fep[res])

        #Check if any existing mutations are added. If so, we need to merge FEP protocols:
        if len(res_fep.keys()) > 1:
            while len(merge_sequence) > 0:
                res = res_order[self.selected_topology.get()][merge_sequence.pop(0)]
                fep_tot = qf.merge_fep_dicts(fep_tot, res_fep[res])

        return fep_tot

    def refresh_residue_list(self):
        """
        Update residues with mutation in listbox for selected topology
        :return:
        """
        self.reslist.delete(0, END)

        for res_nr in sorted(self.topology_mutation[self.selected_topology.get()].keys()):
            resname = self.topology_mutation[self.selected_topology.get()][res_nr][0]
            mutate_to = self.topology_mutation[self.selected_topology.get()][res_nr][1]
            self.reslist.insert(END, u'%s %d \u2192 %s' % (resname, res_nr, mutate_to))

    def update_mutate_from_options(self, mut_from):
        """
        :param mut_from: list with residues to mutate from
        :return: updated optionmenu for mutate from
        """
        menu = self.mutate_from_menu['menu']
        menu.delete(0, END)

        for res in sorted(mut_from):
            menu.add_command(label=res, command=lambda v=res: self.mutate_from.set(v))

    def update_mutate_to_options(self, mutate_to):
        """
        :param mut_to: list with residues to mutate to
        :return: updated optionmenu for mutate to
        """
        menu = self.mutate_to_menu['menu']
        menu.delete(0, END)

        for res in sorted(mutate_to):
            menu.add_command(label=res, command=lambda v=res: self.mutate_to.set(v))


    def get_fep_protocols(self):
        """
        Get all posible mutations from OPLS/FEP and OPLS/.FEP
        :return: {res: [res]}
        """
        fep_paths = ['%s/FF/OPLS/FEP/' % self.app.qgui_path, '%s/FF/OPLS/.FEP/' % self.app.qgui_path]

        res_feps = dict()

        print('Reading FEP protocols:')
        for place in fep_paths:
            for fep in os.listdir(place):
                fepdir = '%s/%s' % (place, fep)
                if os.path.isdir(fepdir):
                    try:
                        res_wt = fep.split('_')[0]
                        res_mut = fep.split('_')[1]
                        if len(fep.split('_')) > 2:
                            res_mut += '+' + '+'.join(fep.split('_')[2:])

                        if not res_wt in res_feps:
                            res_feps[res_wt] = dict()
                        if res_mut in res_feps[res_wt].keys():
                            print('Found duplicate FEP protocol for %s!' % fep)
                            print('Default version of %s will be used' % fep)
                        res_feps[res_wt][res_mut] = fepdir
                        print('%s-->%s' % (res_wt, res_mut))
                    except:
                        print('Failed to read %s' % fepdir)

        if len(res_feps) < 1:
            print('Found no FEP protocols in /QGUI/FF/OPLS/FEP')

        return res_feps


    def update_residue_mutation(self):
        """
        Update mutation from --> to in listbox
        :return:
        """
        #Get selected residue in listbox:
        selection = self.reslist.curselection()
        if len(selection) == 0:
            print('Nothing selected in Residue(s) to mutate!')
            return

        selected = self.reslist.get(selection[0])

        res_nr = int(selected.split()[1])

        self.reslist.delete(selection[0])

        mut_from = self.mutate_from.get()
        mut_to = self.mutate_to.get()
        print(mut_from, mut_to)

        self.topology_mutation[self.selected_topology.get()][res_nr] = [mut_from, mut_to]

        del self.topology_res_fep[self.selected_topology.get()][res_nr]

        #Get FEP files for default mutation
        self.get_fep_files(mut_from, mut_to, res_nr)

        self.mutate_to.set(mut_to)
        self.mutate_from.set(mut_from)

        self.refresh_residue_list()

    def select_missing_atom(self, atomtype):
        """
        Select a atom that will be added to
        :return:
        """
        #Get topology residues
        residues = qf.create_pdb_from_topology(self.topology_paths[self.selected_topology.get()],
                                               self.app.q_settings['library'])

        self.select_res = SelectReturn(self, self.root, elements=residues,
                                       select_title='Select atom %s for FEP' % atomtype, Entry=self.reslist).show()
        self.select_res.configure(bg=self.main_color)
        self.select_res.title('Select atom %s for FEP' % atomtype)
        self.select_res.resizable()

        return self.select_res


    def view_fep(self):
        """

        :return:
        """
        selection = self.feplist.curselection()

        if len(selection) == 0:
            return

        fepfile = self.feplist.get(selection[0])

        nr = int(filter(str.isdigit, fepfile))

        fep = self.topology_fep[self.selected_topology.get()][nr]

        qf.write_fepdict(fep, path=self.app.workdir, printfep=True)

    def duplicate_fepfile(self):
        """

        :return:
        """
        pass

    def add_fepfile(self):
        """

        :return:
        """
        pass

    def delete_fepfile(self):
        """

        :return:
        """
        pass

    def move_fep(self, i=0):
        """
        Moves FEP file in listbox from position j to j + i
        :param value:
        :return:
        """
        pass

    def configure_md(self):
        """
        Opens the MD settings module!
        :return:
        """
        pass

    def write_inputfiles(self):
        """

        :return:
        """
        pass

    def run_fep(self):
        """

        :return:
        """
        pass

    def reslist_event(self, *args):
        """
        :param args:
        :return:
        """
        #Get selected residue and update option menu and FEP files!
        selection = self.reslist.curselection()
        if len(selection) == 0:
            return

        selected = self.reslist.get(selection[0])

        res_wt = selected.split()[0]
        res_nr = int(selected.split()[1])
        res_mut = selected.split()[3]

        print res_wt

    def dialog_window(self):

        self.title('resFEP setup')

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Topology frame
        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        #Frame with residues to mutate and FEP files
        frame2 = Frame(mainframe, bg=self.main_color)
        frame2.grid(row=1, column=0)

        #Frame with configure MD /write/save/close etc.
        frame3 = Frame(mainframe, bg=self.main_color)
        frame3.grid(row=2, column=0)

        #Dropdown menu for topology selection
        self.select_topology = OptionMenu(frame1, self.selected_topology, *self.topology_paths.keys())
        self.select_topology.config(bg=self.main_color, highlightbackground=self.main_color, width=15)
        self.select_topology.grid(row=1, column=0)

        #Load topology button
        load_topology = Button(frame1, text ='Load', highlightbackground=self.main_color, command=self.load_topology)
        load_topology.grid(row=1, column=1)

        #Add a topology button
        add_topology = Button(frame1, text ='Add', highlightbackground=self.main_color, command=self.add_topology)
        add_topology.grid(row=1, column=2)

        #Delete a topology
        del_topology = Button(frame1, text ='Del', highlightbackground=self.main_color, command=self.del_topology)
        del_topology.grid(row=1, column=3)



        #topology name
        self.topology_label = Label(frame1, text='*.top', bg=self.main_color)
        self.topology_label.grid(row=0, column=0, columnspan=3)

        #Residues label
        res_label = Label(frame2, text='Mutate residue(s):', bg=self.main_color)
        res_label.grid(row=0, column=0, columnspan=3, pady=(10,0), sticky='w')

        #Listbox with selected residues to mutate:
        reslist_scroll = Scrollbar(frame2)
        reslist_scroll.grid(row=1, rowspan=3, column=1, sticky='nsw')
        self.reslist = Listbox(frame2, yscrollcommand=reslist_scroll.set, width=17, height=4,
                                 highlightthickness=0, relief=GROOVE, selectmode=SINGLE, exportselection=False)
        reslist_scroll.config(command=self.reslist.yview)
        self.reslist.grid(row=1, rowspan=3, column=0, sticky='nse')
        self.reslist.config(font=tkFont.Font(family="Courier", size=12))
        self.reslist.bind('<<ListboxSelect>>', self.reslist_event)

        #Add residue button
        add_res = Button(frame2, text='Add', highlightbackground=self.main_color, command=self.add_residue)
        add_res.grid(row=1, column=2)

        #Delete residue button
        del_res = Button(frame2, text='Delete', highlightbackground=self.main_color, command=self.del_residue)
        del_res.grid(row=1, column=4)


        #Mutate from dropdown menu:
        self.mutate_from_menu = OptionMenu(frame2, self.mutate_from, 'start')
        self.mutate_from_menu.config(bg=self.main_color, highlightbackground=self.main_color, width=7)
        self.mutate_from_menu.grid(row=2, column=2)

        #right arrow:
        arrow_label = Label(frame2, text='%1s' % u'\u2192', bg=self.main_color)
        arrow_label.grid(row=2, column=3)

        #Mutate to dropdown menu
        self.mutate_to_menu = OptionMenu(frame2, self.mutate_to, 'end')
        self.mutate_to_menu.config(bg=self.main_color, highlightbackground=self.main_color, width=7)
        self.mutate_to_menu.grid(row=2, column=4)

        #Update button
        update_res_mutation = Button(frame2, text='Update', highlightbackground=self.main_color,
                                     command=self.update_residue_mutation)
        update_res_mutation.grid(row=3, column=2, columnspan=3)

        #FEP files label
        feps_label = Label(frame2, text='FEP file(s):', bg=self.main_color)
        feps_label.grid(row=4, column=0, columnspan=3, pady=(10,0), sticky='w')

        #Listbox for FEP files
        feplist_scroll = Scrollbar(frame2)
        feplist_scroll.grid(row=5, rowspan=5, column=1, sticky='nsw')
        self.feplist = Listbox(frame2, yscrollcommand=feplist_scroll.set, width=17, height=8,
                                 highlightthickness=0, relief=GROOVE, selectmode=SINGLE, exportselection=False)
        feplist_scroll.config(command=self.feplist.yview)
        self.feplist.grid(row=5, rowspan=5, column=0, sticky='nse')
        self.feplist.config(font=tkFont.Font(family="Courier", size=12))
        #self.feplist.bind('<<ListboxSelect>>', self.reslist_event)

        #View FEP file button
        view_fep = Button(frame2, text='View', highlightbackground=self.main_color, command=self.view_fep)
        view_fep.grid(row=5, column=2)

        #Duplicate selected FEP file:
        duplicate_fepfile = Button(frame2, text='Duplicate', highlightbackground=self.main_color,
                                   command=self.duplicate_fepfile)
        duplicate_fepfile.grid(row=5, column=3, columnspan=2)

        #Add a new/blank FEP file
        add_fepfile = Button(frame2, text='Add', highlightbackground=self.main_color, command=self.add_fepfile)
        add_fepfile.grid(row=6, column=2)

        #Delete selected FEP file
        del_fepfile = Button(frame2, text='Delete', highlightbackground=self.main_color, command=self.delete_fepfile)
        del_fepfile.grid(row=6, column=3, columnspan=2)

        #Move FEP file up
        move_up = Button(frame2, text=u'\u2191', highlightbackground=self.main_color, command=lambda: self.move_fep(-1))
        move_up.grid(row=7, column=2, sticky='e')

        #Move FEP file down
        move_down = Button(frame2, text=u'\u2193', highlightbackground=self.main_color,
                           command=lambda: self.move_fep(1))
        move_down.grid(row=7, column=3, columnspan=2, sticky='w')

        #Configure MD
        config_md = Button(frame3, text='Configure MD', highlightbackground=self.main_color, command=self.configure_md)
        config_md.grid(row=0, column=0, columnspan=3)

        #Write inputfiles
        write_inp = Button(frame3, text='Write', highlightbackground=self.main_color, command=self.write_inputfiles)
        write_inp.grid(row=1, column=0)

        #Run FEP calculation
        run_fep = Button(frame3, text='Run', highlightbackground=self.main_color, command=self.run_fep)
        run_fep.grid(row=1, column=1)

        #Close resFEP gui
        close_resfep = Button(frame3, text='close', highlightbackground=self.main_color, command=self.destroy)
        close_resfep.grid(row=1, column=2)

