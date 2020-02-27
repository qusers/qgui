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

from tkinter import  Label, Button, Listbox, Scrollbar, EXTENDED, Frame, \
    Toplevel, END, GROOVE, StringVar, OptionMenu, IntVar, Checkbutton

import tkinter.font
from subprocess import call
import numpy as np
import os
from Qplot import Qplot

from tkinter.filedialog import asksaveasfilename, askopenfilenames


class AnalyzeQcalc(Toplevel):
    def __init__(self, app, root, qcalc_what):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root
        #rmsf, rmsd, bond, angle, torsion
        self.qcalc_what = qcalc_what

        self.main_color = self.app.main_color
        self.pdbfile = self.app.pdb_id

        self.residual_check = IntVar()
        self.residual_check.set(0)
        self.residual_input = False

        self.dcd = list()
        self.sorted_dcd = False

        self.title_dict = {'RMSD': 'Selected atoms',
                           'RMSF': 'Selected atoms',
                           'Distance': 'Selected atom pairs',
                           'Angle': 'Selected angle atoms (3)',
                           'Torsion': 'Selected torsion atoms (4)'}

        self.qcalc_nr = {'RMSD': 1,
                         'Fit': 2,
                         'Distance': 3,
                         'Angle': 4,
                         'Torsion': 5,
                         'Entropy': 6,
                         'RMSF': 13}

        self.selected_frame = StringVar()
        self.selected_frame.set('Setup')
        self.selected_frame.trace('w', self.show_var_frame)

        self.dialog_window()

        #Fill atomlist
        self.add_atoms_to_list()

    def add_atoms_to_list(self):
        """
        Add atoms from pdbfile to selection list.
        """
        with open(self.pdbfile) as pdb:
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

    def sel_solute(self):
        """
        Selects solute atoms (all except waters)
        """

        solvent = ['HOH','SPC','WAT','T3P','TIP3', 'H2O']

        self.listbox.select_clear(0, END)

        residues = self.listbox.get(0, END)

        for i in range(len(residues)):
            res = residues[i].split()[2]
            if res not in solvent:
                self.listbox.selection_set(i)

    def sel_sc(self):
        """
        Selects protein sidechains
        """
        bb = ['N', 'CA', 'O', 'C', 'H', 'HA']
        solvent = ['HOH','SPC','WAT','T3P','TIP3', 'H2O']

        self.listbox.select_clear(0, END)

        residues = self.listbox.get(0, END)

        current_res = 0

        bb_atoms = list()
        sc_atoms = list()

        for i in range(len(residues)):
            if current_res != int(residues[i].split()[3]):
                if len(bb_atoms) >= 4:
                    for j in sc_atoms:
                        self.listbox.selection_set(j)
                del bb_atoms[:]
                del sc_atoms[:]
                current_res = int(residues[i].split()[3])
            res = residues[i].split()[2]
            if res not in solvent:
                if residues[i].split()[1] in bb:
                    bb_atoms.append(i)
                else:
                    sc_atoms.append(i)

    def sel_bb(self):
        """
        Selects protein backbone
        """
        bb = ['N', 'CA', 'O', 'C', 'H', 'HA']

        self.listbox.select_clear(0, END)

        residues = self.listbox.get(0, END)

        current_res = 0

        #List containing list indexes for residue backbone atoms
        res_atoms = list()

        for i in range(len(residues)):
            if current_res != int(residues[i].split()[3]):
                if len(res_atoms) >= 4:
                    for j in res_atoms:
                        self.listbox.selection_set(j)
                del res_atoms[:]
                current_res = int(residues[i].split()[3])

            if residues[i].split()[1] in bb:
                res_atoms.append(i)

    def sel_na(self):
        """
        selects entire nucleic acids in listbox
        """

        resnames = ['UR', 'AD', 'GU', 'CY']

        self.listbox.select_clear(0, END)

        residues = self.listbox.get(0, END)

        for i in range(len(residues)):
            if residues[i][13:15] in resnames:
                self.listbox.selection_set(i)

    def clear_sel(self):
        """
        Clears atom selection
        """
        self.listbox.select_clear(0, END)

    def add_sel_atoms(self):
        """
        Adds selected atoms from atomlist to Qcalc mask
        """

        selection = list(map(int, self.listbox.curselection()))

        if len(selection) < 1:
            return
        atom_numbers = list()

        for i in selection:
            atom_numbers.append(self.listbox.get(i).split()[0])

        print(self.qcalc_what)

        if self.qcalc_what == 'RMSD' or self.qcalc_what == 'RMSF':
            self.add_all_atoms(atom_numbers)
        elif self.qcalc_what == 'Distance':
            self.add_atoms(atom_numbers, 2)
        elif self.qcalc_what == 'Angle':
            self.add_atoms(atom_numbers, 3)
        elif self.qcalc_what == 'Torsion':
            self.add_atoms(atom_numbers, 4)

    def add_all_atoms(self, atom_numbers=list()):
        """
        Takes all selected atoms and adds atomnumbers one by one to Qcalc mask
        """

        for atom in atom_numbers:
            self.selected_listbox.insert(END, '%6s' % atom)

    def add_atoms(self, atom_numbers, b=2):
        if len(atom_numbers) != b:
            self.app.errorBox('Warning','Select exactly %d atoms to analyze %s' % (b, self.qcalc_what))
            return
        atom_numbers = ['%6s' % x for x in atom_numbers]

        atoms = ' '.join(atom_numbers)
        self.selected_listbox.insert(END, atoms)

    def del_sel_atoms(self):
        """
        Deletes selected atoms from Qcalc mask
        """
        selection = list(map(int, self.selected_listbox.curselection()))

        if len(selection) < 1:
            return

        continous = True
        j = selection[0]

        #Decide if selection is continous (takes a long time to loop over indeces when selection is big)
        #If selection is large and continous, delete(start, end) works 100000x faster!!
        for i in range(1, len(selection) - 1):
            if j + 1 != i:
                continous = False
                break
            else:
                j = i

        if not continous:
            selection.reverse()
            for i in selection:
                self.selected_listbox.delete(i)
        else:
            self.selected_listbox.delete(selection[0], selection[-1])

    def add_trj(self):
        """
        Select trajectory files for Qcalc
        """

        trjs = askopenfilenames(parent = self, initialdir = self.app.workdir,
                                  filetypes=(("Trajectory", "*.dcd"),("All files","*.*")))

        if trjs:
            self.sorted_dcd = False
            for dcd in trjs:
                self.dcd.append(dcd)
            self.trj_listbox.delete(0, END)

        for dcd in self.dcd:
            self.trj_listbox.insert(END, '/'.join(dcd.split('/')[-2:]))

    def del_trj(self):
        """
        Delete trajectories from listbox
        """
        self.sorted_dcd = False
        selection = list(map(int, self.trj_listbox.curselection()))

        if len(selection) < 1:
            return

        for i in selection:
            dcd = self.trj_listbox.get(i)
            for j in range(len(self.dcd)):
                if '/'.join(self.dcd[j].split('/')[-2:]) == dcd:
                    del self.dcd[j]
                    break

        self.trj_listbox.delete(0, END)
        for dcd in self.dcd:
            self.trj_listbox.insert(END, '/'.join(dcd.split('/')[-2:]))

    def sort_trj(self):
        """
        Sort trajectory files ...
        """
        if not self.sorted_dcd:
            self.dcd.sort()
            self.sorted_dcd = True
        else:
            self.dcd.reverse()

        self.trj_listbox.delete(0, END)
        for dcd in self.dcd:
            self.trj_listbox.insert(END, '/'.join(dcd.split('/')[-2:]))

    def toggle_residual(self):
        """
        Toggle if inputfiles are per residue or in total!
        """
        if self.residual_check.get() == 1:
            self.residual_input = True
            self.app.log(' ', '\nThis will generate one input-file per residue\n')
        else:
            self.residual_input = False
            self.app.log(' ', '\nThis will generate one input-file for all residues\n')

    def import_results(self, filenames=None, run=False):
        """
        reads outputs from Qcalc and inserts into listbox. If run= True, filename is bash script with all
        <input>output files. Else, normal output file is expected
        """
        qcalc_logs = list()
        if not run:
            filenames = askopenfilenames(parent = self, initialdir = self.app.workdir,
                                  filetypes=(("Qcalc", "*.qcalc"),("All files","*.*")))
            if filenames:
                for filename in filenames:
                    qcalc_logs.append(filename)
            else:
                return
        else:
            with open(filenames, 'r') as run_script:
                for line in run_script:
                    if 'Qcalc' in line:
                        qcalc_logs.append(line.split('>')[-1].strip())

        frame = 0
        for output in qcalc_logs:
            insert_results = False
            if os.path.isfile(output):
                with open(output, 'r') as qcalc:
                    for line in qcalc:
                        if insert_results:
                            if 'frame' in line:
                                print(line.strip()[21:])
                                self.results_listbox.insert(END, '%s' % line.strip()[21:])
                            else:
                                frame += 1
                                self.results_listbox.insert(END,'%6d %s' % (frame, line[26:]))

                        if 'Calculation results' in line:
                            insert_results = True


    def run_qcalc(self):
        """
        This will run Qcalc directly (no submission script is used)
        """
        qcalc_input = self.write_qcalc(submit=False)

        tmpfile = open('.tmpfile', 'w')

        call('bash %s' % qcalc_input, shell=True, stdout=tmpfile, stderr=tmpfile)

        self.selected_frame.set('Results')
        self.results_listbox.delete(0, END)
        self.import_results(filenames=qcalc_input, run=True)
        self.update()

    def submit_qcalc(self):
        """
        Submits the job using submission script
        """
        submit_file = self.write_qcalc(submit=True)

        tmpfile = open('.tmpfile', 'w')
        call('bash %s' % submit_file, shell=True, stdout=tmpfile, stderr=tmpfile)

        job_id = open(self.app.workdir + '/.tmpfile','r').readlines()
        self.app.log('info','Submitting %s jobs ...' % self.qcalc_what)
        for line in job_id:
            self.app.main_window.update_txt(line)

        self.app.errorBox('Info','Jobs submitted!')

    def write_qcalc(self, submit=False):
        """
        writes input submission script file for qcalc (.sh) and creatin Qcalc input file(s)
        """
        qcalc = self.app.q_settings[ 'executables' ][3]
        script_name = '%s/%s.sh' % (self.app.workdir, self.qcalc_what)
        sub_script = open(script_name, 'w')
        if submit:
            #Write submission script:
            if os.path.isfile(self.app.settings_path + '/qsubmit'):
                submissionscipt = open(self.app.settings_path + '/qsubmit','r').readlines()
            elif os.path.isfile(self.app.workdir + '/' + 'qsubmit'):
                submissionscipt = open(self.app.workdir + '/' + 'qsubmit','r').readlines()
            else:
                submissionscipt = ['#!/bin/bash\n#Qdyn I/O\n']
                print('submission script not found! Please edit this in settings')
            for line in submissionscipt:
                if '#Qdyn I/O' in line:
                    break
                else:
                    sub_script.write(line)

        if self.residual_input:
            inputfiles = self.write_residual_qcalc()
        else:
            inputfiles = list()
            inputfiles.append(self.write_input_qcalc())

        for inputfile in inputfiles:
            sub_script.write('%s <%s > %s' % (qcalc, inputfile, inputfile.split('.')[0]+'_out.qcalc'))

        if submit:
            #If use submission script, check for end statements (comes after #Qdyn I/O):
            write_end = False
            for k in range(len(submissionscipt)):
                if '#Qdyn I/O' in submissionscipt[k]:
                    end_statements_start = k + 1
                    write_end = True
            if write_end:
                for line in range(end_statements_start, len(submissionscipt)):
                    sub_script.write(submissionscipt[line])

        sub_script.close()

        return script_name

    def find_solute(self):
        resnrs = list()
        solvent = ['HOH','SPC','WAT','TIP3', 'TI3', 'T3P', 'H2O']

        for line in self.listbox.get(0, END):
            res = line.split()[2]
            resnr = int(line.split()[3])

            if res not in solvent:
                if resnr not in resnrs:
                    resnrs.append(resnr)

        start = min(resnrs)
        end = max(resnrs)

        return start, end

    def write_input_qcalc(self, name=None, atoms=None):
        """
        Writes single Qcalc input file
        """
        inputfile = '%s/%s.inp' % (self.app.workdir, self.qcalc_what)
        if name:
            inputfile = name

        input_file = open(inputfile, 'w')

        #Write topology name to input file
        input_file.write('%s\n' % self.app.top_id)

        if not atoms:
            atoms = self.selected_listbox.get(0, END)

        if self.qcalc_what == 'RMSF' or self.qcalc_what == 'RMSD':
            #Write fit proceducre
            start_res, end_res = self.find_solute()
            input_file.write('2\n')
            input_file.write('residue %d %d\n.\n' % (start_res, end_res))

        if len(atoms[0].split()) > 1:
            for i in range(len(atoms)):
                input_file.write('%d\n' % self.qcalc_nr[self.qcalc_what])
                input_file.write('%s\n' % atoms[i])
                input_file.write('.\n')

        else:
            input_file.write('%d\n' % self.qcalc_nr[self.qcalc_what])
            for i in range(len(atoms)):
                input_file.write('%s\n' % atoms[i])
            input_file.write('.\n')

        input_file.write('go\n')
        for dcd in self.dcd:
            input_file.write('%s\n' % dcd)

        input_file.write('.\n')
        input_file.close()

        return inputfile

    def write_residual_qcalc(self):
        """
        Writes n Qcalc inputfile for n residues in selection
        """
        #Go throug atom numbers and create res_nr : [atom_number] dict
        res_atoms = dict()
        pdb = self.listbox.get(0, END)

        for atomnr in self.selected_listbox.get(0, END):
            resnr = None
            for line in pdb:
                if line.split()[0] == atomnr:
                    resnr = line.split()[3]
                    break
            if not resnr:
                self.app.errorBox('Error', 'Could not find atom %s' % atomnr)
                break
            if not resnr in list(res_atoms.keys()):
                res_atoms[resnr] = [atomnr]
            else:
                res_atoms[resnr].append(atomnr)

        inputfiles = list()

        for res in sorted(res_atoms.keys()):
            inputfile = '%s/%s_%04d.inp' % (self.app.workdir, self.qcalc_what, int(res))
            atoms = res_atoms[res]
            inputfiles.append(self.write_input_qcalc(name=inputfile, atoms=atoms))

        return inputfiles

    def plot_it(self):
        y_label = {'RMSD': r'RMSD ($\AA$)',
                         'Distance': r'Distance ($\AA$)',
                         'Angle': 'Angle',
                         'Torsion': 'Torsion',
                         'Entropy': 6,
                         'RMSF': r'RMSF ($\AA$)'}

        x_label = 'Frame'

        y_list = [[]]
        x_list = [[]]

        y_val = list()
        for line in self.results_listbox.get(1, END):
            try:
                y_val.append(float(line.split()[1]))
            except:
                continue

        x_val = np.arange(len(y_val))

        y_list[0].append(y_val)
        x_list[0].append(x_val)

        self.plot_ = Qplot(self, self.root, y_list, x_list, [self.qcalc_what], [[' ']], x_label,
                           y_label[self.qcalc_what])
        self.plot_.resizable()
        self.plot_.config(background=self.main_color)

    def export_table(self, listbox=None):
        """
        Exports data in listbox to file:
        """
        filename = asksaveasfilename(parent=self, title='Export table', initialdir=self.app.workdir,
                                      filetypes=(("Text", "*.txt"), ("All files","*.*")),
                                      initialfile = 'Qdata.txt')

        if not filename:
            return

        if not listbox:
            listbox = self.results_listbox


        data_to_export = listbox.get(1, END)
        output = open(filename, 'w')

        for line in data_to_export:
            output.write(line)
            if not '\n' in line:
                output.write('\n')

        output.close()

        self.app.log('info', '../%s written' % '/'.join(filename.split()[-3:]))

    def show_var_frame(self, *args):
        frames = {'Setup': self.setup_frame,
                  'Results': self.result_frame}

        for i in list(frames.keys()):
            frames[i].grid_forget()
        try:
            frames[self.selected_frame.get()].grid(row=3, column=0, columnspan=2)
        except:
            pass

    def dialog_window(self):
        self.title('Analyze %s' % self.qcalc_what)
        self.mainframe = Frame(self, bg=self.main_color)
        self.mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame vith menu
        topframe = Frame(self.mainframe, bg=self.main_color)
        topframe.grid(row=0, column=0, columnspan=2, pady=(10,10))

        #VAR frames
        self.setup_frame = Frame(self.mainframe, bg=self.main_color)
        self.setup_frame.grid(row=3, column=0, columnspan=2)

        self.result_frame = Frame(self.mainframe, bg=self.main_color)

        save_run = Frame(self.setup_frame, bg=self.main_color)
        save_run.grid(row=14, column=0, columnspan=8)

        #Bottomframe
        bottomframe = Frame(self.mainframe, bg=self.main_color)
        bottomframe.grid(row=4, column=0, columnspan=2)


        #Select frame
        self.view_frame = OptionMenu(topframe, self.selected_frame,
                                   'Setup', 'Results')
        self.view_frame.config(highlightbackground=self.main_color, bg=self.main_color, width=30)
        self.view_frame.grid(row=0, column=0)

        #SETUP FRAME
        #Row 0 reserved for a sync pymol checkbutton... TODO

        listbox_scroll = Scrollbar(self.setup_frame)
        listbox_scroll.grid(row = 1, rowspan = 10, column = 3, sticky = 'nsw')
        self.listbox = Listbox(self.setup_frame, yscrollcommand = listbox_scroll.set, width = 23, height=30,
                               highlightthickness=0, relief=GROOVE, selectmode=EXTENDED)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 1, rowspan = 10, column = 0, columnspan = 3, sticky = 'w')
        self.listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        solute_button = Button(self.setup_frame, text='Solute', highlightbackground=self.main_color,
                               command=self.sel_solute)
        solute_button.grid(row=11, column=0)

        bb_button = Button(self.setup_frame, text='BB', highlightbackground=self.main_color, command=self.sel_bb)
        bb_button.grid(row=11, column=1)

        sc_button = Button(self.setup_frame, text='SC', highlightbackground=self.main_color, command=self.sel_sc)
        sc_button.grid(row=11, column=2)

        na_button = Button(self.setup_frame, text='NA', highlightbackground=self.main_color, command=self.sel_na)
        na_button.grid(row=12, column=0, columnspan=3)

        clear_button = Button(self.setup_frame, text='Clear selection', highlightbackground=self.main_color,
                              command=self.clear_sel)
        clear_button.grid(row=13, column=0, columnspan=3)

        #Column 4 --> right frame reserved
        add_label = Label(self.setup_frame, text=self.title_dict[self.qcalc_what], bg=self.main_color)
        add_label.grid(row=0, column=5, columnspan=4)

        add_atoms = Button(self.setup_frame, text="\u21D2", highlightbackground=self.main_color,
                           command=self.add_sel_atoms)
        add_atoms.grid(row=3, column=4, sticky='s')

        add_atoms = Button(self.setup_frame, text="\u2717", highlightbackground=self.main_color,
                           command=self.del_sel_atoms)
        add_atoms.grid(row=4, column=4, sticky='n')

        selected_scroll = Scrollbar(self.setup_frame)
        selected_scroll.grid(row = 1, rowspan = 6, column = 8, sticky = 'nsw')
        self.selected_listbox = Listbox(self.setup_frame, yscrollcommand = selected_scroll.set, width = 28, height=20,
                               highlightthickness=0, relief=GROOVE, selectmode=EXTENDED)
        selected_scroll.config(command=self.selected_listbox.yview)
        self.selected_listbox.grid(row = 1, rowspan = 6, column = 5, columnspan = 3, sticky = 'w')
        self.selected_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        trj_label = Label(self.setup_frame, text='Trajectory files', bg=self.main_color)
        trj_label.grid(row=7, column=4, columnspan=4)

        trj_scroll = Scrollbar(self.setup_frame)
        trj_scroll.grid(row = 8, rowspan = 4, column = 8, sticky = 'nsw')
        self.trj_listbox = Listbox(self.setup_frame, yscrollcommand = trj_scroll.set, width = 33, height=10,
                               highlightthickness=0, relief=GROOVE, selectmode=EXTENDED)
        trj_scroll.config(command=self.trj_listbox.yview)
        self.trj_listbox.grid(row = 8, rowspan = 4, column = 4, columnspan = 4, sticky = 'e')
        self.trj_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        add_trj = Button(self.setup_frame, text='+', highlightbackground=self.main_color, command=self.add_trj)
        add_trj.grid(row=12, column=4, sticky='e')

        add_trj = Button(self.setup_frame, text='-', highlightbackground=self.main_color, command=self.del_trj)
        add_trj.grid(row=12, column=5, sticky='w')

        sort_trj = Button(self.setup_frame, text='Sort', highlightbackground=self.main_color, command=self.sort_trj)
        sort_trj.grid(row=12, column=6)

        res_label = Label(self.setup_frame, text='Residual:', bg=self.main_color)
        residual_check = Checkbutton(self.setup_frame, bg = self.main_color, variable=self.residual_check,
                                     command=self.toggle_residual)
        if self.qcalc_what == 'RMSF' or self.qcalc_what == 'RMSD':
            res_label.grid(row=13,column=4, columnspan=3)
            residual_check.grid(row=13, column=4, columnspan= 3, sticky='e')

        run_button = Button(save_run, text='Run', highlightbackground=self.main_color, command=self.run_qcalc)
        run_button.grid(row=0, column=0)

        run_button = Button(save_run, text='Submit', highlightbackground=self.main_color, command=self.submit_qcalc)
        run_button.grid(row=0, column=1)

        write_button = Button(save_run, text='Write', highlightbackground=self.main_color,
                              command=self.write_qcalc)
        write_button.grid(row=0, column=2)

        #RESULTS FRAME
        import_button = Button(self.result_frame, text='Import', highlightbackground=self.main_color,
                               command=lambda: self.import_results(filenames=None, run=False))
        import_button.grid(row=0, column=0)

        results_scroll = Scrollbar(self.result_frame)
        results_scroll.grid(row = 1, rowspan = 10, column = 4, sticky = 'nsw')
        self.results_listbox = Listbox(self.result_frame, yscrollcommand = results_scroll.set, width = 60, height=30,
                               highlightthickness=0, relief=GROOVE, selectmode=EXTENDED)
        results_scroll.config(command=self.results_listbox.yview)
        self.results_listbox.grid(row = 1, rowspan = 10, column = 0, columnspan = 3, sticky = 'w')
        self.results_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        export_button = Button(self.result_frame, text='Export table', highlightbackground=self.main_color,
                               command=self.export_table)
        export_button.grid(row=11, column=0)

        plot_button = Button(self.result_frame, text='Plot it!', highlightbackground=self.main_color,
                               command=self.plot_it)
        plot_button.grid(row=11, column=1)


        #Bottom frame Quit/save
        quit_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        quit_button.grid(row=0, column=0)
