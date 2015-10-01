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
    OptionMenu, StringVar, Entry
from subprocess import Popen, PIPE
import tkFont
import time
import os
import signal
import sys
from tkFileDialog import askopenfilenames


class ViewPyMol(Toplevel):
    def __init__(self, app, root, start_cmd=None):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.app.log('info', 'Launching pyMol with %s' % self.app.pdb_id.split('/')[-1])
        self.pdbfile = self.app.pdb_id
        self.start_cmd = start_cmd

        self.pymol_color = StringVar()
        self.pymol_atoms_color = StringVar()
        self.select_mode = StringVar()
        self.pymol_presets = StringVar()
        self.create_bond = StringVar()
        self.bond_atom_name = StringVar()

        self.pymol_color.set('Gray')
        self.pymol_atoms_color.set('C')
        self.select_mode.set('Residues')
        self.pymol_presets.set('Presets')
        self.create_bond.set('H')

        self.select_mode.trace('w', self.select_mode_changed)

        self.dialog_window()

        #QPyMol default settings:
        self.pymol_settings = ['space cmyk', 'set valence, 0.1', 'set sphere_scale, 0.8',
                               'set sphere_quality, 2', 'set stick_radius, 0.17', 'set stick_quality, 10',
                               'set defer_builds_mode, 3', 'set surface_quality, 1',
                               'set spec_power, 250', 'set spec_reflect, 2',
                               'set cartoon_fancy_helices, 1']

        #QpyMol can be launched with additional start commands:
        #self.start_cmd = ['show cartoon', 'hide lines']

        self.add_atoms_to_list()
        self.start_pymol()

        root.protocol('WM_DELETE_WINDOW', self.close_pymol)

    def select_mode_changed(self, *args):
        """
        Change the pymol selectmode
        """
        pymol_value = 1
        if self.select_mode.get() == 'Atoms':
            pymol_value = 0

        self.session.stdin.write('set mouse_selection_mode, %d\n' % pymol_value)

    def start_pymol(self):
        self.app.pymol_running = True
        self.reopen_button.config(state=DISABLED)
        for entry in self.widgetList:
            entry.config(state=NORMAL)

        tmpfile = open(self.app.workdir+'/.tmpfile','wb')
        if 'darwin' in sys.platform:
            try:
                self.session = Popen(["pymol", "-p -x -i", "%s" % self.pdbfile], stdout=tmpfile, stdin=PIPE,
                                     preexec_fn=os.setsid)
            except:
                try:
                    self.session = Popen(["MacPyMol", "-p -x -i", "%s" % self.pdbfile], stdout=tmpfile, stdin=PIPE,
                                     preexec_fn=os.setsid)
                except:
                    print 'No pymol version found'

        else:
            try:
                self.session = Popen(["pymol", "-p", "%s" % self.pdbfile], stdout=tmpfile, stdin=PIPE, preexec_fn=os.setsid)
            except:
                print 'No pymol version found'

        self.update()
        len_log = 35

        self.session.stdin.flush()

        #Load default Q-PyMol settings:
        for settings in self.pymol_settings:
            self.session.stdin.write('%s\n' % settings)



        #Check if other special start-up settings are defined:
        if self.start_cmd:
            for command in self.start_cmd:
                self.session.stdin.write('%s\n' % command)
            self.start_cmd = None

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
                            res = line.split('/')[-2]
                            name = line.split('/')[-1].split('\n')[0]
                            if '`' in res:
                                res = ' '.join(res.split('`'))
                            if '`' in name:
                                name = ''.join(name.split('`'))
                            selected = '%s %s' % (res, name)
                            print selected
                            self.set_selection(selected)

            time.sleep(0.2)

        self.reopen_button.config(state=NORMAL)

        for entry in self.widgetList:
            entry.config(state=DISABLED)

        self.app.pymol_running = False

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

    def close_pymol(self):
        try:
            os.killpg(self.session.pid, signal.SIGTERM)
            self.reopen_button.config(state=NORMAL)
            self.app.pymol_running = False
        except:
            self.app.pymol_running = False
            self.destroy()

    def set_selection(self, atoms='all'):
        """
        Highlights selection from pymol in atomlist.
        """
        if atoms == 'all':
            self.listbox.selection_set(0, END)
        else:
            try:
                resnr = int(atoms.split()[1])
                atomname = atoms.split()[2]
            except:
                print 'Failed to parse selection from PyMol'
                return

            set_selection = []

            #Get list indices to highlight
            tmp_list = self.listbox.get(0, END)
            for index_ in range(len(tmp_list)):
                if resnr == int(tmp_list[index_].split()[-1]):
                    if self.select_mode.get() == 'Atoms':
                        if tmp_list[index_].split()[1].strip() == atomname:
                            set_selection.append(index_)
                            break
                    elif self.select_mode.get() == 'Residues':
                        set_selection.append(index_)

            self.listbox.selection_set(set_selection[0], set_selection[-1])
            self.listbox.yview(set_selection[0])

    def clear_list(self):
        """
        Clears all selection from list
        """
        self.listbox.select_clear(0, END)
        self.pymol_command('select None', '')

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
            self.session.stdin.write('load_traj %s, %s\n' % (dcd, self.pdbfile.split('/')[-1].split('.')[0]))

        self.session.stdin.write('set internal_gui, 1\n')


    def dialog_window(self):
        self.title('Qgui PyMol')
        self.config(background=self.main_color)

        self.widgetList = []

        # Define frames
        frame = Frame(self, bg = self.main_color)
        frame.grid(row=0, column=0, columnspan=2, pady=10, padx=10)

        # Frame with atom list
        frame2 = Frame(self, bg=self.main_color)
        frame2.grid(row=1,column=0, pady=10, padx=(10,0))

        # Frame with pymol commands
        frame3 = Frame(self, bg=self.main_color)
        frame3.grid(row=1, column=1, pady=10, padx=(0,10))

        label_pymol = Label(frame, text='PyMol session (%s)' % self.app.pdb_id.split('/')[-1],
                            bg=self.main_color)
        label_pymol.grid(row=0, column=0, columnspan=2)

        self.close_button = Button(frame, text='Quit', highlightbackground=self.main_color, command=self.close_pymol)
        self.close_button.grid(row=1,column=0)

        self.reopen_button = Button(frame, text='View again', highlightbackground=self.main_color, command=self.start_pymol)
        self.reopen_button.grid(row=1, column=1)

        self.load_trj = Button(frame, text ='Load trajectory', highlightbackground=self.main_color, command=self.load_dcd)
        self.load_trj.grid(row=2, column=0, columnspan=2)
        self.widgetList.append(self.load_trj)

        # fram2 (atom list)
        listbox_scroll = Scrollbar(frame2)
        listbox_scroll.grid(row = 0, rowspan = 10, column = 11, sticky = 'nsw')
        self.listbox = Listbox(frame2, yscrollcommand = listbox_scroll.set, width = 23, height=30,
                               highlightthickness=0, relief=GROOVE, selectmode=EXTENDED)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 0, rowspan = 10, column = 0, columnspan = 10, sticky = 'w')
        self.listbox.config(font=tkFont.Font(family="Courier", size=12))

        sel_all_button = Button(frame2, text='Select all', highlightbackground=self.main_color,
                                command=lambda: self.set_selection('all'))
        sel_all_button.grid(row=11, column=0, columnspan=5)

        clear_button = Button(frame2, text='Clear', highlightbackground=self.main_color, command=self.clear_list)
        clear_button.grid(row=11, column=6, columnspan=5)

        invert_button = Button(frame2, text='Invert selection', highlightbackground=self.main_color,
                               command=self.invert_selection)
        invert_button.grid(row=12, column=0, columnspan=11)

        # frame3 (pymol commands)
        selectmode = Label(frame3, text='Select', bg=self.main_color)
        selectmode.grid(row=0, column=0, sticky='e')

        self.select_box = OptionMenu(frame3, self.select_mode, 'Atoms', 'Residues')
        self.select_box.config(bg=self.main_color, highlightbackground=self.main_color, width=15)
        self.select_box.grid(row=0, column=1, columnspan=2)
        self.widgetList.append(self.select_box)

        self.orient_button = Button(frame3, text='Orient', highlightbackground=self.main_color, command=lambda: self.pymol_command('orient'))
        self.orient_button.grid(row=1, column=1, sticky='e')
        self.widgetList.append(self.orient_button)

        self.zoom_button = Button(frame3, text='Zoom', highlightbackground=self.main_color,
                                  command=lambda: self.pymol_command('zoom'))
        self.zoom_button.grid(row=1, column = 2, sticky='e')
        self.widgetList.append(self.zoom_button)

        self.colorbutton = Button(frame3, text='Color', highlightbackground=self.main_color,
                                  command=lambda: self.pymol_color_atoms_selected(
                                      self.pymol_color.get(), self.pymol_atoms_color.get()))
        self.colorbutton.grid(row=2, column=0, sticky='e')
        self.widgetList.append(self.colorbutton)

        self.color_atoms = OptionMenu(frame3, self.pymol_atoms_color,
                                   'C', 'N', 'O', 'H', 'S', 'P', 'All', 'BG')
        self.color_atoms.config(highlightbackground=self.main_color, bg=self.main_color, width=5)
        self.color_atoms.grid(row=2, column=2, sticky='w')
        self.widgetList.append(self.color_atoms)

        self.colorbox = OptionMenu(frame3, self.pymol_color,
                                   'Brown', 'Black', 'Blue', 'Cyan', 'Orange', 'Gray', 'Green', 'Pink', 'Purple',
                                   'Red', 'Sand', 'Silver', 'White', 'Yellow')
        self.colorbox.config(highlightbackground=self.main_color, bg=self.main_color, width=8)
        self.colorbox.grid(row=2, column=1, sticky='e')
        self.widgetList.append(self.colorbox)

        #Ball and stick:
        show_bs = Label(frame3, text='B & S', bg=self.main_color)
        show_bs.grid(row=3, column=0, sticky='w')
        self.show_bs = Button(frame3, text= 'On', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_special('show b&s'))
        self.show_bs.grid(row=3, column = 1, sticky='e')
        self.widgetList.append(self.show_bs)

        self.hide_bs = Button(frame3, text= 'Off', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_special('hide b&s'))
        self.hide_bs.grid(row=3, column = 2, sticky='w')
        self.widgetList.append(self.hide_bs)


        showcartoon = Label(frame3, text='Cartoon ', bg=self.main_color)
        showcartoon.grid(row=4, column=0, sticky='w')
        self.show_cartoon = Button(frame3, text= 'On', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_command('show cartoon,'))
        self.show_cartoon.grid(row=4, column = 1, sticky='e')
        self.widgetList.append(self.show_cartoon)

        self.hide_cartoon = Button(frame3, text= 'Off', highlightbackground=self.main_color,
                                   command=lambda: self.pymol_command('hide cartoon,'))
        self.hide_cartoon.grid(row=4, column = 2, sticky='w')
        self.widgetList.append(self.hide_cartoon)

        showlines = Label(frame3, text='Lines', bg=self.main_color)
        showlines.grid(row=5, column=0, sticky='w')

        self.show_lines = Button(frame3, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show lines,'))
        self.show_lines.grid(row=5, column=1, sticky='e')

        self.widgetList.append(self.show_lines)
        self.hide_lines = Button(frame3, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide lines,'))
        self.hide_lines.grid(row=5, column=2, sticky='w')
        self.widgetList.append(self.hide_lines)

        showsphere = Label(frame3, text='Spheres', bg=self.main_color)
        showsphere.grid(row=6, column=0, sticky='w')

        self.show_sphere = Button(frame3, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show spheres,'))
        self.show_sphere.grid(row=6, column=1, sticky='e')
        self.widgetList.append(self.show_sphere)

        self.hide_sphere = Button(frame3, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide spheres,'))
        self.hide_sphere.grid(row=6, column=2, sticky='w')
        self.widgetList.append(self.hide_sphere)

        showsticks = Label(frame3, text='Sticks', bg=self.main_color)
        showsticks.grid(row=7, column=0, sticky='w')

        self.show_sticks = Button(frame3, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show sticks,'))
        self.show_sticks.grid(row=7, column=1, sticky='e')
        self.widgetList.append(self.show_sticks)

        self.hide_sticks = Button(frame3, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide sticks,'))
        self.hide_sticks.grid(row=7, column=2, sticky='w')
        self.widgetList.append(self.hide_sticks)

        showlines = Label(frame3, text='Surface', bg=self.main_color)
        showlines.grid(row=8, column=0, sticky='w')

        self.show_surface = Button(frame3, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('show surface,'))
        self.show_surface.grid(row=8, column=1, sticky='e')
        self.widgetList.append(self.show_surface)

        self.hide_surface = Button(frame3, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_command('hide surface,'))
        self.hide_surface.grid(row=8, column=2, sticky='w')
        self.widgetList.append(self.hide_surface)

        self.special_presets = OptionMenu(frame3, self.pymol_presets,
                                   'QuteMol', 'Charged', 'PlastRay', 'Q-pymol')
        self.special_presets.config(highlightbackground=self.main_color, bg=self.main_color, width=10)
        self.special_presets.grid(row=9, column=0, sticky='w')
        self.widgetList.append(self.special_presets)

        self.show_preset = Button(frame3, text='On', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_special('show %s' % self.pymol_presets.get()))
        self.show_preset.grid(row=9, column=1, sticky='e')
        self.widgetList.append(self.show_preset)

        self.hide_preset = Button(frame3, text='Off', highlightbackground=self.main_color,
                                 command=lambda: self.pymol_special('hide %s' % self.pymol_presets.get()))
        self.hide_preset.grid(row=9, column=2, sticky='w')
        self.widgetList.append(self.hide_preset)

        #self.add_h = Button(frame3, text='Add hydrogens', highlightbackground=self.main_color,
        #                    command=lambda: self.pymol_command('h_add '))
        #self.add_h.grid(row=10, column=0, columnspan=3, pady=(20,0))
        #self.widgetList.append(self.add_h)

        #bond_order = Label(frame3, text='Bond order', bg=self.main_color)
        #bond_order.grid(row=11, column=0)

        #self.bond_order_0 = Button(frame3, text='0', bg=self.main_color, highlightbackground=self.main_color,
        #                           command=lambda: self.change_bond_order(0))
        #self.bond_order_0.grid(row=11, column=1, sticky = 'w')
        #self.widgetList.append(self.bond_order_0)

        #self.bond_order_1 = Button(frame3, text='1', bg=self.main_color, highlightbackground=self.main_color,
        #                           command=lambda: self.change_bond_order(1))
        #self.bond_order_1.grid(row=11, column=1, sticky = 'e')
        #self.widgetList.append(self.bond_order_1)

        #self.bond_order_2 = Button(frame3, text='2', bg=self.main_color, highlightbackground=self.main_color,
        #                           command=lambda: self.change_bond_order(2))
        #self.bond_order_2.grid(row=11, column=2, sticky = 'w')
        #self.widgetList.append(self.bond_order_2)

        #self.fix_geom = Button(frame3, text='Fix geometry',  bg=self.main_color, highlightbackground=self.main_color,
        #                       command=self.fix_geometry)
        #self.fix_geom.grid(row=12, column=0,  columnspan=3)
        #self.widgetList.append(self.fix_geom)


