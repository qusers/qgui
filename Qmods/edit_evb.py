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

from Tkinter import Entry, Spinbox, Label, Button,Frame, Toplevel, END, GROOVE, DISABLED, NORMAL, Text, Scrollbar, \
    Listbox, EXTENDED, HORIZONTAL, Radiobutton, IntVar, StringVar

import tkFont


class EditEvbNotes(Toplevel):
    """Edit the note field in EVB/FEP list """

    def __init__(self, app, root, qatom, sel_index, note):
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root
        
        self.main_color = self.app.main_color
        
        self.insert_index = sel_index
        self.note = note
        self.qatom = qatom
        self.sel_index = sel_index

        self.dialog_box()
        self.txt.insert(0, self.note)

    def save_note(self):
        """
        appends note to list
        """
        self.app.q_notes[self.qatom] = self.txt.get().strip()
        self.app.update_q_atoms()

        self.destroy()

    def dialog_box(self):
        """
        Edit EVB notes window
        """

        self.title('Edit EVB/FEP-file notes')
        self.config(background=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        self.txt = Entry(mainframe, width=30, highlightthickness=0)
        self.txt.grid(row=0, column=0, columnspan=2)

        save = Button(mainframe, text='Save', highlightbackground=self.main_color, command=self.save_note)
        save.grid(row=1,column=0, sticky='e')

        cancel = Button(mainframe, text='Cancel', highlightbackground=self.main_color, command=self.destroy)
        cancel.grid(row=1, column=1, sticky='w')


class ImportParameters(Toplevel):
    """Select parameters from list and return to EVB setup"""

    def __init__(self, app, root, insert_listbox, term='[atom_types]'):
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self._root = root
        self.listbox = insert_listbox
        self.term = term

        self.dialog_box()
        self.fill_list()

    def fill_list(self):
        found_atoms = False
        if self.term == '[atom_types]':
            break_term = '[bonds]'
        elif self.term == '[bonds]':
            break_term = '[angles]'
        elif self.term == '[angles]':
            break_term = '[torsions]'
        elif self.term == '[torsions]':
            break_term = '[impropers]'
        elif self.term == '[impropers]':
            break_term = '                   '
        else:
            break_term = '['

        for prmfile in self.app.prms:
            with open(prmfile) as prm:
                for line in prm:
                    if break_term in line:
                        break
                    if found_atoms:
                        self.prm_listbox.insert(END, line)
                    if self.term in line:
                        found_atoms = True

    def get_parameters(self):
        selections = map(int, self.prm_listbox.curselection())
        if len(selections) == 0:
            return
        for selected in selections:
            prm_line = self.prm_listbox.get(selected)
            if self.term == '[atom_types]':
                atomtype = prm_line.split()[0].strip()
                ri = float(prm_line.split()[1].strip())
                ei = float(prm_line.split()[3].strip())
                ci = 70.0
                if atomtype[0] == 'H':
                    ci = 7.0
                ai = 1.6
                ri1_4 = float(prm_line.split()[4].strip())
                ei1_4 = float(prm_line.split()[5].strip())
                mass = float(prm_line.split()[6].strip())
                in_list = map(lambda x: x.split()[0].strip(), self.listbox.get(0,END))
                if atomtype not in in_list:
                    self.listbox.insert(END, '%4s %7.2f %5.2f %5.2f %4.2f %7.2f %5.2f %5.2f' %
                                        (atomtype.ljust(4), ri, ei, ci, ai, ri1_4, ei1_4, mass))
            elif self.term == '[bonds]':
                bond = '%s %s' % (prm_line.split()[0], prm_line.split()[1])

                kb = float(prm_line.split()[2])
                rb = float(prm_line.split()[3])
                alpha = 2.00
                De = (kb / 8.00)
                in_list = map(lambda x: int(x.split()[0]), self.listbox.get(0, END))
                new_nr = max(in_list) + 1
                self.listbox.insert(END, '%2d %6s %6s %5s %7s !%10s' %
                                         (new_nr, str(De), str(alpha), str(rb), str(kb), bond.ljust(10)))
                prm_nr = len(self.app.bond_prm) + 1
                self.app.bond_prm[prm_nr] = [bond, De, alpha, rb, kb]

            elif self.term == '[angles]':
                angle = '%s %s %s' % (prm_line.split()[0], prm_line.split()[1], prm_line.split()[2])
                angle_rev = '%s %s %s' % (prm_line.split()[2], prm_line.split()[1], prm_line.split()[0])
                k = prm_line.split()[3]
                theta = prm_line.split()[4]
                self.app.angle_prm[angle] = [k, theta]
                self.app.angle_prm[angle_rev] = [k, theta]

            elif self.term == '[torsions]':
                t1, t2, t3, t4 = prm_line.split()[0:4]
                torsion = '%s %s %s %s' % (t1, t2, t3, t4)
                torsion_rev = '%s %s %s %s' % (t4, t3, t2, t1)
                kt = float(prm_line.split()[4])
                minima = int(prm_line.split()[5])
                phase = float(prm_line.split()[6])
                paths = float(prm_line.split()[7])
                if torsion not in self.app.torsion_prm.keys():
                    self.app.torsion_prm[torsion] = [[0,0.0,1],[0, 180.0,1],[0, 0.0,1]]
                if torsion_rev not in self.app.torsion_prm.keys():
                    self.app.torsion_prm[torsion_rev] = [[0,0.0,1],[0, 180.0,1],[0, 0.0,1]]

                self.app.torsion_prm[torsion][abs(minima) - 1] = [kt, phase, paths]
                self.app.torsion_prm[torsion_rev][abs(minima) - 1] = [kt, phase, paths]

            elif self.term == '[impropers]':
                t1, t2, t3, t4 = prm_line.split()[0:4]
                improper = '%s %s %s %s' % (t1, t2, t3, t4)
                improper_rev = '%s %s %s %s' % (t4, t3, t2, t1)
                kt = float(prm_line.split()[4])
                phase = float(prm_line.split()[5])

                self.app.improper_prm[improper] = [kt, phase]
                self.app.improper_prm[improper_rev] = [kt, phase]
        if self.term != '[atom_types]':
            self.app.update_all()

    def dialog_box(self):
        self.title('Select parameters to import')
        self.config(bg=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10,10), pady=(10,10))

        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        frame2 = Frame(mainframe, bg=self.main_color)
        frame2.grid(row=1, column=0)

        prm_yscroll = Scrollbar(frame1)
        prm_yscroll.grid(row = 1, rowspan=10, column = 1, sticky = 'nsw', padx=(0,10))
        prm_xscroll = Scrollbar(frame1, orient=HORIZONTAL)
        prm_xscroll.grid(row=11, column=0, sticky='we')

        self.prm_listbox = Listbox(frame1, yscrollcommand = prm_yscroll.set, xscrollcommand = prm_xscroll.set,
                                      width=80, height=30, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        prm_yscroll.config(command=self.prm_listbox.yview)
        prm_xscroll.config(command=self.prm_listbox.xview)
        self.prm_listbox.grid(row=1, rowspan=10, column = 0, sticky = 'e')
        self.prm_listbox.config(font=tkFont.Font(family="Courier", size=12))

        get_button = Button(frame2, text='Get selected', highlightbackground=self.main_color, command=self.get_parameters)
        get_button.grid(row=0, column=0)

        close_button = Button(frame2, text='Exit', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=0, column=1)


class EditParameters(Toplevel):
    """Edit or add Atom parameters in EVB setup"""

    def __init__(self, app, root, insert_listbox, insert_index='END', edit=False):
        Toplevel.__init__(self, root)
        self.app = app
        self._root = root
        self.main_color = self.app.main_color
        self.listbox = insert_listbox
        self.insert_index = insert_index
        self.edit = edit

        if self.edit:
            self.new_title = 'Edit atom type'
        else:
            self.new_title = 'Add atom type'
        self.dialog_box()
        self.fill_entries()

    def fill_entries(self):
        """
        Fill entries when in edit mode
        """
        if self.edit:
            name, ri, ei, ci, ai, r1_4, e1_4, mass = self.listbox.get(self.insert_index).split()[0:]
            self.name.insert(0, name)
            self.ri.insert(0, ri)
            self.ei.insert(0, ei)
            self.ci.insert(0, ci)
            self.ai.insert(0, ai)
            self.ri1_4.insert(0, r1_4)
            self.ei1_4.insert(0, e1_4)
            self.mass.insert(0, mass)

    def save_prm(self):
        """
        Saves parameters to name
        """
        if self.edit:
            #Check if name is changed.
            old_name = self.listbox.get(self.insert_index).split()[0]
            new_name = self.name.get()
            if old_name != new_name:
                #Name is changed, append changes as new entry:
                self.insert_index = END
            else:
                self.listbox.delete(self.insert_index)

        name = self.name.get()
        ri = float(self.ri.get())
        ei = float(self.ei.get())
        ci = float(self.ci.get())
        ai = float(self.ai.get())
        r1_4 = float(self.ri1_4.get())
        e1_4 = float(self.ei1_4.get())
        mass = float(self.mass.get())

        self.app.atomtype_prm[name] = [ri, ei, ci, ai, r1_4, e1_4, mass]
        self.listbox.insert(self.insert_index, '%4s %7.2f %5.2f %5.2f %4.2f %7.2f %5.2f %5.2f' %
                                               (name.ljust(4), ri, ei, ci, ai, r1_4, e1_4, mass))

        self.destroy()

    def dialog_box(self):
        self.title(self.new_title)
        self.config(bg=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10,10), pady=(10,10))

        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        frame2= Frame(mainframe, bg=self.main_color)
        frame2.grid(row=1,column=0)

        #Labels:
        name = Label(frame1, text='Name', bg=self.main_color)
        name.grid(row=0, column=0)

        ri = Label(frame1, text='Ri', bg=self.main_color)
        ri.grid(row=0, column=1)

        ei = Label(frame1, text='Ei', bg=self.main_color)
        ei.grid(row=0, column=2)

        ci = Label(frame1, text='Ci', bg=self.main_color)
        ci.grid(row=0, column=3)

        ai = Label(frame1, text='ai', bg=self.main_color)
        ai.grid(row=0, column=4)

        ri1_4 = Label(frame1, text='Ri(1-4)', bg=self.main_color)
        ri1_4.grid(row=0, column=5)

        ei1_4 = Label(frame1, text='Ei(1-4)', bg=self.main_color)
        ei1_4.grid(row=0, column=6)

        mass = Label(frame1, text='Mass', bg=self.main_color)
        mass.grid(row=0, column=7)

        #Entry fields:
        self.name = Entry(frame1, width=7, bg=self.main_color)
        self.name.grid(row=1, column=0)

        self.ri = Entry(frame1, width=7, bg=self.main_color)
        self.ri.grid(row=1, column=1)

        self.ei = Entry(frame1, width=7, bg=self.main_color)
        self.ei.grid(row=1, column=2)

        self.ci = Entry(frame1, width=7, bg=self.main_color)
        self.ci.grid(row=1, column=3)

        self.ai = Entry(frame1, width=7, bg=self.main_color)
        self.ai.grid(row=1, column=4)

        self.ri1_4 = Entry(frame1, width=7, bg=self.main_color)
        self.ri1_4.grid(row=1, column=5)

        self.ei1_4 = Entry(frame1, width=7, bg=self.main_color)
        self.ei1_4.grid(row=1, column=6)

        self.mass = Entry(frame1, width=7, bg=self.main_color)
        self.mass.grid(row=1, column=7)

        #Save/Cancel
        save = Button(frame2, text='Save', highlightbackground=self.main_color, command=self.save_prm)
        save.grid(row=0,column=0)

        cancel = Button(frame2, text='Cancel', highlightbackground=self.main_color, command=self.destroy)
        cancel.grid(row=0, column=1)


class EditBondParameters(Toplevel):
    """Edit or add Atom parameters in EVB setup"""

    def __init__(self, app, root, insert_listbox, insert_index='END', edit=False):
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self._root = root
        self.listbox = insert_listbox
        self.insert_index = insert_index
        self.edit = edit

        self.kb_var = StringVar()
        self.de_var = StringVar()

        self.kb_var.trace('w', self.kb_changed)

        if self.edit:
            self.new_title = 'Edit bond parameter type'
        else:
            self.new_title = 'Add bond parameter'
        self.dialog_box()

        self.fill_entries()

    def kb_changed(self, *args):
        try:
            alpha = float(self.alpha.get())
        except:
            return

        if alpha == float(2):
            try:
                kb = float(self.kb.get())
                de = (kb / 8.0)
                self.de.delete(0, END)
                self.de.insert(0, '%6.2f' % de)
            except:
                return


    def fill_entries(self):
        if self.edit:
            line = self.app.bondtypes_listbox.get(self.insert_index)
            self.bond_nr, de, alpha, rb = line.split()[0:4]
            try:
                kb = (8.0 * float(de))
            except:
                kb = '??'
            self.bond_nr = 1 + int(self.app.bondtypes_listbox.get(END).split()[0])
            type1, type2 = line.split('!')[1].split()
            self.de.insert(0, de)
            self.alpha.insert(0, alpha)
            self.rb.insert(0, rb)
            self.kb.insert(0, kb)
            self.type1.insert(0, type1.strip())
            self.type2.insert(0, type2.strip())


    def save_prm(self):
        #Add parameters to self.bond_prm.
        #Note if parameter bond exist, it will be overwritten!

        try:
            de = float(self.de.get())
        except:
            try:
                de = float(self.kb.get()) / 8.00
            except:
                de = 0.00
        try:
            alpha = float(self.alpha.get())
        except:
            alpha = 0.00

        try:
            rb = float(self.rb.get())
        except:
            rb = 0.00
        try:
            kb = float(self.kb.get())
        except:
            try:
                kb = de * 8.00
            except:
                kb = 0.00

        state1 = '%s %s' % (self.type1.get(), self.type2.get())
        state1_rev = '%s %s' % (state1.split()[0], state1.split()[1])
        if state1_rev in self.app.bond_prm.keys():
            state1 = state1_rev

        self.app.bond_prm[state1] = [de, alpha, rb, kb]

        print self.app.bond_prm[state1]
        self.app.update_q_bonds()
        self.app.update_status()

    def dialog_box(self):
        self.title(self.new_title)
        self.config(bg=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))


        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=1, column=0)

        bottomframe = Frame(mainframe, bg=self.main_color)
        bottomframe.grid(row=2, column=0)


        #Labels
        de = Label(frame1, text='De', bg = self.main_color)
        de.grid(row=0, column=0)

        alpha = Label(frame1, text=u"\N{GREEK SMALL LETTER ALPHA}", bg=self.main_color)
        alpha.grid(row=0, column=1)

        rb = Label(frame1, text='R0', bg=self.main_color)
        rb.grid(row=0, column=2)

        kb = Label(frame1, text='Kb', bg=self.main_color)
        kb.grid(row=0, column=3)

        type1 = Label(frame1, text='Type 1', bg=self.main_color)
        type1.grid(row=0, column=4)

        type2 = Label(frame1, text='Type 2', bg=self.main_color)
        type2.grid(row=0, column=5)

        self.de = Entry(frame1, widt=7, bg = self.main_color)
        self.de.grid(row=1, column=0)

        self.alpha = Entry(frame1, widt=7, bg=self.main_color)
        self.alpha.grid(row=1, column=1)

        self.rb = Entry(frame1, widt=7, bg=self.main_color)
        self.rb.grid(row=1, column=2)

        self.kb = Entry(frame1, widt=7, bg=self.main_color, textvariable=self.kb_var)
        self.kb.grid(row=1, column=3)

        self.type1 = Entry(frame1, widt=7, bg=self.main_color)
        self.type1.grid(row=1, column=4)

        self.type2 = Entry(frame1, widt=7, bg=self.main_color)
        self.type2.grid(row=1, column=5)

        #Close/Save
        save_button = Button(bottomframe, text='Save', highlightbackground=self.main_color, command=self.save_prm)
        save_button.grid(row=0, column=0, sticky='e')

        close_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=0, column=1, sticky='w')


class EditAngleParameters(Toplevel):
    """Edit or add Atom parameters in EVB setup"""

    def __init__(self, app, root, insert_listbox, insert_index='END', edit=False):
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self._root = root
        self.listbox = insert_listbox
        self.insert_index = insert_index
        self.edit = edit

        if self.edit:
            self.new_title = 'Edit angle parameter type'
        else:
            self.new_title = 'Add angle parameter'
        self.dialog_box()

        self.fill_entries()


    def fill_entries(self):
        if self.edit:
            line = self.app.angletypes_listbox.get(self.insert_index)

            self.angle_nr, k, theta = line.split()[0:3]
            type1, type2, type3 = line.split('!')[1].split()[0:]

            self.force.insert(0, k.strip())
            self.theta.insert(0, theta.strip())
            self.type1.insert(0, type1.strip())
            self.type2.insert(0, type2.strip())
            self.type3.insert(0, type3.strip())

    def save_prm(self):
        #Add parameters to self.bond_prm.
        #Note if parameter bond exist, it will be overwritten!

        force = self.force.get()
        if len(force) < 1:
            force = '??'
        theta = self.theta.get()
        if len(theta) < 1:
            theta = '??'

        angle = '%s %s %s' % (self.type1.get(), self.type2.get(), self.type3.get())
        angle_rev = '%s %s %s' % (self.type3.get(), self.type2.get(), self.type1.get())

        self.app.angle_prm[angle] = [force, theta]
        self.app.angle_prm[angle_rev] = [force, theta]

        self.app.update_angles()
        self.app.update_status()

    def dialog_box(self):
        self.title(self.new_title)
        self.config(bg=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=1, column=0)

        bottomframe = Frame(mainframe, bg=self.main_color)
        bottomframe.grid(row=2, column=0)

        #Labels
        force = Label(frame1, text='K', bg = self.main_color)
        force.grid(row=0, column=0)

        theta = Label(frame1, text=u"\N{GREEK CAPITAL LETTER THETA}", bg=self.main_color)
        theta.grid(row=0, column=1)

        type1 = Label(frame1, text='Type 1', bg=self.main_color)
        type1.grid(row=0, column=2)

        type2 = Label(frame1, text='Type 2', bg=self.main_color)
        type2.grid(row=0, column=3)

        type3 = Label(frame1, text='Type 3', bg=self.main_color)
        type3.grid(row=0, column=4)


        self.force = Entry(frame1, widt=7, bg = self.main_color)
        self.force.grid(row=1, column=0)

        self.theta = Entry(frame1, widt=7, bg=self.main_color)
        self.theta.grid(row=1, column=1)

        self.type1 = Entry(frame1, widt=7, bg=self.main_color)
        self.type1.grid(row=1, column=2)

        self.type2 = Entry(frame1, widt=7, bg=self.main_color)
        self.type2.grid(row=1, column=3)

        self.type3 = Entry(frame1, widt=7, bg=self.main_color)
        self.type3.grid(row=1, column=4)

        #Close/Save
        save_button = Button(bottomframe, text='Save', highlightbackground=self.main_color, command=self.save_prm)
        save_button.grid(row=0, column=0, sticky='e')

        close_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=0, column=1, sticky='w')


class EditTorsionParameters(Toplevel):
    """Edit or add Atom parameters in EVB setup"""

    def __init__(self, app, root, insert_listbox, insert_index='END', edit=False):
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self._root = root
        self.listbox = insert_listbox
        self.insert_index = insert_index
        self.edit = edit

        if self.edit:
            self.new_title = 'Edit torsion parameter type'
        else:
            self.new_title = 'Add torsion parameter'
        self.dialog_box()

        self.fill_entries()

    def fill_entries(self):
        if self.edit:
            line = self.app.torsiontypes_listbox.get(self.insert_index)
            print line
            self.torsion_nr, k, minima, phase, paths = line.split()[0:5]
            type1, type2, type3, type4 = line.split('!')[1].split()[0:]

            self.force.insert(0, k.strip())
            self.minima.insert(0, minima.strip())
            self.paths.insert(0, paths)
            self.phase.insert(0, phase)
            self.type1.insert(0, type1.strip())
            self.type2.insert(0, type2.strip())
            self.type3.insert(0, type3.strip())
            self.type4.insert(0, type4.strip())

    def save_prm(self):
        #Add parameters to self.bond_prm.
        #Note if parameter bond exist, it will be overwritten!

        try:
            force = float(self.force.get())
            minima = int(self.minima.get())
            phase = float(self.phase.get())
            paths = float(self.paths.get())
        except:
            print 'Invalid parameter value(s) given. Parameters not saved!'
            return

        if abs(minima) > 3:
            print 'Minima can not exceed 3!. Parameters not saved!'
            return

        torsion = '%s %s %s %s' % (self.type1.get(), self.type2.get(), self.type3.get(), self.type4.get())
        torsion_rev = '%s %s %s %s' % (self.type4.get(), self.type3.get(), self.type2.get(), self.type1.get())

        if torsion not in self.app.torsion_prm.keys():
            self.app.torsion_prm[torsion] = [[0,0.0,1],[0, 180.0,1],[0, 0.0,1]]
        if torsion_rev not in self.app.torsion_prm.keys():
            self.app.torsion_prm[torsion_rev] = [[0,0.0,1],[0, 180.0,1],[0, 0.0,1]]

        self.app.torsion_prm[torsion][abs(minima) - 1] = [force, phase, paths]
        self.app.torsion_prm[torsion_rev][abs(minima) - 1] = [force, phase, paths]

        self.app.update_torsions()
        self.app.update_status()

    def dialog_box(self):
        self.title(self.new_title)
        self.config(bg=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=1, column=0)

        bottomframe = Frame(mainframe, bg=self.main_color)
        bottomframe.grid(row=2, column=0)

        #Labels
        force = Label(frame1, text='K', bg = self.main_color)
        force.grid(row=0, column=0)

        minima = Label(frame1, text='min.', bg=self.main_color)
        minima.grid(row=0, column=1)

        phase = Label(frame1, text='phase', bg=self.main_color)
        phase.grid(row=0, column=2)

        paths = Label(frame1, text='paths', bg=self.main_color)
        paths.grid(row=0, column=3)

        type1 = Label(frame1, text='Type 1', bg=self.main_color)
        type1.grid(row=0, column=4)

        type2 = Label(frame1, text='Type 2', bg=self.main_color)
        type2.grid(row=0, column=5)

        type3 = Label(frame1, text='Type 3', bg=self.main_color)
        type3.grid(row=0, column=6)

        type4 = Label(frame1, text='Type 4', bg=self.main_color)
        type4.grid(row=0, column=7)

        self.force = Entry(frame1, widt=7, bg = self.main_color)
        self.force.grid(row=1, column=0)

        self.minima = Entry(frame1, widt=7, bg=self.main_color)
        self.minima.grid(row=1, column=1)

        self.phase = Entry(frame1, widt=7, bg=self.main_color)
        self.phase.grid(row=1, column=2)

        self.paths = Entry(frame1, widt=7, bg=self.main_color)
        self.paths.grid(row=1, column=3)

        self.type1 = Entry(frame1, widt=7, bg=self.main_color)
        self.type1.grid(row=1, column=4)

        self.type2 = Entry(frame1, widt=7, bg=self.main_color)
        self.type2.grid(row=1, column=5)

        self.type3 = Entry(frame1, widt=7, bg=self.main_color)
        self.type3.grid(row=1, column=6)

        self.type4 = Entry(frame1, widt=7, bg=self.main_color)
        self.type4.grid(row=1, column=7)

        #Close/Save
        save_button = Button(bottomframe, text='Save', highlightbackground=self.main_color, command=self.save_prm)
        save_button.grid(row=0, column=0, sticky='e')

        close_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=0, column=1, sticky='w')


class EditImproperParameters(Toplevel):
    """Edit or add improper parameters in EVB setup"""

    def __init__(self, app, root, insert_listbox, insert_index='END', edit=False):
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self._root = root
        self.listbox = insert_listbox
        self.insert_index = insert_index
        self.edit = edit

        if self.edit:
            self.new_title = 'Edit improper parameter type'
        else:
            self.new_title = 'Add improper parameter'
        self.dialog_box()

        self.fill_entries()

    def fill_entries(self):
        if self.edit:
            line = self.app.impropertypes_listbox.get(self.insert_index)
            print line
            self.improper_nr, k, phase = line.split()[0:3]
            type1, type2, type3, type4 = line.split('!')[1].split()[0:]

            self.force.insert(0, k.strip())
            self.phase.insert(0, phase)
            self.type1.insert(0, type1.strip())
            self.type2.insert(0, type2.strip())
            self.type3.insert(0, type3.strip())
            self.type4.insert(0, type4.strip())

    def save_prm(self):
        #Add parameters to self.bond_prm.
        #Note if parameter bond exist, it will be overwritten!

        try:
            force = float(self.force.get())
            phase = float(self.phase.get())
        except:
            print 'Invalid parameter value(s) given. Parameters not saved!'
            return

        improper = '%s %s %s %s' % (self.type1.get(), self.type2.get(), self.type3.get(), self.type4.get())
        improper_rev = '%s %s %s %s' % (self.type4.get(), self.type3.get(), self.type2.get(), self.type1.get())

        self.app.improper_prm[improper] = [force, phase]
        self.app.improper_prm[improper_rev] = [force, phase]

        self.app.update_impropers()
        self.app.update_status()

    def dialog_box(self):
        self.title(self.new_title)
        self.config(bg=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=1, column=0)

        bottomframe = Frame(mainframe, bg=self.main_color)
        bottomframe.grid(row=2, column=0)

        #Labels
        force = Label(frame1, text='K', bg = self.main_color)
        force.grid(row=0, column=0)

        phase = Label(frame1, text='phase', bg=self.main_color)
        phase.grid(row=0, column=2)

        type1 = Label(frame1, text='Type 1', bg=self.main_color)
        type1.grid(row=0, column=4)

        type2 = Label(frame1, text='Type 2', bg=self.main_color)
        type2.grid(row=0, column=5)

        type3 = Label(frame1, text='Type 3', bg=self.main_color)
        type3.grid(row=0, column=6)

        type4 = Label(frame1, text='Type 4', bg=self.main_color)
        type4.grid(row=0, column=7)

        self.force = Entry(frame1, widt=7, bg = self.main_color)
        self.force.grid(row=1, column=0)

        self.phase = Entry(frame1, widt=7, bg=self.main_color)
        self.phase.grid(row=1, column=2)

        self.type1 = Entry(frame1, widt=7, bg=self.main_color)
        self.type1.grid(row=1, column=4)

        self.type2 = Entry(frame1, widt=7, bg=self.main_color)
        self.type2.grid(row=1, column=5)

        self.type3 = Entry(frame1, widt=7, bg=self.main_color)
        self.type3.grid(row=1, column=6)

        self.type4 = Entry(frame1, widt=7, bg=self.main_color)
        self.type4.grid(row=1, column=7)

        #Close/Save
        save_button = Button(bottomframe, text='Save', highlightbackground=self.main_color, command=self.save_prm)
        save_button.grid(row=0, column=0, sticky='e')

        close_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=0, column=1, sticky='w')
