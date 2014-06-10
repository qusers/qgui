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

from Tkinter import StringVar, Entry, Spinbox, Text, Label, Button,Frame, Toplevel, Scrollbar, Listbox, Checkbutton, DISABLED, END, GROOVE, LEFT, IntVar
from select_atoms import AtomSelectRange
from select_xyz import AtomSelect
import tkFont

class AddSequenceRestraints(Toplevel):
    """Restraints for MD setup """

    def __init__(self, app, root, pdbfile, newtitle = ' ', restraint_list = []):
        """
        Initialize the add restraint window
        """
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
       
        self.pdbfile = pdbfile
        
        self.restraints = restraint_list
        self.newtitle = newtitle
        
        self.h_restrain = IntVar()
        self.group = IntVar()
        self.h_restrain.set(1)
        self.group.set(1)
        self.temp = []
        
        self.dialog_box()
        self.update_list()
        self.force_entry.delete(0, END)
        self.force_entry.insert(0, 5)

    def update_list(self):
        if len(self.restraints) > 0:
            for restraint in self.restraints:
                self.restraint_list.insert(END, restraint)
                temp = self.restraint_list.get(0,END)
            for entry in temp:
                self.temp.append(entry)
    
    def select_atoms(self):
        """
        opens dialog to select atoms and inserts first into entry 1 and last into entry 2
        """
        entry1 = self.firstatom_entry 
        entry2 = self.lastatom_entry

        self.select_atomrange = AtomSelectRange(self,self.root, self.pdbfile, entry1, entry2)
        self.select_atomrange.configure(bg = self.main_color)
        self.select_atomrange.title('Select atoms')
        self.select_atomrange.resizable()
        
    def add_restraints(self):
        atom_i = self.firstatom_entry.get().strip()
        atom_j = self.lastatom_entry.get().strip()
        force = self.force_entry.get().strip()
        h_rest = self.h_restrain.get()
        group = self.group.get()
        
        self.restraint_list.insert(END,'%6s %6s %8.2f %2s %2s' % (atom_i, atom_j, float(force), h_rest, group))
        self.temp.append('%6s %6s %8.2f %2s %2s' % (atom_i, atom_j, float(force), h_rest, group))
    
    def delete_restraints(self):
        try:
            selected = int(self.restraint_list.curselection()[0])
        except:
            return

        self.restraint_list.delete(selected)
        del self.temp[selected]
        
    
    def save_restraints(self):
        self.app.sequence_restraints = self.temp
        
        self.destroy()

    def dialog_box(self):
        """Defines the outlook of add/edit equilibration window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        #self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=(10,10), padx=(10,10))
        
        atom1_label = Label(left_frame, text = 'Atom i')
        atom1_label.grid(row = 0, column = 0)
        atom1_label.config(bg=self.main_color)
        
        self.firstatom_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.firstatom_entry.grid(row=1,column = 0)
        
        atom2_label = Label(left_frame, text = 'Atom j')
        atom2_label.grid(row = 0, column = 1)
        atom2_label.config(bg=self.main_color)
        
        self.lastatom_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.lastatom_entry.grid(row=1,column = 1)
        
        
        force_label = Label(left_frame, text = 'Force')
        force_label.grid(row = 0, column = 4)
        force_label.config(bg=self.main_color)
        
        self.force_entry = Spinbox(left_frame, width = 6, highlightthickness = 0, relief = GROOVE, from_=0, to=1000)
        self.force_entry.grid(row=1,column = 4)
        
        select_button = Button(left_frame, text = 'select', command = self.select_atoms)
        select_button.grid(row = 1, column = 3, sticky = 'w')
        select_button.config(highlightbackground = self.main_color)
        
        restrain_h = Label(left_frame, text = ' H ')
        restrain_h.grid(row = 0, column = 5)
        restrain_h.config(bg=self.main_color)
        
        group_label = Label(left_frame, text = 'Group')
        group_label.grid(row = 0, column = 6)
        group_label.config(bg=self.main_color)
        
        h_check = Checkbutton(left_frame, variable = self.h_restrain)
        h_check.grid(row=1, column = 5)
        h_check.config(bg=self.main_color)
        
        group_check = Checkbutton(left_frame, variable = self.group)
        group_check.grid(row = 1, column=6) 
        group_check.config(bg=self.main_color)
        
        add_button = Button(left_frame, text = 'Add restraint', command = self.add_restraints)
        add_button.grid(row = 2, column = 0, columnspan=4)
        add_button.config(highlightbackground = self.main_color)

        delete_button = Button(left_frame, text = 'Delete', command = self.delete_restraints)
        delete_button.grid(row = 2, column = 3, columnspan=4)
        delete_button.config(highlightbackground = self.main_color)

        restraint_scroll = Scrollbar(left_frame)
        restraint_scroll.grid(row=3, column = 6, sticky = 'nsw')
        self.restraint_list = Listbox(left_frame, yscrollcommand = restraint_scroll.set, width = 40, height = 6, highlightthickness = 0 , relief = GROOVE)
        restraint_scroll.config(command = self.restraint_list.yview)
        self.restraint_list.grid(row = 3, column = 0, columnspan = 6,padx = (10,0),sticky = 'e')
        self.restraint_list.config(font = tkFont.Font(family="Courier", size=12))
        
        save_button = Button(left_frame, text = 'Save', command = self.save_restraints)
        save_button.grid(row = 4, column = 4) 
        save_button.config(highlightbackground = self.main_color)
        
        cancel_button = Button(left_frame, text = 'Cancel', command=self.destroy)
        cancel_button.grid(row=4, column = 5, columnspan = 2)
        cancel_button.config(highlightbackground = self.main_color)


class AddAtomRestraints(Toplevel):
    """Restraints for MD setup """

    def __init__(self, app, root, pdbfile, newtitle = ' ', restraint_list = []):
        """
        Initialize the add restraint window
        """
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
       
        self.pdbfile = pdbfile
        
        self.restraints = restraint_list
        self.newtitle = newtitle
        
        self.h_restrain = IntVar()
        self.group = IntVar()
        self.h_restrain.set(1)
        self.group.set(1)
        self.temp = []

        self.atom_changed = StringVar()
        self.atom_changed.trace('w', self.atom_written)
        
        self.dialog_box()
        self.update_list()
        self.xforce_entry.delete(0, END)
        self.yforce_entry.delete(0, END)
        self.zforce_entry.delete(0, END)
        self.xforce_entry.insert(0, 5)
        self.yforce_entry.insert(0, 5)
        self.zforce_entry.insert(0, 5)
        
        
    def update_list(self):
        if len(self.restraints) > 0:
            for restraint in self.restraints:
                self.restraint_list.insert(END, restraint)
                temp = self.restraint_list.get(0,END)
            for entry in temp:
                self.temp.append(entry)

    def atom_written(self, *args):
        """
        Traces if a new atom nr. is selected and inserts xyz coordinates.
        """
        x = ''
        y = ''
        z = ''
        atomnr = self.firstatom_entry.get().strip()
        pdbfile = open(self.pdbfile,'r').readlines()
        for line in pdbfile:
            if 'ATOM' in line:
                if line.split()[1] == atomnr:
                    x = line[28:].split()[0]
                    y = line[28:].split()[1]
                    z = line[28:].split()[2]
        self.x_entry.delete(0, END)
        self.y_entry.delete(0, END)
        self.z_entry.delete(0, END)
        self.x_entry.insert(0, x)
        self.y_entry.insert(0, y)
        self.z_entry.insert(0, z)
    
    def select_atoms(self):
        """
        opens dialog to select atoms and inserts first into entry 1 and last into entry 2
        """
        entry1 = self.firstatom_entry 
        entry2 = self.firstatom_entry

        self.select_atomrange = AtomSelectRange(self,self.root, self.pdbfile, entry1, entry2)
        self.select_atomrange.configure(bg = self.main_color)
        self.select_atomrange.title('Select atoms')
        self.select_atomrange.resizable()

    def select_coordinates(self):
        xentry = self.x_entry
        yentry = self.y_entry
        zentry = self.z_entry

        self.getxyz = AtomSelect(self, self.root, self.pdbfile, xentry,yentry,zentry)
        self.getxyz.configure(bg = self.main_color)
        self.getxyz.title('Select coordinates for atom restraint')
        self.getxyz.resizable()


        
    def add_restraints(self):
        atom_i = self.firstatom_entry.get().strip()
        if atom_i.isdigit():
            xforce = float(self.xforce_entry.get().strip())
            yforce = float(self.yforce_entry.get().strip())
            zforce = float(self.zforce_entry.get().strip())

            x = float(self.x_entry.get())
            y = float(self.y_entry.get())
            z = float(self.z_entry.get())

            if self.app.fep:
                lambda_states = self.states.get()
            else:
                lambda_states = '0'
        
            self.restraint_list.insert(END,'%6s %8.3f %8.3f %8.3f %4.1f %4.1f %4.1f %2s' %
                                       (atom_i, x,y,z, xforce, yforce, zforce, lambda_states))
            self.temp.append('%6s %8.3f %8.3f %8.3f %4.1f %4.1f %4.1f %2s' %
                                   (atom_i, x,y,z, xforce, yforce, zforce, lambda_states))
        else:
            return

    def delete_restraints(self):
        try:
            selected = int(self.restraint_list.curselection()[0])
        except:
            return

        self.restraint_list.delete(selected)
        del self.temp[selected]
        
    
    def save_restraints(self):
        self.app.atom_restraints = self.temp
        
        self.destroy()

    def dialog_box(self):
        """Defines the outlook of add/edit equilibration window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        #self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=(10,10), padx=(10,10))
        
        atom1_label = Label(left_frame, text = 'Atom')
        atom1_label.grid(row = 0, column = 0)
        atom1_label.config(bg=self.main_color)
        
        self.firstatom_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE, textvariable=self.atom_changed)
        self.firstatom_entry.grid(row=1,column = 0)
        
        select_button = Button(left_frame, text = 'select', command = self.select_atoms)
        select_button.grid(row = 1, column = 1, sticky = 'w')
        select_button.config(highlightbackground = self.main_color)
        
        x_label = Label(left_frame, text = 'x')
        x_label.grid(row = 0, column = 2)
        x_label.config(bg=self.main_color)

        y_label = Label(left_frame, text = 'y')
        y_label.grid(row = 0, column = 3)
        y_label.config(bg=self.main_color)

        z_label = Label(left_frame, text = 'z')
        z_label.grid(row = 0, column = 4)
        z_label.config(bg=self.main_color)

        
        self.x_entry = Entry(left_frame, width = 8, highlightthickness = 0, relief = GROOVE)
        self.x_entry.grid(row=1,column = 2)

        self.y_entry = Entry(left_frame, width = 8, highlightthickness = 0, relief = GROOVE)
        self.y_entry.grid(row=1,column = 3)

        self.z_entry = Entry(left_frame, width = 8, highlightthickness = 0, relief = GROOVE)
        self.z_entry.grid(row=1,column = 4)

        selectxyz_button = Button(left_frame, text = 'select', command = self.select_coordinates)
        selectxyz_button.grid(row = 1, column = 5, sticky = 'w')
        selectxyz_button.config(highlightbackground = self.main_color)
        
        
        xforce_label = Label(left_frame, text = 'Force x')
        xforce_label.grid(row = 2, column = 2)
        xforce_label.config(bg=self.main_color)

        yforce_label = Label(left_frame, text = 'Force y')
        yforce_label.grid(row = 2, column = 3)
        yforce_label.config(bg=self.main_color)

        zforce_label = Label(left_frame, text = 'Force z')
        zforce_label.grid(row = 2, column = 4)
        zforce_label.config(bg=self.main_color)

        if self.app.fep:
            states_label = Label(left_frame, text= u"\N{GREEK SMALL LETTER LAMDA}-state (0=all)", bg=self.main_color)
            states_label.grid(row=2, column= 0, columnspan=2)

            self.states = Spinbox(left_frame, width=3, highlightthickness=0, relief=GROOVE, from_=0, to=4)
            self.states.grid(row=3, column=0, columnspan=2)

        self.xforce_entry = Spinbox(left_frame, width = 6, highlightthickness = 0, relief = GROOVE, from_=0, to=1000)
        self.xforce_entry.grid(row=3,column = 2)

        self.yforce_entry = Spinbox(left_frame, width = 6, highlightthickness = 0, relief = GROOVE, from_=0, to=1000)
        self.yforce_entry.grid(row=3,column = 3)

        self.zforce_entry = Spinbox(left_frame, width = 6, highlightthickness = 0, relief = GROOVE, from_=0, to=1000)
        self.zforce_entry.grid(row=3,column = 4)

        add_button = Button(left_frame, text = 'Add restraint', command = self.add_restraints)
        add_button.grid(row = 5, column = 0, columnspan=4)
        add_button.config(highlightbackground = self.main_color)

        delete_button = Button(left_frame, text = 'Delete', command = self.delete_restraints)
        delete_button.grid(row = 5, column = 3, columnspan=2)
        delete_button.config(highlightbackground = self.main_color)

        restraint_scroll = Scrollbar(left_frame)
        restraint_scroll.grid(row=6, column = 5, sticky = 'nsw')
        self.restraint_list = Listbox(left_frame, yscrollcommand = restraint_scroll.set, width = 55, height = 6, highlightthickness = 0 , relief = GROOVE)
        restraint_scroll.config(command = self.restraint_list.yview)
        self.restraint_list.grid(row = 6, column = 0, columnspan = 5,padx = (10,0),sticky = 'e')
        self.restraint_list.config(font = tkFont.Font(family="Courier", size=12))
        
        save_button = Button(left_frame, text = 'Save', command = self.save_restraints)
        save_button.grid(row = 7, column = 3)
        save_button.config(highlightbackground = self.main_color)
        
        cancel_button = Button(left_frame, text = 'Cancel', command=self.destroy)
        cancel_button.grid(row=7, column = 4, columnspan = 2)
        cancel_button.config(highlightbackground = self.main_color)


class AddDistanceRestraints(Toplevel):
    """Restraints for MD setup """

    def __init__(self, app, root, pdbfile, newtitle = ' ', restraint_list = []):
        """
        Initialize the add restraint window
        """
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
       
        self.pdbfile = pdbfile
        
        self.restraints = restraint_list
        self.newtitle = newtitle


        self.temp = []
        
        self.dialog_box()
        self.update_list()

        self.rmin_entry.delete(0, END)
        self.rmax_entry.delete(0, END)
        self.force_entry.delete(0, END)
        self.rmin_entry.insert(0, 1.0)
        self.rmax_entry.insert(0, 2.0)
        self.force_entry.insert(0, 5)
        
        
    def update_list(self):
        if len(self.restraints) > 0:
            for restraint in self.restraints:
                self.restraint_list.insert(END, restraint)
                temp = self.restraint_list.get(0,END)
            for entry in temp:
                self.temp.append(entry)
    
    def select_atoms(self):
        """
        opens dialog to select atoms and inserts first into entry 1 and last into entry 2
        """
        entry1 = self.firstatom_entry 
        entry2 = self.lastatom_entry

        self.select_atomrange = AtomSelectRange(self,self.root, self.pdbfile, entry1, entry2, 'MULTIPLE')
        self.select_atomrange.configure(bg = self.main_color)
        self.select_atomrange.title('Select atoms')
        self.select_atomrange.resizable()
        
    def add_restraints(self):
        try:
            atom_i = self.firstatom_entry.get().strip()
            atom_j = self.lastatom_entry.get().strip()
            if len(atom_i) == 0 or len(atom_j) == 0:
                return
            rmin = float(self.rmin_entry.get())
            rmax = float(self.rmax_entry.get())
            force = float(self.force_entry.get())

            if self.app.fep:
                lambda_states = self.states.get()
            else:
                lambda_states = '0'
        
            self.restraint_list.insert(END,'%6s %6s  %3.1f %3.1f %5.1f %2s' %
                                       (atom_i, atom_j, rmin, rmax, force, lambda_states))
            self.temp.append('%6s %6s  %3.1f %3.1f %5.1f %2s' %
                                       (atom_i, atom_j, rmin, rmax, force, lambda_states))
        except:
            return

    def delete_restraints(self):
        try:
            selected = int(self.restraint_list.curselection()[0])
        except:
            return

        self.restraint_list.delete(selected)
        del self.temp[selected]
        
    
    def save_restraints(self):
        self.app.distance_restraints = self.temp
        
        self.destroy()

    def dialog_box(self):
        """Defines the outlook of add/edit equilibration window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        #self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=(10,10), padx=(10,10))
        
        atom1_label = Label(left_frame, text = 'Atom i')
        atom1_label.grid(row = 0, column = 0)
        atom1_label.config(bg=self.main_color)
        
        self.firstatom_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.firstatom_entry.grid(row=1,column = 0)
        
        atom2_label = Label(left_frame, text = 'Atom j')
        atom2_label.grid(row = 0, column = 1)
        atom2_label.config(bg=self.main_color)
        
        self.lastatom_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.lastatom_entry.grid(row=1,column = 1)
        
        
        force_label = Label(left_frame, text = 'Force')
        force_label.grid(row = 0, column = 6)
        force_label.config(bg=self.main_color)
        
        self.force_entry = Spinbox(left_frame, width = 6, highlightthickness = 0, relief = GROOVE, from_=0, to=1000)
        self.force_entry.grid(row=1,column = 6)

        if self.app.fep:
            states_label = Label(left_frame, text= u"\N{GREEK SMALL LETTER LAMDA}-state (0=all)", bg=self.main_color)
            states_label.grid(row=0, column= 7, columnspan=2)

            self.states = Spinbox(left_frame, width=3, highlightthickness=0, relief=GROOVE, from_=0, to=4)
            self.states.grid(row=1, column=7, columnspan=2)
        
        select_button = Button(left_frame, text = 'select', command = self.select_atoms)
        select_button.grid(row = 1, column = 3, sticky = 'w')
        select_button.config(highlightbackground = self.main_color)
        
        rmin_label = Label(left_frame, text = 'r min')
        rmin_label.grid(row = 0, column = 4)
        rmin_label.config(bg = self.main_color)

        rmin_label = Label(left_frame, text = 'r max')
        rmin_label.grid(row = 0, column = 5)
        rmin_label.config(bg = self.main_color)

        self.rmin_entry = Spinbox(left_frame, width = 3, highlightthickness = 0, relief = GROOVE,
                                  from_=0.0, to=99.0, increment=0.1)
        self.rmin_entry.grid(row=1,column = 4)

        self.rmax_entry = Spinbox(left_frame, width = 3, highlightthickness = 0, relief = GROOVE,
                                  from_=0.0, to=99.0, increment=0.1)
        self.rmax_entry.grid(row=1,column = 5)
        
        add_button = Button(left_frame, text = 'Add restraint', command = self.add_restraints)
        add_button.grid(row = 2, column = 0, columnspan=4)
        add_button.config(highlightbackground = self.main_color)

        delete_button = Button(left_frame, text = 'Delete', command = self.delete_restraints)
        delete_button.grid(row = 2, column = 3, columnspan=4)
        delete_button.config(highlightbackground = self.main_color)

        restraint_scroll = Scrollbar(left_frame)
        restraint_scroll.grid(row=3, column = 6, sticky = 'nsw')
        self.restraint_list = Listbox(left_frame, yscrollcommand = restraint_scroll.set, width = 40, height = 6, highlightthickness = 0 , relief = GROOVE)
        restraint_scroll.config(command = self.restraint_list.yview)
        self.restraint_list.grid(row = 3, column = 0, columnspan = 6,padx = (10,0),sticky = 'e')
        self.restraint_list.config(font = tkFont.Font(family="Courier", size=12))
        
        save_button = Button(left_frame, text = 'Save', command = self.save_restraints)
        save_button.grid(row = 4, column = 4) 
        save_button.config(highlightbackground = self.main_color)
        
        cancel_button = Button(left_frame, text = 'Cancel', command=self.destroy)
        cancel_button.grid(row=4, column = 5, columnspan = 2)
        cancel_button.config(highlightbackground = self.main_color)
        

class AddWallRestraints(Toplevel):
    """Restraints for MD setup """

    def __init__(self, app, root, pdbfile, newtitle = ' ', restraint_list = []):
        """
        Initialize the add restraint window
        """
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
       
        self.pdbfile = pdbfile
        
        self.restraints = restraint_list
        self.newtitle = newtitle
        
        self.h_restrain = IntVar()
        self.h_restrain.set(0)
        self.temp = []
        
        self.dialog_box()
        self.update_list()
        self.force_entry.delete(0, END)
        self.force_entry.insert(0, 5)
        self.radius.delete(0, END)
        self.radius.insert(0, 10)
        
        
    def update_list(self):
        if len(self.restraints) > 0:
            for restraint in self.restraints:
                self.restraint_list.insert(END, restraint)
                temp = self.restraint_list.get(0,END)
            for entry in temp:
                self.temp.append(entry)
    
    def select_atoms(self):
        """
        opens dialog to select atoms and inserts first into entry 1 and last into entry 2
        """
        entry1 = self.firstatom_entry 
        entry2 = self.lastatom_entry

        self.select_atomrange = AtomSelectRange(self,self.root, self.pdbfile, entry1, entry2)
        self.select_atomrange.configure(bg = self.main_color)
        self.select_atomrange.title('Select atoms')
        self.select_atomrange.resizable()
        
    def add_restraints(self):
        try:
            force = float(self.force_entry.get())
            morce_dept = float(self.morce_depth.get())
            morce_alpha = float(self.morce_alpha.get())
            radius = float(self.radius.get())
        except:
            return

        atom_i = self.firstatom_entry.get().strip()
        atom_j = self.lastatom_entry.get().strip()
        h_rest = self.h_restrain.get()

        if len(atom_i) > 0 and len(atom_j) > 0:
            self.restraint_list.insert(END,'%6s %6s %5.1f %5.2f %5.2f %5.2f %2s' %
                                       (atom_i, atom_j, radius, force, morce_dept, morce_alpha, h_rest))
            self.temp.append('%6s %6s %5.1f %5.2f %5.2f %5.2f %2s' %
                                       (atom_i, atom_j, radius, force, morce_dept, morce_alpha, h_rest))
    
    def delete_restraints(self):
        try:
            selected = int(self.restraint_list.curselection()[0])
        except:
            return

        self.restraint_list.delete(selected)
        del self.temp[selected]
        
    
    def save_restraints(self):
        self.app.wall_restraints = self.temp
        
        self.destroy()

    def dialog_box(self):
        """Defines the outlook of add/edit equilibration window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        #self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=(10,10), padx=(10,10))
        
        atom1_label = Label(left_frame, text = 'Atom i')
        atom1_label.grid(row = 0, column = 0)
        atom1_label.config(bg=self.main_color)
        
        self.firstatom_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.firstatom_entry.grid(row=1,column = 0)
        
        atom2_label = Label(left_frame, text = 'Atom j')
        atom2_label.grid(row = 0, column = 1)
        atom2_label.config(bg=self.main_color)
        
        self.lastatom_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.lastatom_entry.grid(row=1,column = 1)

        select_button = Button(left_frame, text = 'select', command = self.select_atoms)
        select_button.grid(row = 1, column = 2, sticky = 'w')
        select_button.config(highlightbackground=self.main_color)

        radius = Text(left_frame, width=3, height = 1, bg=self.main_color, borderwidth=0, highlightthickness = 0)
        radius.tag_configure("subscript", offset=-3)
        radius.insert("insert","R","","0","subscript")
        radius.grid(row = 0, column = 3)
        radius.configure(state=DISABLED)

        self.radius = Spinbox(left_frame, width=4,  highlightthickness=0, relief=GROOVE, from_=0, to=99)
        self.radius.grid(row = 1, column=3)
        
        force_label = Label(left_frame, text = 'Force')
        force_label.grid(row = 0, column = 4)
        force_label.config(bg=self.main_color)
        
        self.force_entry = Spinbox(left_frame, width = 6, highlightthickness = 0, relief = GROOVE, from_=0, to=1000)
        self.force_entry.grid(row=1,column = 4)

        morce_depth = Text(left_frame, width=3, height = 1, bg=self.main_color, borderwidth=0, highlightthickness = 0)
        morce_depth.tag_configure("subscript", offset=-3)
        morce_depth.insert("insert","D","","e","subscript")
        morce_depth.grid(row = 0, column = 5)
        morce_depth.configure(state=DISABLED)

        self.morce_depth = Spinbox(left_frame, width=4,  highlightthickness=0, relief=GROOVE, from_=0, to=1000)
        self.morce_depth.grid(row = 1, column=5)

        alpha = Label(left_frame, text = u"\N{GREEK SMALL LETTER ALPHA}", bg=self.main_color)
        alpha.grid(row=0, column=6)

        self.morce_alpha = Spinbox(left_frame, width=4,  highlightthickness=0, relief=GROOVE, from_=0, to=1000)
        self.morce_alpha.grid(row = 1, column=6)

        restrain_h = Label(left_frame, text = ' H ')
        restrain_h.grid(row = 0, column = 7)
        restrain_h.config(bg=self.main_color)
        
        h_check = Checkbutton(left_frame, variable = self.h_restrain)
        h_check.grid(row=1, column = 7)
        h_check.config(bg=self.main_color)

        
        add_button = Button(left_frame, text = 'Add restraint', command = self.add_restraints)
        add_button.grid(row = 2, column = 0, columnspan=4)
        add_button.config(highlightbackground = self.main_color)

        delete_button = Button(left_frame, text = 'Delete', command = self.delete_restraints)
        delete_button.grid(row = 2, column = 3, columnspan=4)
        delete_button.config(highlightbackground = self.main_color)

        restraint_scroll = Scrollbar(left_frame)
        restraint_scroll.grid(row=3, column = 6, sticky = 'nsw')
        self.restraint_list = Listbox(left_frame, yscrollcommand = restraint_scroll.set, width = 50, height = 6, highlightthickness = 0 , relief = GROOVE)
        restraint_scroll.config(command = self.restraint_list.yview)
        self.restraint_list.grid(row = 3, column = 0, columnspan = 6,padx = (10,0),sticky = 'e')
        self.restraint_list.config(font = tkFont.Font(family="Courier", size=12))
        
        save_button = Button(left_frame, text = 'Save', command = self.save_restraints)
        save_button.grid(row = 4, column = 4) 
        save_button.config(highlightbackground = self.main_color)
        
        cancel_button = Button(left_frame, text = 'Cancel', command=self.destroy)
        cancel_button.grid(row=4, column = 5, columnspan = 2)
        cancel_button.config(highlightbackground = self.main_color)
        

