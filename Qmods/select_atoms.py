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

from tkinter import  EXTENDED, MULTIPLE, Button, Frame, Toplevel, Scrollbar, Listbox, END, GROOVE, LEFT

import tkinter.font


class AtomSelectRange(Toplevel):
    """General class to select atoms or atomranges.
    Can input first and last atom in range in entry1 and entry 2.
    If entry1 == entry2, single selection is expected, and atom nr
    is inserted into entry 1 only.
    If writepdb = True, selected atomrange will be written to a new
    pdb file and pdb filename will be inserted into entry 1 (== entry2)"""

    def __init__(self, app, root, pdbfile, entry1, entry2, selectmode='EXTENDED', writepdb=False):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.pdbfile = pdbfile
        if selectmode == 'EXTENDED':
            self._selectmode = EXTENDED
        elif selectmode == 'MULTIPLE':
            self._selectmode = MULTIPLE
        self.writepdb = writepdb
        self.dialog_box()
        self.update_list()

        self.entry1 = entry1
        self.entry2 = entry2

    def update_list(self):
        pdblist = open(self.pdbfile,'r').readlines()
        for line in pdblist:
            if len(line.split()) > 3:
                if 'ATOM' in line.split()[0]:
                    txt = line.split('\n')[0]
                    self.listbox.insert(END, txt)

    def get_selected(self):
        """
        Gets coordinates from selection and returns it to the prepare topology window
        """
        atoms = []
        items = list(map(int, self.listbox.curselection()))

        if self.entry2 != 'qatoms':
            for item in items:
                atom = self.listbox.get(item).split()[1]
                atoms.append(atom)
        #If entry2 = qatoms, EVB setup class:
        elif self.entry2 == 'qatoms':
            for item in items:
                #If atomnumber is not already in list, append it:
                qnr = len(self.app.q_atom_nr)
                if int(self.listbox.get(item).split()[1]) not in list(map(int, list(self.app.q_atom_nr.values()))):
                    qnr += 1
                    self.app.q_atom_nr[qnr] = int(self.listbox.get(item).split()[1])
                    resi = self.listbox.get(item)[17:27].strip()
                    atomname = self.listbox.get(item)[12:17].strip()
                    self.app.q_atom_res[qnr] = resi
                    self.app.q_atom_name[qnr] = atomname
                    self.app.q_notes[qnr] = '%4s %s' % (atomname.ljust(4), resi)

            self.app.update_q_atoms()
            self.app.read_lib()

        if not self.writepdb and self.entry2 != 'qatoms':
            self.entry1.delete(0, END)
            self.entry1.insert(0, atoms[0])
            if self.entry2 == 'qatoms':
                #Q-atom selection. insert atomnr and !atomname:
                for i in range(len(atoms)):
                    self.entry1.insert(END, atoms[i])

            elif self.entry2 != self.entry1 and self.entry2 != 'qatoms':
                self.entry2.delete(0, END)
                self.entry2.insert(0, atoms[-1])

        elif self.writepdb:
            newfile = []
            got_resname = False
            newfilename = self.app.workdir + '/' + 'LIG.pdb'
            pdbfile = open(self.pdbfile,'r').readlines()
            print(atoms)
            for line in pdbfile:
                if 'ATOM' in line:
                    if line.split()[1] in atoms:
                        newfile.append(line)
                        if not got_resname:
                            newfilename = self.app.workdir + '/' + line[17:21].strip() + '.pdb'
                            got_resname = True
            output = open(newfilename,'w')
            for line in newfile:
                output.write(line)
            output.close()
            self.entry1.delete(0, END)
            self.app.ligand_pdb = newfilename
            self.entry1.insert(0, newfilename.split('/')[-1])
            self.app.update_progress()

        self.destroy()

    def cancel(self):
        self.destroy()

    def dialog_box(self):
        """Defines the outlook of Setup MD window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('Select atoms')
        self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=10, padx=(10,0))

        listbox_scroll = Scrollbar(left_frame)
        listbox_scroll.grid(row = 0, rowspan = 10, column = 11, sticky = 'nsw')
        self.listbox = Listbox(left_frame, yscrollcommand = listbox_scroll.set, width = 60, height=30, highlightthickness = 0, relief = GROOVE, selectmode=self._selectmode)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 0, rowspan = 10, column = 0, columnspan = 10, sticky = 'w')
        self.listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        select_button = Button(left_frame, text = 'Select', command = self.get_selected)
        select_button.grid(row = 11, column = 0, columnspan = 6, sticky = 'e')
        select_button.config(highlightbackground = self.main_color)

        cancel_button = Button(left_frame, text = 'Cancel', command = self.cancel)
        cancel_button.grid(row=11, column = 6, columnspan = 6, sticky = 'w' )
        cancel_button.config(highlightbackground = self.main_color)
