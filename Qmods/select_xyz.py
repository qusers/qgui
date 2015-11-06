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

from Tkinter import Button,Frame, Toplevel, Scrollbar, Listbox, END, GROOVE, LEFT
import tkFont

import prepareTopology as pt

class AtomSelect(Toplevel):
    """Implements a dialog-box when Change sim. center is chosen in setup MD"""

    def __init__(self, app, root, pdbfile, xentry,yentry,zentry):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.pdbfile = pdbfile
        self.xentry = xentry
        self.yentry = yentry
        self.zentry = zentry
        self.dialog_box()
        self.update_list()


    def update_list(self):
        pdblist = open(self.pdbfile,'r').readlines()
        for line in pdblist:
            try:
                if 'ATOM' in line.split()[0] or 'HETATM' in line.split()[0]:
                    txt = line.split('\n')[0]
                    self.listbox.insert(END, txt)
            except:
                continue

    def listbox_clicked(self, *args):
        if not self.app.session:
            return

        selected = int( self.listbox.curselection()[0])

        atom_nr = self.listbox.get(selected).split()[1]

        self.app.session.stdin.write('select none \nhide labels\n')
        #self.app.session.stdin.write('select id %s\n' % atom_nr )
        self.app.session.stdin.write('label id %s, name\n' % atom_nr)
        self.app.session.stdin.write('zoom id %s, buffer=%d\n' % (atom_nr, self.app.pymol_zoom))

    def get_selected(self):
        """
        Gets coordinates from selection and returns it to the prepare topology window
        """
        list_index = int(self.listbox.curselection()[0])
        x,y,z = self.listbox.get(list_index).split()[5:8]
        self.xentry.delete(0, END)
        self.yentry.delete(0, END)
        self.zentry.delete(0, END)
        self.xentry.insert(0, x)
        self.yentry.insert(0, y)
        self.zentry.insert(0, z)
        self.destroy()

    def cancel(self):
        self.destroy()

    def dialog_box(self):
        """Defines the outlook of Setup MD window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('Select simulation center')
        self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=10, padx=(10,0))

        listbox_scroll = Scrollbar(left_frame)
        listbox_scroll.grid(row = 0, rowspan = 10, column = 11, sticky = 'nsw')
        self.listbox = Listbox(left_frame, yscrollcommand = listbox_scroll.set, width = 60, height=30, highlightthickness = 0, relief = GROOVE)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 0, rowspan = 10, column = 0, columnspan = 10, sticky = 'w')
        self.listbox.config(font=tkFont.Font(family="Courier", size=12))

        select_button = Button(left_frame, text = 'Select', command = self.get_selected)
        select_button.grid(row = 11, column = 0, columnspan = 6, sticky = 'e')
        select_button.config(highlightbackground = self.main_color)

        cancel_button = Button(left_frame, text = 'Cancel', command = self.cancel)
        cancel_button.grid(row=11, column = 6, columnspan = 6, sticky = 'w' )
        cancel_button.config(highlightbackground = self.main_color)
        self.listbox.bind('<<ListboxSelect>>', self.listbox_clicked)
