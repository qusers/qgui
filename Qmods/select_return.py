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

from Tkinter import EXTENDED, MULTIPLE, Button, Frame, Toplevel, Scrollbar, Listbox, END, GROOVE, LEFT

import tkFont


class SelectReturn(Toplevel):
    """General class to select one element from a list
    and return selection to parent class"""

    def __init__(self, app, root, elements=list(), select_title='select', Entry=None):         #Receives app and root from parent-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.elements = elements
        #Parent variable becomes selected element in list:
        self.parent_entry = Entry

        self.select_title = select_title

        self.dialog_box()
        self.update_list()

    def update_list(self):
        for line in sorted(self.elements):
            self.listbox.insert(END, line)

    def get_selected(self):
        """
        Get selected element and return to parent.
        """

        try:
            selected = map(int, self.listbox.curselection())
        except:
            return

        items_selected = []
        for item in selected:
            items_selected.append(self.listbox.get(item))

        print items_selected

        if 'ligand' in self.select_title:
            #LIE
            for item in items_selected:
                self.app.lig_unique_log.append(item)
        elif 'complex' in self.select_title:
            #LIE
            for item in items_selected:
                self.app.comp_unique_log.append(item)
        elif 'MD logfiles' in self.select_title:
            #ANALYSE TRAJECTORY ENERGIES
            for item in items_selected:
                self.app.unique_logs.append(item)
            self.app.get_energies(self.select_title.split()[0])
        #resFEP
        elif 'Select residue to mutate' in self.select_title:
            for item in items_selected:
                self.app.mutate_residue(item)

        self.destroy()

    def cancel(self):
        if 'ligand' in self.select_title:
            self.app.lig_unique_log.append('cancel')
        elif 'complex' in self.select_title:
            self.app.comp_unique_log.append('cancel')
        self.destroy()

    def dialog_box(self):
        """Defines the outlook of Setup MD window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('Select %s' % self.select_title)
        self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=10, padx=(10,0))

        listbox_scroll = Scrollbar(left_frame)
        listbox_scroll.grid(row = 0, rowspan = 10, column = 11, sticky = 'nsw')
        self.listbox = Listbox(left_frame, yscrollcommand = listbox_scroll.set, width = 30, height=10,
                               highlightthickness=0, relief=GROOVE, selectmode=EXTENDED)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 0, rowspan = 10, column = 0, columnspan = 10, sticky = 'w')
        self.listbox.config(font=tkFont.Font(family="Courier", size=12))

        select_button = Button(left_frame, text = 'Select', command = self.get_selected)
        select_button.grid(row = 11, column = 0, columnspan = 6, sticky = 'e')
        select_button.config(highlightbackground = self.main_color)

        cancel_button = Button(left_frame, text = 'Cancel', command = self.cancel)
        cancel_button.grid(row=11, column = 6, columnspan = 6, sticky = 'w' )
        cancel_button.config(highlightbackground = self.main_color)
