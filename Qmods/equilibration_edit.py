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

from Tkinter import Entry, Spinbox, Label, Button,Frame, Toplevel, END, GROOVE, LEFT


class EditEq(Toplevel):
    """jjj """

    def __init__(self, app, root, nr=1, temp='275', bath=1.0, atoms='All', force=10, ss=1.0, steps=10000):
        """
        :param app:
        :param root:
        :param nr:
        :param temp:
        :param bath:
        :param atoms:
        :param force:
        :param ss:
        :param steps:
        """
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.dialog_box()
        self.nr = nr
        self.temp = temp
        self. bath = bath
        self.atoms = atoms
        self.force = force
        self.ss = ss
        self.steps = steps

        self.nr_entry.insert(0, self.nr)
        #self.nr_entry.config(state=DISABLED)
        self.temp_entry.insert(0,self.temp)
        self.bath_entry.insert(0,self.bath)
        self.restrain_entry.delete(0, END)
        self.restrain_entry.insert(0,self.atoms)
        self.force_entry.insert(0, self.force)
        self.ss_entry.insert(0, self.ss)
        self.steps_entry.insert(0,self.steps)

    def add_eq(self):
        try:
            eq_index = int(self.nr_entry.get()) - 1
        except:
            return
        eq_list = [self.temp_entry.get().strip(),
                   self.bath_entry.get().strip(),
                   self.restrain_entry.get().strip(),
                   self.force_entry.get().strip(),
                   self.ss_entry.get().strip(),
                   self.steps_entry.get().strip()]

        if int(self.nr_entry.get()) <= len(self.app.q_settings[ 'equilibration' ]):
            del self.app.q_settings[ 'equilibration' ][eq_index]
            self.app.q_settings[ 'equilibration' ].insert(eq_index, eq_list)
        else:
            self.app.q_settings[ 'equilibration' ].append(eq_list)

        self.app.updateSettings()
        self.destroy()


    def dialog_box(self):
        """Defines the outlook of add/edit equilibration window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('Add/edit equilibration procedure')
        self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=10, padx=(5,5))

        nr_label = Label(left_frame, text = ' # ')
        nr_label.grid(row = 0, column = 1)
        nr_label.configure(bg=self.main_color)

        self.nr_entry = Entry(left_frame, width = 3, highlightthickness = 0, relief = GROOVE)
        self.nr_entry.grid(row = 1, column = 1, padx=(5,5))

        temp_label = Label(left_frame, text = ' T ')
        temp_label.grid(row = 0, column = 2)
        temp_label.configure(bg=self.main_color)

        self.temp_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.temp_entry.grid(row = 1, column = 2, padx=(5,5))

        bath_label = Label(left_frame, text = ' Bath ')
        bath_label.grid(row = 0, column = 3)
        bath_label.configure(bg=self.main_color)

        self.bath_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.bath_entry.grid(row = 1, column = 3, padx=(5,5))

        restr_label = Label(left_frame, text = 'Restrain')
        restr_label.grid(row = 0, column = 4)
        restr_label.configure(bg=self.main_color)

        self.restrain_entry = Spinbox(left_frame, width = 7, highlightthickness = 0, relief = GROOVE,
                                      values =('All','Solute','Solvent','None'))
        self.restrain_entry.grid(row = 1, column = 4, padx=(5,5))

        force_label = Label(left_frame, text = 'Force')
        force_label.grid(row = 0, column = 5)
        force_label.configure(bg=self.main_color)

        self.force_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.force_entry.grid(row = 1, column = 5, padx=(5,5))

        ss_label = Label(left_frame, text = 'Stepsize')
        ss_label.grid(row = 0, column = 6)
        ss_label.configure(bg=self.main_color)

        self.ss_entry = Entry(left_frame, width = 6, highlightthickness = 0, relief = GROOVE)
        self.ss_entry.grid(row = 1, column = 6, padx=(5,5))

        steps_label = Label(left_frame, text = 'Steps')
        steps_label.grid(row = 0, column = 7)
        steps_label.configure(bg=self.main_color)

        self.steps_entry = Entry(left_frame, width = 12, highlightthickness = 0, relief = GROOVE)
        self.steps_entry.grid(row = 1, column = 7, padx=(5,5))

        add_button = Button(left_frame, text = 'Add', command = self.add_eq)
        add_button.grid(row = 2, column = 6, pady = (15,0))
        add_button.config(highlightbackground = self.main_color)

        cancel_button = Button(left_frame, text = 'Cancel', command= self.destroy)
        cancel_button.grid(row = 2, column = 7, sticky = 'w', pady=(15,0))
        cancel_button.config(highlightbackground = self.main_color)




