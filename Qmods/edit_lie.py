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

from Tkinter import Entry, Spinbox, Label, Button,Frame, Toplevel, END, GROOVE, DISABLED, NORMAL, StringVar


class EditLIE(Toplevel):
    """Class to edit LIE data selected from list """

    def __init__(self, app, root, edit=True, data_to_insert=[], sel_index='END'):
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.insert_index = sel_index
        self.edit = edit
        #Edit or add data?
        if self.edit:
            self.lie_title = 'Edit'
        else:
            self.lie_title = 'Add'

        #Trace changes in parameters to update table/plot
        self.d_el_var = StringVar()
        self.d_vdw_var = StringVar()
        self.alpha_var = StringVar()
        self.beta_var = StringVar()
        self.gamma_var = StringVar()


        self.dialog_box()

        self.gamma_var.trace('w', self.parameters_changed)
        self.beta_var.trace('w', self.parameters_changed)
        self.alpha_var.trace('w', self.parameters_changed)
        self.d_vdw_var.trace('w', self.parameters_changed)
        self.d_el_var.trace('w', self.parameters_changed)
        #Insert values:
        for entry in self.entries:
            entry.delete(0, END)

        for i in range(len(data_to_insert)):
            self.entries[i].insert(0, data_to_insert[i])

        #If edit, disable entries:
        if self.edit:
            for entry in self.edit_disabled:
                entry.config(state=DISABLED)

    def parameters_changed(self, *args):
        """
        Update dG when parameters are toggled:
        """
        self.dG.config(state=NORMAL)
        self.d_el.config(state=NORMAL)
        self.d_vdw.config(state=NORMAL)
        #Get data from entries:
        lie = []
        for i in range(1, 6):
            try:
                lie.append(float(self.entries[i].get()))
            except:

                return

        dG = lie[0]*lie[1] + (lie[2]*lie[3]) + lie[4]

        self.dG.delete(0, END)
        self.dG.insert(0, '%8.2f' % dG)

        if self.edit:
            self.dG.config(state=DISABLED)
            self.d_el.config(state=DISABLED)
            self.d_vdw.config(state=DISABLED)

    def save_lie(self):
        """
        Replaces or adds new date to Fit Lie list (parent)
        """
        #Enable entries:
        for entry in self.edit_disabled:
            entry.config(state=NORMAL)

        #Get data:
        lie_data = []
        i = 0
        for entry in self.entries:
            if i > 0:
                lie_data.append(float(entry.get()))
            else:
                lie_data.append('#%s' % entry.get().split('#')[-1])
            i += 1

        lie_data.append(0.0)

        print lie_data
        if self.edit:
            self.app.lie_data[self.insert_index] = lie_data
        else:
            self.app.lie_data.append(lie_data)

        print self.app.lie_data

        self.app.update_table()
        self.app.update_plot()
        self.destroy()


    def dialog_box(self):
        """Defines the outlook of add/edit equilibration window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('%s LIE data' % self.lie_title)
        self.config(background=self.main_color)

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        frame2 = Frame(mainframe, bg=self.main_color)
        frame2.grid(row=1, column=0)

        name = Label(frame1, text='#Name     ', bg=self.main_color)
        name.grid(row=0, column=0)

        beta = Label(frame1, text=u"\N{GREEK SMALL LETTER BETA}", bg=self.main_color)
        beta.grid(row=0, column=1)

        d_el = Label(frame1, text=u"\N{GREEK CAPITAL LETTER DELTA}E(el)", bg=self.main_color)
        d_el.grid(row=0, column=2)

        alpha = Label(frame1, text=u"\N{GREEK SMALL LETTER ALPHA}", bg=self.main_color)
        alpha.grid(row=0, column=3)

        d_vdw = Label(frame1, text=u"\N{GREEK CAPITAL LETTER DELTA}E(vdW)", bg=self.main_color)
        d_vdw.grid(row=0, column=4)

        gamma = Label(frame1, text=u"\N{GREEK SMALL LETTER GAMMA}", bg=self.main_color)
        gamma.grid(row=0, column=5)

        dG = Label(frame1, text=u"\N{GREEK CAPITAL LETTER DELTA}G", bg=self.main_color)
        dG.grid(row=0, column=6)

        dG_exp = Label(frame1, text=u"\N{GREEK CAPITAL LETTER DELTA}G(exp)", bg=self.main_color)
        dG_exp.grid(row=0, column=7)

        self.name = Entry(frame1, width=10, highlightthickness=0)
        self.name.grid(row=1, column=0)

        self.beta = Spinbox(frame1, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-1.00, to=1.00, increment=0.01, textvariable=self.beta_var)
        self.beta.grid(row=1, column=1)

        self.d_el = Entry(frame1, width=7, highlightthickness=0, textvariable=self.d_el_var)
        self.d_el.grid(row=1, column=2)

        self.alpha = Spinbox(frame1, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-1.00, to=1.00, increment=0.01, textvariable=self.alpha_var)
        self.alpha.grid(row=1, column=3)

        self.d_vdw = Entry(frame1, width=7, highlightthickness=0, textvariable=self.d_vdw_var)
        self.d_vdw.grid(row=1, column=4)

        self.gamma = Spinbox(frame1, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-99.99, to=99.99, increment=0.5, textvariable=self.gamma_var)
        self.gamma.grid(row=1, column=5)

        self.dG = Entry(frame1, width=7, highlightthickness=0)
        self.dG.grid(row=1, column=6)

        self.dG_exp = Entry(frame1, width=7, highlightthickness=0)
        self.dG_exp.grid(row=1, column=7)

        save_button = Button(frame2, text='Save', highlightbackground=self.main_color, command=self.save_lie)
        save_button.grid(row=0, column=0)

        close_button = Button(frame2, text='Close', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=0, column=1)

        self.entries = [self.name, self.beta, self.d_el, self.alpha, self.d_vdw, self.gamma, self.dG, self.dG_exp]
        self.edit_disabled = [self.d_el, self.d_vdw, self.dG]
