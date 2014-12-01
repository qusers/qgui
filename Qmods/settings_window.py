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

from Tkinter import X, TOP, Entry, Label, Button, Spinbox, Frame, Toplevel, Scrollbar, Listbox, Checkbutton, DISABLED, NORMAL, END, GROOVE, IntVar
import cPickle
from tkFileDialog import askopenfilename
from tkFileDialog import askdirectory
import os
import tkFont
from equilibration_edit import EditEq
from edit_file import FileEdit


class QguiSettings(Toplevel):
    """
    self.q_settings = [0' workdir',1[prm],2[lib],3[equilibration],4[use sub script (1/0), command], 5[executables],6'schrodinger path']
    """

    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.q_settings = cPickle.load(open(self.app.settings_path + '/Qsettings','rb'))
        self.submit_check = IntVar()
        self.dialog_window()
        self.updateSettings()
        self.toggle_submission()
        self.list_font = tkFont.Font(family="Monaco", size=12)

    def close(self):
        self.submitcommand.config(state=NORMAL)
        self.q_settings[4][1] = self.submitcommand.get().strip()
        self.q_settings[5][0] = self.qprep.get().strip()
        self.q_settings[5][1] = self.qdyn.get().strip()
        self.q_settings[5][2] = self.qfep.get().strip()
        self.q_settings[5][3] = self.qcalc.get().strip()
        self.schrodinger_path.config(state=NORMAL)
        self.q_settings[6] = self.schrodinger_path.get().strip()


        cPickle.dump(self.q_settings, open(self.app.settings_path + '/Qsettings','wb'))
        self.app.getSettings()
        self.destroy()

    def open_file(self):
        """
        opens the submit file for editing
        """

        if not os.path.isfile(self.app.settings_path + '/qsubmit'):
            newfile = open(self.app.settings_path + '/qsubmit','w')
            newfile.write('#!/bin/bash\n#\n#Qdyn I/O\n')
            newfile.close()

        self.fileEdit = FileEdit(self, self.app.settings_path + '/qsubmit')
        self.fileEdit.config(bg=self.main_color)
        self.fileEdit.title('Edit submission script')
        self.fileEdit.resizable()

    def edit_md_settings(self):
        if not os.path.isfile(self.app.settings_path + '/qmd_settings'):
            newfile = open(self.app.settings_path + '/qmd_settings','w')
            newfile.write('0.5   #simulation time\n1.0   #stepsize\n1     #inputfiles\n'
                          '300   #temperature\n10    #bath coupling\n'
                          '1     #shake solvent\n0     #shake solute\n0     #shake hydrogens\n1     #lrf\n'
                          '10    #solute-solute cutoff\n10    #solvent-solvent cutoff\n'
                          '10    #solute-solvent cutoff\n99    #Q-atoms cutoff\n99    #lrf cutoff\n'
                          '10.0  #shell force\n60.0  #radial force\n'
                          '1     #polarization restraint\n20.0  #polarization force\n25    #update non-bonded list\n'
                          '5     #update energy summary\n10    #write energy file\n100   #write trajectory\n'
                          'all   #trajectory atoms\n')
            newfile.close()

        self.fileEdit = FileEdit(self, self.app.settings_path + '/qmd_settings')
        self.fileEdit.config(bg=self.main_color)
        self.fileEdit.title('Edit default MD settings')
        self.fileEdit.resizable()

    def reset_md_settings(self):
        newfile = open(self.app.settings_path + '/qmd_settings','w')
        newfile.write('0.5   #simulation time\n1.0   #stepsize\n1     #inputfiles\n'
                        '300   #temperature\n10    #bath coupling\n'
                        '1     #shake solvent\n0     #shake solute\n0     #shake hydrogens\n1     #lrf\n'
                        '10    #solute-solute cutoff\n10    #solvent-solvent cutoff\n'
                        '10    #solute-solvent cutoff\n99    #Q-atoms cutoff\n99    #lrf cutoff\n'
                        '10.0  #shell force\n60.0  #radial force\n'
                        '1     #polarization restraint\n20.0  #polarization force\n25    #update non-bonded list\n'
                        '5     #update energy summary\n10    #write energy file\n100   #write trajectory\n'
                        'all   #trajectory atoms\n')
        newfile.close()

    def updateSettings(self):
        """
        Read settings from .settings and update lists
        self.qsettings = ['workdir',[prm],[lib]]
        """

        if self.q_settings[0] == 'default':
            self.workdir = os.path.dirname(os.path.abspath(__file__))
        else:
            self.workdir = self.q_settings[0]

        #Get parameter file(s) for force field and update prm list:
        self.prm_list.delete(0,END)
        for prmfile in self.q_settings[1]:
            self.prm_list.insert(END, '.../%s/%s' % (prmfile.split('/')[-2], prmfile.split('/')[-1]))

        #Get Libaray file(s) and update lib list:
        self.lib_list.delete(0,END)
        for libraryfile in self.q_settings[2]:
            self.lib_list.insert(END,'.../%s/%s' % (libraryfile.split('/')[-2], libraryfile.split('/')[-1]))

        #Get equilibration procedure and update list:
        self.eq_list.delete(0,END)
        self.eq_list.insert(END,'%3s %4s %4s %8s %4s %10s %8s' % ('# ',' T ','Bath', 'Restr.',' F ','st. size', 'Steps'))

        for i in range(len(self.q_settings[3])):
            temp, bath, atoms, force, stepsize, steps = self.q_settings[3][i][0:]
            self.eq_list.insert(END, '%3s %4s %4s %7s %6s %8s %10s' % ((i+1), temp, bath, atoms, force, stepsize, steps))

        self.submitcommand.config(state=NORMAL)
        self.submitcommand.delete(0, END)
        self.submitcommand.insert(0, self.q_settings[4][1])

        submit = int(self.q_settings[4][0])
        self.submit_check.set(submit)

        self.qprep.delete(0, END)
        self.qdyn.delete(0, END)
        self.qfep.delete(0, END)
        self.qcalc.delete(0, END)

        self.qprep.insert(0, self.q_settings[5][0])
        self.qdyn.insert(0, self.q_settings[5][1])
        self.qfep.insert(0, self.q_settings[5][2])
        self.qcalc.insert(0, self.q_settings[5][3])

        #Get schrodinger path if it exists:
        if self.q_settings[6]:
            path = self.q_settings[6]
        else:
            path = 'NA'
        self.schrodinger_path.config(state=NORMAL)
        self.schrodinger_path.delete(0,END)
        self.schrodinger_path.insert(0, path)
        #self.schrodinger_path.config(state=DISABLED)



    def toggle_submission(self):
        """
        Use submission script or not..
        """
        submit = self.submit_check.get()
        self.q_settings[4][0] = submit

        if submit == 0:
            self.edit_file.config(state=DISABLED)
            self.submitcommand.config(state=DISABLED)
        else:
            self.edit_file.config(state=NORMAL)
            self.submitcommand.config(state=NORMAL)


    def remove_from_prm(self):
        """
        Removes selected item in prm list from q_settings list and saves .q_settings file
        """
        if len(self.q_settings[1]) > 0:
            try:
                selected_index = int(self.prm_list.curselection()[0])
            except:
                return
            self.prm_list.delete(selected_index)
            del self.q_settings[1][selected_index]


    def remove_from_lib(self):
        """
        Removes selected item in prm list from q_settings list and saves Qsettings file
        """
        if len(self.q_settings[2]) > 0:
            try:
                selected_index = int(self.lib_list.curselection()[0])
            except:
                return
            self.lib_list.delete(selected_index)
            del self.q_settings[2][selected_index]



    def add_lib(self):
        """
        Opens file dialog (.lib) and selected lib is added to Qsettings.
        """
        libfile = askopenfilename(parent = self, initialdir = self.workdir,
                                  filetypes=(("Library", "*.lib"),("All files","*.*")))

        if libfile not in self.q_settings[2]:
            if len(libfile) > 4:
                self.q_settings[2].append(libfile)
        else:
            self.app.errorBox('Warning','A library with the same name already exist.')

        self.updateSettings()

    def add_prm(self):
        """
        Opens file dialog (.prm) and selected prm is added to Qsettings.
        """
        prmfile = askopenfilename(parent = self, initialdir = self.workdir,
                                  filetypes=(("Parameters", "*.prm"),("All files","*.*")))

        if prmfile not in self.q_settings[1]:
            if len(prmfile) > 4:
                self.q_settings[1].append(prmfile)
        else:
            self.app.errorBox('Warning','A parameter file with the same name already exist.')

        self.updateSettings()

    def edit_eq_button(self):
        """
        jjj
        """
        try:
            list_index = int(self.eq_list.curselection()[0])
        except:
            return
        if list_index != 0:
            nr, temp, bath, atoms, force, ss, steps = self.eq_list.get(list_index).split()
            self.edit_eq(nr, temp, bath, atoms, force, ss, steps)

    def add_eq_button(self):
        try:
            list_index = int(self.eq_list.curselection()[0])
        except:
            list_index = len(self.q_settings[3])
        if list_index != 0:
            nr = list_index + 1
            self.edit_eq(nr)

    def edit_eq(self, nr=1, temp='275', bath=1.0, atoms='All', force=10, ss=1.0, steps=10000):
        """
        opens dialog to edit/add equilibration steps
        """

        self.equilibration = EditEq(self,self.root, nr, temp, bath, atoms, force, ss, steps)
        self.equilibration.configure(bg = self.main_color)
        self.equilibration.title('Equilibration settings')
        self.equilibration.resizable()

    def remove_eq(self):
        """
        Removes selected equilibration step
        """

        if len(self.q_settings[3]) > 0:
            try:
                selected_index = int(self.eq_list.curselection()[0])
            except:
                return
            if selected_index != 0:
                self.eq_list.delete(selected_index)
                del self.q_settings[3][selected_index - 1]

    def select_path(self):
        """
        select schrodinger path
        """
        schrodinger = askdirectory(parent=self, mustexist=False, title='Set Schrodinger path', initialdir = self.app.workdir)
        if schrodinger != '':
            self.q_settings[6] = schrodinger
        else:
            if self.q_settings[6]:
                schrodinger = self.q_settings[6]
            else:
                schrodinger = 'NA'

        self.schrodinger_path.config(state=NORMAL)
        self.schrodinger_path.delete(0, END)
        self.schrodinger_path.insert(0, '%s' % schrodinger)
        #self.schrodinger_path.config(state=DISABLED)


    def dialog_window(self):
        """Defines the outlook of Setup MD window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('Settings')
        self.config(background=self.main_color)

        frame1 = Frame(self, bg = self.main_color)
        frame1.pack(side = TOP, pady=10, padx=(10,10), fill = X, expand = 1)

        frame2 = Frame(self, bg = self.main_color)
        frame2.pack(side = TOP, padx=(100,10), fill = X)


        prm_label = Label(frame1, text = 'Force Field Parameters')
        prm_label.grid(row = 0, column = 0, columnspan = 4)
        prm_label.configure(bg=self.main_color)

        prm_scroll = Scrollbar(frame1)
        prm_scroll.grid(row=1,column = 3, sticky = 'nsw')
        self.prm_list = Listbox(frame1, yscrollcommand = prm_scroll.set, width = 30, height = 3, highlightthickness = 0 , relief = GROOVE)
        prm_scroll.config(command = self.prm_list.yview)
        self.prm_list.grid(row = 1, column = 1, columnspan = 2, sticky = 'e')

        prm_add = Button(frame1, text = 'Add', command=self.add_prm)
        prm_add.grid(row = 2, column = 1)
        prm_add.config(highlightbackground = self.main_color)

        prm_delete = Button(frame1, text = 'Remove', command=self.remove_from_prm)
        prm_delete.grid(row = 2, column = 2)
        prm_delete.config(highlightbackground = self.main_color)

        lib_label = Label(frame1, text = 'Force Field Libraries')
        lib_label.grid(row = 3, column = 0, columnspan = 4, pady = (10,0))
        lib_label.configure(bg=self.main_color)

        lib_scroll = Scrollbar(frame1)
        lib_scroll.grid(row=4,column = 3, sticky = 'nsw')
        self.lib_list = Listbox(frame1, yscrollcommand = lib_scroll.set, width = 30, height = 3, highlightthickness = 0 , relief = GROOVE)
        lib_scroll.config(command = self.prm_list.yview)
        self.lib_list.grid(row = 4, column = 1, columnspan = 2, sticky = 'e')

        lib_add = Button(frame1, text = 'Add', command=self.add_lib)
        lib_add.grid(row = 5, column = 1)
        lib_add.config(highlightbackground = self.main_color)

        lib_delete = Button(frame1, text = 'Remove', command=self.remove_from_lib)
        lib_delete.grid(row = 5, column = 2)
        lib_delete.config(highlightbackground = self.main_color)

        #Equilibration procedure:
        eq_label = Label(frame1, text = 'Default equilibration procedure')
        eq_label.grid(row = 0, column = 5, columnspan = 4, padx = (10,0))
        eq_label.configure(bg=self.main_color)

        eq_scroll = Scrollbar(frame1)
        eq_scroll.grid(row=1, rowspan = 4, column = 9, sticky = 'nsw')
        self.eq_list = Listbox(frame1, yscrollcommand = eq_scroll.set, width = 50, height = 13, highlightthickness = 0 , relief = GROOVE)
        eq_scroll.config(command = self.eq_list.yview)
        self.eq_list.grid(row = 1, rowspan = 4, column = 5, columnspan = 4, padx = (10,0),sticky = 'e')
        self.eq_list.config(font = tkFont.Font(family="Courier", size=12))

        eq_add = Button(frame1, text = 'Add', command=self.add_eq_button)
        eq_add.grid(row = 5, column = 6)
        eq_add.config(highlightbackground = self.main_color)

        eq_edit = Button(frame1, text = 'Edit', command=self.edit_eq_button)
        eq_edit.grid(row = 5, column = 7)
        eq_edit.config(highlightbackground = self.main_color)

        eq_delete = Button(frame1, text = 'Remove', command=self.remove_eq)
        eq_delete.grid(row = 5, column = 8)
        eq_delete.config(highlightbackground = self.main_color)

        submit_label = Label(frame1, text = 'Use submission script:', bg=self.main_color)
        submit_label.grid(row = 6, column = 0, columnspan = 3, sticky='e')

        submit_check = Checkbutton(frame1, bg = self.main_color, variable=self.submit_check, command=self.toggle_submission)
        submit_check.grid(row=6, column=4, sticky='w')

        self.edit_file = Button(frame1, text='Edit', command=self.open_file)
        self.edit_file.grid(row=6, column=5)
        self.edit_file.config(highlightbackground=self.main_color)

        subcommand = Label(frame1, text='Command:', bg=self.main_color)
        subcommand.grid(row=6, column = 6)

        self.submitcommand = Spinbox(frame1, width =7, highlightthickness=0, relief=GROOVE, values=('./', 'qsub'))
        self.submitcommand.grid(row=6, column=7, columnspan=1)

        md_settings = Label(frame1, text='Default MD settings:', bg=self.main_color)
        md_settings.grid(row=7, column=0, columnspan=3, sticky='e')

        edit_md_settings = Button(frame1, text='Edit', highlightbackground=self.main_color, command=self.edit_md_settings)
        edit_md_settings.grid(row=7, column=5)

        reset_md_settings = Button(frame1, text='Reset', highlightbackground=self.main_color, command=self.reset_md_settings)
        reset_md_settings.grid(row=7, column=6)

        schrod_path = Label(frame1, text='Schrodinger path:', bg=self.main_color)
        schrod_path.grid(row=8, column=0, columnspan=3, sticky='e')

        self.schrodinger_path = Entry(frame1, width=40, highlightthickness=0, relief=GROOVE)
        self.schrodinger_path.grid(row=8, column=3, columnspan=5, sticky='w')

        select_path = Button(frame1, text='Select', highlightbackgroun=self.main_color, command = self.select_path)
        select_path.grid(row=8, column=8)

        executables = Label(frame1, text = 'Executables:', bg=self.main_color)
        executables.grid(row=9, column=0, columnspan=3, sticky='e')

        self.qprep = Spinbox(frame1, width=7, highlightthickness=0, relief=GROOVE, values='Qprep5')
        self.qprep.grid(row=9, column=3, columnspan=2)

        self.qdyn = Spinbox(frame1, width=7, highlightthickness=0, relief=GROOVE, values=('Qdyn5','Qdyn5p'))
        self.qdyn.grid(row=9, column=5, columnspan=2)

        self.qfep = Spinbox(frame1, width=7, highlightthickness=0, relief=GROOVE, values='Qfep5')
        self.qfep.grid(row=9, column=7, columnspan=1)

        self.qcalc = Spinbox(frame1, width=7, highlightthickness=0, relief=GROOVE, values='Qcalc5')
        self.qcalc.grid(row=9, column=8, columnspan=1)


        #Save/cancel
        save_button = Button(frame2, text = 'Save', command=self.close)
        save_button.grid(row=10,column=6, columnspan = 1, sticky = 'e')
        save_button.config(highlightbackground = self.main_color)

        cancel_button = Button(frame2, text = 'Close', command = self.destroy)
        cancel_button.grid(row=10, column = 7, columnspan = 1, sticky = 'e')
        cancel_button.config(highlightbackground = self.main_color)
