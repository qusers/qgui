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
    OptionMenu, StringVar, Spinbox, SINGLE, LabelFrame, MULTIPLE
from tkFileDialog import askopenfilename
import qgui_functions as qf
from select_return import SelectReturn
from subprocess import Popen, PIPE
import tkFont
import time
import os
import signal
import sys
import numpy as np
import shutil
from copy import deepcopy


class ResFEP(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root

        self.topology_start = None
        self.topology_end = None


        self.selected_topology = StringVar()
        self.selected_topology.set('Topology start')

        self.mutate_to = StringVar()
        self.mutate_from = StringVar()

        self.mutate_from.set('start')
        self.mutate_to.set('end')



        #Check if a topology is loaded in main window
        if self.app.top_id:
            self.topolgy_start = self.app.top_id

        #Trace stuff
        self.selected_topology.trace('w', self.topology_changed)

        self.dialog_window()

    def insert_topology_name(self):
        """
        Write topology name on screen
        """
        topname = '*.top'

        if self.selected_topology.get() == 'Topology end':
            topology = self.topology_end
        else:
            topology = self.topology_start

        if topology:
            topname = topology.split('/')[-1]

        self.topology_label.config(text=topname)

    def topology_changed(self, *args):
        """
        Update values in window when topology is toggeld
        """
        self.insert_topology_name()

    def load_topology(self):
        """
        :return:
        """
        filename = askopenfilename(parent=self, initialdir=self.app.workdir,
                                   filetypes=(("TOP", "*.top"), ("All files", '*.*')))
        if filename != '':
            if self.selected_topology.get() == 'Topology start':
                self.topology_start = filename
            else:
                self.topology_end = filename

            self.insert_topology_name()

    def add_residue(self):
        """
        Opens a new window to select residue to mutate

        """

        self.select_res = SelectReturn(self, self.root, elements=list(), select_title='select residue',
                                       Entry=self.reslist)
        self.select_res.configure(bg=self.main_color)
        self.select_res.title('Select residue to mutate')
        self.select_res.resizable()

    def mutate_residue(self, residue):
        """
        :param residue: RES 123
        Checks if RES exist in FEP library and loads FEP files. Abort if not existing!
        """
        resname = residue.split()[0]
        #TODO

    def del_residue(self):
        """
        Delete selected residue from mutation list
        :return:
        """
        pass

    def update_residue_mutation(self):
        """
        Update mutation from --> to in listbox
        :return:
        """
        self.reslist.insert(END,'TYR193 > ALA193')

    def view_fep(self):
        """

        :return:
        """
        pass

    def duplicate_fepfile(self):
        """

        :return:
        """
        pass

    def add_fepfile(self):
        """

        :return:
        """
        pass

    def delete_fepfile(self):
        """

        :return:
        """
        pass

    def move_fep(self, i=0):
        """
        Moves FEP file in listbox from position j to j + i
        :param value:
        :return:
        """
        pass

    def configure_md(self):
        """
        Opens the MD settings module!
        :return:
        """
        pass

    def write_inputfiles(self):
        """

        :return:
        """
        pass

    def run_fep(self):
        """

        :return:
        """
        pass

    def reslist_event(self, *args):
        """
        :param args:
        :return:
        """
        pass


    def dialog_window(self):

        self.title('resFEP setup')

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Topology frame
        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        #Frame with residues to mutate and FEP files
        frame2 = Frame(mainframe, bg=self.main_color)
        frame2.grid(row=1, column=0)

        #Frame with configure MD /write/save/close etc.
        frame3 = Frame(mainframe, bg=self.main_color)
        frame3.grid(row=2, column=0)

        #Dropdown menu for topology selection
        select_topology = OptionMenu(frame1, self.selected_topology, 'Topology start', 'Topology end')
        select_topology.config(bg=self.main_color, highlightbackground=self.main_color, width=20)
        select_topology.grid(row=1, column=0)

        #Load topology button
        load_topology = Button(frame1, text ='Load', highlightbackground=self.main_color, command=self.load_topology)
        load_topology.grid(row=1, column=1)

        #topology name
        self.topology_label = Label(frame1, text='*.top', bg=self.main_color)
        self.topology_label.grid(row=0, column=0, columnspan=2)

        #Residues label
        res_label = Label(frame2, text='Mutate residue(s):', bg=self.main_color)
        res_label.grid(row=0, column=0, columnspan=3, pady=(10,0), sticky='w')

        #Listbox with selected residues to mutate:
        reslist_scroll = Scrollbar(frame2)
        reslist_scroll.grid(row=1, rowspan=3, column=1, sticky='nsw')
        self.reslist = Listbox(frame2, yscrollcommand=reslist_scroll.set, width=17, height=4,
                                 highlightthickness=0, relief=GROOVE, selectmode=SINGLE, exportselection=False)
        reslist_scroll.config(command=self.reslist.yview)
        self.reslist.grid(row=1, rowspan=3, column=0, sticky='nse')
        self.reslist.config(font=tkFont.Font(family="Courier", size=12))
        self.reslist.bind('<<ListboxSelect>>', self.reslist_event)

        #Add residue button
        add_res = Button(frame2, text='Add', highlightbackground=self.main_color, command=self.add_residue)
        add_res.grid(row=1, column=2)

        #Delete residue button
        del_res = Button(frame2, text='Delete', highlightbackground=self.main_color, command=self.del_residue)
        del_res.grid(row=1, column=4)


        #Mutate from dropdown menu:
        self.mutate_from_menu = OptionMenu(frame2, self.mutate_from, 'start')
        self.mutate_from_menu.config(bg=self.main_color, highlightbackground=self.main_color, width=7)
        self.mutate_from_menu.grid(row=2, column=2)

        #right arrow:
        arrow_label = Label(frame2, text='%1s' % u'\u2192', bg=self.main_color)
        arrow_label.grid(row=2, column=3)

        #Mutate to dropdown menu
        self.mutate_to_menu = OptionMenu(frame2, self.mutate_to, 'end')
        self.mutate_to_menu.config(bg=self.main_color, highlightbackground=self.main_color, width=7)
        self.mutate_to_menu.grid(row=2, column=4)

        #Update button
        update_res_mutation = Button(frame2, text='Update', highlightbackground=self.main_color,
                                     command=self.update_residue_mutation)
        update_res_mutation.grid(row=3, column=2, columnspan=3)

        #FEP files label
        feps_label = Label(frame2, text='FEP file(s):', bg=self.main_color)
        feps_label.grid(row=4, column=0, columnspan=3, pady=(10,0), sticky='w')

        #Listbox for FEP files
        feplist_scroll = Scrollbar(frame2)
        feplist_scroll.grid(row=5, rowspan=5, column=1, sticky='nsw')
        self.feplist = Listbox(frame2, yscrollcommand=feplist_scroll.set, width=17, height=8,
                                 highlightthickness=0, relief=GROOVE, selectmode=SINGLE, exportselection=False)
        feplist_scroll.config(command=self.feplist.yview)
        self.feplist.grid(row=5, rowspan=5, column=0, sticky='nse')
        self.feplist.config(font=tkFont.Font(family="Courier", size=12))
        #self.feplist.bind('<<ListboxSelect>>', self.reslist_event)

        #View FEP file button
        view_fep = Button(frame2, text='View', highlightbackground=self.main_color, command=self.view_fep)
        view_fep.grid(row=5, column=2)

        #Duplicate selected FEP file:
        duplicate_fepfile = Button(frame2, text='Duplicate', highlightbackground=self.main_color,
                                   command=self.duplicate_fepfile)
        duplicate_fepfile.grid(row=5, column=3, columnspan=2)

        #Add a new/blank FEP file
        add_fepfile = Button(frame2, text='Add', highlightbackground=self.main_color, command=self.add_fepfile)
        add_fepfile.grid(row=6, column=2)

        #Delete selected FEP file
        del_fepfile = Button(frame2, text='Delete', highlightbackground=self.main_color, command=self.delete_fepfile)
        del_fepfile.grid(row=6, column=3, columnspan=2)

        #Move FEP file up
        move_up = Button(frame2, text=u'\u2191', highlightbackground=self.main_color, command=lambda: self.move_fep(-1))
        move_up.grid(row=7, column=2, sticky='e')

        #Move FEP file down
        move_down = Button(frame2, text=u'\u2193', highlightbackground=self.main_color,
                           command=lambda: self.move_fep(1))
        move_down.grid(row=7, column=3, columnspan=2, sticky='w')

        #Configure MD
        config_md = Button(frame3, text='Configure MD', highlightbackground=self.main_color, command=self.configure_md)
        config_md.grid(row=0, column=0, columnspan=3)

        #Write inputfiles
        write_inp = Button(frame3, text='Write', highlightbackground=self.main_color, command=self.write_inputfiles)
        write_inp.grid(row=1, column=0)

        #Run FEP calculation
        run_fep = Button(frame3, text='Run', highlightbackground=self.main_color, command=self.run_fep)
        run_fep.grid(row=1, column=1)

        #Close resFEP gui
        close_resfep = Button(frame3, text='close', highlightbackground=self.main_color, command=self.destroy)
        close_resfep.grid(row=1, column=2)

