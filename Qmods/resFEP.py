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


        #Check if a topology is loaded in main window
        if self.app.top.id:
            self.topolgy_start = self.app.top_id

        #Trace stuff
        self.selected_topology.trace('w', self.topology_changed)

    def insert_topology_name(self):
        """
        Write topology name on screen
        """
        topology = self.topolgy_start
        if self.selected_topology.get() == 'Topology end':
            topology = self.topology_end


        topname = topology.split('/')[-1]
        self.topology_label.config(text=topname)

    def topology_changed(self, *args):
        """
        Update values in window when topology is toggeld
        """
        pass

    def load_topology(self):
        """
        :return:
        """
        filename = askopenfilename(parent=self, initialdir=self.app.workdir,
                                   filetypes=(("TOP", "*.top"), ("All files", '*.*')))
        if filename != '':
            #entry.delete(0,END)
            #entry.insert(0, filename.split('/')[-1])

            if self.selected_topology.get() == 'Topology start':
                self.topolgy_start = filename
            else:
                self.ligand_end = filename

            self.insert_topology_name()

    def dialog_window(self):

        self.title('resFEP setup')

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Topology frame
        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        #Dropdown menu for topology selection
        select_topology = OptionMenu(frame1, self.selected_topology, 'Topology start', 'Topology end')
        select_topology.config(bg=self.main_color, highlightbackground=self.main_color, width=15)
        select_topology.grid(row=0, column=0)

        #Load topology button
        load_topology = Button(frame1, text ='Load', highlightbackground=self.main_color, command=self.load_topolgy)
        load_topology.grid(row=0, column=1)

        #topology name
        self.topology_label = Label(frame1, text='*.top', bg=self.main_color)
        self.topology_label.grid(row=0, column=3)


