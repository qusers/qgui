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

from Tkinter import Text, Label, Button, Frame, Scrollbar, YES, Menu, DISABLED, NORMAL, GROOVE, END, TOP, PhotoImage,\
    HORIZONTAL, NONE
from tkFileDialog import askdirectory
import sys


TITLE = 'Q-gui'

class MainWindow(Frame):
    """Makes the mainwindow. Mainwindow inherits the Frame-class. """

    def __init__(self, app, root): #Receives app and root from Qgui class.
        Frame.__init__(self, root) # app = self from qgui and root = self.root from qgui
        self.root = root
        self.app = app #This is the reference to Qgui-class
        self.write = True
        
        self.main_color = self.app.main_color

        self.make_grid()
        self.menubar()

    def set_entryfield(self, pdb_id):
        """Deletes characters from entry_structure field when a new pdb_id 
        is retrived from pdb.org. Gets pdb_id from Qgui-class
        - send_pdb_id_to_entryfield() method - when a correct pdb_id is submitted to
        pdb.org. Inserts this pdb_id in to entry_structure.

        The deletion happens up to string length 10. """
        self.entry_structure.config(state = NORMAL)
        self.entry_structure.delete(0.0, END)
        self.entry_structure.insert(0.0, pdb_id)
        self.entry_structure.config(state = DISABLED)
        


    def set_topology_entryfield(self, top_id):
        """Deletes characters from topology_entry field when a new pdb_id 
        is retrived from pdb.org. Gets pdb_id from Qgui-class
        - send_pdb_id_to_entryfield() method - when a correct pdb_id is submitted to
        pdb.org. Inserts this pdb_id in to entry_topology field.
        """
        self.entry_topology.config(state = NORMAL)
        self.entry_topology.delete(0.0, END)
        self.entry_topology.insert(0.0, top_id)
        self.entry_topology.config(state = DISABLED)

    def select_workdir(self):
        """
        Opens a dialog where user can browse to desired working directory.
        """
        new_workdir = askdirectory(parent=self.root, mustexist=False, title='Change workdir',
                                   initialdir = self.app.workdir)

        if new_workdir == '':
            new_workdir = self.app.q_settings[ 'workdir' ]
        self.app.workdir = new_workdir
        self.app.q_settings[ 'workdir' ] = self.app.workdir
        self.app.saveSettings()

        print self.app.workdir

    def make_grid(self):
        """
        Makes to Frames, upper_frame and bottom_frame, and uses grid to organize
        widgets in to these two Frames."""

        self.root.title('Qgui') #TODO

        xmax = self.winfo_screenwidth()
        ymax = self.winfo_screenheight()

        x0 = self.x0 = xmax/2 - 589/2
        y0 = self.y0 = ymax/2 - 508/2

        self.root.config(bg=self.main_color)
        if 'darwin' in sys.platform:
            self.root.geometry("599x509+%d+%d" % (x0, y0))
            self.root.maxsize(589,508)
            self.root.minsize(589,508)
        else:
            self.root.geometry("599x509+%d+%d" % (x0, y0))
            self.root.maxsize(589,518)
            self.root.minsize(589,518)

        upper_frame = Frame(self, bg = self.main_color)
        upper_frame.pack(side = TOP)

        bottom_frame = Frame(self, bg = self.main_color)
        bottom_frame.pack(side = TOP)

        #Makes the labels for 'Structure' and 'Topology'

        label_structure = Label(upper_frame, text = 'Structure:')
        label_structure.grid(row = 0, column = 1, padx = 5, pady = 5)
        label_structure.config(background = self.main_color)
        label_topology = Label(upper_frame, text = 'Topology:')
        label_topology.grid(row = 1, column = 1)
        label_topology.config(background = self.main_color)

        #Makes entry-space for 'Structure' and 'Topology'.

        self.entry_structure = Text(upper_frame,width=20, height=1)
        self.entry_structure.grid(row = 0, column = 2)
        self.entry_structure.config(highlightthickness = 1, relief = GROOVE, highlightbackground='Blue')
        self.entry_structure.insert(0.0, '*.pdb') #pdb_id is set here from set_structure_entryfield()
        self.entry_structure.config(state = DISABLED)

        self.entry_topology = Text(upper_frame, widt=20, height=1)
        self.entry_topology.grid(row = 1, column = 2)
        self.entry_topology.config(highlightthickness = 1, relief = GROOVE, highlightbackground='Blue')
        self.entry_topology.insert(0.0, '*.top') #Get topology_id here!
        self.entry_topology.config(state = DISABLED)

        #Makes the 'Load' and 'View in pymol' buttons

        button_structure = Button(upper_frame, text = 'Load', command = self.app.load_dialog) #Here is a reference to Qgui-class, this enables to access load_structure()
        button_structure.grid(row = 0, column = 3)
        button_structure.config(highlightbackground = self.main_color, relief = GROOVE)
        button_view = Button(upper_frame, text = 'View in PyMOL', command=self.app.view_pymol)
        button_view.grid(row = 0, column = 4)
        button_view.config(highlightbackground = self.main_color, relief = GROOVE)
        button_topology = Button(upper_frame, text = 'Load', command = self.app.get_top_file)
        button_topology.grid(row = 1, column = 3)
        button_topology.config(highlightbackground = self.main_color, relief = GROOVE)

        button_clear = Button(upper_frame, text = 'Clear', command = self.clear_button_pressed)
        button_clear.grid(row = 0, column = 5)
        button_clear.config(highlightbackground = self.main_color, relief = GROOVE)

        photo = PhotoImage(file=self.app.qgui_path + "/Qmods/logo_qgui.gif")
        w = Label(upper_frame, image=photo, bg=self.main_color)
        w.photo = photo
        w.grid(row=0, rowspan=2, column=0)


        button_topo_clear = Button(upper_frame, text = 'Clear', command = self.clear_button_topo_pressed)
        button_topo_clear.grid(row = 1, column = 5)
        button_topo_clear.config(highlightbackground = self.main_color, relief = GROOVE)

        #Makes the text-window where the logfile can be displayed
        scrollb = Scrollbar(bottom_frame)
        xscroll = Scrollbar(bottom_frame, orient=HORIZONTAL)
        scrollb.grid(row = 3, column = 4, sticky = 'ns', padx = (0,5), pady = 5)
        xscroll.grid(row=4, column=0, columnspan=4, sticky='we', padx=(5,0))
        self.txt = Text(bottom_frame, xscrollcommand=xscroll.set, yscrollcommand=scrollb.set, wrap=NONE)
        #self.txt.config(undo= True)
        self.txt.grid(row = 3, columnspan = 4, sticky="nsew", padx=(5, 0), pady=(5,0))
        self.txt.config(highlightthickness=1, relief=GROOVE, highlightbackground='Blue')
        scrollb.config(command=self.txt.yview)
        xscroll.config(comman=self.txt.xview)

        close_button = Button(bottom_frame, text = 'Quit', command = self.app.exit)
        close_button.grid(row = 5, column = 3, sticky = 'e')
        close_button.config(highlightbackground = self.main_color, relief = GROOVE)

        self.pack(expand=YES)

    def clear_button_pressed(self):
        """Deletes any entry in entry_structure-field. Enables to download a new pdb file.
        Calls Qgui controllers reset_pdb() method. """
        self.entry_structure.config(state = NORMAL)
        self.entry_structure.delete(0.0, END)
        self.app.reset_pdb(None)
        self.entry_structure.config(state = DISABLED)

    def clear_button_topo_pressed(self):
        """Deletes any entry in entry_topology-field. Enables to download a new pdb file.
        Calls Qgui controllers reset_pdb() method. """
        self.entry_topology.config(state = NORMAL)
        self.entry_topology.delete(0.0, END)
        self.app.reset_topology(None)
        self.entry_topology.config(state = DISABLED)

    def update_txt(self,msg):
        """Receives msg from log() method from Qgui-class.
        Makes the txt-field inmutable. """
        self.txt.config(state = NORMAL)
        self.txt.insert(END,msg)
        self.txt.yview(END)
        self.txt.config(state = DISABLED)

    def menubar(self):
        """Creates a menubar on top of the screen."""

        menubar = Menu(self.root)
        self.root.config(menu = menubar)

        #Defines the File-menu.
        convert_menu = Menu(menubar, tearoff=0)
        filemenu = Menu(menubar, tearoff = 0)
        filemenu.add_command(label = 'Change workdir', command = self.select_workdir)
        filemenu.add_command(label = 'Import structure', command = self.app.load_dialog)
        filemenu.add_command(label = 'Import Topology', command = self.app.get_top_file)
        filemenu.add_cascade(label='Convert', menu=convert_menu, underline=0)
        filemenu.add_command(label = 'Settings', command = self.app.load_settings)
        filemenu.add_separator()
        filemenu.add_command(label = 'Close', command = self.app.exit) #Here we use the reference to Qgui-class
        menubar.add_cascade(label = 'File', menu = filemenu)

        #Convert files submenu:
        convert_menu.add_command(label=(u'*.top \u2192 *.pdb'), command=self.app.convert_top_pdb)

        #Defines the Prepare-menu.

        preparemenu = Menu(menubar, tearoff = 0)
        preparemenu.add_command(label = 'PDB', command = self.app.prepare_pdb)
        preparemenu.add_command(label = 'Topology', command = self.app.prepare_topo)
        preparemenu.add_command(label = 'Parameters', command = self.app.create_parameters)
        preparemenu.add_command(label = 'Trjmask', command = self.app.open_trjmask)
        menubar.add_cascade(label = 'Prepare', menu = preparemenu)

        #Defines the Setup-menu.

        setupmenu = Menu(menubar, tearoff = 0)
        setupmenu.add_command(label = 'MD', command = self.app.setup_md)
        setupmenu.add_command(label = 'LIE', command = self.app.setup_lie)
        setupmenu.add_command(label = 'FEP', command = self.app.setup_fep)
        setupmenu.add_command(label = 'resFEP', command=self.app.setup_resfep)
        setupmenu.add_command(label = 'EVB', command = self.app.setup_evb)
        menubar.add_cascade(label = 'Setup', menu = setupmenu)

        #Defines the analyze-menu.
        subanalyze_trajectory = Menu(menubar, tearoff = 0)
        analyzemenu = Menu(menubar, tearoff = 0)
        analyzemenu.add_cascade(label='Trajectory', menu = subanalyze_trajectory, underline = 0)
        analyzemenu.add_command(label='LIE', command = self.app.analyze_lie)
        analyzemenu.add_command(label='FEP', command = self.app.analyze_fep)
        analyzemenu.add_command(label='resFEP', command=self.app.analyze_resFEP)
        menubar.add_cascade(label='Analyze', menu=analyzemenu)

        #Defines the submenu of Trajectory in Analyze-menu.

        #subanalyze_trajectory = Menu(menubar, tearoff = 0)
        subanalyze_trajectory.add_command(label = 'Energies', command = self.app.analyze_energies)
        subanalyze_trajectory.add_command(label = 'RMSF', command = lambda: self.app.analyze_qcalc('RMSF'))
        subanalyze_trajectory.add_command(label = 'RMSD', command = lambda: self.app.analyze_qcalc('RMSD'))
        subanalyze_trajectory.add_command(label = 'Distance', command = lambda: self.app.analyze_qcalc('Distance'))
        subanalyze_trajectory.add_command(label = 'Angle', command = lambda: self.app.analyze_qcalc('Angle'))
        subanalyze_trajectory.add_command(label = 'Torsion', command = lambda: self.app.analyze_qcalc('Torsion'))
        #subanalyze_trajectory.add_command(label = 'Entropy', command = self.button_command)
        #Defines the submenu of EVB in Analyze-menu.

        subanalyze_EVB = Menu(analyzemenu)
        subanalyze_EVB.add_command(label = 'Reference reaction', command = self.app.evb_calibration)
        subanalyze_EVB.add_command(label = 'Reaction energies', command = self.app.analyze_evb_reactions)
        subanalyze_EVB.add_command(label = 'Thermodynamic parameters', command = self.app.evb_arrhenius)
        analyzemenu.add_cascade(label = 'EVB', menu = subanalyze_EVB, underline = 0)

        #Defines the Help-menu.

        helpmenu = Menu(menubar, tearoff = 0)
        helpmenu.add_command(label = 'User manual', command = self.button_command)
        helpmenu.add_command(label = 'About', command = self.app.aboutQ)
        menubar.add_cascade(label = 'Help', menu = helpmenu)


    def button_command(self):
        pass
