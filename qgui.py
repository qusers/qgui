#!/usr/bin/env python3

# -*- coding: utf-8 -*-
__author__ = "Geir Villy Isaksen"
__copyright__ = "Copyright (C) 2017 University of Tromso / Geir Villy Isaksen"
__credits__ = ["Bjorn Olav Brandsdal", "Tor Arne Heim Andberg", "Laura Liikanen", "Johan Aqvist", "Christoffer Lind",
               "Masoud Kazemi", "Jaka Socan"]
__license__ = "GPL v3"
__version__ = "2.00"
__maintainer__ = "Geir Isaksen"
__email__ = "geir.isaksen@uit.no"
__status__ = "Production"

# Copyright 2016 Geir Villy Isaksen/University of Tromso, all rights reserved

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

from tkinter import Tk, Entry, Label, Button, Toplevel, GROOVE, messagebox
from tkinter.filedialog import askopenfilename, asksaveasfilename
from urllib.request import urlopen
import os
import sys
import datetime
import getpass

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/Qmods')

from splash_qgui import SplashScreen
from pdbpreparewindow import PdbFilePrepare
from mainwindow import MainWindow
from topoprep import TopologyPrepare
from trjmask import TrjMask
from setup_md import SetupMd
from setup_lie import SetupLie
from settings_window import QguiSettings
from prmlib import CreatePrmLib
from analyze_LIE import AnalyzeLie
from setup_evb import SetupEVB
from setup_FEP import SetupFEP
from resFEP import ResFEP
from evb_rrc import EvbCalibration
from analyze_evb_re import EvbReactions
from aboutQ import AboutQ
from qpymol import ViewPyMol
from evb_arrhenius import EvbArrhenius
from analyze_FEP import AnalyzeFEP
from analyse_resFEP import Analyse_resFEP
from analyze_energies import AnalyzeEnergies
from qcalc import AnalyzeQcalc
import prepareTopology as pt
import qgui_functions as qf
import pickle
import time

class QGui(object):
    """
    This is the main Qgui controller-class.

     """

    def __init__(self, root): #receives root from main-method
        self.root = root
        #Path to user settings
        self.settings_path = os.path.expanduser('~') + '/.qgui'

        #Check if settings directory exists, if not create it
        if not os.path.exists(self.settings_path):
            os.makedirs(self.settings_path)

        #path to Qgui code
        self.qgui_path = os.path.dirname(os.path.realpath(__file__))
        self.root.withdraw()
        self.show_splash()

        #or 'LightCyan', 'Ivory', 'Gainsboro', 'AliceBlue' ?
        self.main_color = 'AliceBlue'

        self.main_window = MainWindow(self, self.root) #sends self and self.root to MainWindow to make a reference
        self.main_window.config(bg=self.main_color)

        self.pymol_running = False

        #On mac the following is needed to make window pop to front when opened:
        #if 'darwin' in sys.platform:
        #    try:
        #        os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python3" to true' ''')
        #    except:
        #        os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

        self.structure = None
        self.top_id = None
        self.pdb_id = None
        self.log('info', 'Welcome to Qgui v2.0')
        self.log('info', 'The Free Energy package for Q.')
        self.log(' ', '\n')

        self.getSettings()

        self.checkUpdates()

        #Check for workdir key upon launce (-p)
        try:
            this_dir = sys.argv[1]
            if this_dir == '-p':
                new_workdir = os.getcwd()
                print('\nWorkdir: %s\n' % new_workdir)
                self.q_settings[ 'workdir' ] = new_workdir
                self.workdir = new_workdir
                self.saveSettings()

        except:
            print('\nUsing workdir defined in settings')
            print('TIPS: Launch Qgui -p to automatically set current directory as workdir.\n')

    def show_splash(self):
        SplashScreen(self, self.root, self.qgui_path + "/Qmods/qgui_box.gif")

        time.sleep(2)
        self.root.deiconify()

    def checkUpdates(self):
        if os.path.isfile(self.settings_path + '/qupdate'):
            update_msg = open(self.settings_path + '/qupdate').readlines()
            update_msg = '\n'.join(update_msg)
            self.errorBox('Info', update_msg)

            os.system('rm %s' % self.settings_path + '/qupdate')

    def set_default_qsettings(self):
        """
        Sets q_settings to default
        :return: q_settings dictionary
        """
        q_settings = {'workdir': 'default',
                      'parameter': [],
                      'library': [],
                      'equilibration': [
                          [1, 0.2,'All',10.0, 0.1, 10000 ],
                          [50, 1.0,'All',10.0, 1.0, 10000 ],
                          [150, 1.0,'All',5.0, 1.0, 10000 ],
                          [275, 1.0,'All',5.0, 1.0, 10000 ],
                          ['End', 10.0,'None',0, 1.0, 100000]
                      ],
                      'subscript': [1, './'],
                      'executables': ['qprep','qdyn','qfep','qcalc'],
                      'schrodinger path':  None}

        return q_settings


    def getSettings(self):
        """
        Reads settings and returns workdir, [prm] and [lib]

        self.q_settings = dict()
        """
        try:
            self.q_settings = pickle.load(open(self.settings_path + '/Qsettings','rb'))

            if type(self.q_settings) == list:
                Keys = [ 'workdir', 'parameter', 'library', 'equilibration', 'subscript', 'executables', 'schrodinger path' ]
                new_settings = {}

                map(lambda K: new_settings.update({Keys[K[0]]: K[1]}), enumerate(self.q_settings))

                self.q_settings = new_settings

                print('Converted old Q-settings format to new format.')
                print(self.q_settings)
                print('Please see too that you settings are they way you set it up to be.')

                pickle.dump(self.q_settings, open(self.settings_path + '/Qsettings','wb'))

        except:

            self.q_settings = self.set_default_qsettings()

            pickle.dump(self.q_settings, open(self.settings_path + '/Qsettings','wb'))

        #Get workdirectory:
        if self.q_settings[ 'workdir' ] == 'default':
            #self.workdir = os.path.dirname(os.path.realpath(__file__))
            if 'darwin' in sys.platform:
                self.workdir = '/Users/' + getpass.getuser()
            else:
                self.workdir = '/home/' + getpass.getuser()
        else:
            self.workdir = self.q_settings[ 'workdir' ]

        #Get parameter file(s) for force field:
        self.prms = self.q_settings[ 'parameter' ]

        #Get Libaray file(s):
        self.libs = self.q_settings[ 'library' ]

    def saveSettings(self):
        """
        saves q_settings to file Qsettings
        """
        pickle.dump(self.q_settings, open(self.settings_path + '/Qsettings','wb'))

    def log(self, useDate = 'nope', msg='text'):
        if useDate == 'info':
            now = datetime.datetime.now()
            self.main_window.update_txt('%s: %s\n' % (now.strftime('%Y-%m-%d %H:%M:%S'),msg))
        else:
            self.main_window.update_txt(msg)

    def load_dialog(self):
        """Opens the dialogwindow when the Load button is pressed for Structure. """
        self.dialog = DialogWindow(self, self.root) #Sends self and self.root to DialogWindow
        self.dialog.configure(background = self.main_color)

    def prepare_pdb(self):
        """Starts PDB prepare if a legitimed pdb file is downloaded. """
        if not self.pdb_id:
            self.errorBox('Error','Pease load a PDB file before starting PDB prepare.')
            return
        else:
            self.prepare = PdbFilePrepare(self, self.root, self.pdb_id, self.workdir)
            self.prepare.configure(background = self.main_color)
            self.prepare.resizable()

    def create_parameters(self):
        """
        Opens prepare parameters (requires Impact from schrodinger!)
        """
        if not self.pdb_id:
            self.errorBox('Info','No pdb loaded.')
            self.log('info','Please load and preapare a pdb to generate parameters.')
            return
        else:
            self.log('info','Generate Q parameter session started')
            if not pt.checkPdb(self.pdb_id):
                pdb_id_q = pt.convertPdb(self.pdb_id, self.workdir, 0, 0)
                self.log('info','%s converted to Q format pdb (%s).' % (self.pdb_id.split('/')[-1], pdb_id_q.split('/')[-1]))
                self.pdb_id = self.workdir + '/' + pdb_id_q
                self.main_window.set_entryfield(self.pdb_id.split('/')[-1])

        self.makeparameters = CreatePrmLib(self, self.root)
        self.makeparameters.configure(bg=self.main_color)
        self.makeparameters.resizable()

    def open_trjmask(self):
        self.trjmask = TrjMask(self, self.root)
        self.trjmask.configure(background=self.main_color)

    def convert_top_pdb(self):
        """
        Convert a topology file to a pdb file
        :return:
        """
        self.log('info', 'Select topology file\n')
        topfile = askopenfilename(parent=self.root, initialdir=self.workdir, title='Select topology file',
                                  filetypes=(("top", "*.top"), ("All files", "*.*")))
        if not topfile:
            self.log(' ', 'No topology file loaded.\n')
            return

        pdb = qf.create_pdb_from_topology(topfile, self.q_settings['library'])

        if not pdb:
            self.log(' ', 'Could not convert topology to pdb.\n')
            return

        pdb_file = asksaveasfilename(parent=self.root, initialdir=self.workdir, title='Save pdb file as',
                                     filetypes=(("top", "*.top"), ("All files", "*.*")))

        if not pdb_file:
            self.log(' ', 'pdb file not saved\n')
            return

        pdbout = open(pdb_file, 'w')
        for line in pdb:
            pdbout.write(line)

        pdbout.close()
        self.log('info', 'Converted %s to %s' % (topfile.split('/')[-1], pdb_file.split('/')[-1]))

    def evb_calibration(self):
        """
        open evb calibration
        """
        self.evb_rcc = EvbCalibration(self, self.root)
        self.evb_rcc.configure(bg=self.main_color)
        self.log('info', 'EVB RRC session started.')

    def evb_arrhenius(self):
        """
        Open EVB arrhenius class to analyze thermodynamic activation parameters
        """
        self.evbarrhenius = EvbArrhenius(self, self.root)
        self.evbarrhenius.configure(bg=self.main_color)
        self.log('info', 'EVB Arrhenius session started.')

    def errorBox(self, tit ='Error', msg=''):
        if tit == 'Error':
            messagebox.showerror(tit, msg)
        elif tit == 'Info':
            messagebox.showinfo(tit, msg)
        else:
            messagebox.showwarning(tit, msg)

    def setup_md(self):
        """Opens the Setup MD for Qdyn when Setup --> MD is selected from the menubar
        """
        if not self.top_id:
            self.errorBox('Error','Please load topolgy before setting up MD')
            self.log('info', 'Please load topology before setting up MD.')
            return

        if not self.pdb_id:
            pdb_path = '%s/.tmp' % self.workdir
            if not os.path.isdir(pdb_path):
                os.makedirs(pdb_path)
            pdb_name = '%s.pdb' % self.top_id.split('/')[-1].split('.')[0]
            qf.write_top_pdb(self.top_id, pdb_name, pdb_path, self.q_settings['library'])
            self.pdb_id = '%s/%s' % (pdb_path, pdb_name)

        self.setup_md = SetupMd(self,self.root, self.pdb_id, self.top_id)
        self.setup_md.configure(background = self.main_color)
        self.setup_md.title('Setup MD')
        self.setup_md.resizable()

        self.log('info', 'MD session started.')

    def setup_lie(self):
        """
        Opens window to configure LIE run
        """
        self.lie_setup = SetupLie(self,self.root)
        self.lie_setup.configure(background = self.main_color)
        self.lie_setup.title('Setup LIE')
        self.lie_setup.resizable()

        self.log('info', 'LIE session started.')

    def setup_evb(self):
        """
        Opens window to configure MD/FEP/EVB simuations
        """
        if self.top_id:
            self.evb_setup = SetupEVB(self, self.root)
            self.evb_setup.configure(background=self.main_color)
            self.evb_setup.title('Setup EVB')
            self.evb_setup.resizable()
            self.log('info','Setup EVB session started.')
        else:
            self.log('info', 'Please load topology before setting up EVB.')
            self.errorBox('Error','Please load topolgy before setting up EVB')

    def setup_fep(self):
        """
        Opens window to prepare and setup FEP
        """
        if self.top_id:
            self.fep_setup = SetupFEP(self, self.root)
            self.fep_setup.configure(background=self.main_color)
            self.fep_setup.title('Setup FEP')
            self.fep_setup.resizable()
            self.log('info','Setup FEP session started.')
        else:
            self.log('info', 'Please load topology before setting up FEP.')
            self.errorBox('Error','Please load topolgy before setting up FEP')

    def setup_resfep(self):
        """
        Open the FEP residue mutation window (resFEP) that utilizes the predefined FEP protocol
        :return:
        """
        self.resfep_setup = ResFEP(self, self.root)
        self.resfep_setup.configure(background=self.main_color)
        self.log('info', 'Setup residue FEP session started')

    def analyze_lie(self):
        """
        Opens window to analyze LIE simulations
        """
        self.lie_analyze = AnalyzeLie(self,self.root)
        self.lie_analyze.configure(background= self.main_color)
        self.lie_analyze.title('Analyze LIE')
        self.lie_analyze.resizable()

    def analyze_evb_reactions(self):
        """
        Opens window to analyze EVB reaction profiles
        """
        self.evb_reactions = EvbReactions(self, self.root)
        self.evb_reactions.configure(background=self.main_color)
        self.evb_reactions.title('EVB Reaction Energies')

    def analyze_fep(self):
        """
        Opens window to analyze the free energy perturbation (part 1) from Qfep
        """
        self.fep_analyze = AnalyzeFEP(self, self.root)
        self.fep_analyze.configure(background=self.main_color)
        self.fep_analyze.title('Free Energy Perturbation Analyze')

    def analyze_resFEP(self):
        """
        Opens window to analyse resFEP (FEP protocol runs)
        :return:
        """
        self.resFEP_analyze = Analyse_resFEP(self, self.root)
        self.resFEP_analyze.configure(background=self.main_color)
        self.resFEP_analyze.title('Analyse resFEP')

    def analyze_energies(self):
        """
        Opens window to analyze md.log file energies
        """
        self.energy_analyze = AnalyzeEnergies(self, self.root)
        self.energy_analyze.configure(background=self.main_color)
        self.energy_analyze.title('Analyse trajectory energies')

    def analyze_qcalc(self, what):
        if self.pdb_id and self.top_id:
            self.qcalc_analyze = AnalyzeQcalc(self, self.root, what)
            self.qcalc_analyze.configure(background=self.main_color)
        else:
            self.log('info', 'Please load structure and topology files')
            self.errorBox('Error','Please load structure and topology files')



    def load_settings(self):
        """
        Opens the settings window
        """
        self.load_settings = QguiSettings(self, self.root)
        self.load_settings.configure(bg = self.main_color)
        self.load_settings.title('Settings')


    def prepare_topo(self):
        """Opens Topology Prepare when Prepare->Topology from the menubar
        is chosen. Checks if the file used is in the correct from,
        if not then converts it to pdb_q format.
        check_pdb_id has two outcomes, either it is True = the file is in Q format
        or False = the file is in standard pdb format.

        If check_pdb_id == False the convertPdb()" method from prepareTopology file
        is called."""

        if not self.pdb_id:
            self.errorBox('Error','Please load a structure (pdb) to generate topology from.')
        else:
            if not pt.checkPdb(self.pdb_id):
                pdb_id_q = pt.convertPdb(self.pdb_id, self.workdir, 0, 0)
                self.log('info','%s converted to Q format pdb (%s).' % (self.pdb_id.split('/')[-1], pdb_id_q.split('/')[-1]))
                self.pdb_id = self.workdir + '/' + pdb_id_q
                self.main_window.set_entryfield(self.pdb_id.split('/')[-1])
                self.topo_prepare = TopologyPrepare(self, self.root, self.pdb_id, self.prms, self.libs)
                self.topo_prepare.configure(background = self.main_color)
                self.topo_prepare.resizable()


            else:
                self.topo_prepare = TopologyPrepare(self, self.root, self.pdb_id, self.prms, self.libs)
                self.topo_prepare.configure(background = self.main_color)
                self.topo_prepare.resizable()


    def get_pdb_entry(self, pdb_id): #Fetch pdb file from website.
        """This method sends pdb_id to pdb.org server and
        get back a pdb-file. When an incorrect pdb_id is sent the response
        from pdb.org is (currently) a short file consisting of one sentence.
        The Exception is called when the length of the pdbfile is shorter than 20 lines.
        If the length of the response - for wrong pdb_id - from pdb.org is
        changed it can possibly cause failure in the feedback of the mentioned method."""

        if not '.pdb' in pdb_id:
            pdb_id = pdb_id + '.pdb'

        try:
            response = urlopen('%s/%s' % ('http://pdb.org/pdb/files', pdb_id))
            pdbfile = response.readlines()
            #The following can cause trouble if the length of the response from
            # pdb.org is changed when the pdb_id is faulty.
            if len(pdbfile) < 20:
                #PdbEntryNotFoundError() method is raised when the pdb_id is faulty.
                self.errorBox('Error', 'Error retrieving PDB file from pdb.org. Make sure to use the correct PDB id.')
                self.log('info', 'Error retrieving PDB file from pdb.org. Make sure to use the correct PDB id.')
                self.main_window.clear_button_pressed()

            else:
                with open(self.workdir + '/' + pdb_id, 'wb') as f:
                    for line in pdbfile:
                        f.write(line)
                self.pdb_id = '%s/%s' % (self.workdir, pdb_id)
                self.log('info', '%s successfully downloaded from pdb.org.' % self.pdb_id.split('/')[-1])
                self.dialog.dialog_destroy()
                self.main_window.set_entryfield(pdb_id.split('/')[-1])


        except IOError:
            self.errorBox('Error','Problems connecting to pdb.org. Make sure to use correct PDB id.')
            self.log('info', 'Could not connect to pdb.')

    def get_local_file(self):
        """Fetches a local file.""" #TODO
        filename = askopenfilename(parent = self.root, initialdir = self.workdir,
                                   filetypes=(("pdb", "*.pdb"),("All files","*.*")))
        pdb_id = filename
        if pdb_id != '':
            self.pdb_id = pdb_id
            self.log('info', '%s loaded' % self.pdb_id.split('/')[-1])
            self.main_window.set_entryfield(pdb_id.split('/')[-1])
        self.dialog.dialog_destroy()

    def get_top_file(self):
        filename = askopenfilename(parent = self.root, initialdir = self.workdir,
                                   filetypes=(("Topology", "*.top"),("All files","*.*")))
        top_id = filename
        if top_id != '':
            self.top_id = top_id
            self.log('info', '%s loaded' % self.top_id.split('/')[-1])
            self.main_window.set_topology_entryfield(top_id.split('/')[-1])

            #Check header in topology file for lib and prm files!
            print(self.libs)

            with open(self.top_id, 'r') as top:
                for line in top:
                    if 'LIB_FILES' in line:
                        for lib in line.split()[1].split(';'):
                            if lib not in self.libs:
                                print(lib)
                                if os.path.isfile(lib):
                                    print('From topologyfile, %s was added to library files (not in settings).')
                                    self.libs.append(lib)
                                else:
                                    print("Warning: Lib file %s in topology not found!" % lib)
                        break


    def update_pdb_id_entryfield(self, pdb_id=None):
        """
        Sends pdb_id to MainWindow-class method set_entryfield(). """
        if pdb_id != None:
            self.main_window.set_entryfield(pdb_id.split('/')[-1])
        else:
            self.main_window.set_entryfield(self.pdb_id.split('/')[-1])

    def reset_pdb(self, new_pdb):
        """Resets self.pdb_id to empty.
        This method is called from MainWindow-class, when the structure-entryfield is emptied.
        Does not delete already downloaded pdbfile from the directory. """
        self.pdb_id = new_pdb

    def reset_topology(self, new_top):
        """Resets self.pdb_id to empty.
        This method is called from MainWindow-class, when the structure-entryfield is emptied.
        Does not delete already downloaded pdbfile from the directory. """
        self.top_id = new_top

    def view_pymol(self):

        self.pymol_session = ViewPyMol(self, self.root)
        self.pymol_session.configure(background=self.main_color)
        self.pymol_session.resizable()

        self.pymol_running = True

    def exit(self):
        """This method is called from the mainwindow menubar->filemenu->Close is chosen."""
        if self.pymol_running:
            self.errorBox('Info','Please close Q-pymol first!')
            print('Pymol is running!')
        self.root.quit()

    def aboutQ(self):
        """
        Opens up the about Qgui window
        """
        self.about_q = AboutQ(self,self.root)
        self.about_q.configure(background = self.main_color)


class DialogWindow(Toplevel):
    """Implements the dialog-box that opens upp when 'Load' button for structure is
    pressed. DialogWindow-class har got a reference to MainWindow-class and calls
    to methods from controller-class: get_local_file and import_button_clicked. """

    def __init__(self, app, root): #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root

        self.main_color = self.app.main_color

        self.dialog_box()

    def dialog_box(self):
        """Defines a dialog-box when the Load-button is pressed from MainWindow."""

        self.title('Load/import')
        self.geometry("250x100+400+100")
        self.maxsize(300,150)
        self.minsize(250,100)
        browse = Button(self, text = 'Browse', command = self.app.get_local_file) #Reference to Qgui-class
        browse.grid(row = 0, column = 1, sticky = 'e')
        browse.config(highlightbackground = self.main_color, relief = GROOVE)
        browse_label = Label(self, text = 'Browse for local file: *.pdb')
        browse_label.grid(row = 0, column = 0, sticky = 'w')
        browse_label.config(background = self.main_color)

        self.pdb_entryfield = Entry(self)
        self.pdb_entryfield.grid(row = 2, column = 0)
        self.pdb_entryfield.config(highlightthickness = 0, relief = GROOVE)
        import_pdb = Button(self, text = 'Import', command = self.import_button_clicked)
        import_pdb.grid(row = 2, column = 1)
        import_pdb.config(highlightbackground = self.main_color, relief = GROOVE)
        import_label = Label(self, text = 'Import structure from PDB:')
        import_label.grid(row = 1, column = 0)
        import_label.config(background = self.main_color)

    def import_button_clicked(self):
        #Get pdb id from textfield
        pdb_id = self.pdb_entryfield.get()
        #Call get_pdb_entry in controller(Qgui-class)
        self.app.get_pdb_entry(pdb_id)

    def dialog_destroy(self):

        self.destroy()


if __name__ == '__main__':
    root = Tk()
    app = QGui(root)
    root.mainloop()
