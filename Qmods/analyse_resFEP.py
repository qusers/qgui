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
    OptionMenu, StringVar, Spinbox, SINGLE, LabelFrame, MULTIPLE, Entry
from tkFileDialog import askopenfilename, askdirectory
import qgui_functions as qf
import numpy as np
from tkSimpleDialog import askstring
import tkFont
import os


class Analyse_resFEP(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root


        #Store paths for resFEP
        self.feps_paths = dict()

        #store resFep energies:
        self.feps = dict()

        self.dialog_window()

    def add_fep_run(self):
        """

        :return:
        """
        #Select main path for FEP1..FEPN
        self.app.log(' ', '\nSelect resFEP main directory (inputfiles, FEP1, FEP2...FEPN)\n')
        title_dialog = 'Select resFEP main directory'
        path = self.app.workdir
        rundir = askdirectory(parent=self, mustexist=True, title=title_dialog, initialdir=path)

        if not rundir:
            return
        print(rundir)

        #Ask for a suitable name for this FEP run from user
        fep_title = askstring('Add title', 'resFEP title:', parent=self)

        if not fep_title:
            fep_title = rundir.split('/')[-1]

        if not fep_title in self.feps.keys():
            self.feps[fep_title] = dict()
            self.feps_paths[fep_title] = rundir
            self.feplist.insert(END, fep_title)
        else:
            self.app.errorBox('Warning', 'Title: %s already exists!' % fep_title)
            self.app.log(' ', 'Title: %s already exist. Give it a different name.\n' % fep_title)
            return

        #Time to run Qfep
        self.calculate_fep(fep_title)

    def calculate_fep(self, fep_title, recalcFEP=False):
        """

        :param fep_title:
        :return:
        """
        fepdirs = self.get_fep_dirs(self.feps_paths[fep_title])
        print fepdirs

        if len(fepdirs) < 1:
            self.app.errorBox('Warning', 'Found no FEP directories!')
            return

        fepout = dict()

        #Go through all FEP directories
        for i in sorted(fepdirs.keys()):
            fepdir = fepdirs[i]
            fep = fepdir.split('/')[-1]
            if not fep in fepout.keys():
                fepout[fep] = dict()

            os.chdir(fepdir)
            print('Analysing %s' % fep)

            temps = list()
            #Find temperature directories
            for temp in sorted(filter(os.path.isdir, os.listdir(os.getcwd()))):
                if temp[0].isdigit:
                    temps.append(temp)

            #Go through temperature directories and analyse FEP runs
            for temp in temps:
                if temp not in fepout[fep].keys():
                    fepout[fep][temp] = dict()
                    fepout[fep][temp]['dG'] = list()
                    fepout[fep][temp]['dGf'] = list()
                    fepout[fep][temp]['dGr'] = list()
                print('Analysing FEP at temperature %s' % temp)
                tdir = '%s/%s' % (fepdir, temp)
                os.chdir(tdir)

                #Go through parallel runs
                for i in sorted(filter(os.path.isdir, os.listdir(os.getcwd()))):
                    calcFEP = False
                    rundir = '%s/%s' % (tdir, i)
                    os.chdir(rundir)

                    #check if qfep.out exist
                    if not os.path.isfile('qfep.out'):
                        calcFEP = True
                    elif recalcFEP:
                        calcFEP = True

                    #Run Qfep
                    if calcFEP:
                        runQfep = qf.write_qfep_input(qfep_path=rundir, states=2, temp=float(temp))

                        if runQfep:
                            qf.run_Qfep(qpath=rundir, qfep_inp='qfep.inp', qfep=self.app.q_settings['executables'][2])
                        else:
                            print('Found no energy files in %s' % rundir)

                    #Read qfep.out and collect FEP summary
                    FEPdata = qf.get_qfep_part1(qpath=rundir, qfep='qfep.out')

                    if FEPdata:
                        fepout[fep][temp]['dGf'].append(FEPdata['sum_dGf'][-1])
                        fepout[fep][temp]['dGr'].append(FEPdata['sum_dGr'][0])
                        fepout[fep][temp]['dG'].append(FEPdata['dG'][-1])

        #Calculate average FEP and update self.feps
        self.calc_ave_fep(fep_title, fepout)

    def calc_ave_fep(self, fep_title, fepdata):
        """
        Calculates average FEP values and updates self.feps
        :param fepdata: dict() generated by calculate_fep
        :return:
        """

        for fep in sorted(fepdata.keys()):
            for temp in sorted(fepdata[fep].keys()):
                if temp not in self.feps[fep_title].keys():

                    self.feps[fep_title][temp] = dict()
                    self.feps[fep_title][temp]['dG'] = [0, 0]
                    self.feps[fep_title][temp]['dGf'] = [0, 0]
                    self.feps[fep_title][temp]['dGr'] = [0, 0]
                    self.feps[fep_title][temp][fep] = dict()

                ave_frwd = np.average(fepdata[fep][temp]['frwd'])
                frwd_sem = np.std(fepdata[fep][temp]['frwd'])/np.sqrt(len(fepdata[fep][temp]['frwd']))

                ave_rev = np.average(fepdata[fep][temp]['rev'])
                rev_sem = np.std(fepdata[fep][temp]['rev'])/np.sqrt(len(fepdata[fep][temp]['rev']))

                ave_dg = np.average(fepdata[fep][temp]['ave'])
                dg_sem = np.std(fepdata[fep][temp]['ave'])/np.sqrt(len(fepdata[fep][temp]['ave']))

                self.feps[fep_title][temp][fep]['dGf'] = [ave_frwd, frwd_sem]
                self.feps[fep_title][temp][fep]['dGr'] = [ave_rev, rev_sem]
                self.feps[fep_title][temp][fep]['dG'] = [ave_dg, dg_sem]

                self.feps[fep_title][temp]['dG'][0] += ave_dg
                self.feps[fep_title][temp]['dG'][1] = np.sqrt(self.feps[fep_title][temp]['dG'][1] + dg_sem)

                self.feps[fep_title][temp]['dGf'][0] += ave_frwd
                self.feps[fep_title][temp]['dGf'][1] = np.sqrt(self.feps[fep_title][temp]['dGf'][1] + frwd_sem)

                self.feps[fep_title][temp]['dGr'][0] += ave_rev
                self.feps[fep_title][temp]['dGr'][1] = np.sqrt(self.feps[fep_title][temp]['dGr'][1] + rev_sem)

    def get_fep_dirs(self, rundir):
        """
        Collects all paths FEP$i ($i = integer)
        :param rundir:
        :return: fepdirs
        """
        fepdirs = dict()

        os.chdir(rundir)
        for f in sorted(filter(os.path.isdir, os.listdir(os.getcwd()))):
            print f
            if f.startswith('FEP'):
                if f[3:].isdigit():
                    nr = int(f[3:])
                    fepdirs[nr] = '%s/%s' % (rundir, f)

        os.chdir(self.app.workdir)

        return fepdirs

    def del_fep_run(self):
        """

        :return:
        """
        pass

    def add_fep_combine(self, sign):
        """
        adds the selected FEP run to the combine FEP runs with + or - (thermodynamic cycle)
        :param sign: + or -
        :return:
        """
        selected = self.feplist.curselection()

        if len(selected) < 1:
            self.change_signs_combine(sign)
            return

        for i in selected:
            fep = self.feplist.get(i)
            self.comb_list.insert(END, '%s %s' % (sign, fep))

    def change_signs_combine(self, sign):
        """
        Changes the sign of selected entry in combine list
        :param sign:
        :return:
        """
        selected = self.comb_list.curselection()

        if len(selected) < 1:
            return

        for i in selected:
            fep = self.comb_list.get(i).split()[-1]
            self.comb_list.delete(i)
            self.comb_list.insert(i, '%s %s' % (sign, fep))

    def recalc_fep(self):
        """
        Takes selected FEP title en recalculates all FEPs.
        :return:
        """
        selected = self.feplist.curselection()

        if len(selected) < 1:
            return

        for i in selected:
            fep_title = self.feplist.get(i)
            self.app.log('info', 'Recalculating FEP for %s. Please wait..\n' % fep_title)
            self.calculate_fep(fep_title=fep_title, recalcFEP=True)

    def combine_feps(self):
        """
        Combines FEP dGs as defined in Combine list (add or substract)
        :return:
        """
        comb = dict()
        fep_temps = set()
        selected = self.comb_list.curselection()

        if len(selected) < 1:
            return

        #How to combine?
        for i in selected:
            sign = self.comb_list.get(i).split()[0]
            fep = self.comb_list.get(i).split()[1]
            if sign == '+':
                comb[fep] = 1.
            else:
                comb[fep] = -1.

            #Check that fep exist in self.feps, if not abort!
            if not fep in self.feps.keys():
                print('Could not find FEP title %s' % fep)
                return
            else:
                if len(self.feps[fep].keys()) < 1:
                    self.app.errorBox('Warning', 'Found no FEP data for %s.' % fep)
                    return

            #Add temperature(s) to set and remove non-identical temperatures
            tmp = set()
            for temp in self.feps[fep].keys():
                tmp.add(temp)
            if len(fep_temps) == 0:
                fep_temps = tmp
            else:
                fep_temps = tmp & fep_temps

        #Create new entry in self.feps
        fep_title = self.combine_title.get()
        self.feps[fep_title] = dict()

        #Calculate average FEP per temperature:
        for temp in fep_temps:
            if not temp in self.feps[fep_title].keys():
                self.feps[fep_title][temp] = dict()
                self.feps[fep_title][temp]['dG'] = [0, 0]
                self.feps[fep_title][temp]['dGf'] = [0, 0]
                self.feps[fep_title][temp]['dGr'] = [0, 0]

            for i in comb.keys():
                dG = comb[i] * self.feps[i][temp]['dG'][0]
                dG_sem = self.feps[i][temp]['dG'][1]

                dGf = comb[i] * self.feps[i][temp]['dGf'][0]
                dGf_sem = self.feps[i][temp]['dGf'][1]

                dGr = comb[i] * self.feps[i][temp]['dGr'][0]
                dGr_sem = self.feps[i][temp]['dG'][1]

                self.feps[fep_title][temp]['dG'][0] += dG
                self.feps[fep_title][temp]['dG'][1] = np.sqrt(self.feps[fep_title][temp]['dG'][1] + dG_sem)

                self.feps[fep_title][temp]['dGf'][0] += dGf
                self.feps[fep_title][temp]['dGf'][1] = np.sqrt(self.feps[fep_title][temp]['dGf'][1] + dGf_sem)

                self.feps[fep_title][temp]['dGr'][0] += dGr
                self.feps[fep_title][temp]['dGr'][1] = np.sqrt(self.feps[fep_title][temp]['dGr'][1] + dGr_sem)

    def clear_combined(self):
        """
        Deletes entries in combine FEP listbox
        :return:
        """

        self.comb_list.delete(0, END)

    def export_table(self):
        """

        :return:
        """
        pass

    def save_session(self):
        """

        :return:
        """
        pass

    def load_session(self):
        """

        :return:
        """
        pass

    def dialog_window(self):

        self.title('Analyse resFEP')

        mainframe = Frame(self, bg=self.main_color)
        mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame with FEP runs and combine list
        frame1 = Frame(mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0)

        #Frame with listbox of FEP summary/dG values with sem
        frame2 = Frame(mainframe, bg=self.main_color)
        frame2.grid(row=1, column=0, pady=(10,0))

        #Frame with Quit / save / load
        frame3 = Frame(mainframe, bg=self.main_color)
        frame3.grid(row=2, column=0, pady=(10,0))

        ##### FRAME 1
        #Label for FEP listbox
        fepruns = Label(frame1, text='FEP runs:', bg= self.main_color)
        fepruns.grid(row=0, column=0, columnspan=2)

        #listbox with FEP run entries
        feplist_scroll = Scrollbar(frame1)
        feplist_scroll.grid(row=1, rowspan=4, column=1, sticky='nsw')
        self.feplist = Listbox(frame1, yscrollcommand=feplist_scroll.set, width=25, height=4,
                                 highlightthickness=0, relief=GROOVE, selectmode=EXTENDED, exportselection=True)
        feplist_scroll.config(command=self.feplist.yview)
        self.feplist.grid(row=1, rowspan=4, column=0, sticky='nse')
        self.feplist.config(font=tkFont.Font(family="Courier", size=12))
        #self.feplist.bind('<<ListboxSelect>>', self.feplist_event)

        #Add fep runs
        add_fep = Button(frame1, text='Add FEP runs', width=10, highlightbackground=self.main_color,
                         command=self.add_fep_run)
        add_fep.grid(row=5, column=0, columnspan=2)

        #Delete fep runs
        del_fep = Button(frame1, text='Delete', width=10, highlightbackground=self.main_color, command=self.del_fep_run)
        del_fep.grid(row=6, column=0, columnspan=2)


        #Recalculate FEP (run Qfep again)
        recalc_fep = Button(frame1, text='Recalc FEP', width=10, highlightbackground=self.main_color,
                            command=self.recalc_fep)
        recalc_fep.grid(row=7, column=0, columnspan=2)

        #Right arrow to combine box label
        right_arrow = Label(frame1, text=u'\u2192', bg=self.main_color)
        right_arrow.grid(row=2, rowspan=2, column=2)

        #+ value to selected FEP run
        plus_dg = Button(frame1, text ='+', highlightbackground=self.main_color,
                         command=lambda: self.add_fep_combine('+'))
        plus_dg.grid(row=2, column=3)

        #- value to selected fep run
        minus_dg = Button(frame1, text ='-', highlightbackground=self.main_color,
                          command=lambda: self.add_fep_combine('-'))
        minus_dg.grid(row=3, column=3)

        #Label for FEP listbox
        combine = Label(frame1, text='Combine:', bg= self.main_color)
        combine.grid(row=0, column=4, columnspan=2)

        #listbox with added FEP runs to be combined (added/substracted)
        comb_list_scroll = Scrollbar(frame1)
        comb_list_scroll.grid(row=1, rowspan=4, column=5, sticky='nsw')
        self.comb_list = Listbox(frame1, yscrollcommand=comb_list_scroll.set, width=25, height=4,
                                 highlightthickness=0, relief=GROOVE, selectmode=SINGLE, exportselection=True)
        comb_list_scroll.config(command=self.comb_list.yview)
        self.comb_list.grid(row=1, rowspan=4, column=4, sticky='nse')
        self.comb_list.config(font=tkFont.Font(family="Courier", size=12))
        #self.feplist.bind('<<ListboxSelect>>', self.feplist_event)

        #Clear comb_list
        clear_comblist = Button(frame1, text='Clear', width=10, highlightbackground=self.main_color,
                                command=self.clear_combined)
        clear_comblist.grid(row=7, column=4, columnspan=2)

        #Entry field for user defined title for the combine FEP
        self.combine_title = Entry(frame1, width=17, highlightthickness=0, relief=GROOVE)
        self.combine_title.grid(row=5, column=4, columnspan=2)
        self.combine_title.insert(0, 'ddG title')

        #Button for combining the FEPS defined in comb_list
        calc_combine = Button(frame1, text='Combine', width=10, highlightbackground=self.main_color,
                              command=self.combine_feps)
        calc_combine.grid(row=6, column=4, columnspan=2)

        ##### FRAME 2
        ave_ene = Label(frame2, text='Average energies:', bg=self.main_color)
        ave_ene.grid(row=0, column=0, columnspan=2)

        #List with sumarised dG values
        dg_list_scroll = Scrollbar(frame2)
        dg_list_scroll.grid(row=1, rowspan=4, column=1, sticky='nsw')
        self.dg_list = Listbox(frame2, yscrollcommand=dg_list_scroll.set, width=55, height=10,
                                 highlightthickness=0, relief=GROOVE, selectmode=SINGLE, exportselection=True)
        dg_list_scroll.config(command=self.dg_list.yview)
        self.dg_list.grid(row=1, rowspan=4, column=0, sticky='nse')
        self.dg_list.config(font=tkFont.Font(family="Courier", size=12))
        #self.feplist.bind('<<ListboxSelect>>', self.feplist_event)

        #Export table
        export_table = Button(frame2, text='Export table', width=10, highlightbackground=self.main_color,
                              command=self.export_table)
        export_table.grid(row=5, column=0, columnspan=2)


        ##### FRAME 3
        #SAVE
        save_session = Button(frame3, text='Save', width=10, highlightbackground=self.main_color,
                              command=self.save_session)
        save_session.grid(row=0, column=0)
        #LOAD
        load_session = Button(frame3, text='Load', width=10, highlightbackground=self.main_color,
                              command=self.load_session)
        load_session.grid(row=0, column=1)

        #QUIT
        quit_session = Button(frame3, text='Quit', width=10, highlightbackground=self.main_color,
                              command=self.destroy)
        quit_session.grid(row=0, column=2)