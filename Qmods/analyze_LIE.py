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

from tkinter import Label, Button, Listbox, EXTENDED, Spinbox, Entry, LabelFrame, Text, Frame, \
    Toplevel, DISABLED, END, GROOVE, NORMAL
import tkinter.font
from tkinter.filedialog import askdirectory, asksaveasfilename
import os
import time
import numpy as np
from select_return import SelectReturn
from Qplot import Qplot
import mdlog_energies as mdle
from fit_qlie import FitQlie


class AnalyzeLie(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root
        self.app.log('info', 'Analyze LIE session started.')

        self.main_color = self.app.main_color

        self.dialog_window()

        #Insert default LIE parameter values:
        self.alpha_entry.delete(0, END)
        self.alpha_entry.insert(0, '0.18')
        self.beta_entry.delete(0,END)
        self.beta_entry.insert(0, '0.50')
        self.gamma_entry.delete(0,END)
        self.gamma_entry.insert(0, '0.00')


        #Use Qfep or log files to get Q-energies:
        self.use_qfep = False

        #Initialize ligand and complex lists
        self.lig_qene = []
        self.comp_qene = []
        self.lig_log = []
        self.comp_log = []

        #Check if ligand and complex directories exist in workdir:
        #if os.path.isdir(self.app.workdir + '/ligand'):
        #    self.get_lie_files('ligand', False)
        #if os.path.isdir(self.app.workdir + '/complex'):
        #    self.get_lie_files('complex', False)


    def get_lie_files(self, qsystem='complex', askforfile=False):
        """
        Gets directory to ligand/complex and collects ene files for Qfep
        and log files for interaction energy analyzis.
        """

        if qsystem =='ligand':
            self.lig_ene = []
            self.lig_unique_log = []
            self.lig_log = []
            self.lig_qene = None
            enefiles = self.lig_ene
            unique_logs = self.lig_unique_log
            md_logs = self.lig_log
            entry = self.ligand_dir
            runs_entry = self.ligand_runs_entry
            qrow=0
        elif qsystem == 'complex':
            self.comp_ene = []
            self.comp_unique_log = []
            self.comp_log = []
            self.comp_qene = None
            enefiles = self.comp_ene
            unique_logs = self.comp_unique_log
            md_logs = self.comp_log
            entry = self.complex_dir
            runs_entry = self.complex_runs_entry
            qrow=1

        #Get directory:
        if askforfile:
            qdir = askdirectory(parent=self, mustexist=False, title='Select LIE directory for %s' % qsystem, initialdir= self.app.workdir)
            if qdir == '':
                runs_entry.config(state=NORMAL)
                runs_entry.delete(0.0, END)
                runs_entry.config(state=DISABLED)
                return
            else:
                entry.config(state=NORMAL)
                entry.delete(0.0, END)
                entry.insert(0.0, '../%s' % qdir.split('/')[-1])


        else:
            if qsystem == 'ligand':
                qdir = self.app.workdir + '/ligand'
            elif qsystem == 'complex':
                qdir = self.app.workdir + '/complex'

        #Fill all Qdyn logs to logfiles firs:
        logfiles = []

        dir_nr = 0
        multiple_runs = False
        #Try to look for md.log files in qdir/1
        if os.path.isdir('%s/1' % qdir):
            self.app.log('info','Collecting MD files from subdirectories...')
            logline = len(self.app.main_window.txt.get(0.0, END).split('\n'))
            logline = str(logline) + '.0'
            self.update()
            multiple_runs = True
            done = False
            while not done:
                self.update()
                if os.path.isdir('%s/%d' % (qdir, (dir_nr + 1))):
                    dir_nr += 1
                    self.app.main_window.txt.config(state=NORMAL)
                    self.app.main_window.txt.insert(END, '\n../%s/%d:\n' % (qdir.split('/')[-1], dir_nr))
                    logfiles_in_dir = self.find_files(['.log'], '%s/%d' % (qdir, dir_nr))
                    for files in logfiles_in_dir:
                        if mdle.is_md_log(files):
                            logfiles.append(files)
                            print(files)
                            self.app.main_window.txt.insert(END, '%s\n' % files.split('/')[-1])
                            self.app.main_window.txt.yview(END)
                            self.update()
                    self.app.main_window.txt.delete(logline, END)
                else:
                    done = True

        #If 1 does not exist, look in current folder
        if not multiple_runs:
            self.app.log('info','Collecting MD files from directory...')
            logline = len(self.app.main_window.txt.get(0.0, END).split('\n'))
            logline = str(logline) + '.0'
            self.update()
            logfiles_in_dir = self.find_files(['.log'], qdir)
            self.app.main_window.txt.config(state=NORMAL)
            self.app.main_window.txt.insert(END, '\n \n')
            for files in logfiles_in_dir:
                if mdle.is_md_log(files):
                    logfiles.append(files)
                    print(files)
                    self.app.main_window.txt.insert(END, '%s\n' % files.split('/')[-1])
                    self.app.main_window.txt.yview(END)
                    self.update()
            self.app.main_window.txt.delete(logline, END)

        runs_entry.grid(in_=self.general_label, row=qrow, column=3)
        runs_entry.config(state=NORMAL)
        runs_entry.delete(0.0, END)
        runs_entry.config(state=NORMAL)
        if dir_nr > 1 or dir_nr == 0:
            sub = 'subdirs.'
        elif dir_nr == 1:
            sub = 'subdir.'
        if len(logfiles) > 1 or len(logfiles) == 0:
            fil = 'files'
        elif len(logfiles) == 1:
            fil = 'file'
        runs_entry.insert(0.0, '%d %s / %d MD %s' % (dir_nr, sub, len(logfiles), fil))
        runs_entry.config(state=DISABLED)


        #Get main name of MD runs from dirs, and ask user to select:
        tmptitle = '%s MD file(s) to analyze' % qsystem
        select_logs = []
        for logs in logfiles:
            if logs.split('/')[-1] not in select_logs:
                select_logs.append(logs.split('/')[-1])
        if len(select_logs) > 1:
            self.app.log('info', 'Found %d unique MD files. Select the ones to analyze.' % len(select_logs))
            self.select_name = SelectReturn(self, self.root, select_logs, tmptitle, None)
            self.select_name.configure(background=self.main_color)
            self.select_name.resizable()
            #Wait for user to select MD files to analyze:
            self.app.main_window.txt.config(state=NORMAL)
            self.app.main_window.txt.insert(END, '\nPlease select %s files to analyze\n' % qsystem)
            logline = len(self.app.main_window.txt.get(0.0, END).split('\n')) - 2
            logline = str(logline) + '.0'
            dotcount = 0
            while len(unique_logs) == 0:
                self.app.main_window.txt.delete(logline, END)
                self.app.main_window.txt.insert(END, '\nPlease select %s files to analyze'
                                                     % qsystem + (dotcount * '.') + '\n')
                self.update()
                time.sleep(0.2)
                dotcount += 1
                if dotcount == 40:
                    dotcount = 0
            if 'cancel' not in unique_logs:
                self.app.log('info', 'LIE %s MD files specified' % qsystem)
            else:
                unique_logs = []
                self.app.log('info', 'LIE MD files not selected. Please specify this to run the analyze.')
                entry.delete(0.0, END)
                runs_entry.config(state=NORMAL)
                runs_entry.delete(0.0, END)
                runs_entry.config(state=DISABLED)

        for logs in logfiles:
            if logs.split('/')[-1] in unique_logs:
                md_logs.append(logs)

        print(md_logs)

        if len(self.comp_log) != 0 and len(self.lig_log) != 0:
            self.compute_button.config(state=NORMAL)

    def find_files(self, infile, qdir):
        """
        looks for files in directory. Contents of filename can be several,
        and must be specified as a list. Last entry in list is the file ends with
        line.
        Returns a list of all files.
        """
        list_of_files = []
        for files in os.listdir(qdir):
            matches = 0
            if files.endswith(infile[-1]):
                for term in range(len(infile)):
                    if infile[term] not in files:
                        break
                    elif infile[term] in files:
                        matches += 1
                    if matches >= len(infile):
                        list_of_files.append('%s/%s' % (qdir, files))

        return list_of_files

    def do_qfep(self):
        """
        Collects all .en files and runs Qfep.
        Q-energies is then collected from qfep.out
        """
        pass

    def extract_energies(self, system='complex'):
        if system == 'complex':
            logfiles = self.comp_log
        else:
            logfiles = self.lig_log

        self.app.log(' ', '\nExtracting %s energies ...\n (this may take a while)\n' % system)
        self.app.main_window.txt.yview(END)
        self.update()

        qene, ave, stderr = mdle.get_q_energies(logfiles)

        self.app.log(' ', '\nDone extracting %s energies\n' % system)
        self.app.main_window.txt.yview(END)

        return qene, ave, stderr

    def compute_lie(self):
        """
        Gets Q-energies (from log files or qfep.out) for ligand
        and complex and computes dG with specified parameters
        defined in the General section of the Analyze LIE window

        ['Q-Q', 'Q-prot', 'Q-wat', 'Q-surr.', 'Q-any']
        """
        if len(self.lig_log) == 0 and len(self.comp_log) == 0:
            self.app.errorBox('Warning', 'MD files missing for ligand and complex!')
            return
        elif len(self.lig_log) == 0 and len(self.comp_log) > 0:
            self.app.errorBox('Warning','MD files missing for ligand.')
            return
        elif len(self.lig_log) > 0 and len(self.comp_log) == 0:
            self.app.errorBox('Warning', 'MD files missing for complex')
            return
        else:
            pass

        if self.use_qfep:
            #Get Q-energies from qfep.out (use self.do_qfep)
            self.app.errorBox('Info', 'Sorry, Qfep is not implemented here yet.')
            return
        if not self.use_qfep:
            #Get Q-energies from MD log files if they do not exist:
            if not self.lig_qene:
                self.lig_qene, self.lig_ave, self.lig_stderr = self.extract_energies('ligand')
                self.update()

            if not self.comp_qene:
                self.comp_qene, self.comp_ave, self.comp_stderr = self.extract_energies('complex')
                self.update()

            lig_ave = self.lig_ave
            lig_stderr = self.lig_stderr
            comp_ave = self.comp_ave
            comp_stderr = self.comp_stderr


        alpha = float(self.alpha_entry.get())
        beta = float(self.beta_entry.get())
        gamma = float(self.gamma_entry.get())

        dg = (beta * ((comp_ave[2][0] + comp_ave[1][0]) - lig_ave[2][0]) +
              (alpha * ((comp_ave[2][1] + comp_ave[1][1]) - lig_ave[2][1])) + gamma)

        dg_stderr = np.sqrt(beta**2 * (lig_stderr[2][0]**2 + comp_stderr[2][0]**2 + comp_stderr[1][0]**2) +
                            (alpha**2 * (lig_stderr[2][1]**2 + comp_stderr[2][1]**2 + comp_stderr[1][1]**2)))

        insert_in_entries = [lig_ave[2][0], lig_stderr[2][0],
                             lig_ave[2][1], lig_stderr[2][1],
                             comp_ave[2][0], comp_stderr[2][0],
                             comp_ave[2][1], comp_stderr[2][1],
                             comp_ave[1][0], comp_stderr[1][0],
                             comp_ave[1][1], comp_stderr[1][1],
                             dg, dg_stderr]

        #Collect all entries:
        entries = [self.ligand_el_entry, self.ligand_el_stderr_entry,
                   self.ligand_vdw_entry, self.ligand_vdw_stderr_entry,
                   self.complex_el_w_entry, self.complex_el_w_stderr_entry,
                   self.complex_vdw_w_entry, self.complex_vdw_w_stderr_entry,
                   self.complex_el_p_entry, self.complex_el_p_stderr_entry,
                   self.complex_vdw_p_entry, self.complex_vdw_p_stderr_entry,
                   self.dg_entry, self.result_stderr_entry]
        #Config all entry boxes:
        for entry_ in entries:
            entry_.config(state=NORMAL)
            entry_.config(font=tkinter.font.Font(family="Courier", size=12))

        #Insert results:
        for result in range(len(entries)):
            entries[result].delete(0, END)
            entries[result].insert(0, '%10.2f' % insert_in_entries[result])

        self.app.log('info', 'LIE computation completed.')
        self.save_button.config(state=NORMAL)
        self.ligname.config(state=NORMAL)

    def plot_selection(self):
        """
        Takes system selection (ligand/complex or both) and interacion selection(s)
        and opens a new window with the corrresponding plot (matplotlib)
        """
        vdw_ylist = []
        el_ylist = []

        vdw_xlist = []
        el_xlist = []

        vdw_titles = []
        el_titles = []

        qtitles = []
        plottitles = []

        #Get selection from list:
        try:
            system_index = list(map(int, self.system_list.curselection()))
            type_index = list(map(int, self.interaction_list.curselection()))
            print(system_index)
            if len(type_index) == 0:
                return
            if len(system_index) == 0:
                return

            for index_ in system_index:
                for indextype in type_index:
                    system_name = self.system_list.get(index_).strip()
                    if system_name == 'COMPLEX':
                        title_system = '(comp.)'
                        mdlogs = self.comp_log
                        qenergies = self.comp_qene
                        try:
                            base_name = self.comp_log[0].split('.log')[0] + '.inp'
                        except:
                            base_name = None
                    elif system_name == 'LIGAND':
                        title_system = '(lig)'
                        mdlogs = self.lig_log
                        qenergies = self.lig_qene
                        try:
                            base_name = self.lig_log[0].split('.log')[0] + '.inp'
                        except:
                            base_name = None
                    if not qenergies:
                        if mdlogs:
                            if system_name == 'COMPLEX':
                                self.comp_qene, self.comp_ave, self.comp_stderr = self.extract_energies('complex')
                                qenergies = self.comp_qene
                            elif system_name == 'LIGAND':
                                self.lig_qene, self.lig_ave, self.lig_stderr = self.extract_energies('ligand')
                                qenergies = self.lig_qene
                        else:
                            self.app.errorBox('Error', 'No %s dir specified' % system_name)
                            return
                    title_interactions = self.interaction_list.get(indextype)
                    if title_interactions.split()[1] == 'el':
                        list_to_fill = el_ylist
                        x_to_fill = el_xlist
                        title_list = el_titles
                        qindex_2 = 0
                        if not 'El' in plottitles:
                            plottitles.append('El')
                    else:
                        list_to_fill = vdw_ylist
                        x_to_fill = vdw_xlist
                        title_list = vdw_titles
                        qindex_2 = 1
                        if not 'vdW' in plottitles:
                            plottitles.append('vdW')
                    title_list.append('%s %s' % (title_interactions.split()[0], title_system))
                    if title_interactions.split()[0] == 'Q-Q':
                        qindex_1 = 0
                    elif title_interactions.split()[0] == 'Q-prot':
                        qindex_1 = 1
                    elif title_interactions.split()[0] == 'Q-wat':
                        qindex_1 = 2
                    elif title_interactions.split()[0] == 'Q-surr':
                        qindex_1 = 3
                    list_to_fill.append(qenergies[qindex_1][qindex_2])
                    tmplist = []

                    #Find stepsize and output intervall to fill in X-values (ns)
                    stepsize = None
                    interval = None

                    if base_name:
                        inputfile = open(base_name,'r').readlines()
                        for line in inputfile:
                            if 'stepsize ' in line:
                                stepsize = float(line.split()[1])
                            if 'output ' in line:
                                interval = float(line.split()[1])
                    else:
                        self.app.errorBox('Warning','Could not find stepsize and interval.'
                                                    ' Using 1 ps stepsize and interval 5.')
                        stepsize = 1.00
                        interval = 5.00
                    for j in range(0, len(qenergies[qindex_1][qindex_2])):
                        tmplist.append((float(j) * stepsize * interval) / 1000000.00)
                    x_to_fill.append(tmplist)

        except:
            print('No selections')
            return

        qtitles = []
        qenergies = []
        xlist = []
        if len(el_titles) > 0:
            qtitles.append(el_titles)
            qenergies.append(el_ylist)
            xlist.append(el_xlist)
        if len(vdw_titles) > 0:
            qtitles.append(vdw_titles)
            qenergies.append(vdw_ylist)
            xlist.append(vdw_xlist)

        #qtitles = [el_titles,vdw_titles]
        #qenergies = [el_ylist, vdw_ylist]
        #xlist = [el_xlist, vdw_xlist]

        print(plottitles)
        print(qtitles)

        self.plot_ = Qplot(self, self.root, qenergies, xlist, plottitles, qtitles, 'Time (ns)', 'Energy (kcal/mol)')
        self.plot_.resizable()
        self.plot_.config(background=self.main_color)

    def save_results(self):
        """
        Save results to new or existing file
        """
        savefile = None
        savefile = asksaveasfilename(parent=self, title='Save LIE results', initialdir=self.app.workdir,
                                      filetypes=(("Q LIE", "*.qlie"), ("All files","*.*")),
                                      initialfile = 'results.qlie')

        if savefile:
            if not savefile.endswith('.qlie'):
                savefile += '.qlie'
                print(savefile)
            if os.path.isfile(savefile):
                oldfile = open(savefile, 'r').readlines()
                oldfile.append('\n')
            else:
                oldfile = []

            ligname = self.ligname.get().strip()
            alpha = float(self.alpha_entry.get())
            beta = float(self.beta_entry.get())
            gamma = float(self.gamma_entry.get())
            dg = float(self.dg_entry.get())
            dg_stderr = float(self.result_stderr_entry.get())
            dG_exp = float(self.dg_exp_entry.get())
            heading = '#%10s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n' % \
                          (ligname.ljust(10), 'el_q','vdW_q', 'el_w', 'vdW_w', 'el_p', 'vdW_p',
                           'el_s', 'vdW_s', 'alpha', 'beta','gamma', 'dG', 'dG_exp')
            oldfile.append(heading)

            lig_ene = '%11s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' % ('ligand'.ljust(11),
                            self.lig_ave[0][0], self.lig_ave[0][1], self.lig_ave[2][0], self.lig_ave[2][1], self.lig_ave[1][0],
                            self.lig_ave[1][1], self.lig_ave[3][0], self.lig_ave[3][1], alpha, beta, gamma, dg, dG_exp)
            oldfile.append(lig_ene)

            lig_stderr = '%11s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' % ('std.err'.ljust(11),
                            self.lig_stderr[0][0], self.lig_stderr[0][1], self.lig_stderr[2][0], self.lig_stderr[2][1],
                            self.lig_stderr[1][0], self.lig_stderr[1][1], self.lig_stderr[3][0], self.lig_stderr[3][1],
                            0.0, 0.0, 0.0, dg_stderr, 0.00)
            oldfile.append(lig_stderr)

            comp_ene = '%11s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' % ('complex'.ljust(11),
                            self.comp_ave[0][0], self.comp_ave[0][1], self.comp_ave[2][0], self.comp_ave[2][1], self.comp_ave[1][0],
                            self.comp_ave[1][1],self.comp_ave[3][0], self.comp_ave[3][1], alpha, beta, gamma, dg, dG_exp)
            oldfile.append(comp_ene)

            comp_stderr = '%11s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' % ('std.err'.ljust(11),
                            self.comp_stderr[0][0], self.comp_stderr[0][1], self.comp_stderr[2][0], self.comp_stderr[2][1],
                            self.comp_stderr[1][0],self.comp_stderr[1][1],
                            self.comp_stderr[3][0], self.comp_stderr[3][1], 0.0, 0.0, 0.0, dg_stderr, 0.00)
            oldfile.append(comp_stderr)

            output = open(savefile, 'w')
            for line in oldfile:
                output.write(line)
            output.close()

            self.app.log('info', 'LIE results saved to ../%s' % savefile.split('/')[-1])

        else:
            return

    def fitlie(self):
        """
        Opens up window to fit beta and gamma LIE parameters
        """
        self.fit_lie = FitQlie(self, self.root)
        self.fit_lie.resizable()
        self.fit_lie.configure(bg=self.main_color)

    def dialog_window(self):
        """
        Plot window
        """
        self.title('Analyze Lie')
        self.mainframe=Frame(self, bg=self.main_color)
        self.mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame with dirs and compute button (frame)
        frame = Frame(self.mainframe, bg=self.main_color)
        frame.grid(row=0, column=0, columnspan=2)


        #Frame with ligand results (frame2)
        frame2 = Frame(self.mainframe, bg=self.main_color)
        frame2.grid(row=1, column=0, sticky='nswe', pady=(5, 10))


        #Frame with complex results (frame3)
        frame3 = Frame(self.mainframe, bg=self.main_color)
        frame3.grid(row=2, column=0)


        #Frame with dG result and Close button (frame4)
        frame4 = Frame(self.mainframe, bg=self.main_color)
        frame4.grid(row=3, column=0, pady=10)


        #Fame with plotting functions (frame5)
        frame5 = Frame(self.mainframe, bg=self.main_color)
        frame5.grid(row=1, rowspan=3, column=1, padx=(10, 0))

        #Frame with save/close button
        frame6 = Frame(self.mainframe, bg=self.main_color)
        frame6.grid(row=4, column=0, columnspan=2)


        #General frame (frame)
        self.general_label = LabelFrame(frame, text='General', bg=self.main_color)
        self.general_label.grid()

        ligand_dir = Label(frame, text='Ligand dir.:', bg=self.main_color)
        ligand_dir.grid(in_=self.general_label, row=0, column=0, sticky='we')

        self.ligand_dir = Text(frame, width=15, height=1)
        self.ligand_dir.grid(in_=self.general_label, row=0, column=1)
        self.ligand_dir.config(highlightthickness=1)

        select_ligand_dir = Button(frame, text='Select', command=lambda: self.get_lie_files('ligand', True),
                                   highlightbackground=self.main_color)
        select_ligand_dir.grid(in_=self.general_label, row=0, column=2)

        self.ligand_runs_entry = Text(frame, width=30, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)


        complex_dir = Label(frame, text='Complex dir.:', bg=self.main_color)
        complex_dir.grid(in_=self.general_label, row=1, column=0, sticky='e')

        self.complex_dir = Text(frame, width=15, height=1)
        self.complex_dir.grid(in_=self.general_label, row=1, column=1)
        self.complex_dir.config(highlightthickness=1)

        select_complex_dir = Button(frame, text='Select', command=lambda: self.get_lie_files('complex', True),
                                    highlightbackground=self.main_color)
        select_complex_dir.grid(in_=self.general_label, row=1, column=2)

        self.complex_runs_entry = Text(frame, width=30, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)

        alpha_label = Label(frame, text= "\N{GREEK SMALL LETTER ALPHA}", bg=self.main_color)
        alpha_label.grid(in_=self.general_label, row=2, column=0, sticky='e')

        self.alpha_entry = Spinbox(frame,width=7, highlightthickness=0, relief=GROOVE,
                                   from_=0.00, to=1.50, increment=0.01)
        self.alpha_entry.grid(in_=self.general_label, row=2, column=1, sticky='w')

        beta_label = Label(frame, text= "\N{GREEK SMALL LETTER BETA}", bg=self.main_color)
        beta_label.grid(in_=self.general_label, row=3, column=0, sticky='e')

        self.beta_entry = Spinbox(frame,width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-1.00, to=1.00, increment=0.01)
        self.beta_entry.grid(in_=self.general_label, row=3, column=1, sticky='w')

        gamma_label = Label(frame, text= "\N{GREEK SMALL LETTER GAMMA}", bg=self.main_color)
        gamma_label.grid(in_=self.general_label, row=4, column=0, sticky='e')

        self.gamma_entry = Spinbox(frame,width=7, highlightthickness=0, relief=GROOVE,
                                   from_=-99.99, to=99.99, increment=0.01)
        self.gamma_entry.grid(in_=self.general_label, row=4, column=1, sticky='w')

        fit_button = Button(frame, text='Fit', highlightbackground=self.main_color, command=self.fitlie)
        fit_button.grid(in_=self.general_label, row=2, rowspan=3, column=2)

        self.compute_button = Button(frame, text='Compute LIE', highlightbackground=self.main_color, command=self.compute_lie)
        self.compute_button.grid(row=5, column=0, columnspan=3)
        self.compute_button.config(state=DISABLED)


        #Ligand frame (frame2)
        ligand_label = LabelFrame(frame2, text='Ligand', padx=10, bg=self.main_color)
        ligand_label.grid()

        ligand_el = Text(frame2, width=6, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        ligand_el.tag_configure("subscript", offset=-3)
        ligand_el.insert("insert",' <el>',"",'w','subscript')
        ligand_el.grid(in_=ligand_label, row=0, column=0)
        ligand_el.config(state=DISABLED)

        self.ligand_el_entry = Entry(frame2, width=10, highlightthickness=0)
        self.ligand_el_entry.grid(in_=ligand_label, row=0, column=1)
        self.ligand_el_entry.config(state=DISABLED)

        ligand_el_stderr = Label(frame2, text='stderr = ', bg=self.main_color)
        ligand_el_stderr.grid(in_=ligand_label, row=0, column=2)

        self.ligand_el_stderr_entry = Entry(frame2, width=10, highlightthickness=0)
        self.ligand_el_stderr_entry.grid(in_=ligand_label, row=0, column=3)
        self.ligand_el_stderr_entry.config(state=DISABLED)

        ligand_vdw = Text(frame2, width=6, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        ligand_vdw.tag_configure("subscript", offset=-3)
        ligand_vdw.insert("insert",'<vdW>',"",'w','subscript')
        ligand_vdw.grid(in_=ligand_label, row=1, column=0)
        ligand_vdw.config(state=DISABLED)

        self.ligand_vdw_entry = Entry(frame2, width=10, highlightthickness=0)
        self.ligand_vdw_entry.grid(in_=ligand_label, row=1, column=1)
        self.ligand_vdw_entry.config(state=DISABLED)

        ligand_vdw_stderr = Label(frame2, text='stderr = ', bg=self.main_color)
        ligand_vdw_stderr.grid(in_=ligand_label, row=1, column=2)

        self.ligand_vdw_stderr_entry = Entry(frame2, width=10, highlightthickness=0)
        self.ligand_vdw_stderr_entry.grid(in_=ligand_label, row=1, column=3)
        self.ligand_vdw_stderr_entry.config(state=DISABLED)

        #complex frame (frame3)
        complex_label = LabelFrame(frame3, text='Complex', padx=10, bg=self.main_color)
        complex_label.grid()

        complex_el_w = Text(frame3, width=6, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        complex_el_w.tag_configure("subscript", offset=-3)
        complex_el_w.insert("insert",' <el>',"",'w','subscript')
        complex_el_w.grid(in_=complex_label, row=0, column=0)
        complex_el_w.config(state=DISABLED)

        self.complex_el_w_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_el_w_entry.grid(in_=complex_label, row=0, column=1)
        self.complex_el_w_entry.config(state=DISABLED)

        complex_el_w_stderr = Label(frame3, text='stderr = ', bg=self.main_color)
        complex_el_w_stderr.grid(in_=complex_label, row=0, column=2)

        self.complex_el_w_stderr_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_el_w_stderr_entry.grid(in_=complex_label, row=0, column=3)
        self.complex_el_w_stderr_entry.config(state=DISABLED)

        complex_vdw_w = Text(frame3, width=6, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        complex_vdw_w.tag_configure("subscript", offset=-3)
        complex_vdw_w.insert("insert",'<vdW>',"",'w','subscript')
        complex_vdw_w.grid(in_=complex_label, row=1, column=0)
        complex_vdw_w.config(state=DISABLED)

        self.complex_vdw_w_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_vdw_w_entry.grid(in_=complex_label, row=1, column=1)
        self.complex_vdw_w_entry.config(state=DISABLED)

        complex_vdw_w_stderr = Label(frame3, text='stderr = ', bg=self.main_color)
        complex_vdw_w_stderr.grid(in_=complex_label, row=1, column=2)

        self.complex_vdw_w_stderr_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_vdw_w_stderr_entry.grid(in_=complex_label, row=1, column=3)
        self.complex_vdw_w_stderr_entry.config(state=DISABLED)

        complex_el_p = Text(frame3, width=6, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        complex_el_p.tag_configure("subscript", offset=-1)
        complex_el_p.insert("insert",' <el>',"",'p','subscript')
        complex_el_p.grid(in_=complex_label, row=2, column=0)
        complex_el_p.config(state=DISABLED)

        self.complex_el_p_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_el_p_entry.grid(in_=complex_label, row=2, column=1)
        self.complex_el_p_entry.config(state=DISABLED)

        complex_el_p_stderr = Label(frame3, text='stderr = ', bg=self.main_color)
        complex_el_p_stderr.grid(in_=complex_label, row=2, column=2)

        self.complex_el_p_stderr_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_el_p_stderr_entry.grid(in_=complex_label, row=2, column=3)
        self.complex_el_p_stderr_entry.config(state=DISABLED)

        complex_vdw_p = Text(frame3, width=6, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        complex_vdw_p.tag_configure("subscript", offset=-1)
        complex_vdw_p.insert("insert",'<vdW>',"",'p','subscript')
        complex_vdw_p.grid(in_=complex_label, row=3, column=0)
        complex_vdw_p.config(state=DISABLED)

        self.complex_vdw_p_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_vdw_p_entry.grid(in_=complex_label, row=3, column=1)
        self.complex_vdw_p_entry.config(state=DISABLED)

        complex_vdw_p_stderr = Label(frame3, text='stderr = ', bg=self.main_color)
        complex_vdw_p_stderr.grid(in_=complex_label, row=3, column=2)

        self.complex_vdw_p_stderr_entry = Entry(frame3, width=10, highlightthickness=0)
        self.complex_vdw_p_stderr_entry.grid(in_=complex_label, row=3, column=3)
        self.complex_vdw_p_stderr_entry.config(state=DISABLED)

        #Result of LIE (frame4)
        result_label = LabelFrame(frame4, text='Binding Free Energy', padx=10, bg=self.main_color)
        result_label.grid()

        dg_label = Label(frame4, text="\N{GREEK CAPITAL LETTER DELTA}G = ", bg=self.main_color)
        dg_label.grid(in_=result_label, row=0, column=0)

        self.dg_entry = Entry(frame4, width=10, highlightthickness=0)
        self.dg_entry.grid(in_=result_label, row=0, column=1)
        self.dg_entry.config(state=DISABLED)

        result_stderr = Label(frame4, text='stderr = ', bg=self.main_color)
        result_stderr.grid(in_=result_label, row=0, column=2)

        self.result_stderr_entry = Entry(frame4, width=10, highlightthickness=0)
        self.result_stderr_entry.grid(in_=result_label, row=0, column=3)
        self.result_stderr_entry.config(state=DISABLED)

        #Plot frame (frame5)
        plot_label = LabelFrame(frame5, text='Analyze', padx=10, bg=self.main_color)
        plot_label.grid()

        system_label = Label(frame5, text='Select system(s):', bg=self.main_color)
        system_label.grid(in_=plot_label, row=0, column=0, columnspan=2)

        self.system_list = Listbox(frame5, width=22, height=2, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        self.system_list.grid(in_=plot_label, row=1, column=0, columnspan=2, sticky='w')
        self.system_list.insert(END, 'LIGAND', 'COMPLEX')
        self.system_list.config(font=tkinter.font.Font(family="Courier", size=12))

        interactions_label = Label(frame5, text='Select interaction(s):', bg=self.main_color)
        interactions_label.grid(in_=plot_label, row=3, column=0, columnspan=2)

        self.interaction_list = Listbox(frame5, width=22, height=8, highlightthickness=0, relief=GROOVE,
                                        selectmode=EXTENDED, exportselection=False)
        self.interaction_list.grid(in_=plot_label, row=4, column=0, columnspan=2, sticky='w')
        interactions = ['Q-Q    el', 'Q-Q    vdW','Q-prot el', 'Q-prot vdW', 'Q-wat  el',
                                     'Q-wat  vdW', 'Q-surr el', 'Q-surr vdW']
        for qtype in interactions:
            self.interaction_list.insert(END,'%10s%10s' % (qtype.split()[0].ljust(10),qtype.split()[1]))
        self.interaction_list.config(font=tkinter.font.Font(family="Courier", size=12))

        self.plot_button = Button(frame5, text='Plot', command=self.plot_selection, highlightbackground=self.main_color)
        self.plot_button.grid(in_=plot_label, row=10, column=0, columnspan=2)

        #Frame 6: Save/Close
        ligname = Label(frame6, text='#Ligand name:', bg=self.main_color)
        ligname.grid(row=0, column=0)
        self.ligname = Entry(frame6, width=10, highlightthickness=0)
        self.ligname.grid(row=0, column=1)
        self.ligname.insert(0, 'LIG')

        dg_exp = Label(frame6, text="\N{GREEK CAPITAL LETTER DELTA}G(exp):", bg=self.main_color)
        dg_exp.grid(row=1, column=0, sticky='e')
        self.dg_exp_entry = Entry(frame6, width=10, highlightthickness=0)
        self.dg_exp_entry.grid(row=1, column=1)
        self.dg_exp_entry.insert(0, '0.00')

        self.save_button = Button(frame6, text='Save', highlightbackground=self.main_color, command=self.save_results)
        self.save_button.grid(row=0, rowspan=2, column=3, sticky='w')
        self.save_button.config(state=DISABLED)
        close_button = Button(frame6, text='Close', highlightbackground=self.main_color, command=self.destroy)
        close_button.grid(row=3, column=4, sticky='e',pady=(10,0), padx=(20,0))
