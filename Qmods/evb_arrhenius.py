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

from tkinter import Label, TOP, Button, Listbox, Scrollbar, EXTENDED, Spinbox, Entry, Text, Frame, \
    Toplevel, DISABLED, END, GROOVE, NORMAL, BOTH, IntVar, StringVar, Checkbutton, OptionMenu, HORIZONTAL
from tkinter.filedialog import askopenfilename
from tkinter.simpledialog import askstring
import tkinter.font
from subprocess import call
from cycler import cycler
import numpy as np
import os
import shutil
import matplotlib
matplotlib.use('TkAgg')
#Implement default mpl key bindings
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from tkinter.filedialog import asksaveasfilename, askdirectory
matplotlib.rcParams['text.usetex'] = True

class EvbArrhenius(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root

        #This can be used to make global designs, ratter than locals:
        #---> self.root.option_add('*background', self.main_color)

        self.dg_plot = None

        self.selected_frame = StringVar()

        #{'name': {T1: {1: path, 2:path ....}}
        self.titles = dict()

        #{name: T1: {1: [dG_act, dG_react]}, 2...}
        self.titles_dG = dict()

        #{name: {T1:[<dG>, SEM]}}
        self.titles_ave_act = dict()

        #{name: {T1:[<dG>, SEM]}}
        self.titles_ave_rxn = dict()

        #{name: [Temperature]}
        self.title_exclude = dict()

        #{name: {dH, dS, COD}}
        self.titles_parameters = dict()

        #{name: {temp: {act/rxn: [upper, lower]}}}
        self.dg_upper_lower = dict()

        self.dialog_window()

        self.selected_frame.set('Show ...')
        self.selected_frame.trace('w', self.show_var_frame)

    def add_title(self):
        """
        add a new project title (main name for plot)
        """

        title = askstring('Add title', 'EVB project title:', parent=self)
        if not title:
            return
        if title in self.titles_listbox.get(0, END):
            self.app.errorBox('Warning', 'Title already exist in Titles.')
            return

        self.titles_listbox.insert(END, title)

        for i in [self.titles, self.titles_dG, self.titles_parameters, self.titles_ave_act, self.titles_ave_rxn,
                  self.dg_upper_lower]:
            i[title] = dict()

        #Highlight latest title in listbox
        self.titles_listbox.select_clear(0, END)
        titles = self.titles_listbox.get(0, END)
        for ind_ in range(len(titles)):
            if titles[ind_] == title:
                self.titles_listbox.selection_set(ind_)

    def del_title(self):
        selections = list(map(int, self.titles_listbox.curselection()))
        for selected in selections:
            title = self.titles_listbox.get(selected)
            self.titles_listbox.delete(selected)
            del self.titles[title]
            for i in [self.titles_dG, self.titles_parameters, self.titles_ave_act, self.titles_ave_rxn,
                      self.dg_upper_lower]:
                if title in list(i.keys()):
                    del i[title]

        self.update_tables()

    def add_runs(self):
        """
        Add Temperature one-by-one or select temperature directory to append all runs
        """
        try:
            title = list(map(int, self.titles_listbox.curselection()))
        except:
            self.app.errorBox('Warning', 'Select a project title to append temperature(s) to.')
            return

        if len(title) != 1:
            self.app.errorBox('Warning', 'Select exactly one project title to append temperature(s) to.')
            return

        title = self.titles_listbox.get(title[0])

        #Start appending directories to collect energy files from
        temperatures = list()
        title_dialog = 'Select temperature directory'
        path = self.app.workdir
        rundir = askdirectory(parent=self, mustexist=True, title=title_dialog, initialdir= path)
        if not rundir:
            return

        #Check of a single temperature is selected, or directory with several temperetures
        try:
            temp = float(rundir.split('/')[-1])
            if temp < 270:
                if len(os.listdir(rundir)) > 3:
                    for t in os.listdir(rundir):
                        try:
                            temp = float(t.split('/')[-1])
                            if temp > 270.0:
                                temperatures.append('%s/%s' % (rundir, t))
                        except:
                            continue
                else:
                    print(('Directory ../%s not recognized as temperature directory' % '/'.join(rundir.split('/')[-3:])))
                    return
            temperatures.append(rundir)

        except:
            for t in os.listdir(rundir):
                try:
                    temp = float(t.split('/')[-1])
                    if temp > 270.0:
                        temperatures.append('%s/%s' % (rundir, t))
                except:
                    continue

        if len(temperatures) == 0:
            self.app.errorBox('Error', 'Did not recognize selected temperature directories.')
            return

        print(temperatures)

        #Go throught temperature directories and collect subdir paths
        for t in temperatures:
            #Collect run dirs and break
            temp = float(t.split('/')[-1])


            nr_dir = 1
            while True:
                subdir = '%s/%d' % (t, nr_dir)
                if os.path.isdir(subdir):
                    enefiles = False
                    for ene in os.listdir(subdir):
                        if ene.endswith('.en'):
                            if os.path.getsize('%s/%s' % (subdir, ene)) != 0:
                                enefiles = True
                                break
                            else:
                                break
                    if enefiles:
                        if not temp in list(self.titles[title].keys()):
                            self.titles[title][temp] = dict()
                            self.titles_dG[title][temp] = dict()
                            self.titles_ave_act[title][temp] = dict()
                            self.titles_ave_rxn[title][temp] = dict()
                            self.dg_upper_lower[title][temp] = dict()
                            self.dg_upper_lower[title][temp]['activation'] = [-999999, 999999]
                            self.dg_upper_lower[title][temp]['reaction'] = [-999999, 9999999]

                        self.titles[title][temp][nr_dir] = subdir
                        self.app.log(' ','...../%s added\n' % '/'.join(subdir.split('/')[-2:]))
                        self.update()
                        print(subdir)
                    nr_dir += 1
                else:
                    break

        #Get all dG(act) dG(react) and normalize
        self.get_dG(title)

        #Genereate themodynamic parameters
        self.getParameters(title)

        #update tables
        self.update_tables()

        self.list_titles_event()

    def del_runs(self):
        title_sel = list(map(int, self.titles_listbox.curselection()))
        if len(title_sel) != 1:
            self.app.log(' ', '\nSelect exactly one title to delete temperature from!\n')
            return

        selections = list(map(int, self.runs_listbox.curselection()))
        if len(selections) == 0:
            self.app.log(' ', '\nNo temperature selected for deletion!\n')
            return

        title = self.titles_listbox.get(title_sel[0])

        for selected in reversed(selections):
            temp = self.runs_listbox.get(selected)
            if '.' in temp:
                temp = float(temp)
            else:
                temp = int(temp)
            self.runs_listbox.delete(selected)

            del self.titles[title][temp]
            for i in [self.titles_dG, self.titles_parameters, self.titles_ave_act, self.titles_ave_rxn]:
                if temp in list(i[title].keys()):
                    del i[title][temp]

        #self.get_dG(title)

        self.getParameters(title)

        self.update_tables()
        self.list_titles_event()

    def add_exluded_temp(self):
        """
        Add temperatures to exclude from Arrhenius calculations for title
        """
        sel_titles = list(map(int, self.titles_listbox.curselection()))
        if len(sel_titles) > 1:
            print('Select exactly one title to exclude temperatures')
            return
        title = self.titles_listbox.get(sel_titles[0])

        sel_temp = list(map(int, self.runs_listbox.curselection()))
        if len(sel_temp) == 0:
            'No temperatures selected for %s' % title
            return

        for t in sel_temp:
            temp = self.runs_listbox.get(t)
            if '.' in temp:
                temp = float(temp)
            else:
                temp = int(temp)

            if title not in list(self.title_exclude.keys()):
                self.title_exclude[title] = list()

            if temp not in self.title_exclude[title]:
                self.title_exclude[title].append(temp)
                self.ae_nb_listbox.insert(END, temp)

        #Genereate themodynamic parameters
        self.getParameters(title)

        #update tables
        self.update_tables()

        self.list_titles_event()

    def del_excluded_temp(self):
        title_sel = list(map(int, self.titles_listbox.curselection()))
        if len(title_sel) != 1:
            self.app.log(' ', '\nSelect exactly one title to delete excluded temperature from!\n')
            return

        title = self.titles_listbox.get(title_sel[0])

        selection = list(map(int, self.exclude_listbox.curselection()))
        for i in selection:
            temp = self.exclude_listbox.get(i)
            if '.' in temp:
                temp = float(temp)
            else:
                temp = int(temp)

            del self.title_exclude[title][self.title_exclude[title].index(temp)]

        #Genereate themodynamic parameters
        self.getParameters(title)

        #update tables
        self.update_tables()

        self.list_titles_event()

    def write_qfep_inp(self, path):
        """
        Write qfep.inp in path for md*.en
        """
        found_ene_files = True
        #Collect all md*.en files:
        enefiles = list()
        for qfile in os.listdir(path):
            if qfile.endswith('.en'):
                if qfile.startswith('md'):
                    enefiles.append(qfile)

        print(('Found %d MD energy files in ../%s' % (len(enefiles), '/'.join(path.split('/')[-3:]))))

        if len(enefiles) == 0:
            found_ene_files = False
            return found_ene_files

        if not self.kT.get()[0].isdigit():
            try:
                T = float(path.split('/')[-2])
            except:
                print('Could not recognize temperature, using 300 K')
                T = 300.00
        else:
            T = float(self.kT.get())

        kT = 0.001987209 * T

        #Sort energy files (if not Qfep will give messed up free energy profile!)
        enefiles = sorted(enefiles, key=lambda x: x.split('/')[-1])

        #make input for Qfep
        inputfile = open('%s/qfep.inp' % path,'w')
        inputfile.write('%d\n' % len(enefiles))
        inputfile.write('2  0\n')
        inputfile.write('%.4f  %s\n' % (kT, self.skip_points.get()))
        inputfile.write('%s\n' % self.bins.get())
        inputfile.write('%s\n' % self.binpoints_min.get())
        inputfile.write('%s\n' % self.alpha_entry.get())
        inputfile.write('1\n')
        inputfile.write('1  2  %s  0.0 0 0.000\n' % self.hij_entry.get())
        inputfile.write('%s\n' % self.linear_comb.get())

        for ene in enefiles:
            inputfile.write('%s\n' % ene.split('/')[-1])

        inputfile.close()

        return found_ene_files

    def run_qfep(self, path, qfep_inp):
        """
        Calls Qfep with qfep_inp and produces qfep.out
        """
        #Get Qfep executable
        qfep = self.app.q_settings[ 'executables' ][2]
        #Move to path and run Qfep
        self.app.log(' ','Running Qfep in ../%s\n' % '/'.join(path.split('/')[-3:]))
        os.chdir(path)

        tmpfile = open('qfep.out', 'w')
        call('%s <%s' % (qfep, qfep_inp), shell=True, stdout=tmpfile, stderr=tmpfile)

        tmpfile.close()
        os.chdir(self.app.workdir)

    def get_enegaps_dg(self, path):
        """
        Read qfepout #part3 and returns 2 lists: Energy gaps and delta G.
        """
        #Get enegap and dG
        part3 = False
        count = 0
        enegaps = list()
        dG = list()
        qfepout = '%s/qfep.out' % path

        #Make sure that qfep.out exist in dir. If not, run Qfep
        if not os.path.isfile(qfepout):
            if not os.path.isfile('%s/qfep.inp' % path):
                self.app.log(' ', '\n#####\nqfep.inp does not exist in %s. New input generated!\n####\n'
                                  % '/'.join(path.split('/')[-3:]))
                self.write_qfep_inp(path)
            self.run_qfep(path, 'qfep.inp')

        if os.path.isfile(qfepout):
            with open(qfepout, 'r') as qfepout:
                for line in qfepout:
                    if '# Part 3: Bin-averaged summary' in line:
                        part3 = True
                    if part3:
                        if line == '' or "# Part 4:" in line:
                            break
                        if count > 1:
                            try:
                                enegaps.append(float(line.split()[1]))
                                dG.append(float(line.split()[3]))
                            except:
                                continue
                        else:
                            count += 1
        else:
            self.app.log(' ', '\n\nWARNING! Could not create qfep.out in %s\n\n' % '/'.join(path.split('/')[-3:]))

        return enegaps, dG

    def import_qfep(self):
        """
        Import EVB parameters from qfep file
        """
        filename = None
        filename = askopenfilename(parent = self, initialdir = self.app.workdir,
                                   filetypes=(("Qfep", "*.qfep"),("All files","*.*")))
        if not filename:
            return

        entries = [[None], [None], [self.kT, self.skip_points], [self.bins], [self.binpoints_min], [self.alpha_entry],
                           [None], [None, None, self.hij_entry], [self.linear_comb]]

        qfep = open(filename, 'r').readlines()
        try:
            for i in range(len(qfep)):
                if 2 < i < 8:
                    for j in range(len(entries[i])):
                        entry = entries[i][j]
                        if entry:
                            val = qfep[i].split()[j]
                            entry.delete(0, END)
                            entry.insert(0, val)
                if i == 8:
                    entry = entries[i][0]
                    val = qfep[i].split('\n')[0].strip()
                    entry.delete(0, END)
                    entry.insert(0, val)
                if i > 8:
                    break
        except:
            self.app.errorBox('Error', 'Could not read file!')
            return

    def import_project(self):
        """
        Import saved project
        """
        self.app.errorBox('Info', 'This function is not available at the moment')
        return

    def on_key_event(self, event):
        print(('you pressed %s' % event.key))
        key_press_handler(event, self.canvas, self.toolbar)

    def update_plot(self, t_inv, dg_t, sem, parameters, titles, t_exc=list(), dg_exc=list(), sem_exc=list()):
        if len(titles) < 1:
            return

        if self.dg_plot:
            self.dg_plot.clear()

        else:
            #Set color cycle for plots:
            #matplotlib.rcParams['axes.color_cycle'] = ['k', 'b', 'g', 'r', 'm', 'y', 'c', 'brown',
            #                                           'burlyWood', 'cadetBlue', 'DarkGreen', 'DarkBlue',
            #                                           'DarkMagenta', 'DarkSalmon', 'DimGray', 'Gold']
            matplotlib.rcParams['axes.prop_cycle'] = cycler(color=['k','b', 'g', 'r', 'm', 'y', 'c', 'brown',
                                                       'burlyWood', 'cadetBlue', 'DarkGreen', 'DarkBlue',
                                                       'DarkMagenta', 'DarkSalmon', 'DimGray', 'Gold'])


            #Create subplot
            self.dg_plot = self.plot_window.add_subplot(111, facecolor='white')
            self.plot_window.subplots_adjust(hspace=0.5)

            #X/Y labels
            self.dg_plot.set_xlabel(r'$1/T$')
            self.dg_plot.set_ylabel(r'$\Delta G^{\ddag}/T$')

            #Move label box outside plot region
            box = self.dg_plot.get_position()
            self.dg_plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            #Fit subplot to figure/canvas
            #rect=(left,bottom,top,right)
            self.plot_window.tight_layout(rect=(0.005, 0, 0.8, 1))

        temp_labels = list()
        for i in range(len(titles)):
            cycle = matplotlib.rcParams['axes.prop_cycle']
            line_color = list(matplotlib.rcParams['axes.prop_cycle'])[i%len(cycle)]['color']
            print(line_color)
            title = titles[i]
            dg = np.array(dg_t[i])
            temp = np.array(t_inv[i])
            error = np.array(sem[i])
            dh = parameters[i][0]
            ds = parameters[i][1]

            #genereate temperature labels
            all_temps = temp
            if len(t_exc[i]) > 0:
                all_temps = np.array(sorted(t_inv[i] + t_exc[i], reverse=True))

            for inv in all_temps:
                t_real = 1. / inv
                temp_labels.append(r'$%.0f^{-1}$' % t_real)

            #plt.xticks([list of tick locations], [list of tick lables])
            self.dg_plot.errorbar(temp, dg, error, fmt='o', color=line_color, linestyle='None')
            if len(t_exc) > 0:
                if len(t_exc[i]) > 0:
                    exc_t = np.array(t_exc[i])
                    exc_dg = np.array(dg_exc[i])
                    exc_error = np.array(sem_exc[i])
                    self.dg_plot.errorbar(exc_t, exc_dg, exc_error, fmt='x', color=line_color, linestyle='None')

            self.dg_plot.plot(all_temps, (dh * all_temps) - ds, '-', linewidth=2.0, color=line_color, label=title)
            self.dg_plot.autoscale(enable=True)

            #Set new xtick labels
            self.dg_plot.axes.set_xticks(all_temps)
            self.dg_plot.axes.set_xticklabels(temp_labels)

            self.dg_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})

        self.canvas.draw()

    def recomp_evb(self):
        """
        Recomputes EVB with current settings (alpha, Hij etc.) from main window.
        """
        selections = list(map(int, self.titles_listbox.curselection()))
        if len(selections) == 0:
            self.app.errorBox('Warning', 'Select titles to recompute EVB.')
            return
        self.titles_listbox.select_clear(0, END)

        titles = list()
        for selected in selections:
            titles.append(self.titles_listbox.get(selected))

        for title in titles:
            ##{name: T1: {1: [dG_act, dG_react]}, 2...}
            for temp in sorted(self.titles[title].keys()):
                for i in sorted(self.titles[title][temp].keys()):
                    path = self.titles[title][temp][i]
                    self.write_qfep_inp(path)
                    self.run_qfep(path, 'qfep.inp')

            self.get_dG(title)
            self.getParameters(title)

        self.update_tables()
        self.list_titles_event()

    def get_dG(self, title):
        """
        Collects activation and reaction energies from self.title_ene ({title: {nr:[[enegaps],[dG]]}})
        """

        for temp in list(self.titles[title].keys()):
            act = list()
            rxn = list()
            for nr in list(self.titles[title][temp].keys()):
                self.app.log(' ', 'Collecting Reaction Free Energies for %s: %.2f/%d\n' % (title, temp, nr))
                self.update()
                dG = self.get_enegaps_dg(self.titles[title][temp][nr])[1]
                rsE, tsE, rE = self.find_dg(dG)
                if tsE != 'na' and rE != 'na':
                    self.titles_dG[title][temp][nr] = [tsE, rE]
                    act.append(tsE)
                    rxn.append(rE)

                    #Update upper_lower:
                    if tsE > self.dg_upper_lower[title][temp]['activation'][0]:
                        self.dg_upper_lower[title][temp]['activation'][0] = tsE
                    if tsE < self.dg_upper_lower[title][temp]['activation'][1]:
                        self.dg_upper_lower[title][temp]['activation'][1] = tsE

                    if rE > self.dg_upper_lower[title][temp]['reaction'][0]:
                        self.dg_upper_lower[title][temp]['reaction'][0] = rE
                    if rE < self.dg_upper_lower[title][temp]['reaction'][1]:
                        self.dg_upper_lower[title][temp]['reaction'][1] = rE

            self.titles_ave_act[title][temp] = [np.average(act), np.std(act) / np.sqrt(len(act))]
            self.titles_ave_rxn[title][temp] = [np.average(rxn), np.std(rxn) / np.sqrt(len(rxn))]

    def compute_ave_dg(self, title, temp):
        """
        This is called whenever the user has modified the upper and lower limits of dG.
        """
        act = list()
        rxn = list()

        act_upper = self.dg_upper_lower[title][temp]['activation'][0]
        act_lower = self.dg_upper_lower[title][temp]['activation'][1]
        rxn_upper = self.dg_upper_lower[title][temp]['reaction'][0]
        rxn_lower = self.dg_upper_lower[title][temp]['reaction'][1]

        for nr in sorted(self.titles_dG[title][temp].keys()):
            tsE, rE = self.titles_dG[title][temp][nr][:]

            #Will only add energies if both activation and reaction energies are withing thresholds
            if (tsE < act_upper) and (tsE > act_lower):
                if (rE < rxn_upper) and (rE > rxn_lower):
                    act.append(tsE)
                    rxn.append(rE)

        self.titles_ave_act[title][temp] = [np.average(act), np.std(act) / np.sqrt(len(act))]
        self.titles_ave_rxn[title][temp] = [np.average(rxn), np.std(rxn) / np.sqrt(len(rxn))]

        #Genereate themodynamic parameters
        self.getParameters(title)

        #update tables
        self.update_tables()

        self.list_titles_event()

    def set_upper_lower(self):
        """
        Takes upper and lower values for dG(act) and dG(rxn) and recomputes parameters.
        """
        selections = list(map(int, self.dg_listbox.curselection()))

        if len(selections) != 1:
            return

        try:
            #title = self.dg_listbox.get(selections[0]).split()[0]
            #temp = int(self.dg_listbox.get(selections[0]).split()[1])
            title = self.titles_listbox.get(int(self.titles_listbox.curselection()[0]))
            temp = int(self.dg_listbox.get(selections[0]).split()[1 + (len(title.split()) - 1)])
        except:
            print('Returned error in "set_upper_lower" function.')
            return

        try:
            act_upper = float(self.dg_act_upper.get())
            act_lower = float(self.dg_act_lower.get())
            rxn_upper = float(self.dg_rxn_upper.get())
            rxn_lower = float(self.dg_rxn_lower.get())
        except:
            print('Invalid value encountered. Unable to update thermodynamic parameters!')
            return

        self.dg_upper_lower[title][temp]['activation'][0] = act_upper
        self.dg_upper_lower[title][temp]['activation'][1] = act_lower
        self.dg_upper_lower[title][temp]['reaction'][0] = rxn_upper
        self.dg_upper_lower[title][temp]['reaction'][1] = rxn_lower
        print('New upper and lower values set!')

        self.compute_ave_dg(title, temp)


    def find_dg(self,dG):
        """
        Finds dG(act) and dG(react) in in self.dG_ene[title][nr]
        """

        if len(dG) < 1:
            return 0, 'na', 'na'
        rsE = dG[0]
        tsE = 'na'
        rE = 'na'

        found_TS = False
        found_RE = False

        for i in range(1, len(dG) - 1):
            prev = float(dG[i - 1])
            mid = float(dG[i])
            nxt = float(dG[i + 1])

            if mid <= prev and mid <= nxt:
                if not found_TS:
                    if mid < rsE:
                        rsE = mid
                elif found_TS:
                    if not found_RE:
                        rE = mid
                        found_RE = True
                    elif found_RE:
                        if mid < rE:
                            rE = mid

            elif mid >= prev and mid >= nxt:
                if not found_TS:
                    if abs(mid) > 0:
                        tsE = mid
                        found_TS = True
                elif found_TS:
                    if mid > tsE:
                        tsE = mid
                        found_RE = False
                        rE = dG[-1]

        try:
            tsE = (float(tsE) - float(rsE))
        except:
            tsE = 'na'

        try:
            rE = (float(rE) - float(rsE))

        except:
            rE = 'na'

        return rsE, tsE, rE

    def getParameters(self, title):
        """
        Returns the thermodynamic activation parameters (dH and dS)
        from linear regresion of dG(1/T) = dH*1/T - dS

        """
        dG = list()
        t = list()
        dG_T = list()
        T_inv = list()

        for temp in sorted(self.titles_ave_act[title].keys()):
            excluded = False
            if title in list(self.title_exclude.keys()):
                if temp in self.title_exclude[title]:
                    excluded = True
            if not excluded:
                dG_T.append(self.titles_ave_act[title][temp][0]/ float(temp))
                T_inv.append(1.00 / temp)
                dG.append(self.titles_ave_act[title][temp][0])
                t.append(temp)

        if len(t) < 2:
            print('Need more than 1 temperature to extract parameters ..')
            self.titles_parameters[title]['dH'] = 0
            self.titles_parameters[title]['dS'] = 0
            self.titles_parameters[title]['COD'] = 0
            self.titles_parameters[title]['s_dg'] = 0
            self.titles_parameters[title]['s_dh'] = 0
            self.titles_parameters[title]['s_ds'] = 0
            self.titles_parameters[title]['s_r'] = 0
            self.titles_parameters[title]['SEM_dg'] = 0
            self.titles_parameters[title]['SEM_dh'] = 0
            self.titles_parameters[title]['SEM_ds'] = 0
            self.titles_parameters[title]['sample'] = len(dG)
            return

        _y = np.array(dG_T)
        _x = np.array(T_inv)
        _A = np.vstack([_x,np.ones(len(_x))]).T        #Vector matrix for linear algebra least sq.

        model, resid = np.linalg.lstsq(_A, _y)[:2]     #generate model and residuals
        r2 = 1 - resid / (_y.size * _y.var())         # R^2

        #Remove this:
        print(model)
        print(r2)

        dH = model[0]
        dS = -model[1]

        try:
            r2 = r2[0]
        except:
            pass


        s2_dg = 0
        s2_dg2 = 0
        sst = 0
        sst2 = 0
        t_ave = np.average(t)
        #Get residuals squared
        #s2_dg is SSR (sum squared residuals)
        for i in range(len(dG)):
            temp = T_inv[i] #t[i]
            dG_i = dG_T[i] #dG[i]
            dG_m = (dH * temp) - dS
            s2_dg += (dG_i - dG_m) ** 2

            #This is just to compute the actual standard error for dG = dH - TdS instead of (dG/T = dH/T - dS)
            dG_i2 = dG[i]
            dG_m2 = dH - (t[i] * dS)
            s2_dg2 += (dG_i2 - dG_m2) ** 2

            #print 'dG_i = %.5f    dG_m = %.5f' % (dG_i, dG_m)
            print(('dG_i = %.5f    dG_m = %.5f' % (dG_i2, dG_m2)))
            #Sum of squares for the temperature
            sst += (temp - np.average(T_inv)) ** 2
            sst2 += (t[i] - np.average(t)) ** 2
            # sst += (temp - t_ave) ** 2

        #standard error squared of the model (remove 2 degrees of freedom <-- linear regression parameters)
        deg_freedom = 2
        if len(dG) < 3:
            print('Warning: linear regression with only 2 points')
            print('--> not possible to remove 2 degrees of freedom!')
            s2_dg = 0
            s2_dg2 = 0

            #Standard error of the model
            s_dg = 0
            s_dg2 = 0

            #Standard error of the mean for the model
            SEM_dg = np.sqrt(s_dg2) / np.sqrt(len(dG))

            #Standard error for the slope (dH)
            s_dh = 0

            #Standar error of the mean for the slope (dH)
            SEM_dh = 0

            #Standard error of the intercept (dS)
            s_ds = 0
            SEM_ds = 0

            #Standard error for the correlation coefficient
            s_r = 0
            r2 = 1

        else:
            s2_dg /= (float(len(dG)) - deg_freedom)
            s2_dg2 /= (float(len(dG)) - deg_freedom)

            #Standard error of the model
            s_dg = np.sqrt(s2_dg)
            s_dg2 = np.sqrt(s2_dg2)

            #Standard error of the mean for the model
            SEM_dg = np.sqrt(s_dg2) / np.sqrt(len(dG))

            #Standard error for the slope (dH)
            s_dh = (dH / np.sqrt(len(dG) - deg_freedom)) * np.sqrt((1.0 / r2) - 1.0)

            #Standar error of the mean for the slope (dH)
            SEM_dh = s_dh / np.sqrt(len(dG))

            #Standard error of the intercept (dS)
            s_ds = s_dg * np.sqrt((1.0 / len(dG)) + (np.average(T_inv) ** 2 / sst))
            SEM_ds = s_ds / np.sqrt(len(dG))

            #Standard error for the correlation coefficient
            s_r = np.sqrt((1.0 - r2) / (len(dG) - deg_freedom))
            r = np.sqrt(r2)

        self.titles_parameters[title]['dH'] = dH
        self.titles_parameters[title]['dS'] = dS
        self.titles_parameters[title]['COD'] = r2
        self.titles_parameters[title]['s_dg'] = s_dg2
        self.titles_parameters[title]['s_dh'] = s_dh
        self.titles_parameters[title]['s_ds'] = s_ds
        self.titles_parameters[title]['s_r'] = s_r
        self.titles_parameters[title]['SEM_dg'] = SEM_dg
        self.titles_parameters[title]['SEM_dh'] = SEM_dh
        self.titles_parameters[title]['SEM_ds'] = SEM_ds
        self.titles_parameters[title]['sample'] = len(dG)

        print((self.titles_parameters))

    def update_tables(self):
        """
        Updates all listboxes
        """
        self.update_dg_listbox()
        self.update_re_listbox()
        self.update_model_computed_listbox()

    def update_dg_listbox(self):
        """
        Updates self.dg_listbox with dG(act), dG(rxn) and corresponding SEM
        """
        self.dg_listbox.delete(0, END)

        self.dg_listbox.insert(END, "Title       Temp "
                                       " \N{GREEK CAPITAL LETTER DELTA}G(act)"
                                       "   +/-"
                                       "   \N{GREEK CAPITAL LETTER DELTA}G(rxn)"
                                       "   +/-")

        for title in list(self.titles_dG.keys()):
            for temp in sorted(self.titles_dG[title].keys()):
                act, act_sem = self.titles_ave_act[title][temp][0:]
                rxn, rxn_sem = self.titles_ave_rxn[title][temp][0:]
                self.dg_listbox.insert(END, '%10s %5.0f %7.2f %7.2f %7.2f %7.2f' % (title.ljust(10), temp,
                                                                                    act, act_sem,
                                                                                    rxn, rxn_sem))

    def update_re_listbox(self):
        self.re_tot_listbox.delete(0, END)

        self.re_tot_listbox.insert(END, "Title      "
                                       "\N{GREEK CAPITAL LETTER DELTA}H(act)"
                                       "  \N{GREEK CAPITAL LETTER DELTA}S(act)   COD"
                                       "  SE(\N{GREEK CAPITAL LETTER DELTA}G)"
                                       "  SE(\N{GREEK CAPITAL LETTER DELTA}H)"
                                       "  SE(\N{GREEK CAPITAL LETTER DELTA}S)"
                                       "  s(R)  #pts")

        for title in list(self.titles_parameters.keys()):
            try:
                dH = self.titles_parameters[title]['dH']
                dS = self.titles_parameters[title]['dS']
                COD = self.titles_parameters[title]['COD']
                s_dg = self.titles_parameters[title]['s_dg']
                s_dh = self.titles_parameters[title]['s_dh']
                s_ds = self.titles_parameters[title]['s_ds']
                s_dr = self.titles_parameters[title]['s_r']
                pts = self.titles_parameters[title]['sample']
                self.re_tot_listbox.insert(END, '%10s %6.2f  %8.5f  %5.2f %5.2f   %5.2f %8.5f %5.2f %3d' %
                                                (title.ljust(10), dH, dS, COD, s_dg, s_dh, s_ds, s_dr, pts))
            except:
                continue

    def update_model_computed_listbox(self):
        """
        Computed dG from the arrhenius parameters with SEM for each temperature and compare to average dG from EVB
        """
        self.ae_nb_listbox.delete(0, END)
        self.ae_nb_listbox.insert(END, "Title         T    <\N{GREEK CAPITAL LETTER DELTA}G>    +/- "
                                       "    \N{GREEK CAPITAL LETTER DELTA}H    +/-    "
                                       "T\N{GREEK CAPITAL LETTER DELTA}S "
                                       "   +/-    \N{GREEK CAPITAL LETTER DELTA}G     +/-   "
                                       " \N{GREEK CAPITAL LETTER DELTA}\N{GREEK CAPITAL LETTER DELTA}G")
        for title in list(self.titles_ave_act.keys()):
            dH = self.titles_parameters[title]['dH']
            dS = self.titles_parameters[title]['dS']
            exc_T = 0
            if title in list(self.title_exclude.keys()):
                exc_T = len(self.title_exclude[title])

            sem_ds = self.titles_parameters[title]['s_ds']/ np.sqrt(len(self.titles_ave_act[title]) - exc_T)
            sem_dh = self.titles_parameters[title]['s_dh']/ np.sqrt(len(self.titles_ave_act[title]) - exc_T)
            for temp in sorted(self.titles_ave_act[title].keys()):
                dg_arr = dH - (temp * dS)
                #s_(x+y)^2 = s_x^2 + s_y^2 --> s_(x+y) = sqrt(s_x^2 + s_y^2)
                sem_dg_arr = float(np.sqrt(sem_dh**2 + (temp * sem_ds)**2))
                dg, sem_dg = self.titles_ave_act[title][temp][0:]
                ddg = dg - dg_arr
                self.ae_nb_listbox.insert(END, '%10s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f' %
                                               (title.ljust(10), temp, dg,sem_dg, dH, sem_dh, temp * dS, temp * sem_ds,
                                                dg_arr, sem_dg_arr, ddg))

    def get_qfep_part0(self, path, filename):
        """
        Reads qfep.out and returns Part 0 as a list (part0 = [[state1],[state2]])
        """
        part0 = list()
        found_part0 = False

        qfepout = '%s/%s' % (path, filename)
        tmp = list()
        with open(qfepout, 'r') as qfepout:
            for line in qfepout:
                if found_part0:
                    if len(line.split()) < 2 or '# Part 1:' in line:
                        break
                    try:
                        if int(line.split()[1]) == 1 or int(line.split()[1]) == 2:
                            l = line[32:40]
                            eqtot = float(line[40:48])
                            eqbnd = float(line[48:56])
                            eqang = float(line[56:64])
                            eqtor = float(line[64:72])
                            eqimp = float(line[72:80])
                            eqel = float(line[80:88])
                            eqvdw = float(line[88:96])
                            qqel = float(line[96:104])
                            qqvdw = float(line[104:112])
                            qpel = float(line[112:120])
                            qpvdw = float(line[120:128])
                            qwel = float(line[128:136])
                            qwvdw = float(line[136:144])

                            tmp.append([l, eqtot, eqbnd, eqang, eqtor, eqimp, eqel, eqvdw, qqel, qqvdw, qpel, qpvdw,
                                        qwel, qwvdw])

                            if int(line.split()[1]) == 2:
                                part0.append(tmp)
                                tmp = list()
                    except:
                        continue

                if '# Part 0:' in line:
                    found_part0 = True

        return part0

    def get_qfep_part3(self, path, filename):
        """
        Reads qfep.out in path, and returns part 3 as list
        """
        part3 = list()
        found_part3 = False
        qfepout = '%s/%s' % (path, filename)
        with open(qfepout, 'r') as qfepout:
            for line in qfepout:
                if "# Part 4:" in line:
                    break
                if found_part3:
                    if not line.startswith('#') and len(line.split()) > 4:
                        part3.append(line)
                    if len(line.split()) < 4:
                        break
                if '# Part 3:' in line:
                    found_part3 = True

        return part3

    def export_table(self, listbox):
        """
        Exports data in listbox to file:
        """
        filename = asksaveasfilename(parent=self, title='Export table', initialdir=self.app.workdir,
                                      filetypes=(("Text", "*.txt"), ("All files","*.*")),
                                      initialfile = 'Qdata.txt')

        if not filename:
            return



        data_to_export = listbox.get(1, END)
        output = open(filename, 'w')

        for line in data_to_export:
            output.write(line)
            if not '\n' in line:
                output.write('\n')

        output.close()

        self.app.log('info', '../%s written' % '/'.join(filename.split()[-3:]))

    def list_titles_event(self, *args):

        selections = list(map(int, self.titles_listbox.curselection()))
        if len(selections) == 0:
            return

        plt_titles = list()
        plt_temps = list()
        plt_dg = list()
        plt_errors = list()
        plt_parameters = list()
        plt_dg_exc = list()
        plt_sem_exc = list()
        plt_t_exc = list()

        self.runs_listbox.delete(0, END)
        self.exclude_listbox.delete(0, END)
        for selected in selections:
            title = self.titles_listbox.get(selected)
            plt_titles.append(title)
            plt_parameters.append([self.titles_parameters[title]['dH'], self.titles_parameters[title]['dS']])
            tmp_temps = list()
            tmp_errors = list()
            tmp_dg = list()
            tmp_t_exc = list()
            tmp_dg_exc = list()
            tmp_sem_exc = list()
            if title in list(self.title_exclude.keys()):
                for t_ex in self.title_exclude[title]:
                    self.exclude_listbox.insert(END, '%5s' % t_ex)
                if len(self.title_exclude[title]) > 0:
                    tmp_t_exc = self.title_exclude[title]
            for temp in sorted(self.titles[title]):
                self.runs_listbox.insert(END, '%5s' % temp)
                if temp not in tmp_t_exc:
                    tmp_temps.append(1.00/temp)
                    tmp_dg.append(self.titles_ave_act[title][temp][0]/ temp)
                    tmp_errors.append(self.titles_ave_act[title][temp][1]/ temp)
                else:
                    tmp_dg_exc.append(self.titles_ave_act[title][temp][0]/ temp)
                    tmp_sem_exc.append(self.titles_ave_act[title][temp][1]/ temp)
            if len(tmp_t_exc) > 0:
                tmp_t_exc = [1.0 / x for x in tmp_t_exc]

            plt_temps.append(tmp_temps)
            plt_dg.append(tmp_dg)
            plt_errors.append(tmp_errors)
            plt_t_exc.append(tmp_t_exc)
            plt_dg_exc.append(tmp_dg_exc)
            plt_sem_exc.append(tmp_sem_exc)

        if len(plt_dg) > 0:
            self.update_plot(plt_temps, plt_dg, plt_errors, plt_parameters, plt_titles, plt_t_exc, plt_dg_exc,
                             plt_sem_exc)

    def list_runs_event(self, *args):
        selections = list(map(int, self.runs_listbox.curselection()))
        if len(selections) == 0:
            return
        pass

    def list_dg_events(self, *args):

        #delete existing values
        self.dg_act_upper.delete(0, END)
        self.dg_act_lower.delete(0, END)
        self.dg_rxn_upper.delete(0, END)
        self.dg_rxn_lower.delete(0, END)

        selections = list(map(int, self.dg_listbox.curselection()))

        if len(selections) != 1:
            return

        try:
            #title = self.dg_listbox.get(selections[0]).split()[0]
            title = self.titles_listbox.get(int(self.titles_listbox.curselection()[0]))
            temp = int(self.dg_listbox.get(selections[0]).split()[1 + (len(title.split()) - 1)])

            act_upper = self.dg_upper_lower[title][temp]['activation'][0]
            act_lower = self.dg_upper_lower[title][temp]['activation'][1]
            rxn_upper = self.dg_upper_lower[title][temp]['reaction'][0]
            rxn_lower = self.dg_upper_lower[title][temp]['reaction'][1]

            self.dg_act_upper.insert(0, act_upper)
            self.dg_act_lower.insert(0, act_lower)
            self.dg_rxn_upper.insert(0, rxn_upper)
            self.dg_rxn_lower.insert(0, rxn_lower)

        except:
            pass


    def show_var_frame(self, *args):
        frames = {'Arrhenius Plot': self.plot_frame,
                  'Reaction Free Energies': self.dg_frame,
                  'Regression Parameters': self.re_tot_frame,
                  'Model VS Computed': self.ae_nb_frame}

        for i in list(frames.keys()):
            frames[i].grid_forget()
        try:
            frames[self.selected_frame.get()].grid(row=3, column=0, columnspan=2)
        except:
            pass

    def dialog_window(self):
        self.title('EVB Thermodynamic Parameters')
        self.mainframe = Frame(self, bg=self.main_color)
        self.mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame with import button etc
        topframe = Frame(self.mainframe, bg=self.main_color)
        topframe.grid(row=0, column=0, pady=(10,10))

        #plot frame
        self.plot_frame = Frame(self.mainframe, bg=self.main_color)
        self.plot_frame.grid(row=3, column=0)
        frame_plot= Frame(self.plot_frame, bg=self.main_color)
        frame_plot.pack(side=TOP, padx=(10, 10), pady=(10, 0), fill=BOTH)

        #Frame selector frame:
        sel_frame = Frame(self.mainframe, bg=self.main_color)
        sel_frame.grid(row=2, column=0)

        #Bottomframe
        bottomframe = Frame(self.mainframe, bg=self.main_color)
        bottomframe.grid(row=4, column=0)

        #Variable frames:

        #Project entries and energy files
        self.proj_frame = Frame(self.mainframe, bg=self.main_color)
        self.proj_frame.grid(row=1, column=0, pady=(0,10))

        self.re_tot_frame =  Frame(self.mainframe, bg=self.main_color)

        self.ae_nb_frame = Frame(self.mainframe, bg=self.main_color)

        self.dg_frame = Frame(self.mainframe, bg=self.main_color)

        #self.rxn_energies = Frame(self.mainframe, bg=self.main_color)


        #topframe content
        import_label = Label(topframe, text='Import: ', bg=self.main_color)
        import_label.grid(row=0, column=0)

        qfep_import = Button(topframe, text='EVB parameters', highlightbackground=self.main_color, command=self.import_qfep)
        qfep_import.grid(row=0, column=1, columnspan=2)

        alpha = Text(topframe, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        alpha.tag_configure("subscript", offset=-1)
        alpha.insert("insert","\N{GREEK SMALL LETTER ALPHA}","",'ij','subscript')
        alpha.grid(row=1, column=0, sticky='e')
        alpha.config(state=DISABLED)

        self.alpha_entry = Spinbox(topframe, width=7, highlightthickness=0, relief=GROOVE, from_=-500.00, to=500.00,
                                   increment=1)
        self.alpha_entry.grid(row=1, column=1)
        self.alpha_entry.delete(0, END)
        self.alpha_entry.insert(0, '0')

        hij = Text(topframe, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        hij.tag_configure("subscript", offset=-1)
        hij.insert("insert",'H',"",'ij','subscript')
        hij.grid(row=1, column=2, sticky='e')
        hij.config(state=DISABLED)

        self.hij_entry = Spinbox(topframe, width=7, highlightthickness=0, relief=GROOVE, from_=0, to=200, increment=1.0)
        self.hij_entry.grid(row=1, column=3)

        kt = Label(topframe, text='    kT     ', bg=self.main_color)
        kt.grid(row=0, column=4)

        bins = Label(topframe, text='   Bins    ', bg=self.main_color)
        bins.grid(row=0, column=5)

        binpoints = Label(topframe, text='Min. pts', bg=self.main_color)
        binpoints.grid(row=0, column=6)

        skip_points = Label(topframe, text='   Skip    ', bg=self.main_color)
        skip_points.grid(row=0, column=7)

        linear_combination = Label(topframe, text='Lin. comb.', bg=self.main_color)
        linear_combination.grid(row=0, column=8)

        self.kT = Entry(topframe, width=5, highlightthickness=0)
        self.kT.grid(row=1, column=4)
        self.kT.delete(0,END)
        self.kT.insert(0, 'auto')

        self.bins = Entry(topframe, width=5, highlightthickness=0)
        self.bins.grid(row=1, column=5)
        self.bins.insert(0, '50')

        self.binpoints_min = Entry(topframe, width=5, highlightthickness=0)
        self.binpoints_min.grid(row=1, column=6)
        self.binpoints_min.insert(0, '30')

        self.skip_points = Entry(topframe, width=5, highlightthickness=0)
        self.skip_points.grid(row=1, column=7)
        self.skip_points.insert(0, '100')

        self.linear_comb = Spinbox(topframe, width=6, highlightthickness=0, relief=GROOVE,
                                   values=(' 1    -1 ', '  -1    1'))
        self.linear_comb.grid(row=1, column=8)



        #Frame3 (plot)
        self.plot_window = Figure(figsize=(5.5,3), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.plot_window, master=frame_plot)
        self.plot_window.patch.set_facecolor('white')

        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame_plot)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        #Select frame
        self.view_frame = OptionMenu(sel_frame, self.selected_frame,
                                   'Arrhenius Plot', 'Regression Parameters',
                                   'Reaction Free Energies', 'Model VS Computed')
        self.view_frame.config(highlightbackground=self.main_color, bg=self.main_color, width=30)
        self.view_frame.grid(row=4, column=0)

        recomp_evb = Button(sel_frame, text='Re-compute EVB', highlightbackground=self.main_color,
                            command=self.recomp_evb)
        recomp_evb.grid(row=4, column=2, columnspan=2)

        #Project frame
        entry_label = Label(self.proj_frame, text='Titles', bg=self.main_color)
        entry_label.grid(row=0, column=0)

        add_title = Button(self.proj_frame, text='+', highlightbackground=self.main_color, command=self.add_title)
        add_title.grid(row=0, column=1)

        del_title = Button(self.proj_frame, text='-', highlightbackground=self.main_color, command=self.del_title)
        del_title.grid(row=0, column=2)

        runs_label = Label(self.proj_frame, text=' T ', bg=self.main_color)
        runs_label.grid(row=0, column=4)

        add_runs = Button(self.proj_frame, text='+', highlightbackground=self.main_color, command=self.add_runs)
        add_runs.grid(row=0, column=5)

        del_title = Button(self.proj_frame, text='-', highlightbackground=self.main_color, command=self.del_runs)
        del_title.grid(row=0, column=6)

        excl_label = Label(self.proj_frame, text=' X ', bg=self.main_color)
        excl_label.grid(row=0, column=8)

        add_excl = Button(self.proj_frame, text='+', highlightbackground=self.main_color, command=self.add_exluded_temp)
        add_excl.grid(row=0, column=9)

        del_excl = Button(self.proj_frame, text='-', highlightbackground=self.main_color, command=self.del_excluded_temp)
        del_excl.grid(row=0, column=10)


        titles_yscroll = Scrollbar(self.proj_frame)
        titles_yscroll.grid(row = 1, rowspan=10, column = 3, sticky = 'nsw', padx=(0,10))
        titles_xscroll = Scrollbar(self.proj_frame, orient=HORIZONTAL)
        titles_xscroll.grid(row=11, column=0, columnspan=3, sticky='we')

        self.titles_listbox = Listbox(self.proj_frame, yscrollcommand = titles_yscroll.set,
                                      xscrollcommand=titles_xscroll.set,
                                      width=18, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        titles_yscroll.config(command=self.titles_listbox.yview)
        titles_xscroll.config(command=self.titles_listbox.xview)
        self.titles_listbox.grid(row=1, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.titles_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        self.titles_listbox.bind('<<ListboxSelect>>', self.list_titles_event)

        runs_yscroll = Scrollbar(self.proj_frame)
        runs_yscroll.grid(row = 1, rowspan=10, column = 7, sticky = 'nsw', padx=(0,10))
        runs_xscroll = Scrollbar(self.proj_frame, orient=HORIZONTAL)
        runs_xscroll.grid(row=11, column=4, columnspan=3, sticky='we')

        self.runs_listbox = Listbox(self.proj_frame, yscrollcommand = runs_yscroll.set,
                                      xscrollcommand=runs_xscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        runs_yscroll.config(command=self.runs_listbox.yview)
        runs_xscroll.config(command=self.runs_listbox.xview)
        self.runs_listbox.grid(row=1, rowspan=10, column = 4, columnspan=3, sticky='e')
        self.runs_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        self.runs_listbox.bind('<<ListboxSelect>>', self.list_runs_event)

        exclude_yscroll = Scrollbar(self.proj_frame)
        exclude_yscroll.grid(row = 1, rowspan=10, column = 11, sticky = 'nsw', padx=(0,10))
        exclude_xscroll = Scrollbar(self.proj_frame, orient=HORIZONTAL)
        exclude_xscroll.grid(row=11, column=8, columnspan=3, sticky='we')

        self.exclude_listbox = Listbox(self.proj_frame, yscrollcommand = exclude_yscroll.set,
                                      xscrollcommand=exclude_xscroll.set,
                                      width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        exclude_yscroll.config(command=self.exclude_listbox.yview)
        exclude_xscroll.config(command=self.exclude_listbox.xview)
        self.exclude_listbox.grid(row=1, rowspan=10, column = 8, columnspan=3, sticky='e')
        self.exclude_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        #self.exclude_listbox.bind('<<ListboxSelect>>', self.list_runs_event)


        #Regression stats
        re_tot_yscroll = Scrollbar(self.re_tot_frame)
        re_tot_yscroll.grid(row = 1, rowspan=16, column = 3, sticky = 'nsw', padx=(0,10))

        self.re_tot_listbox = Listbox(self.re_tot_frame, yscrollcommand = re_tot_yscroll.set,
                                      width=82, height=16, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        re_tot_yscroll.config(command=self.re_tot_listbox.yview)
        self.re_tot_listbox.grid(row=1, rowspan=16, column = 0, columnspan=3, sticky='e')
        self.re_tot_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        #self.ae_tot_listbox.bind('<<ListboxSelect>>', self.list_runs_event)

        export_tot = Button(self.re_tot_frame, text='Export table', highlightbackground=self.main_color,
                           command=lambda: self.export_table(self.re_tot_listbox))
        export_tot.grid(row=17, column=0, columnspan=3, sticky='e')

        #Model VS computed
        ae_nb_yscroll = Scrollbar(self.ae_nb_frame)
        ae_nb_yscroll.grid(row = 1, rowspan=16, column = 3, sticky = 'nsw', padx=(0,10))

        self.ae_nb_listbox = Listbox(self.ae_nb_frame, yscrollcommand = ae_nb_yscroll.set,
                                      width=82, height=16, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        ae_nb_yscroll.config(command=self.ae_nb_listbox.yview)
        self.ae_nb_listbox.grid(row=1, rowspan=16, column = 0, columnspan=3, sticky='e')
        self.ae_nb_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        #self.ae_tot_listbox.bind('<<ListboxSelect>>', self.list_runs_event)

        export_nb = Button(self.ae_nb_frame, text='Export table', highlightbackground=self.main_color,
                           command=lambda: self.export_table(self.ae_nb_listbox))
        export_nb.grid(row=17, column=0, columnspan=3, sticky='e')

        #Reaction Energies
        dg_yscroll = Scrollbar(self.dg_frame)
        dg_yscroll.grid(row = 1, rowspan=16, column = 3, sticky = 'nsw', padx=(0,10))

        self.dg_listbox = Listbox(self.dg_frame, yscrollcommand = dg_yscroll.set,
                                      width=57, height=16, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        dg_yscroll.config(command=self.dg_listbox.yview)
        self.dg_listbox.grid(row=1, rowspan=16, column = 0, columnspan=3, sticky='w')
        self.dg_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        self.dg_listbox.bind('<<ListboxSelect>>', self.list_dg_events)

        #### Remover outliers ####
        upper_label = Label(self.dg_frame, text='upper', bg=self.main_color)
        upper_label.grid(row=1, column=5)

        lower_label = Label(self.dg_frame, text='lower', bg=self.main_color)
        lower_label.grid(row=1, column=6)

        dg_act_label = Label(self.dg_frame, text='dG(act)', bg=self.main_color)
        dg_act_label.grid(row=2, column=4)

        dg_rxn_label = Label(self.dg_frame, text='dG(rxn)', bg=self.main_color)
        dg_rxn_label.grid(row=3, column=4)

        self.dg_act_upper = Entry(self.dg_frame, width=7, highlightthickness=0)
        self.dg_act_upper.grid(row=2, column=5)

        self.dg_act_lower = Entry(self.dg_frame, widt=7, highlightthickness=0)
        self.dg_act_lower.grid(row=2,column=6)

        self.dg_rxn_upper = Entry(self.dg_frame, width=7, highlightthickness=0)
        self.dg_rxn_upper.grid(row=3, column=5)

        self.dg_rxn_lower = Entry(self.dg_frame, widt=7, highlightthickness=0)
        self.dg_rxn_lower.grid(row=3,column=6)

        apply_tresh = Button(self.dg_frame, text='Apply', highlightbackground=self.main_color,
                             command=self.set_upper_lower)
        apply_tresh.grid(row=4, column=5, columnspan=2)

        ############################
        export_dg = Button(self.dg_frame, text='Export table', highlightbackground=self.main_color,
                           command=lambda: self.export_table(self.dg_listbox))
        export_dg.grid(row=17, column=0, columnspan=3, sticky='e')


        #Bottom frame Quit/save
        quit_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        quit_button.grid(row=0, column=0)
