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

#Changes for matplotlib v2.2.3 on lines 525, 530, 564, 724, 1689, 1809, 1812, 633
# NavigationToolbar2TkAgg->NavigationToolbar2Tk

from tkinter import Label, TOP, Button, Listbox, Scrollbar, EXTENDED, Spinbox, Entry, Text, Frame, \
    Toplevel, DISABLED, END, GROOVE, NORMAL, BOTH, IntVar, StringVar, Checkbutton, OptionMenu, HORIZONTAL
from tkinter.filedialog import askopenfilename
from tkinter.simpledialog import askstring
import tkinter.font
from subprocess import call
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
from cycler import cycler

class EvbReactions(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root

        self.main_color = self.app.main_color
        #This can be used to make global designs, ratter than locals:
        #---> self.root.option_add('*background', self.main_color)

        self.dg_plot = None

        #Toggle between plotting diabatic states (reorganization energy calculations) and reaction free energy profiles
        self.diabatic_plot = False

        #Plot options in reorganization energy window
        self.plot_hii = IntVar()
        self.plot_fit = IntVar()
        self.plot_norm_fit = IntVar()
        self.plot_hij = IntVar()
        self.plot_c = IntVar()
        self.plot_mix = IntVar()

        self.plot_hii.set(1)

        #Selector when changing windows from dropdown menu
        self.selected_frame = StringVar()

        #{'name': {1:path, 2:path ....}}
        self.titles = dict()

        #{'name':{1: [[d_eps],[d_G]], 2:...}}
        self.titles_ene = dict()

        #{name:[[d_eps], [<d_G>]]}
        self.titles_ave = dict()

        #{name: {1: [dG_act, dG_react]}, 2...}
        self.titles_dG = dict()

        # self.titles[title] = [dU, sem] where dU and stdm have the following keys:
        #'EQtot', 'EQbnd', 'EQang', 'EQtor', 'EQimp' , 'EQel', 'EQvdw', 'qq_el', 'qq_vdw', 'qp_el', 'qp_vdw',
        # 'qw_el', 'qw_vdw', 'alpha', 'Hij'}
        self.titles_dU = dict()

        #{title: {run:{energy_gap: [dg1, dg2, c1, c2, Hij]}}}
        self.title_reorge = dict()


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

        self.titles[title] = dict()
        self.titles_ene[title] = dict()

        self.titles_dG[title] = dict()

        #Highlight latest title in listbox
        self.titles_listbox.select_clear(0, END)
        titles = self.titles_listbox.get(0, END)
        for ind_ in range(len(titles)):
            if titles[ind_] == title:
                self.titles_listbox.selection_set(ind_)
                #self.list_titles_event()

        #self.update_plot([[0]],[[0]],[title])

    def del_title(self):
        selections = list(map(int, self.titles_listbox.curselection()))
        for selected in selections:
            title = self.titles_listbox.get(selected)
            self.titles_listbox.delete(selected)
            del self.titles[title]
            for i in [self.titles_dG, self.titles_ave, self.titles_dU, self.titles_ene]:
                if title in list(i.keys()):
                    del i[title]

        self.update_plot([[0]],[[0]],[' '])
        self.update_tables()

    def add_runs(self):
        """
        Add subrun one-by-one or select temperature directory to append all runs
        """
        try:
            prj_entry = list(map(int, self.titles_listbox.curselection()))
        except:
            self.app.errorBox('Warning', 'Select a project title to append runs to.')
            return

        if len(prj_entry) != 1:
            self.app.errorBox('Warning', 'Select exactly one project title to append runs to.')
            return

        temp_selected = False
        prj_title = self.titles_listbox.get(prj_entry[0])

        #Start appending directories to collect energy files from
        dirs = list()
        title = 'Select EVB temperature or run directory'
        path = self.app.workdir
        while True:
            rundir = askdirectory(parent=self, mustexist=True, title=title, initialdir= path)
            if not rundir:
                break
            #Check if the first directory is a temperature directory
            if len(dirs) == 0:
                try:
                    print((int(float(rundir.split('/')[-1]))))
                    if int(float(rundir.split('/')[-1])) > 270:
                        self.app.log(' ', '%s Added. Collecting runs for temperature ...\n\n'
                                          % '/'.join(rundir.split('/')[-3:]))

                        temp_selected = True
                except:
                    pass

            if not temp_selected:
                path = '%s' % '/'.join(rundir.split('/')[:-1])

                #Check if en files exist in dir, and if size above 0:
                enefiles = False
                for ene in sorted(os.listdir(rundir)):
                    if ene.endswith('.en'):
                        if os.path.getsize('%s/%s' % (rundir, ene)) != 0:
                            enefiles = True
                        else:
                            enefiles = False
                            break
                if enefiles:
                    dirs.append(rundir)
                    title = '%s Added. Select next dir or cancel' % dirs[-1].split('/')[-1]
                    self.app.log(' ', '%s Added. Select next dir or cancel to finish\n'
                                      % '/'.join(dirs[-1].split('/')[-3:]))
                    try:
                        self.titles[prj_title][int(rundir.split('/')[-1])] = rundir
                    except:
                        self.titles[prj_title][rundir.split('/')[-1]] = rundir
            else:
                break

        #If temperature direcory is recognized, run thorugh all subfolders and collect qfep out, or run Qfep
        if temp_selected:
            #Temperature directory. Collect run dirs and break
            nr_dir = 1
            while True:
                subdir = '%s/%d' % (rundir, nr_dir)

                if os.path.isdir(subdir):
                    enefiles = False
                    for ene in os.listdir(subdir):
                        if ene.endswith('.en'):
                            if os.path.getsize('%s/%s' % (subdir, ene)) != 0:
                                enefiles = True
                            else:
                                enefiles = False
                                break
                    if enefiles:
                        dirs.append(subdir)
                        self.titles[prj_title][nr_dir] = subdir
                        self.titles_ene[prj_title][nr_dir] = list()
                        self.titles_dG[prj_title][nr_dir] = list()
                        self.app.log(' ','...../%s added\n' % '/'.join(subdir.split('/')[-2:]))
                        print(subdir)
                    nr_dir += 1
                else:
                    break

        #Collect Bin-Ave energies
        self.get_binave_out(prj_title, False)

        #Get all dG(act) dG(react) and normalize
        self.get_dG(prj_title)

        #Compute activation energies
        self.compute_ave_activation_ene(prj_title)

        self.list_titles_event()

    def copy_all_enefiles(self, prj_title, qfep_inp=False, delete_existing=False):
        """
        Copies all energy files belonging to prj_title to the all directory. If qfep_inp is False, a qfep.inp
        file will be copied (if existing from the subdirectories). If delete_existing=True, all existing files
        will be removed from all_dir before copying
        """

        all_dir = '%s/%s_all' % (self.app.workdir, prj_title)
        if not os.path.isdir(all_dir):
            os.mkdir(all_dir)
            qfep_inp = False

        if delete_existing:
            for qene in os.listdir(all_dir):
                if qene.endswith('.en'):
                    os.remove('%s/%s' % (all_dir, qene))

        enefiles = list()
        for qdir in list(self.titles[prj_title].keys()):
            if qdir != 'all':
                path = self.titles[prj_title][qdir]
                qfiles = os.listdir(path)
                for qfile in qfiles:
                    if qfile.endswith('.en'):
                        if qfile.startswith('md'):
                            enefiles.append('%s/%s' % (path, qfile))

        enefiles = sorted(enefiles, key=lambda x: x.split('/')[-1])
        count = 0
        for enefile in enefiles:
            if os.path.getsize(enefile) != 0:
                count += 1
                shutil.copy2(enefile, '%s/md%05d.en' % (self.titles[prj_title]['all'], count))
                print(('../%s ---> md%05d.en' % ('/'.join(enefile.split('/')[-3:]), count)))

        if not qfep_inp:
            self.write_qfep_inp(all_dir)
        self.run_qfep(all_dir, 'qfep.inp')
        self.get_binave_out(prj_title, False)

    def del_runs(self):
        title_sel = list(map(int, self.titles_listbox.curselection()))
        if len(title_sel) != 1:
            self.app.log(' ', '\nSelect exactly one title to delete runs from!\n')
            return

        selections = list(map(int, self.runs_listbox.curselection()))
        if len(selections) == 0:
            self.app.log(' ', '\nNo runs selected for deletion!\n')
            return

        title = self.titles_listbox.get(title_sel[0])

        titles_to_recopy = list()

        for selected in reversed(selections):
            search = '/'.join(self.runs_listbox.get(selected).split('/')[-3:])
            self.runs_listbox.delete(selected)

            if 'all' in list(self.titles[title].keys()):
                if title not in titles_to_recopy:
                    titles_to_recopy.append(title)
            for nr in list(self.titles[title].keys()):
                if '/'.join(self.titles[title][nr].split('/')[-3:]) == search:
                    if nr != 'all':
                        for i in [self.titles, self.titles_dG, self.titles_ene]:
                            if title in list(i.keys()):
                                if nr in list(i[title].keys()):
                                    del i[title][nr]

        for title in titles_to_recopy:
            if 'all' in list(self.titles_ene[title].keys()):
                self.titles_ene[title]['all'] = list()
                self.titles_dG[title]['all'] = list()
                self.copy_all_enefiles(title, qfep_inp=False, delete_existing=True)

        #Get all dG(act) dG(react) and normalize
        self.get_dG(title)

        #Compute activation energies
        self.compute_ave_activation_ene(title)

        self.update_tables()
        self.list_titles_event()

    def make_evb_average(self):
        """
        Creates a new directory and copies all energy files (renamed and sorted) to directory, and
        runs Qfep on all energy files.
        """
        try:
            prj_entry = list(map(int, self.titles_listbox.curselection()))
        except:
            self.app.errorBox('Warning', 'Select a project title to append runs to.')
            return

        self.app.log(' ', '\n')
        self.app.log('info','Generating EVB average profile ...')


        for i in range(len(prj_entry)):
            prj_title = self.titles_listbox.get(prj_entry[i])

            #Create entries for averages
            #Add directory for average of all runs
            self.titles[prj_title]['all'] = '%s/%s_all' % (self.app.workdir, prj_title)
            self.titles_ene[prj_title]['all'] = list()
            self.titles_dG[prj_title]['all'] = list()

            #Copy all energy files to average directory:
            self.copy_all_enefiles(prj_title, qfep_inp=False, delete_existing=True)

            #Collect Bin-Ave energies
            self.get_binave_out(prj_title, False)

            #Get all dG(act) dG(react) and normalize
            self.get_dG(prj_title)

            #Update listboxes to latest added
            list_titles = self.titles_listbox.get(0, END)
            for ind_ in range(len(list_titles)):
                if self.titles_listbox.get(ind_) == prj_title:
                    self.titles_listbox.selection_set(ind_)
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

        #Sort energy files (if not Qfep will give messed up free energy profile!)
        enefiles = sorted(enefiles, key=lambda x: x.split('/')[-1])

        #make input for Qfep
        inputfile = open('%s/qfep.inp' % path,'w')
        inputfile.write('%d\n' % len(enefiles))
        inputfile.write('2  0\n')
        inputfile.write('%s  %s\n' % (self.kT.get(), self.skip_points.get()))
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

    def get_binave_out(self, prj_title, run_qfep = False):
        """
        Goes through all subdirectories for prj_title and reads qfep.out bin-ave section and appends to
        self.titles_ene. If qfep.out does not exist, Qfep is run. If run_qfep == True, Qfep will be run in all
        direcotories
        """

        self.app.log('info', 'Collecting Reaction Free Energies ...')

        for nr in sorted([x for x in self.titles[prj_title].keys() if type(x)==int]) + [x for x in self.titles[prj_title].keys() if type(x)!=int]:
            path = self.titles[prj_title][nr]
            qfepfile = '%s/qfep.out' % self.titles[prj_title][nr]
            if not os.path.isfile(qfepfile):
                if not os.path.isdir('%s.inp' % qfepfile.split('.')[-2]):
                    #Write qfep.inp
                    qfep_inp = self.write_qfep_inp(path)
                    self.app.log(' ', '\nGenerating new qfep.inp: ../%s\n' % self.titles[prj_title][nr].split('/')[-3:])
                    if qfep_inp:
                        self.run_qfep(path, 'qfep.inp')
                    else:
                        self.app.log(' ', '\nFailed to create qfep.inp in ../%s\n' % self.titles[prj_title][nr].split('/')[-3:])
            elif run_qfep:
                #Write new qfep.inp
                self.app.log(' ', '\nGenerating new qfep.inp: ../%s\n' % self.titles[prj_title][nr].split('/')[-3:])
                qfep_inp = self.write_qfep_inp(path)
                if qfep_inp:
                    self.run_qfep(path, 'qfep.inp')
                else:
                    self.app.log(' ', '\nFailed to create qfep.inp in ../%s\n' % self.titles[prj_title][nr].split('/')[-3:])
            enegaps, dG = self.get_enegaps_dg(qfepfile)
            self.titles_ene[prj_title][nr] = [enegaps, dG]

    def get_enegaps_dg(self, qfepout):
        """
        Read qfepout #part3 and returns 2 lists: Energy gaps and delta G.
        """
        #Get enegap and dG
        part3 = False
        count = 0
        enegaps = list()
        dG = list()
        if os.path.isfile(qfepout):
            with open(qfepout, 'r') as qfepout:
                for line in qfepout:
                    if '# Part 3: Bin-averaged summary' in line:
                        part3 = True
                    if part3:
                        if line == '':
                            break
                        if count > 1:
                            try:
                                enegaps.append(float(line.split()[1]))
                                dG.append(float(line.split()[3]))
                            except:
                                continue
                        else:
                            count += 1

        return enegaps, dG

    def import_qfep(self):
        """
        Import EVB parameters from qfep file
        """
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

    def update_plot(self, d_eps, d_G, titles ):
        if len(titles) < 1:
            return

        if self.dg_plot:
            self.dg_plot.clear()

        #Set color cycle for plots:
        #Made a change here
        matplotlib.rcParams['axes.prop_cycle'] = cycler('color', ['k', 'b', 'g', 'r', 'm', 'y', 'c', 'brown',
                                                   'burlyWood', 'cadetBlue', 'DarkGreen', 'DarkBlue',
                                                   'DarkMagenta', 'DarkSalmon', 'DimGray', 'Gold'])

        #Create subplot
        self.dg_plot = self.plot_window.add_subplot(111, facecolor='white')
        self.plot_window.subplots_adjust(hspace=0.5)

        #X/Y labels
        self.dg_plot.set_xlabel(r'$\Delta \epsilon$ (kcal/mol)')
        self.dg_plot.set_ylabel(r'$\Delta G$ (kcal/mol)')

        #Move label box outside plot region
        box = self.dg_plot.get_position()
        self.dg_plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        #Fit subplot to figure/canvas
        #rect=(left,bottom,top,right)
        self.plot_window.tight_layout(rect=(0.005, 0, 0.8, 1))

        for i in range(len(titles)):
            title = str(titles[i])
            nr = 'ave'

            if len(title.split('_')) > 1:
                if title.split('_')[-1].isdigit():
                    nr = title.split('_')[-1]
                    if len(title.split('_')) > 2:
                        title = '\_'.join(title.split('_')[0:-1])
                    else:
                        title = title.split('_')[0]
                else:
                    title = '_'.join(title.split('\_'))

            dG = d_G[i]
            deps = d_eps[i]
            self.dg_plot.plot(deps, dG, '-', linewidth=2.0, label=r'%s$_{%s}$' % (title, nr))
            self.dg_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})

        self.canvas.draw()

    def update_diabatic_plot(self):
        """
        Plots diabatic free energy functions and updates reorganization energies
        {title: {run: 'diabatic': {energy_gap: [dg1, dg2, c1, c2, Hij], 'marcus1': [[e], [dg]], 'marcus2': {[e]:[dg]},
         'norm1':[[e], [dg]], 'norm2': [[e], [dg]], 'reorg': value}}}
        """
        selections = list(map(int, self.titles_listbox.curselection()))
        if len(selections) == 0:
            return

        if self.dg_plot:
            self.dg_plot.clear()

        reorgs = list()

        title = None
        #Get reorganization energy avarage for selected title(s)
        for i in selections:
            title = self.titles_listbox.get(i)
            if not title in list(self.title_reorge.keys()):
                return
            for run in list(self.title_reorge[title].keys()):
                reorgs.append(self.title_reorge[title][run]['reorg'])

        reorg_ave = np.average(reorgs)
        reorg_se = np.std(reorgs) / np.sqrt(len(reorgs))

        self.reorg_avg.config(state=NORMAL)
        self.reorg_se.config(state=NORMAL)
        self.reorg_avg.delete(0.0, END)
        self.reorg_se.delete(0.0, END)

        self.reorg_avg.insert(0.0, '%.1f' % reorg_ave)
        self.reorg_se.insert(0.0, '%.1f' % reorg_se)

        self.reorg_avg.config(state=DISABLED)
        self.reorg_se.config(state=DISABLED)

        #Plot selected run
        runs = list(map(int, self.runs_listbox.curselection()))
        if len(runs) > 1:
            print('Select exactly one title with one run to plot diabatic energy functions')
            return

        if len(runs) == 0:
            return

        run = self.runs_listbox.get(runs[0]).split('/')[-1]
        path = None

        for i in list(self.title_reorge[title].keys()):
            if i.split('/')[-1] == run:
                path = i

        if not path:
            return

        self.reorg_curr.config(state=NORMAL)
        self.reorg_curr.delete(0.0, END)

        self.reorg_curr.insert(0.0, '%.1f' % self.title_reorge[title][path]['reorg'])

        self.reorg_curr.config(state=DISABLED)


        #Create subplot
        self.dg_plot = self.plot_window.add_subplot(111, facecolor='white')
        self.plot_window.subplots_adjust(hspace=0.5)

        #X/Y labels
        self.dg_plot.set_xlabel(r'$\Delta \epsilon$ (kcal/mol)')
        self.dg_plot.set_ylabel(r'$\Delta g$ (kcal/mol)')

        #Move label box outside plot region
        box = self.dg_plot.get_position()
        self.dg_plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        #Fit subplot to figure/canvas
        #rect=(left,bottom,top,right)
        self.plot_window.tight_layout(rect=(0.005, 0, 0.8, 1))

        deps = list()
        dg1_list = list()
        dg2_list = list()
        dg1_c1_list = list()
        dg2_c2_list = list()
        hij_list = list()
        dg_mix = list()
        dg_c_mix = list()

        for enegap in sorted(self.title_reorge[title][path]['diabatic'].keys()):
            deps.append(enegap)
            dg1, dg2, c1, c2, hij = self.title_reorge[title][path]['diabatic'][enegap][0:]
            dg1_list.append(dg1)
            dg2_list.append(dg2)
            dg1_c1_list.append(dg1 * c1)
            dg2_c2_list.append(dg2 * c2)
            hij_list.append(2 * hij * np.sqrt(c1*c2))
            dg_mix.append(dg1 + dg2)
            dg_c_mix.append((dg1*c1) + (dg2*c2))

        #Check if something was plotted
        plot_it = False

        #Plot diabatic free energy functions?
        if self.plot_hii.get() == 1:
            plot_it = True
            if self.plot_mix.get() == 0:
                if self.plot_c.get() == 1:
                    self.dg_plot.plot(deps, dg1_c1_list, 'bo', linewidth=1.0, label = r'$\Delta g_1$')
                    self.dg_plot.plot(deps, dg2_c2_list, 'ro', linewidth=1.0, label = r'$\Delta g_2$')
                else:
                    self.dg_plot.plot(deps, dg1_list, 'bo', linewidth=1.0, label = r'$\Delta g_1$')
                    self.dg_plot.plot(deps, dg2_list, 'ro', linewidth=1.0, label = r'$\Delta g_2$')
            else:
                self.dg_plot.plot(deps, dg1_list, 'bo', linewidth=1.0, label = r'$\Delta g_1$')
                self.dg_plot.plot(deps, dg2_list, 'ro', linewidth=1.0, label = r'$\Delta g_2$')

        #Plot fit?
        if self.plot_fit.get() == 1:
            plot_it = True
            eps1 = self.title_reorge[title][path]['marcus1'][0]
            marcus1 = self.title_reorge[title][path]['marcus1'][1]
            eps2 = self.title_reorge[title][path]['marcus2'][0]
            marcus2 = self.title_reorge[title][path]['marcus2'][1]

            self.dg_plot.plot(eps1, marcus1, 'k', linewidth=2.0, label='Fit')
            self.dg_plot.plot(eps2, marcus2, 'k', linewidth=2.0)

        #Plot normalized fit?
        if self.plot_norm_fit.get() == 1:
            plot_it = True
            eps1 = self.title_reorge[title][path]['norm1'][0]
            norm1 = self.title_reorge[title][path]['norm1'][1]
            eps2 = self.title_reorge[title][path]['norm2'][0]
            norm2 = self.title_reorge[title][path]['norm2'][1]

            self.dg_plot.plot(eps1, norm1, 'y', linewidth=2.0, label='Fit*')
            self.dg_plot.plot(eps2, norm2, 'y', linewidth=2.0)

        #Plot adiabatic mixing?
        if self.plot_mix.get() == 1:
            plot_it = True

            if self.plot_c.get() == 1:
                if self.plot_hij.get() == 1:
                    self.dg_plot.plot(deps, [a - b for a, b in zip(dg_c_mix, hij_list)], 'kx',
                                      label=r'$\Delta g_1 + \Delta g_2$')
                else:
                    self.dg_plot.plot(deps, dg_c_mix, 'kx', label=r'$\Delta g_1 + \Delta g_2$')

            else:
                self.dg_plot.plot(deps, dg_mix, 'kx', label=r'$\Delta g_1 + \Delta g_2$')



        if plot_it:
            self.dg_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
        self.canvas.draw()

    def recomp_evb(self):
        """
        Recomputes EVB with current settings (alpha, Hij etc.) from main window.
        """
        selections = list(map(int, self.titles_listbox.curselection()))
        if len(selections) == 0:
            self.app.errorBox('Warning', 'Select titles to recompute EVB.')
        self.titles_listbox.select_clear(0, END)
        existing = self.titles_listbox.get(0, END)

        titles = list()
        for selected in selections:
            titles.append(self.titles_listbox.get(selected))

        for title in titles:
            self.get_binave_out(title, True)


        for title in titles:
            self.titles_listbox.selection_set(existing.index(title))
            self.get_dG(title)
            self.compute_ave_activation_ene(title)

        self.update_tables()
        self.list_titles_event()

    def get_dG(self, title):
        """
        Collects activation and reaction energies from self.title_ene ({title: {nr:[[enegaps],[dG]]}})
        """

        dG_act_all = list()
        dG_rxn_all = list()

        for nr in self.titles_ene[title]:
            rsE, tsE, rE = self.find_dg(title, nr)
            self.titles_dG[title][nr] = [tsE, rE]
            if nr != 'all':
                if tsE != 'na':
                    dG_act_all.append(tsE)
                if rE != 'na':
                    dG_rxn_all.append(rE)

        self.titles_dG[title]['ave'] = [np.average(dG_act_all), np.average(dG_rxn_all)]

    def find_dg(self,title, nr):
        """
        Finds dG(act) and dG(react) in in self.dG_ene[title][nr]
        """
        dG = self.titles_ene[title][nr][1]
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

            #Normalize reaction profile if rE == 0
            if int(rsE) != 0:
                self.titles_ene[title][nr][1] = [-rsE + x for x in self.titles_ene[title][nr][1]]

        except:
            rE = 'na'

        return rsE, tsE, rE

    def update_tables(self):
        """
        Updates all listboxes
        """
        self.update_total_activation_ene_listbox()
        self.update_nonbond_activation_ene_listbox()
        self.update_dg_listbox()

    def update_total_activation_ene_listbox(self):
        """
        Updates listbox with total activation energies
        """

        self.ae_tot_listbox.delete(0, END)
        self.ae_tot_listbox.insert(END, "Title       \N{GREEK CAPITAL LETTER DELTA}U(el) "
                                        "\N{GREEK CAPITAL LETTER DELTA}U(vdW)"
                                        "  \N{GREEK CAPITAL LETTER DELTA}U(bnd)"
                                        "  \N{GREEK CAPITAL LETTER DELTA}U(ang)"
                                        "  \N{GREEK CAPITAL LETTER DELTA}U(tor)"
                                        "  \N{GREEK CAPITAL LETTER DELTA}U(imp)"
                                        " \N{GREEK CAPITAL LETTER DELTA}U(tot) "
                                        "\N{GREEK CAPITAL LETTER DELTA}U(rr+rs)")
        for title in list(self.titles_dU.keys()):
            self.ae_tot_listbox.insert(END, '%10s%8.2f%8.2f%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' % (title.ljust(10),
                                    self.titles_dU[title][0]['EQel'], self.titles_dU[title][0]['EQvdw'],
                                    self.titles_dU[title][0]['EQbnd'], self.titles_dU[title][0]['EQang'],
                                    self.titles_dU[title][0]['EQtor'], self.titles_dU[title][0]['EQimp'],
                                    self.titles_dU[title][0]['EQtot'], self.titles_dU[title][0]['EQtot'] +
                                                                    self.titles_dU[title][0]['alpha'] +
                                                                    self.titles_dU[title][0]['Hij']))

            u_rr_rs_sem = float(np.sqrt(self.titles_dU[title][1]['EQtot']**2 + self.titles_dU[title][1]['alpha']**2 +
                                  self.titles_dU[title][1]['Hij']**2))
            self.ae_tot_listbox.insert(END, '%10s%8.2f%8.2f%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' % (' +/- SEM'.ljust(10),
                                    self.titles_dU[title][1]['EQel'], self.titles_dU[title][1]['EQvdw'],
                                    self.titles_dU[title][1]['EQbnd'], self.titles_dU[title][1]['EQang'],
                                    self.titles_dU[title][1]['EQtor'], self.titles_dU[title][1]['EQimp'],
                                    self.titles_dU[title][1]['EQtot'], u_rr_rs_sem))
            self.ae_tot_listbox.insert(END, '')

    def update_nonbond_activation_ene_listbox(self):
        """
        Updates listbox with total activation energies
        """

        self.ae_nb_listbox.delete(0, END)
        self.ae_nb_listbox.insert(END, "Title       \N{GREEK CAPITAL LETTER DELTA}U(qq_el) "
                                        "\N{GREEK CAPITAL LETTER DELTA}U(qq_vdW)"
                                        " \N{GREEK CAPITAL LETTER DELTA}U(qp_el)"
                                        " \N{GREEK CAPITAL LETTER DELTA}U(qp_vdw)"
                                        " \N{GREEK CAPITAL LETTER DELTA}U(qw_el)"
                                        " \N{GREEK CAPITAL LETTER DELTA}U(qw_vdw)")
        for title in list(self.titles_dU.keys()):
            self.ae_nb_listbox.insert(END, '%10s%9.2f %10.2f%10.2f %10.2f %10.2f %10.2f' % (title.ljust(10),
                                    self.titles_dU[title][0]['qq_el'], self.titles_dU[title][0]['qq_vdw'],
                                    self.titles_dU[title][0]['qp_el'], self.titles_dU[title][0]['qp_vdw'],
                                    self.titles_dU[title][0]['qw_el'], self.titles_dU[title][0]['qw_vdw']))


            self.ae_nb_listbox.insert(END, '%10s%9.2f %10.2f%10.2f %10.2f %10.2f %10.2f' % (' +/- SEM'.ljust(10),
                                    self.titles_dU[title][1]['qq_el'], self.titles_dU[title][1]['qq_vdw'],
                                    self.titles_dU[title][1]['qp_el'], self.titles_dU[title][1]['qp_vdw'],
                                    self.titles_dU[title][1]['qw_el'], self.titles_dU[title][1]['qw_vdw']))
            self.ae_nb_listbox.insert(END, '')

    def update_dg_listbox(self):
        """
        Updates self.dg_listbox with dG(act), dG(rxn) and corresponding SEM
        """
        self.dg_listbox.delete(0, END)

        self.dg_listbox.insert(END, "Title       Temp "
                                       " \N{GREEK CAPITAL LETTER DELTA}G(act)"
                                       "   +/-"
                                       "   \N{GREEK CAPITAL LETTER DELTA}G(rxn)"
                                       "   +/-"
                                       "    pts")
        temp = 300.00
        for title in list(self.titles_dG.keys()):
            got_temp = False
            dg_act = list()
            dg_rxn = list()
            for nr in list(self.titles_dG[title].keys()):
                if nr != 'ave':
                    if not got_temp:
                        with open(self.titles[title][nr]+'/qfep.out', 'r') as qfepout:
                            for line in qfepout:
                                if '--> Give kT & no, of pts to skip: # kT' in line:
                                    temp = round(float(line.split('=')[1]) / 0.001987209, 0)
                                    print(temp)
                                    got_temp = True
                                    break
                    try:
                        dg_act.append(float(self.titles_dG[title][nr][0]))
                        dg_rxn.append(float(self.titles_dG[title][nr][1]))
                    except:
                        continue
            if len(dg_act) > 0:
                act_sem = np.std(dg_act)/np.sqrt(len(dg_act))
                rxn_sem = np.std(dg_rxn)/np.sqrt(len(dg_rxn))
            else:
                act_sem = 0
                rxn_sem = 0
            self.dg_listbox.insert(END, '%10s %5.0f %7.2f %7.2f %7.2f %7.2f %5d' % (title.ljust(10), temp,
                                                                                    np.average(dg_act), act_sem,
                                                                                    np.average(dg_rxn), rxn_sem,
                                                                                    len(dg_act)))

    def compute_ave_activation_ene(self, title):
        """
        Collects all weighted and eigenvector coefficent scaled Q EVB energies from paths and updates global dicts
        """
        self.titles_dU[title] = [dict(), dict()]
        all_dQ = {  'EQtot': list(),
                    'EQbnd':list(),
                    'EQang': list(),
                    'EQtor': list(),
                    'EQimp': list(),
                    'EQel': list(),
                    'EQvdw': list(),
                    'qq_el': list(),
                    'qq_vdw': list(),
                    'qp_el': list(),
                    'qp_vdw': list(),
                    'qw_el': list(),
                    'qw_vdw': list(),
                    'alpha': list(),
                    'Hij': list()}

        runs = 0

        #Collect all activetion potential energy terms
        #dQ = {'EQtot', 'EQbnd', 'EQang', 'EQtor', 'EQimp' , 'EQel', 'EQvdw', 'qq_el', 'qq_vdw',
        #'qp_el', 'qp_vdw', 'qw_el', 'qw_vdw', 'alpha', 'Hij'}
        for nr in list(self.titles[title].keys()):
            if self.titles_dG[title][nr][0] != 'na':
                runs += 1
                dQ = self.compute_activation_ene(title, nr)
                for term in list(dQ.keys()):
                    all_dQ[term].append(dQ[term])

        #Compute averages and standard deviation of the mean
        for term in list(all_dQ.keys()):
            self.titles_dU[title][0][term] = np.average(all_dQ[term])
            self.titles_dU[title][1][term] = np.std(all_dQ[term]) / np.sqrt(runs)

        self.update_tables()

    def compute_activation_ene(self, title, nr):
        """
        Computes dU (TS - RS) from rs_energies and ts_energies derived from self.get_rs_ts_energies. All these energies
        are scaled by the corresponding eigenvector coefficients (c1/c2) and weighted by the lambda contribution to the
        TS and RS bin.

        rs/ts_energies = {lambda_i: {w_i, EQtot, EQbnd, EQang, EQtor, EQimp, EQel, EQvdw, qq_el, qq_vdw, qp_el, qp_vdw,
        qw_el, qw_vdw, alpha, Hij}}
        """

        rs_energies, ts_energies = self.get_rs_ts_energies(title, nr)

        #Compute averages for all lambdas
        #ave_ene = { 'lambda', 'w', 'EQtot', 'EQbnd', 'EQang', 'EQtor', 'EQimp' , 'EQel', 'EQvdw', 'qq_el', 'qq_vdw',
        #'qp_el', 'qp_vdw', 'qw_el', 'qw_vdw', 'alpha', 'Hij'}
        rs_ave = self.compute_avg_lambda_energies(rs_energies)
        ts_ave = self.compute_avg_lambda_energies(ts_energies)

        #Go through dictionaries and
        dQ = dict()
        for term in list(ts_ave.keys()):
            if term != 'lambda' and term != 'w':
                dQ[term] = ts_ave[term] - rs_ave[term]

        return dQ

    def compute_avg_lambda_energies(self, energies):
        """
        Goes trhough dictionary on the form:
        energies = {lambda_i: {w_i, EQtot, EQbnd, EQang, EQtor, EQimp, EQel, EQvdw, qq_el, qq_vdw, qp_el, qp_vdw,
        qw_el, qw_vdw, alpha, Hij}}
        and computes the average for all lambdas in dictionary. Returns dictionary with a single entry.
        Note that average here means the sum of the weighted and eigenvector coefficient scaled energies, since
        sum(w_i) = 1

        ave_ene = { 'lambda', 'w', 'EQtot', 'EQbnd', 'EQang', 'EQtor', 'EQimp' , 'EQel', 'EQvdw', 'qq_el', 'qq_vdw',
        #'qp_el', 'qp_vdw', 'qw_el', 'qw_vdw', 'alpha', 'Hij'}
        """
        ave_ene = { 'lambda': 0,
                    'w': 0,
                    'EQtot': 0,
                    'EQbnd':0,
                    'EQang': 0,
                    'EQtor': 0,
                    'EQimp': 0,
                    'EQel': 0,
                    'EQvdw': 0,
                    'qq_el': 0,
                    'qq_vdw': 0,
                    'qp_el': 0,
                    'qp_vdw': 0,
                    'qw_el': 0,
                    'qw_vdw': 0,
                    'alpha': 0,
                    'Hij': 0}

        for l in list(energies.keys()):
            w = energies[l]['w']
            ave_ene['lambda'] += (w * float(l))
            for qterm in list(energies[l].keys()):
                ave_ene[qterm] += energies[l][qterm]

        return ave_ene

    def get_rs_ts_energies(self, title, nr):
        """
        Goes through part 0 of qfep.out and extracts average EVB potential energies (state 1 + state 2) scaled with
        corresponding eigenvector coefficient and weighted by the lambda contribution (w_i) to RS/TS bin.

        U(lambda_i) = w_i * (c1**2 * u1(lambda_i) + c2**2 * u2(1 - lambda_i) + C2**2 alpha - 2|c1c2|Hij

        returns rs_energies and ts_energies
        rs/ts_energies = {lambda_i: {w_i, EQtot, EQbnd, EQang, EQtor, EQimp, EQel, EQvdw, qq_el, qq_vdw, qp_el, qp_vdw,
        qw_el, qw_vdw, alpha, Hij}}
        """
        pass
        alpha = None
        hij = None

        #Get alpha and Hij from qfep.inp
        qfepinp = '%s/qfep.inp' % self.titles[title][nr]
        if not os.path.isfile(qfepinp):
            self.write_qfep_inp(self.titles[title][nr])
            print(('No qfep.inp in %s. Creating new...' % self.titles[title][nr]))

        linecount = 0
        with open(qfepinp) as qfepinp:
            for line in qfepinp:
                linecount += 1
                if linecount == 6:
                    alpha = float(line.split()[0])
                if linecount == 8:
                    hij = float(line.split()[2])
                if linecount == 9:
                    break

        #rs ==> {lambda_i : [w_i, c1, c2]}
        rs_lambdas, ts_lambdas = self.get_bin_lambdas(title, nr)

        rs_energies = dict()
        ts_energies = dict()
        part0 = self.get_qfep_part0(self.titles[title][nr], 'qfep.out')
        #Go through part 0 and collect weighted Q energies for relevant lambda (rs_lambda/ts_lambda)
        #part0 = [ [[state1], [state2]], [[state1], [state2]]]
        # state1/2 = [l, eqtot, eqbnd, eqang, eqtor, eqimp, eqel, eqvdw, qqel, qqvdw, qpel, qpvdw, qwel, qwvdw]
        for i in range(len(part0)):
            l = part0[i][0][0]
            found_lambda = False
            if l in list(rs_lambdas.keys()):
                energies = rs_energies
                lambdas = rs_lambdas
                found_lambda = True
            if l in list(ts_lambdas.keys()):
                energies = ts_energies
                lambdas = ts_lambdas
                found_lambda = True

            if found_lambda:
                w_i, c1, c2 = lambdas[l][0:]

                energies[l] = {'w': w_i,
                               'EQtot': w_i * (c1 * part0[i][0][1] + c2 * part0[i][1][1]),
                               'EQbnd': w_i * (c1 * part0[i][0][2] + c2 * part0[i][1][2]),
                               'EQang': w_i * (c1 * part0[i][0][3] + c2 * part0[i][1][3]),
                               'EQtor': w_i * (c1 * part0[i][0][4] + c2 * part0[i][1][4]),
                               'EQimp': w_i * (c1 * part0[i][0][5] + c2 * part0[i][1][5]),
                               'EQel': w_i * (c1 * part0[i][0][6] + c2 * part0[i][1][6]),
                               'EQvdw': w_i * (c1 * part0[i][0][7] + c2 * part0[i][1][7]),
                               'qq_el': w_i * (c1 * part0[i][0][8] + c2 * part0[i][1][8]),
                               'qq_vdw': w_i * (c1 * part0[i][0][9] + c2 * part0[i][1][9]),
                               'qp_el': w_i * (c1 * part0[i][0][10] + c2 * part0[i][1][10]),
                               'qp_vdw': w_i * (c1 * part0[i][0][11] + c2 * part0[i][1][11]),
                               'qw_el': w_i * (c1 * part0[i][0][12] + c2 * part0[i][1][12]),
                               'qw_vdw': w_i * (c1 * part0[i][0][13] + c2 * part0[i][1][13]),
                               'alpha': w_i * (c2 * alpha),
                               'Hij': w_i * (-2.0 * np.sqrt(c1) * np.sqrt(c2) * hij)}


        return rs_energies, ts_energies

    def get_bin_lambdas(self, title, nr):
        """
        Goes through part 2 of qfep.out and collects all lambda values and eigenvector coefficients that contributes
        to the transition and reactant state bins. returns a RS and TS dictionary {lambda_i: [w_i, c1, c2]} where
        w_i is the weight of the lambda-value (points from lambda to bin / sum points in bin) and c1/c2 are the
        eignevector coefficients.
        """
        bin_rs, bin_ts = self.find_rs_ts_bins(title, nr)

        #rs ==> {lambda_i : [pts_i, c1, c2]}
        rs = dict()
        rs_pts = 0
        ts = dict()
        ts_pts = 0

        qfepout = '%s/qfep.out' % self.titles[title][nr]

        found_part2 = False
        bin_found = False
        with open(qfepout,'r') as qfepout:
            for line in qfepout:
                if found_part2:
                    if not line.startswith('#') and len(line.split()) > 7:
                        lamda = line.split()[0]
                        pts = float(line.split()[6])
                        c1 = float(line.split()[7])
                        c2 = float(line.split()[8])
                        if line.split()[1] == bin_rs:
                            rs[lamda] = [pts, c1, c2]
                            rs_pts += pts

                        elif line.split()[1] == bin_ts:
                            ts[lamda] = [pts, c1, c2]
                            ts_pts += pts

                    if len(line.split()) < 6 or '# Part 3:' in line:
                        break
                if '# Part 2' in line:
                    found_part2 = True

        #rs ==> {lambda_i : [w_i, c1, c2]}
        for state in (rs, ts):
            if state == rs:
                scale = float(rs_pts)
            else:
                scale = float(ts_pts)
            for lamda in list(state.keys()):
                state[lamda][0] /= scale

        return rs, ts

    def find_rs_ts_bins(self, title, nr):
        """
        Finds the bins corresponding to the reactant and transition state (qfep.out).
        Returns bin(RS) and bin(TS), alpha and Hij
        """
        qfepinp = '%s/qfep.inp' % self.titles[title][nr]
        if not os.path.isfile(qfepinp):
            self.write_qfep_inp(self.titles[title][nr])
            print(('No qfep.inp in %s. Creating new...' % self.titles[title][nr]))

        qfepout = '%s/qfep.out' % self.titles[title][nr]
        if not os.path.isfile(qfepout):
            print(('No qfep.out in %s. Running Qfep to generate file...' % self.titles[title][nr]))
            self.run_qfep(self.titles[title][nr], 'qfep.inp')

        #Go through part3 from qfep.out and find bin corresponding to RS and TS
        part3 = self.get_qfep_part3(self.titles[title][nr], 'qfep.out')

        bin_rs = None
        bin_ts = None
        found_TS = False
        rsE = float(part3[0].split()[3])
        tsE = 'na'

        for i in range(1, len(part3) - 1):
            prev = float(part3[i - 1].split()[3])
            mid = float(part3[i].split()[3])
            nxt = float(part3[i + 1].split()[3])

            if mid <= prev and mid <= nxt:
                if not found_TS:
                    if mid < rsE:
                        rsE = mid
                        bin_rs = part3[i].split()[0]

            elif mid >= prev and mid >= nxt:
                if not found_TS:
                    if abs(mid) > 0:
                        tsE = mid
                        bin_ts = part3[i].split()[0]
                        found_TS = True
                elif found_TS:
                    if mid > tsE:
                        tsE = mid
                        bin_ts = part3[i].split()[0]
        if not bin_ts or not bin_rs:
            print('Problems locating reactant and transition state ... ')

        return bin_rs, bin_ts

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
                            #OLD QFEP FORMAT:
                            if len(line.split()) < 17:
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


                            #New and nicer Qfep format
                            else:
                                l = line.split()[3]
                                eqtot = float(line.split()[4])
                                eqbnd = float(line.split()[5])
                                eqang = float(line.split()[6])
                                eqtor = float(line.split()[7])
                                eqimp = float(line.split()[8])
                                eqel = float(line.split()[9])
                                eqvdw = float(line.split()[10])
                                qqel = float(line.split()[11])
                                qqvdw = float(line.split()[12])
                                qpel = float(line.split()[13])
                                qpvdw = float(line.split()[14])
                                qwel = float(line.split()[15])
                                qwvdw = float(line.split()[16])

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

    def get_qfep_part2(self, path, filename):
        """
        Returns dictionary with energy gaps as keys and for each key a list on the form
        [dg1, dg2, c1**2, c2**2]
        A key with alpha and H12 will also be generated
        """
        de_dg = dict()
        found_p2 = False
        qfepout = '%s/%s' % (path, filename)

        hij = 0
        dg1 = 3
        dg2 = 4
        c1 = 7
        c2 = 8

        with open(qfepout, 'r') as p2:
            for i in p2:
                if '--> i, j, A_ij, mu_ij, eta_ij, r_xy0: #' in i:
                    hij = float(i.split()[10])

                if '# Linear combination co-efficients' in i:
                    try:
                        if float(i.split('=')[-1].split()[0]) == -1.0:
                            dg1 = 4
                            dg2 = 3
                            c1 = 8
                            c2 = 7
                    except:
                        continue
                if '# Part 3: Bin-averaged summary:' in i:
                    break
                if found_p2:
                    if len(i.split()) > 6 and '#' not in i:
                        de_dg[float(i.split()[2])] = [float(i.split()[dg1]), float(i.split()[dg2]),
                                                      float(i.split()[c1]), float(i.split()[c2]), hij]
                if '# Part 2: Reaction free energy summary:' in i:
                    found_p2 = True

        return de_dg

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
                    found_part3 = False
                    break
                if found_part3:
                    if not line.startswith('#') and len(line.split()) > 4:
                        part3.append(line)
                    if len(line.split()) < 4:
                        break
                if '# Part 3:' in line:
                    found_part3 = True

        return part3

    def get_parab_data(self, de_dg12):
        """
        Sorts energy gaps from smalles to biggest and generate lists with points to crossing point.
        These data points are used to make the fitted marcus parabolas.
        """
        e1_list = list()
        dg1_list = list()
        e2_list = list()
        dg2_list = list()

        e1_cont = list()
        dg1_cont = list()
        e2_cont = list()
        dg2_cont = list()

        for i in sorted(de_dg12.keys()):
            dg1 = de_dg12[i][0]
            dg2 = de_dg12[i][1]

            if dg1 <= dg2:
                e1_list.append(i)
                dg1_list.append(dg1)
                e2_cont.append(i)
                dg2_cont.append(dg2)
            else:
                e2_list.append(i)
                dg2_list.append(dg2)
                e1_cont.append(i)
                dg1_cont.append(dg1)

        return e1_list, dg1_list, e2_list, dg2_list, e1_cont, dg1_cont, e2_cont, dg2_cont

    def f_parabel(self, prm, x):
        """
        2nd order polynomial function where prm is a list with the three constants a,b,c (ax^2 + bx + c).
        Returns parabel value for x.
        """
        return (prm[0] * x**2) + (prm[1] * x) + prm[2]

    def reorgE_intrinsic(self, norm1, norm2):
        """
        Finds energy gap when the two normalized Marcus parabalas
        are equal (lambda/4) and returns the reorganization energy.
        """
        a = norm1[0] - norm2[0]
        b = norm1[1] - norm2[1]
        c = norm1[2] - norm2[2]
        de = (-b + np.sqrt(b**2 - (4*a*c)))/ (2*a)
        de2 = (-b - np.sqrt(b**2 - (4*a*c)))/ (2*a)

        l_4 = (norm1[0] * de**2) + (norm1[1] * de) + norm1[2]
        l_4_2 = (norm1[0] * de2**2) + (norm1[1] * de2) + norm1[2]

        if l_4_2 < l_4:
            l_4 = l_4_2

        return 4.*l_4

    def normalize_parabel(self, prm):
        """
        Shifts a 2nd-order polynom to f(x) = 0 when f'(x) = 0. Returns new parameters.
        f(x) = ax**2 + bx + c
        prm = [a,b,c]
        """
        a,b,c = prm[0:]

        f_shift = (b**2 - (2 * b**2)) / (4 * a) + c

        return [a, b, (c - f_shift)]

    def compute_reorganization(self):
        """
        Extractes the diabatic free energy functions, with corresponding eigenvector coefficients and computes
        the (intrinsic) reaorganization energy by fitting marcus parabolas to the free energy functions. This is done
        for all runs in the selected title.

        {title: {run: 'diabatic': {energy_gap: [dg1, dg2, c1, c2, Hij], 'marcus1': [[e], [dg]], 'marcus2': {[e]:[dg]},
         'norm1':[[e], [dg]], 'norm2': [[e], [dg]], 'reorg': value}}}
        """

        #Get selected title(s)
        selections = list(map(int, self.titles_listbox.curselection()))
        if len(selections) == 0:
            return

        #Collect dibatic energies for all runs in title(s)
        for selected in selections:
            title = self.titles_listbox.get(selected)
            self.app.log(' ', '\n\nComputing reorganization energies for %s\n' % title)
            if title not in list(self.title_reorge.keys()):
                self.title_reorge[title] = dict()
            for i in list(self.titles[title].keys()):
                path = self.titles[title][i]
                self.app.log(' ', '%s\n' % path)
                print(('%s' % path))
                self.title_reorge[title][path] = dict()
                de_dg = self.get_qfep_part2(path, 'qfep.out')
                self.title_reorge[title][path]['diabatic'] = de_dg

                e1, dg1, e2, dg2, e1_2, dg1_2, e2_2, dg2_2 = self.get_parab_data(de_dg)

                #Generate energy gaps for parapbel function
                e1_mod = np.arange(min(e1) * 1.6, max(e1) + 60, 10)
                e2_mod = np.arange(min(e2) - 60, 1.6 * max(e2), 10)

                #Genere fit parameters (a,b,c)
                dg1_prm = np.polyfit(e1, dg1, 2)
                dg2_prm = np.polyfit(e2, dg2, 2)

                marc1 = self.f_parabel(dg1_prm, e1_mod)
                marc2 = self.f_parabel(dg2_prm, e2_mod)

                self.title_reorge[title][path]['marcus1'] = [e1_mod, marc1]
                self.title_reorge[title][path]['marcus2'] = [e2_mod, marc2]

                #Generate marcus parabolas with equal bottom point
                #dg1_norm = np.polyfit(e1, map(lambda x: x - min(dg1 + dg1_2), dg1), 2)
                #dg2_norm = np.polyfit(e2, map(lambda x: x - min(dg2 + dg2_2), dg2), 2)

                dg1_norm = self.normalize_parabel(dg1_prm)
                dg2_norm = self.normalize_parabel(dg2_prm)

                norm1 = self.f_parabel(dg1_norm, e1_mod)
                norm2 = self.f_parabel(dg2_norm, e2_mod)

                self.title_reorge[title][path]['norm1'] = [e1_mod, norm1]
                self.title_reorge[title][path]['norm2'] = [e2_mod, norm2]

                reorge = self.reorgE_intrinsic(dg1_norm, dg2_norm)

                self.title_reorge[title][path]['reorg'] = reorge

            self.app.log(' ', '\n\nCompleted %s\n' % title)

        self.update_diabatic_plot()

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

        titles = list()
        plot_titles = list()
        d_G = list()
        d_eps = list()
        dG_act = list()
        dG_rxn = list()

        try:
            for selected in selections:
                title = self.titles_listbox.get(selected)
                titles.append(title)
                dG_act.append(float(self.titles_dG[title]['ave'][0]))
                dG_rxn.append(float(self.titles_dG[title]['ave'][1]))
                print(dG_act)
                try:
                    tmp_dg = self.titles_ene[title]['all'][1]
                    tmp_eps = self.titles_ene[title]['all'][0]
                    d_G.append(tmp_dg)
                    d_eps.append(tmp_eps)
                    plot_titles.append(title)
                except:
                    continue

        except:
            return
        print(titles)
        #Show runs for titles:
        self.runs_listbox.delete(0, END)
        for title in titles:
            for run in sorted([x for x in self.titles[title].keys() if type(x)==int]) + [x for x in self.titles[title].keys() if type(x)!=int]:
                if run != 'all':
                    self.runs_listbox.insert(END,'../%s' % '/'.join(self.titles[title][run].split('/')[-3:]))

        if len(d_eps) == 0:
            for i in range(len(titles)):
                d_G.append([0])
                d_eps.append([0])

        self.ave_dg_act.config(state=NORMAL)
        self.ave_dg_act.delete(0.0, END)
        self.ave_dg_rxn.config(state=NORMAL)
        self.ave_dg_rxn.delete(0.0, END)

        if len(dG_act) != 0:
            #Update average dG values:
            self.ave_dg_act.insert(0.0, '%.2f' % np.average(dG_act))
            self.ave_dg_rxn.insert(0.0, '%.2f' % np.average(dG_rxn))
        else:
            self.ave_dg_act.insert(0.0, '--')
            self.ave_dg_rxn.insert(0.0, '--')

        self.ave_dg_act.config(state=DISABLED)
        self.ave_dg_rxn.config(state=DISABLED)

        #Update current dG values, if exactly one is selected:
        self.dg_act.config(state=NORMAL)
        self.dg_act.delete(0.0, END)
        self.dg_react.config(state=NORMAL)
        self.dg_react.delete(0.0, END)

        if len(titles) == 1:
            if 'all' in list(self.titles_dG[title].keys()):
                self.dg_act.insert(0.0, '%.2f' % self.titles_dG[title]['all'][0])
                self.dg_react.insert(0.0, '%.2f' % self.titles_dG[title]['all'][1])
            else:
                self.dg_act.insert(0.0, '--')
                self.dg_react.insert(0.0, '--')
        else:
            self.dg_act.insert(0.0, '--')
            self.dg_react.insert(0.0, '--')

        self.dg_react.config(state=DISABLED)
        self.dg_act.config(state=DISABLED)

        #plot average Qfep for titles:
        if not self.diabatic_plot:
            self.update_plot(d_eps, d_G, plot_titles)
        else:
            self.update_diabatic_plot()

    def list_runs_event(self, *args):
        selections = list(map(int, self.runs_listbox.curselection()))
        if len(selections) == 0:
            return

        #Show average profile or not?
        show_average = True

        titles = list()
        d_G = list()
        d_eps = list()
        dg_act = list()
        dg_rxn = list()

        dg_act_title = list()
        dg_rxn_title = list()

        for selected in selections:
            search = '/'.join(self.runs_listbox.get(selected).split('/')[-3:])
            for title in list(self.titles.keys()):
                for nr in self.titles[title]:
                    if '/'.join(self.titles[title][nr].split('/')[-3:]) == search:
                        if show_average:
                            if title not in titles:
                                if 'all' in list(self.titles_ene[title].keys()):
                                    titles.append(title)
                                    d_G.append(self.titles_ene[title]['all'][1])
                                    d_eps.append(self.titles_ene[title]['all'][0])
                                if self.titles_dG[title]['ave'][0]:
                                    dg_act_title.append(float(self.titles_dG[title]['ave'][0]))
                                    dg_rxn_title.append(float(self.titles_dG[title]['ave'][1]))
                        titles.append('%s_%s' % (title, nr))
                        d_G.append(self.titles_ene[title][nr][1])
                        d_eps.append(self.titles_ene[title][nr][0])
                        print((self.titles_dG[title][nr][0]))
                        if self.titles_dG[title][nr][0] != 'na':
                            dg_act.append(self.titles_dG[title][nr][0])
                            dg_rxn.append(self.titles_dG[title][nr][1])

                        break
        if not self.diabatic_plot:
            self.update_plot(d_eps, d_G, titles)
        else:
            self.update_diabatic_plot()


        #Update dG values
        self.ave_dg_act.config(state=NORMAL)
        self.ave_dg_rxn.config(state=NORMAL)
        self.dg_act.config(state=NORMAL)
        self.dg_react.config(state=NORMAL)
        self.dg_act.delete(0.0, END)
        self.dg_react.delete(0.0, END)
        self.ave_dg_act.delete(0.0, END)
        self.ave_dg_rxn.delete(0.0, END)

        if len(selections) > 1:
            self.ave_dg_act.delete(0.0, END)
            self.ave_dg_rxn.delete(0.0, END)
            try:
                self.ave_dg_act.insert(0.0, '%.2f' % np.average(dg_act))
                self.ave_dg_rxn.insert(0.0, '%.2f' % np.average(dg_rxn))
            except:
                self.ave_dg_act.insert(0.0, '--')
                self.ave_dg_rxn.insert(0.0, '--')
            self.dg_act.insert(0.0, '--')
            self.dg_react.insert(0.0, '--')
        elif len(selections) == 1:
            try:
                self.dg_act.insert(0.0, dg_act[0])
                self.dg_react.insert(0.0, dg_rxn[0])
            except:
                self.dg_act.insert(0.0, '--')
                self.dg_react.insert(0.0, '--')
            self.ave_dg_act.insert(0.0, '%.2f' % np.average(dg_act_title))
            self.ave_dg_rxn.insert(0.0, '%.2f' % np.average(dg_rxn_title))

        self.ave_dg_act.config(state=DISABLED)
        self.ave_dg_rxn.config(state=DISABLED)
        self.dg_act.config(state=DISABLED)
        self.dg_react.config(state=DISABLED)

    def show_var_frame(self, *args):
        frames = {'Reaction Free Energy Profiles': self.plot_frame,
                  'Reaction Free energies': self.dg_frame,
                  'Activation Energies Total': self.ae_tot_frame,
                  'Activation Energies Non Bonded': self.ae_nb_frame,
                  'Reorganization Energies': self.plot_frame}

        for i in list(frames.keys()):
            frames[i].grid_forget()
        try:
            frames[self.selected_frame.get()].grid(row=3, column=0, columnspan=2)
            if self.selected_frame.get() == 'Reorganization Energies':
                self.diabatic_plot = True
                self.reorg_frame.grid(row=4, column=0, columnspan=2)
                if self.dg_plot:
                    self.dg_plot.clear()
                    self.canvas.draw()
            else:
                self.diabatic_plot = False
                self.reorg_frame.grid_forget()
        except:
            pass

    def dialog_window(self):

        self.title('EVB Reaction Energies')
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
        bottomframe.grid(row=5, column=0)

        #Variable frames:

        #Project entries and energy files
        self.proj_frame = Frame(self.mainframe, bg=self.main_color)
        self.proj_frame.grid(row=1, column=0, pady=(0,10))

        self.ae_tot_frame = Frame(self.mainframe, bg=self.main_color)

        self.ae_nb_frame = Frame(self.mainframe, bg=self.main_color)

        self.dg_frame = Frame(self.mainframe, bg=self.main_color)

        self.reorg_frame = Frame(self.mainframe, bg=self.main_color)


        #topframe content
        import_label = Label(topframe, text='Import: ', bg=self.main_color)
        import_label.grid(row=0, column=0)

        qfep_import = Button(topframe, text='EVB parameters', highlightbackground=self.main_color, command=self.import_qfep)
        qfep_import.grid(row=0, column=1, columnspan=2)

        proj_import = Button(topframe, text='Project', highlightbackground=self.main_color, command=self.import_project)
        proj_import.grid(row=0, column=3)

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
        self.kT.insert(0, '0.596')

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
                                   'Reaction Free Energy Profiles', 'Reaction Free energies',
                                   'Activation Energies Total', 'Activation Energies Non Bonded',
                                   'Reorganization Energies')
        self.view_frame.config(highlightbackground=self.main_color, bg=self.main_color, width=30)
        self.view_frame.grid(row=4, column=0)

        #Project frame
        entry_label = Label(self.proj_frame, text='Titles', bg=self.main_color)
        entry_label.grid(row=0, column=0)

        add_title = Button(self.proj_frame, text='+', highlightbackground=self.main_color, command=self.add_title)
        add_title.grid(row=0, column=1)

        del_title = Button(self.proj_frame, text='-', highlightbackground=self.main_color, command=self.del_title)
        del_title.grid(row=0, column=2)

        runs_label = Label(self.proj_frame, text='Runs', bg=self.main_color)
        runs_label.grid(row=0, column=4)

        add_runs = Button(self.proj_frame, text='+', highlightbackground=self.main_color, command=self.add_runs)
        add_runs.grid(row=0, column=5)

        del_title = Button(self.proj_frame, text='-', highlightbackground=self.main_color, command=self.del_runs)
        del_title.grid(row=0, column=6)


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
                                      width=25, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        runs_yscroll.config(command=self.runs_listbox.yview)
        runs_xscroll.config(command=self.runs_listbox.xview)
        self.runs_listbox.grid(row=1, rowspan=10, column = 4, columnspan=3, sticky='e')
        self.runs_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        self.runs_listbox.bind('<<ListboxSelect>>', self.list_runs_event)


        ave_dG = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        ave_dG.tag_configure("superscript", offset=3)
        ave_dG.insert("insert","<\N{GREEK CAPITAL LETTER DELTA}G","","\\u2021",'superscript')
        ave_dG.insert(END,'>')
        ave_dG.grid(row=1, column=8)
        ave_dG.config(state=DISABLED)

        self.ave_dg_act = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.ave_dg_act.grid(row=1, column=9)
        self.ave_dg_act.insert(0.0, '--')

        ave_dG_rxn = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        ave_dG_rxn.tag_configure("superscript", offset=-3)
        ave_dG_rxn.insert("insert","<\N{GREEK CAPITAL LETTER DELTA}G","","o",'superscript')
        ave_dG_rxn.insert(END,'>')
        ave_dG_rxn.grid(row=3, column=8)
        ave_dG_rxn.config(state=DISABLED)

        self.ave_dg_rxn = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.ave_dg_rxn.grid(row=3, column=9)
        self.ave_dg_rxn.insert(0.0, '--')


        real_dG = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        real_dG.tag_configure("superscript", offset=3)
        real_dG.insert("insert","\N{GREEK CAPITAL LETTER DELTA}G","","\\u2021",'superscript')
        real_dG.grid(row=5, column=8)
        real_dG.config(state=DISABLED)

        self.dg_act = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.dg_act.grid(row=5, column=9)
        self.dg_act.insert(0.0, '--')

        react_dG = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        react_dG.tag_configure("superscript", offset=-3)
        react_dG.insert("insert","\N{GREEK CAPITAL LETTER DELTA}G","","o",'superscript')
        react_dG.grid(row=7, column=8)
        react_dG.config(state=DISABLED)

        self.dg_react = Text(self.proj_frame, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.dg_react.grid(row=7, column=9)
        self.dg_react.insert(0.0, '--')


        recomp_evb = Button(self.proj_frame, text='Re-compute EVB', highlightbackground=self.main_color,
                            command=self.recomp_evb)
        recomp_evb.grid(row=0, column=8, columnspan=2)

        make_avg = Button(self.proj_frame, text='Make average profile', highlightbackground=self.main_color,
                          command=self.make_evb_average)
        make_avg.grid(row=10, column=8, columnspan=2)

        ae_tot_yscroll = Scrollbar(self.ae_tot_frame)
        ae_tot_yscroll.grid(row = 1, rowspan=16, column = 3, sticky = 'nsw', padx=(0,10))

        self.ae_tot_listbox = Listbox(self.ae_tot_frame, yscrollcommand = ae_tot_yscroll.set,
                                      width=82, height=16, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        ae_tot_yscroll.config(command=self.ae_tot_listbox.yview)
        self.ae_tot_listbox.grid(row=1, rowspan=16, column = 0, columnspan=3, sticky='e')
        self.ae_tot_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        #self.ae_tot_listbox.bind('<<ListboxSelect>>', self.list_runs_event)

        export_tot = Button(self.ae_tot_frame, text='Export table', highlightbackground=self.main_color,
                           command=lambda: self.export_table(self.ae_tot_listbox))
        export_tot.grid(row=17, column=0, columnspan=3, sticky='e')

        #Activation energies non bonded
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
                                      width=82, height=16, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        dg_yscroll.config(command=self.dg_listbox.yview)
        self.dg_listbox.grid(row=1, rowspan=16, column = 0, columnspan=3, sticky='e')
        self.dg_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        #self.ae_tot_listbox.bind('<<ListboxSelect>>', self.list_runs_event)

        export_dg = Button(self.dg_frame, text='Export table', highlightbackground=self.main_color,
                           command=lambda: self.export_table(self.dg_listbox))
        export_dg.grid(row=17, column=0, columnspan=3, sticky='e')

        #Reorganization Energies
        hii_plot_check = Checkbutton(self.reorg_frame, bg=self.main_color, variable=self.plot_hii,
                                          command=self.update_diabatic_plot)
        hii_plot_check.grid(row=0, column=0, sticky='e')

        hii = Text(self.reorg_frame, width=11, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        hii.tag_configure("subscript", offset=-1)
        hii.insert("insert","\N{GREEK CAPITAL LETTER DELTA}g","",'i,j','subscript')
        hii.grid(row=0, column=1, sticky='w')
        hii.config(state=DISABLED)

        fit_plot_check = Checkbutton(self.reorg_frame, bg=self.main_color, variable=self.plot_fit,
                                          command=self.update_diabatic_plot)
        fit_plot_check.grid(row=0, column=2, sticky='e')

        fit = Text(self.reorg_frame, width=11, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        fit.insert("insert", 'Fit')
        fit.grid(row=0, column=3, sticky='w')
        fit.config(state=DISABLED)

        fitnorm_plot_check = Checkbutton(self.reorg_frame, bg=self.main_color, variable=self.plot_norm_fit,
                                          command=self.update_diabatic_plot)
        fitnorm_plot_check.grid(row=0, column=4, sticky='e')

        fit_norm = Text(self.reorg_frame, width=11, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        fit_norm.insert("insert", 'Fit shifted')
        fit_norm.grid(row=0, column=5, sticky='w')
        fit_norm.config(state=DISABLED)

        cj_plot_check = Checkbutton(self.reorg_frame, bg=self.main_color, variable=self.plot_mix,
                                          command=self.update_diabatic_plot)
        cj_plot_check.grid(row=1, column=0, sticky='e')

        mix = Text(self.reorg_frame, width=11, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        mix.insert("insert", "\N{GREEK CAPITAL LETTER DELTA}gi + \N{GREEK CAPITAL LETTER DELTA}gj")
        mix.grid(row=1, column=1, sticky='w')
        mix.config(state=DISABLED)

        ci_plot_check = Checkbutton(self.reorg_frame, bg=self.main_color, variable=self.plot_c,
                                          command=self.update_diabatic_plot)
        ci_plot_check.grid(row=1, column=2, sticky='e')

        ci = Text(self.reorg_frame, width=11, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        ci.tag_configure("subscript", offset=-1)
        ci.insert("insert",'C',"",'i,j','subscript')
        ci.grid(row=1, column=3, sticky='w')
        ci.config(state=DISABLED)

        hij_plot_check = Checkbutton(self.reorg_frame, bg=self.main_color, variable=self.plot_hij,
                                          command=self.update_diabatic_plot)
        hij_plot_check.grid(row=1, column=4, sticky='e')

        hij = Text(self.reorg_frame, width=11, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        hij.tag_configure("subscript", offset=-1)
        hij.insert("insert",'H',"",'ij','subscript')
        hij.grid(row=1, column=5, sticky='w')
        hij.config(state=DISABLED)





        avg_lamdbda = Text(self.reorg_frame, width=4, height=1, bg=self.main_color, borderwidth=0,
                           highlightthickness=0)
        avg_lamdbda.insert('insert', "<\N{GREEK small LETTER LAMDA}>")
        avg_lamdbda.grid(row=0, column=6, sticky='e', padx=(20,5))
        avg_lamdbda.config(state=DISABLED)

        curr_lamdbda = Text(self.reorg_frame, width=4, height=1, bg=self.main_color, borderwidth=0,
                           highlightthickness=0)
        curr_lamdbda.insert('insert', " \N{GREEK small LETTER LAMDA} ")
        curr_lamdbda.grid(row=1, column=6, sticky='e', padx=(20,5))
        curr_lamdbda.config(state=DISABLED)

        self.reorg_avg = Text(self.reorg_frame, width=8, height=1, bg=self.main_color, borderwidth=0,
                              highlightthickness=0)
        self.reorg_avg.insert('insert', '--')
        self.reorg_avg.grid(row=0, column=7)
        self.reorg_avg.config(state=DISABLED)

        reorg_se = Text(self.reorg_frame, width=4, height=1, bg=self.main_color, borderwidth=0,
                           highlightthickness=0)
        reorg_se.insert('insert', "+/-")
        reorg_se.grid(row=0, column=8, sticky='e', padx=(5,5))
        reorg_se.config(state=DISABLED)

        self.reorg_se = Text(self.reorg_frame, width=8, height=1, bg=self.main_color, borderwidth=0,
                              highlightthickness=0)
        self.reorg_se.insert('insert', '--')
        self.reorg_se.grid(row=0, column=9)
        self.reorg_se.config(state=DISABLED)

        self.reorg_curr = Text(self.reorg_frame, width=8, height=1, bg=self.main_color, borderwidth=0,
                               highlightthickness=0)
        self.reorg_curr.insert('insert', '--')
        self.reorg_curr.grid(row=1, column=7)
        self.reorg_curr.config(state=DISABLED)

        compute_reorg = Button(self.reorg_frame, text='Compute', highlightbackground=self.main_color,
                               command=self.compute_reorganization)
        compute_reorg.grid(row=1, column=8, columnspan=2)



        #Bottom frame Quit/save
        quit_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        quit_button.grid(row=0, column=0)
