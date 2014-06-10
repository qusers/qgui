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

from Tkinter import Label, Button, Listbox, Scrollbar, EXTENDED, Spinbox, Entry, Text, Frame, \
    Toplevel, END, GROOVE, StringVar, OptionMenu, HORIZONTAL

from tkSimpleDialog import askstring
import tkFont
from subprocess import call
import numpy as np
import os
from Qplot import Qplot

from tkFileDialog import askdirectory


class AnalyzeFEP(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root

        self.main_color = self.app.main_color

        #Global dicts
        #{title : {1: path, 2: path...}}
        self.titles = dict()

        #{title: 1: {lambda: [], dgf: [], dgr: [], dg:[]}}
        self.titles_dg = dict()

        self.states_var = StringVar()

        self.selected_frame = StringVar()

        self.selected_frame.set('Show ...')
        self.selected_frame.trace('w', self.show_var_frame)


        self.dialog_window()
        self.states_var.trace('w', self.states_changed)

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
        self.titles_dg[title] = dict()

        #Highlight latest title in listbox
        self.titles_listbox.select_clear(0, END)
        titles = self.titles_listbox.get(0, END)
        for ind_ in range(len(titles)):
            if titles[ind_] == title:
                self.titles_listbox.selection_set(ind_)
                #self.list_titles_event()

        #self.update_plot([[0]],[[0]],[title])

    def del_title(self):
        selections = map(int, self.titles_listbox.curselection())
        for selected in selections:
            title = self.titles_listbox.get(selected)
            self.titles_listbox.delete(selected)
            del self.titles[title]
            for i in [self.titles_dg]:
                if title in i.keys():
                    del i[title]

        self.update_tables()

    def add_runs(self):
        """
        Add subrun one-by-one or select temperature directory to append all runs
        """
        try:
            prj_entry = map(int, self.titles_listbox.curselection())
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
        title = 'Select FEP temperature or run directory'
        path = self.app.workdir
        while True:
            rundir = askdirectory(parent=self, mustexist=True, title=title, initialdir= path)
            if not rundir:
                break
            #Check if the first directory is a temperature directory
            if len(dirs) == 0:
                try:
                    print int(float(rundir.split('/')[-1]))
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
                for ene in os.listdir(rundir):
                    if ene.endswith('.en'):
                        if os.path.getsize('%s/%s' % (rundir, ene)) != 0:
                            enefiles = True
                            break
                        else:
                            break
                if enefiles:
                    dirs.append(rundir)
                    title = '%s Added. Select next dir or cancel' % dirs[-1].split('/')[-1]
                    self.app.log(' ', '%s Added. Select next dir or cancel to finish\n'
                                      % '/'.join(dirs[-1].split('/')[-3:]))
                    try:
                        self.titles[prj_title][int(rundir.split('/')[-1])] = rundir
                        self.titles_dg[prj_title][int(rundir.split('/')[-1])] = dict()
                    except:
                        self.titles[prj_title][rundir.split('/')[-1]] = rundir
                        self.titles_dg[prj_title][rundir.split('/')[-1]] = dict()
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
                                break
                            else:
                                break
                    if enefiles:
                        dirs.append(subdir)
                        self.titles[prj_title][nr_dir] = subdir
                        self.titles_dg[prj_title][nr_dir] = dict()
                        self.app.log(' ','...../%s added\n' % '/'.join(subdir.split('/')[-2:]))
                        print subdir
                    nr_dir += 1
                else:
                    break

        #Get all dG(act) dG(react) and normalize
        self.get_dg(prj_title)

        self.list_titles_event()


    def del_runs(self):
        title_sel = map(int, self.titles_listbox.curselection())
        if len(title_sel) != 1:
            self.app.log(' ', '\nSelect exactly one title to delete runs from!\n')
            return

        selections = map(int, self.runs_listbox.curselection())
        if len(selections) == 0:
            self.app.log(' ', '\nNo runs selected for deletion!\n')
            return

        title = self.titles_listbox.get(title_sel[0])

        for selected in reversed(selections):
            search = '/'.join(self.runs_listbox.get(selected).split('/')[-3:])
            self.runs_listbox.delete(selected)

            for nr in self.titles[title].keys():
                if '/'.join(self.titles[title][nr].split('/')[-3:]) == search:
                     for i in [self.titles, self.titles_dg]:
                        if title in i.keys():
                            if nr in i[title].keys():
                                del i[title][nr]


        self.update_tables()
        self.list_titles_event()

    def write_qfep_inp(self, path):
        """
        Write qfep.inp in path for md*.en
        """

        states = int(self.states_entry.get())

        found_ene_files = True
        #Collect all md*.en files:
        enefiles = list()
        for qfile in os.listdir(path):
            if qfile.endswith('.en'):
                if qfile.startswith('md'):
                    enefiles.append(qfile)

        print 'Found %d MD energy files in ../%s' % (len(enefiles), '/'.join(path.split('/')[-3:]))

        if len(enefiles) == 0:
            found_ene_files = False
            return found_ene_files

        #Sort energy files (if not Qfep will give messed up free energy profile!)
        enefiles = sorted(enefiles, key=lambda x: x.split('/')[-1])

        #make input for Qfep
        inputfile = open('%s/qfep.inp' % path,'w')
        inputfile.write('%d\n' % len(enefiles))
        inputfile.write('%d  0\n' % states)
        inputfile.write('%s  %s\n' % (self.kT.get(), self.skip_points.get()))
        inputfile.write('%s\n' % self.bins.get())
        inputfile.write('%s\n' % self.binpoints_min.get())
        #Write states - 1 alpha values = 0
        for i in range(states - 1):
            inputfile.write('0\n')
        #Number of off-diagonal elements
        inputfile.write('0\n')
        inputfile.write('%s\n' % self.linear_comb.get())

        for ene in enefiles:
            inputfile.write('%s\n' % ene.split('/')[-1])

        inputfile.close()

        return found_ene_files

    def run_qfep(self, path, qfep_inp):
        """
        Calls Qfep with qfep_inp and produces qfep.out
        """
        if not os.path.isfile('%s/%s' % (path, qfep_inp)):
            self.write_qfep_inp(path)

        #Get Qfep executable
        qfep = self.app.q_settings[5][2]
        #Move to path and run Qfep
        self.app.log(' ','Running Qfep in ../%s\n' % '/'.join(path.split('/')[-3:]))
        os.chdir(path)

        tmpfile = open('qfep.out', 'w')
        call('%s <%s' % (qfep, qfep_inp), shell=True, stdout=tmpfile, stderr=tmpfile)

        tmpfile.close()
        os.chdir(self.app.workdir)

    def get_dg(self, prj_title):
        """
        Collects FEP summary for all entries in title
        {title: {entry1:{lambda, dg_f, dg_r, dg}}
        """

        for entry in sorted(self.titles[prj_title]):
            if len(self.titles_dg[prj_title][entry]) == 0:
                self.titles_dg[prj_title][entry] = self.get_part1(self.titles[prj_title][entry])

        self.update_tables()

    def get_part1(self, path, qfep = 'qfep.out'):

        part1 = {'lambda': list(), 'sum_dGf': list(), 'sum_dGr': list(), 'dG': list() }

        if not os.path.isfile('%s/%s' % (path, qfep)):
            self.run_qfep(path, 'qfep.inp')

        with open('%s/%s' % (path, qfep), 'r') as qfep_out:
            found_part1 = False
            for line in qfep_out:
                if found_part1:
                    if len(line.split()) < 6:
                        break
                    if not line.startswith('#'):
                        part1['lambda'].append(float(line.split()[0]))
                        part1['sum_dGf'].append(float(line.split()[2]))
                        part1['sum_dGr'].append(float(line.split()[4]))
                        part1['dG'].append(float(line.split()[5]))
                if '# Part 1' in line:
                    found_part1 = True

        return part1

    def update_tables(self):
        self.fep_listbox.delete(0, END)
        self.fep_listbox.insert(END, 'TITLE      sum(dGf)   +/-  sum(dGr)   +/-     <dG>    +/-')

        for title in self.titles_dg.keys():
            dgf = list()
            dgr = list()
            dg = list()
            for nr in sorted(self.titles_dg[title].keys()):
                dgf.append(self.titles_dg[title][nr]['sum_dGf'][-1])
                dgr.append(self.titles_dg[title][nr]['sum_dGr'][0])
                dg.append(self.titles_dg[title][nr]['dG'][-1])
            self.fep_listbox.insert(END, '%10s %8.2f %6.2f %8.2f %6.2f %8.2f %6.2f' % (title.ljust(10), np.average(dgf),
                                                                                       np.std(dgf) / np.sqrt(len(dgf)),
                                                                                       np.average(dgr),
                                                                                       np.std(dgr) / np.sqrt(len(dgr)),
                                                                                       np.average(dg),
                                                                                       np.std(dg) / np.sqrt(len(dg))))



    def list_titles_event(self, *args):

        self.runs_listbox.delete(0, END)
        selected = map(int, self.titles_listbox.curselection())

        if len(selected) < 1:
            return

        for i in selected:
            title = self.titles_listbox.get(i)
            for entry in sorted(self.titles[title].keys()):
                self.runs_listbox.insert(END, '%s/%s' % (title, entry))


    def list_runs_event(self, *args):
        pass


    def states_changed(self, *args):
        """
        Just trace when numbers of states is changed ==> update linear combination field with correct number of entries.
        """
        lin_comp = ['1 0', '1 0 0', '1 0 0 0']
        try:
            states = int(self.states_entry.get())
            new_lin = lin_comp[states - 2]
        except:
            return

        self.linear_comb.delete(0, END)
        self.linear_comb.insert(0, new_lin)

    def compute_fep(self):
        pass

    def plot_fep(self):
        """
        plot selected parts from FEP list
        """
        dict_translated = {'FEP summary': self.titles_dg}
        term_translated = {'Sum dGf': 'sum_dGf', 'Sum dGr': 'sum_dGr', '<dG>': 'dG'}

        #Get selected part
        part = map(int, self.plot_listbox.curselection())
        if len(part) != 1:
            return

        part = self.plot_listbox.get(part[0]).strip()
        energies = dict_translated[part]

        #Get selected dG to plot
        sel_terms = map(int, self.term_listbox.curselection())
        if len(sel_terms) < 1:
            return

        #Get title(s) to plot
        sel_titles = map(int, self.titles_listbox.curselection())
        if len(sel_titles) < 1:
            print 'No titles selected for plot!'
            return

        #Create dictionary for each title to contain runs to plot
        titles = dict()
        for tit in sel_titles:
            titles[self.titles_listbox.get(tit)] = list()

        #Check if runs are selected, if not, average all from title
        sel_runs = map(int, self.runs_listbox.curselection())
        if len(sel_runs) == 0:
            for tit in titles.keys():
                titles[tit] = sorted(self.titles[tit].keys())
        else:
            for i in map(int, self.runs_listbox.curselection()):
                run = self.runs_listbox.get(i)
                titles[run.split('/')[0]].append(int(run.split('/')[1]))


        ylist = [[]]
        xlist = [[]]
        plot_titles = [[]]


        for tit in titles.keys():
            for t in sel_terms:
                y_dict = dict()
                term = self.term_listbox.get(t)
                plot_titles[0].append('%s %s' % (tit, term))
                term = term_translated[self.term_listbox.get(t).strip()]
                for i in titles[tit]:
                    nr = 0
                    for entry in energies[tit][i][term]:
                        nr += 1
                        if nr not in y_dict:
                            y_dict[nr] = list()
                        y_dict[nr].append(entry)

                tmp_y = list()

                for nr in y_dict.keys():
                    tmp_y.append(np.average(y_dict[nr]))
                ylist[0].append(tmp_y)

                tmp_x = np.arange(len(tmp_y))
                xlist[0].append(tmp_x)

        main_title = ['FEP']

        self.plot_ = Qplot(self, self.root, ylist, xlist, main_title, plot_titles, ' ', r'$\Delta$G (kcal/mol)')
        self.plot_.resizable()
        self.plot_.config(background=self.main_color)



    def show_var_frame(self, *args):
        frames = {'FEP summary': self.fep_frame,
                  'Plot': self.plot_frame}

        for i in frames.keys():
            frames[i].grid_forget()
        try:
            frames[self.selected_frame.get()].grid(row=3, column=0, columnspan=2)
        except:
            pass

    def dialog_window(self):
        self.title('Analyze FEP')
        self.mainframe = Frame(self, bg=self.main_color)
        self.mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame with import button etc
        topframe = Frame(self.mainframe, bg=self.main_color)
        topframe.grid(row=0, column=0, pady=(10,10))

        #plot frame
        self.plot_frame = Frame(self.mainframe, bg=self.main_color)

        #FEP frame
        self.fep_frame = Frame(self.mainframe, bg=self.main_color)


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

        compute_fep = Button(topframe, text='Compute FEP', highlightbackground=self.main_color, command=self.compute_fep)
        compute_fep.grid(row=1, column=0, columnspan=2)

        states_label = Label(topframe, text='States', bg=self.main_color)
        states_label.grid(row=0, column=2)

        self.states_entry = Spinbox(topframe, width=3, highlightthickness=0, relief=GROOVE, from_=2, to=4, increment=1,
                                    textvariable=self.states_var)
        self.states_entry.grid(row=1, column=2)

        kt = Label(topframe, text='    kT     ', bg=self.main_color)
        kt.grid(row=0, column=3)

        bins = Label(topframe, text='   Bins    ', bg=self.main_color)
        bins.grid(row=0, column=4)

        binpoints = Label(topframe, text='Min. pts', bg=self.main_color)
        binpoints.grid(row=0, column=5)

        skip_points = Label(topframe, text='   Skip    ', bg=self.main_color)
        skip_points.grid(row=0, column=6)

        linear_combination = Label(topframe, text='Lin. comb.', bg=self.main_color)
        linear_combination.grid(row=0, column=7)

        self.kT = Entry(topframe, width=5, highlightthickness=0)
        self.kT.grid(row=1, column=3)
        self.kT.delete(0,END)
        self.kT.insert(0, '0.596')

        self.bins = Entry(topframe, width=5, highlightthickness=0)
        self.bins.grid(row=1, column=4)
        self.bins.insert(0, '50')

        self.binpoints_min = Entry(topframe, width=5, highlightthickness=0)
        self.binpoints_min.grid(row=1, column=5)
        self.binpoints_min.insert(0, '30')

        self.skip_points = Entry(topframe, width=5, highlightthickness=0)
        self.skip_points.grid(row=1, column=6)
        self.skip_points.insert(0, '100')

        self.linear_comb = Entry(topframe, width=8, highlightthickness=0)
        self.linear_comb.grid(row=1, column=7)
        self.linear_comb.insert(0, '1 0')

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
        self.titles_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.titles_listbox.bind('<<ListboxSelect>>', self.list_titles_event)

        runs_yscroll = Scrollbar(self.proj_frame)
        runs_yscroll.grid(row = 1, rowspan=10, column = 7, sticky = 'nsw', padx=(0,10))
        runs_xscroll = Scrollbar(self.proj_frame, orient=HORIZONTAL)
        runs_xscroll.grid(row=11, column=4, columnspan=3, sticky='we')

        self.runs_listbox = Listbox(self.proj_frame, yscrollcommand = runs_yscroll.set,
                                      xscrollcommand=runs_xscroll.set,
                                      width=20, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        runs_yscroll.config(command=self.runs_listbox.yview)
        runs_xscroll.config(command=self.runs_listbox.xview)
        self.runs_listbox.grid(row=1, rowspan=10, column = 4, columnspan=3, sticky='e')
        self.runs_listbox.config(font=tkFont.Font(family="Courier", size=12))
        self.runs_listbox.bind('<<ListboxSelect>>', self.list_runs_event)

        #Select frame
        self.view_frame = OptionMenu(sel_frame, self.selected_frame,
                                   'FEP summary', 'Plot')
        self.view_frame.config(highlightbackground=self.main_color, bg=self.main_color, width=30)
        self.view_frame.grid(row=4, column=0)

        #FEP FRAME
        fep_yscroll = Scrollbar(self.fep_frame)
        fep_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.fep_listbox = Listbox(self.fep_frame, yscrollcommand=fep_yscroll.set,
                                   width=60, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        fep_yscroll.config(command=self.fep_listbox.yview)
        self.fep_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.fep_listbox.config(font=tkFont.Font(family="Courier", size=12))

        #PLOT FRAME
        #This can potentially be made more detailed later...
        self.plot_parts = {'FEP summary': ['Sum dGf', 'Sum dGr', '<dG>']}

        plot_yscroll = Scrollbar(self.plot_frame)
        plot_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.plot_listbox = Listbox(self.plot_frame, yscrollcommand=plot_yscroll.set,
                                   width=20, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        plot_yscroll.config(command=self.plot_listbox.yview)
        self.plot_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.plot_listbox.config(font=tkFont.Font(family="Courier", size=12))

        term_yscroll = Scrollbar(self.plot_frame)
        term_yscroll.grid(row=0, rowspan=10, column=7, sticky='nsw', padx=(0,10))
        self.term_listbox = Listbox(self.plot_frame, yscrollcommand=term_yscroll.set,
                                   width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        term_yscroll.config(command=self.plot_listbox.yview)
        self.term_listbox.grid(row=0, rowspan=10, column = 4, columnspan=3, sticky='e')
        self.term_listbox.config(font=tkFont.Font(family="Courier", size=12))

        for part in self.plot_parts.keys():
            self.plot_listbox.insert(END, part)
            for term in self.plot_parts[part]:
                self.term_listbox.insert(END, term)

        self.plot_listbox.selection_set(0)

        plotit = Button(self.plot_frame, text='Plot it!', highlightbackground=self.main_color, command=self.plot_fep)
        plotit.grid(row=0, rowspan=10, column=8)

        #Bottom frame Quit/save
        quit_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        quit_button.grid(row=0, column=0)
