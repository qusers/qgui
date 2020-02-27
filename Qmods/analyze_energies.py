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

from tkinter import Label, Button, Listbox, Scrollbar, EXTENDED, Frame, \
    Toplevel, END, GROOVE, StringVar, OptionMenu, HORIZONTAL
from tkinter.simpledialog import askstring
import tkinter.font
import numpy as np
import os
from select_return import SelectReturn
from Qplot import Qplot

from tkinter.filedialog import askdirectory, asksaveasfilename


class AnalyzeEnergies(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root

        self.main_color = self.app.main_color

        self.selected_frame = StringVar()

        self.selected_frame.set('Show ...')
        self.selected_frame.trace('w', self.show_var_frame)

        #Global dicts
        #{title : {1: path, 2: path...}}
        self.titles = dict()
        self.titles_ene = dict()
        self.titles_files = dict()
        self.unique_logs = list()
        self.dialog_window()

    def add_title(self):
        """
        add a new project title (main name for plot)
        """

        title = askstring('Add title', 'Project title:', parent=self)
        if not title:
            return
        if title in self.titles_listbox.get(0, END):
            self.app.errorBox('Warning', 'Title already exist in Titles.')
            return

        self.titles_listbox.insert(END, title)

        self.titles[title] = dict()
        self.titles_files[title] = dict()
        self.titles_ene[title] = dict()

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
            del self.titles_files[title]
            del self.titles_ene[title]
            #for i in [self.titles_ene]:
            #    if title in i.keys():
            #        del i[title]

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

        #Collect unique MD logfiles
        logfiles = list()
        del self.unique_logs[0:len(self.unique_logs)]

        #Start appending directories to collect energy files from
        dirs = list()
        title = 'Select temperature or run directory'
        path = self.app.workdir
        while True:
            rundir = askdirectory(parent=self, mustexist=True, title=title, initialdir= path)
            if not rundir:
                break
            #Check if the first directory is a temperature/Master directory
            if len(dirs) == 0:
                try:
                    if os.path.isdir('%s/1' % rundir):
                        self.app.log(' ', '%s Added. Collecting logfiles from subdirectories...\n\n'
                                          % '/'.join(rundir.split('/')[-3:]))

                        temp_selected = True
                        break
                except:
                    pass

            if not temp_selected:
                path = '%s' % '/'.join(rundir.split('/')[:-1])

                #Check if en files exist in dir, and if size above 0:
                enefiles = False
                for ene in os.listdir(rundir):
                    if ene.endswith('.log'):
                        if 'md' in ene:
                            if os.path.getsize('%s/%s' % (rundir, ene)) != 0:
                                enefiles = True
                                nr_path = int(rundir.split('/')[-1])
                                if nr_path not in list(self.titles_files[prj_title].keys()):
                                    self.titles_files[prj_title][nr_path] = list()
                                if ene not in self.titles_files[prj_title][nr_path]:
                                    self.titles_files[prj_title][nr_path].append(ene)
                                if ene not in logfiles:
                                    logfiles.append(ene)

                if enefiles:
                    dirs.append(rundir)
                    title = '%s Added. Select next dir or cancel' % dirs[-1].split('/')[-1]
                    self.app.log(' ', '%s Added. Select next dir or cancel to finish\n'
                                      % '/'.join(dirs[-1].split('/')[-3:]))
                    try:
                        self.titles[prj_title][int(rundir.split('/')[-1])] = rundir
                        self.titles_ene[prj_title][int(rundir.split('/')[-1])] = dict()
                    except:
                        self.titles[prj_title][rundir.split('/')[-1]] = rundir
                        self.titles_ene[prj_title][rundir.split('/')[-1]] = dict()
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
                        if ene.endswith('.log'):
                            if 'md' in ene:
                                if os.path.getsize('%s/%s' % (subdir, ene)) != 0:
                                    enefiles = True

                                    if nr_dir not in list(self.titles_files[prj_title].keys()):
                                        self.titles_files[prj_title][nr_dir] = list()
                                    if ene not in self.titles_files[prj_title][nr_dir]:
                                        self.titles_files[prj_title][nr_dir].append(ene)
                                    if ene not in logfiles:
                                        logfiles.append(ene)

                    if enefiles:
                        dirs.append(subdir)
                        self.titles[prj_title][nr_dir] = subdir
                        self.titles_ene[prj_title][nr_dir] = dict()
                        self.app.log(' ','../%s added\n' % '/'.join(subdir.split('/')[-2:]))
                        print(subdir)
                    nr_dir += 1
                else:
                    break

        #Select log files, and extract energies:
        self.app.log('info', 'Found %d unique MD files. Select the ones to analyze.' % len(logfiles))
        self.select_name = SelectReturn(self, self.root, logfiles, '%s MD logfiles' % prj_title, None)
        self.select_name.configure(background=self.main_color)
        self.select_name.resizable()

        self.list_titles_event()

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

        for selected in reversed(selections):
            search = '/'.join(self.runs_listbox.get(selected).split('/')[-3:])
            self.runs_listbox.delete(selected)

            for nr in list(self.titles[title].keys()):
                if '/'.join(self.titles[title][nr].split('/')[-3:]) == search:
                    for i in [self.titles, self.titles_ene]:
                        if title in list(i.keys()):
                            if nr in list(i[title].keys()):
                                del i[title][nr]


        #self.update_tables()
        self.list_titles_event()

    def get_energies(self, title):
        """
        Extract energies for project title and logfiles specified self.unique_logs
        """
        #Delete logfiles that is not selected
        for entry in sorted(self.titles_files[title].keys()):
            del_log = list()
            for log in sorted(self.titles_files[title][entry]):
                print(log)
                if log not in self.unique_logs:
                    del_log.append(log)
            for log_del in del_log:
                del self.titles_files[title][entry][self.titles_files[title][entry].index(log_del)]

        self.list_titles_event()

        self.app.log(' ', '\nCollecting energies, please wait!\n')
        self.update()

        for entry in sorted(self.titles_files[title].keys()):
            if entry not in list(self.titles_ene[title].keys()):
                self.titles_ene[title][entry] = dict()
            for log in sorted(self.titles_files[title][entry]):
                self.titles_ene[title][entry][log] = dict()
                self.titles_ene[title][entry][log]['solute'] = {'el': list(),
                                                                'vdW': list(),
                                                                'bond': list(),
                                                                'angle': list(),
                                                                'torsion': list(),
                                                                'improper': list(),
                                                                'total': list()}
                self.titles_ene[title][entry][log]['solvent'] = {'el': list(),
                                                                'vdW': list(),
                                                                'bond': list(),
                                                                'angle': list(),
                                                                'torsion': list(),
                                                                'improper': list(),
                                                                'total': list()}
                self.titles_ene[title][entry][log]['solute-solvent'] = {'el': list(),
                                                                        'vdW': list(),
                                                                        'total': list()}

                self.titles_ene[title][entry][log]['LRF'] = {'el': list()}

                self.titles_ene[title][entry][log]['Q-atom'] = {'el': list(),
                                                                'vdW': list(),
                                                                'bond': list(),
                                                                'angle': list(),
                                                                'torsion': list(),
                                                                'improper': list(),
                                                                'total': list()}

                self.titles_ene[title][entry][log]['restraints'] = {'total': list(),
                                                                    'exc. atoms': list(),
                                                                    'solv. rad.': list(),
                                                                    'solv. pol.': list(),
                                                                    'shell atoms': list(),
                                                                    'solute': list()}

                self.titles_ene[title][entry][log]['SUM'] = {'total': list(),
                                                             'potential': list(),
                                                             'kinetic': list()}

                self.titles_ene[title][entry][log]['temperature'] = {'step': list(),
                                                                     'T (K)': list()}

                logfile = '%s/%s' % (self.titles[title][entry], log)
                if os.path.isfile(logfile):
                    with open(logfile, 'r') as mdlog:
                        print('Collecting energies from')
                        print(('   %s' % logfile))
                        for line in mdlog:
                            if 'Stepsize' in line:
                                self.titles_ene[title][entry][log]['stepsize'] = float(line.split('=')[-1])
                            elif 'Energy summary print-out interval' in line:
                                self.titles_ene[title][entry][log]['interval'] = float(line.split('=')[-1])
                            elif 'solute        ' in line:
                                el, vdw, bnd, ang, tor, imp = list(map(float, line.split()[1:7]))
                                total = (el + vdw + bnd + ang + tor + imp)
                                self.titles_ene[title][entry][log]['solute']['el'].append(el)
                                self.titles_ene[title][entry][log]['solute']['vdW'].append(vdw)
                                self.titles_ene[title][entry][log]['solute']['bond'].append(bnd)
                                self.titles_ene[title][entry][log]['solute']['angle'].append(ang)
                                self.titles_ene[title][entry][log]['solute']['torsion'].append(tor)
                                self.titles_ene[title][entry][log]['solute']['improper'].append(imp)
                                self.titles_ene[title][entry][log]['solute']['total'].append(total)

                            elif line.startswith('solvent       '):
                                el, vdw, bnd, ang, tor, imp = list(map(float, line.split()[1:7]))
                                total = (el + vdw + ang + tor + imp)
                                self.titles_ene[title][entry][log]['solvent']['el'].append(el)
                                self.titles_ene[title][entry][log]['solvent']['vdW'].append(vdw)
                                self.titles_ene[title][entry][log]['solvent']['bond'].append(bnd)
                                self.titles_ene[title][entry][log]['solvent']['angle'].append(ang)
                                self.titles_ene[title][entry][log]['solvent']['torsion'].append(tor)
                                self.titles_ene[title][entry][log]['solvent']['improper'].append(imp)
                                self.titles_ene[title][entry][log]['solvent']['total'].append(total)

                            elif 'solute-solvent' in line:
                                el, vdw = list(map(float, line.split()[1:3]))
                                total = (el + vdw)
                                self.titles_ene[title][entry][log]['solute-solvent']['el'].append(el)
                                self.titles_ene[title][entry][log]['solute-solvent']['vdW'].append(vdw)
                                self.titles_ene[title][entry][log]['solute-solvent']['total'].append(total)

                            elif 'LRF           ' in line:
                                el = float(line.split()[1])
                                self.titles_ene[title][entry][log]['LRF']['el'].append(el)

                            elif 'Q-atom        ' in line:
                                el, vdw, bnd, ang, tor, imp = list(map(float, line.split()[1:7]))
                                total = (el + vdw + bnd + ang + tor + imp)
                                self.titles_ene[title][entry][log]['Q-atom']['el'].append(el)
                                self.titles_ene[title][entry][log]['Q-atom']['vdW'].append(vdw)
                                self.titles_ene[title][entry][log]['Q-atom']['bond'].append(bnd)
                                self.titles_ene[title][entry][log]['Q-atom']['angle'].append(ang)
                                self.titles_ene[title][entry][log]['Q-atom']['torsion'].append(tor)
                                self.titles_ene[title][entry][log]['Q-atom']['improper'].append(imp)
                                self.titles_ene[title][entry][log]['Q-atom']['total'].append(total)

                            elif line.startswith('restraints    '):
                                total, fix, slvnt_rad, slvnt_pol, shell, solute = list(map(float, line.split()[1:7]))
                                self.titles_ene[title][entry][log]['restraints']['total'].append(total)
                                self.titles_ene[title][entry][log]['restraints']['exc. atoms'].append(fix)
                                self.titles_ene[title][entry][log]['restraints']['solv. rad.'].append(slvnt_rad)
                                self.titles_ene[title][entry][log]['restraints']['solv. pol.'].append(slvnt_pol)
                                self.titles_ene[title][entry][log]['restraints']['shell atoms'].append(shell)
                                self.titles_ene[title][entry][log]['restraints']['solute'].append(solute)

                            elif 'SUM           ' in line and len(line.split()) > 3:
                                total, potential, kinetic = list(map(float, line.split()[1:4]))
                                self.titles_ene[title][entry][log]['SUM']['total'].append(total)
                                self.titles_ene[title][entry][log]['SUM']['potential'].append(potential)
                                self.titles_ene[title][entry][log]['SUM']['kinetic'].append(kinetic)

                            elif 'Temperature at step' in line and 'T_free' in line:
                                step = int(line.split(':')[0].split()[-1])
                                t_free = float(line.split()[-1])
                                self.titles_ene[title][entry][log]['temperature']['step'].append(step)
                                self.titles_ene[title][entry][log]['temperature']['T (K)'].append(t_free)

        self.app.log(' ', '\nEnergy collection completed!\n')

        self.update_tables()

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

    def update_tables(self):
        tables = [self.fep_listbox, self.fep2_listbox, self.fep3_listbox, self.restraints_listbox, self.lrf_listbox,
                  self.sum_listbox]

        terms = ['solute', 'solute-solvent', 'solvent', 'restraints', 'LRF', 'SUM']

        heads = [['Title', '  total', '     el', '    vdW', '    bond', '   angle', '  torsion', ' improper' ],
                 ['Title', '  total', '     el', '    vdW'],
                 ['Title', '  total', '     el', '    vdW', '    bond', '   angle', '  torsion', ' improper' ],
                 ['Title', '  total', 'exc. atoms', 'solv. rad.', 'solv. pol.', 'shell atoms', '  solute'],
                 ['Title', '     el'],
                 ['Title', 'potential', ' kinetic']]

        for i in range(len(tables)):
            table = tables[i]
            table.delete(0, END)
            head = ' '.join(['%10s' % x.ljust(10) for x in heads[i]])
            table.insert(END, head)

            term = terms[i]

            for title in list(self.titles_ene.keys()):
                ene = list()
                se = list()

                for j in range(1, len(heads[i])):
                    tmp_ene = list()
                    ene_type = heads[i][j].strip()

                    for entry in list(self.titles_ene[title].keys()):
                        for log in self.titles_ene[title][entry]:
                            tmp_ene += self.titles_ene[title][entry][log][term][ene_type]
                    ene.append(round(np.average(tmp_ene), 2))
                    se.append(round(float(np.std(tmp_ene)), 2))
                ene_str = ' '.join(['%10.2f' % x for x in ene])
                se_str = ' '.join(['%10.2f' % x for x in se])

                tit = '%10s' % title.ljust(10)
                se_tit = 'stdev     '

                ene_str = tit + ene_str
                se_str = se_tit + se_str

                table.insert(END, ene_str)
                table.insert(END, se_str)
                table.insert(END, ' ')

    def plot_it(self):
        """
        Open Qplot and plots selected titles and corresponding energy terms
        """
        #{title: {entry: [log]}}

        tit_entry = dict()

        titles_sel = list(map(int, self.titles_listbox.curselection()))
        if len(titles_sel) < 1:
            return

        runs_sel = list(map(int, self.runs_listbox.curselection()))

        #Collect titles, runs and specific md files to plot
        if len(runs_sel) > 0:
            for i in runs_sel:
                run = self.runs_listbox.get(i)
                title = run.split('/')[0]
                entry = int(run.split('/')[1])
                logfile = run.split('/')[-1]
                if title not in list(tit_entry.keys()):
                    tit_entry[title] = dict()
                if entry not in list(tit_entry[title].keys()):
                    tit_entry[title][entry] = list()

                tit_entry[title][entry].append(logfile)

        else:
            for i in titles_sel:
                title = self.titles_listbox.get(i)
                tit_entry[title] = self.titles_files[title]

        print(tit_entry)


        #Collect terms to plot (terms becomes the subplot titles):
        main_title = list()
        plot_titles = list()
        ylist = list()
        xlist = list()
        sel_terms = list(map(int, self.term_listbox.curselection()))
        if len(sel_terms) < 1:
            print('Please select terms to plot (el, vdw...)')
            return
        for i in sel_terms:
            main_title.append(self.term_listbox.get(i))
            plot_titles.append([])
            xlist.append([])
            ylist.append([])

        print(main_title)

        plot_types = list()
        #Collect types to plot
        sel_types = list(map(int, self.plot_listbox.curselection()))
        if len(sel_types) < 1:
            return
        for i in sel_types:
            plot_types.append(self.plot_listbox.get(i))

        #Collect energies and append to xlist, ylist and plot_titles
        for i_term in range(len(main_title)):
            subtitle = main_title[i_term]

            for title in list(tit_entry.keys()):
                for plot_type in plot_types:
                    line_title = '%s %s' % (title, plot_type)
                    plot_titles[i_term].append(line_title)
                    tmp_y = list()
                    tmp_x = list()
                    time = 0.0

                    for entry in sorted(tit_entry[title].keys()):
                        for log in sorted(tit_entry[title][entry]):
                            interval = float(self.titles_ene[title][entry][log]['interval'])
                            stepsize = float(self.titles_ene[title][entry][log]['stepsize'])
                            for i in self.titles_ene[title][entry][log][plot_type][subtitle]:
                                tmp_y.append(float(i))
                                tmp_x.append(time)
                                time += (interval * stepsize)/ 1000000.00
                    ylist[i_term].append(tmp_y)
                    xlist[i_term].append(tmp_x)


        #PLOT IT
        self.plot_ = Qplot(self, self.root, ylist, xlist, main_title, plot_titles, 'Time (ns)', 'E (kcal/mol)')
        self.plot_.resizable()
        self.plot_.config(background=self.main_color)

    def list_titles_event(self, *args):

        self.runs_listbox.delete(0, END)
        selected = list(map(int, self.titles_listbox.curselection()))

        if len(selected) < 1:
            return

        for i in selected:
            title = self.titles_listbox.get(i)
            for entry in sorted(self.titles[title].keys()):
                for log in sorted(self.titles_files[title][entry]):
                    self.runs_listbox.insert(END, '%s/%s/%s' % (title, entry,log))


    def list_runs_event(self, *args):
        pass

    def list_plot_event(self, *args):
        self.term_listbox.delete(0, END)

        selections = list(map(int, self.plot_listbox.curselection()))

        if len(selections) == 0:
            return
        terms = list()
        for i in selections:
            plot = self.plot_listbox.get(i).strip()
            for term in self.plot_parts[plot]:
                if term not in terms:
                    terms.append(term)
                    self.term_listbox.insert(END, term)

    def show_var_frame(self, *args):
        #'Solute/solvent', 'Restraints', 'LRF', 'SUM', 'Plot'
        frames = {'Solute': self.ss_frame,
                  'Solute-solvent': self.ss2_frame,
                  'Solvent': self.ss3_frame,
                  'Restraints': self.restraints_frame,
                  'LRF': self.lrf_frame,
                  'SUM': self.sum_frame,
                  'Plot': self.plot_frame}

        for i in list(frames.keys()):
            frames[i].grid_forget()
        try:
            frames[self.selected_frame.get()].grid(row=3, column=0, columnspan=2)
        except:
            pass

    def dialog_window(self):
        self.title('Analyze trajectory energies')
        self.mainframe = Frame(self, bg=self.main_color)
        self.mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame with import button etc
        topframe = Frame(self.mainframe, bg=self.main_color)
        topframe.grid(row=0, column=0, pady=(10,10))

        #plot frame
        self.plot_frame = Frame(self.mainframe, bg=self.main_color)

        #solute frame
        self.ss_frame = Frame(self.mainframe, bg=self.main_color)

        #solute-solvent frame
        self.ss2_frame = Frame(self.mainframe, bg=self.main_color)

        #solvent frame
        self.ss3_frame = Frame(self.mainframe, bg=self.main_color)

        #Restraints frame
        self.restraints_frame = Frame(self.mainframe, bg=self.main_color)

        #Locar reaction field (LRF) frame
        self.lrf_frame = Frame(self.mainframe, bg=self.main_color)

        #SUM/ Total kinetic/potential energies
        self.sum_frame = Frame(self.mainframe, bg=self.main_color)

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
        runs_yscroll.grid(row = 1, rowspan=10, column = 11, sticky = 'nsw', padx=(0,10))
        runs_xscroll = Scrollbar(self.proj_frame, orient=HORIZONTAL)
        runs_xscroll.grid(row=11, column=4, columnspan=6, sticky='we')

        self.runs_listbox = Listbox(self.proj_frame, yscrollcommand = runs_yscroll.set,
                                      xscrollcommand=runs_xscroll.set,
                                      width=40, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                      exportselection=False)
        runs_yscroll.config(command=self.runs_listbox.yview)
        runs_xscroll.config(command=self.runs_listbox.xview)
        self.runs_listbox.grid(row=1, rowspan=10, column = 4, columnspan=6, sticky='e')
        self.runs_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        self.runs_listbox.bind('<<ListboxSelect>>', self.list_runs_event)

        #Select frame
        self.view_frame = OptionMenu(sel_frame, self.selected_frame,
                                   'Solute', 'Solute-solvent', 'Solvent', 'Restraints', 'LRF', 'SUM', 'Plot')
        self.view_frame.config(highlightbackground=self.main_color, bg=self.main_color, width=30)
        self.view_frame.grid(row=4, column=0)

        #Solute energies
        fep_yscroll = Scrollbar(self.ss_frame)
        fep_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.fep_listbox = Listbox(self.ss_frame, yscrollcommand=fep_yscroll.set,
                                   width=90, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        fep_yscroll.config(command=self.fep_listbox.yview)
        self.fep_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.fep_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        export_solute = Button(self.ss_frame, text='Export table', highlightbackground=self.main_color,
                               command=lambda: self.export_table(self.fep_listbox))
        export_solute.grid(row=11, column=0, columnspan=3)

        #Soute-solvent energies
        fep2_yscroll = Scrollbar(self.ss2_frame)
        fep2_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.fep2_listbox = Listbox(self.ss2_frame, yscrollcommand=fep2_yscroll.set,
                                   width=50, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        fep2_yscroll.config(command=self.fep2_listbox.yview)
        self.fep2_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.fep2_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        export_solutesolvent = Button(self.ss2_frame, text='Export table', highlightbackground=self.main_color,
                               command=lambda: self.export_table(self.fep2_listbox))
        export_solutesolvent.grid(row=11, column=0, columnspan=3)

        #solvent energies
        fep3_yscroll = Scrollbar(self.ss3_frame)
        fep3_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.fep3_listbox = Listbox(self.ss3_frame, yscrollcommand=fep3_yscroll.set,
                                   width=90, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        fep3_yscroll.config(command=self.fep3_listbox.yview)
        self.fep3_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.fep3_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        export_solvent = Button(self.ss3_frame, text='Export table', highlightbackground=self.main_color,
                               command=lambda: self.export_table(self.fep3_listbox))
        export_solvent.grid(row=11, column=0, columnspan=3)

        #Restraints energies
        restraints_yscroll = Scrollbar(self.restraints_frame)
        restraints_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.restraints_listbox = Listbox(self.restraints_frame, yscrollcommand=restraints_yscroll.set,
                                   width=80, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        restraints_yscroll.config(command=self.restraints_listbox.yview)
        self.restraints_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.restraints_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        export_restraints = Button(self.restraints_frame, text='Export table', highlightbackground=self.main_color,
                               command=lambda: self.export_table(self.restraints_listbox))
        export_restraints.grid(row=11, column=0, columnspan=3)

        #LRF energies
        lrf_yscroll = Scrollbar(self.lrf_frame)
        lrf_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.lrf_listbox = Listbox(self.lrf_frame, yscrollcommand=lrf_yscroll.set,
                                   width=25, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        lrf_yscroll.config(command=self.lrf_listbox.yview)
        self.lrf_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.lrf_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        export_lrf = Button(self.lrf_frame, text='Export table', highlightbackground=self.main_color,
                               command=lambda: self.export_table(self.lrf_listbox))
        export_lrf.grid(row=11, column=0, columnspan=3)

        #SUM potential / kinetic energies
        sum_yscroll = Scrollbar(self.sum_frame)
        sum_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.sum_listbox = Listbox(self.sum_frame, yscrollcommand=sum_yscroll.set,
                                   width=40, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        sum_yscroll.config(command=self.sum_listbox.yview)
        self.sum_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.sum_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        export_sum = Button(self.sum_frame, text='Export table', highlightbackground=self.main_color,
                               command=lambda: self.export_table(self.sum_listbox))
        export_sum.grid(row=11, column=0, columnspan=3)

        #PLOT FRAME
        #This can potentially be made more detailed later...
        self.plot_parts = {'solute': ['el', 'vdW', 'bond', 'angle', 'torsion', 'improper', 'total'],
                           'solute-solvent': ['el', 'vdW', 'total'],
                           'solvent': ['el', 'vdW', 'total'],
                           'Q-atom': ['el', 'vdW', 'bond', 'angle', 'torsion', 'improper', 'total'],
                           'LRF': ['el'],
                           'restraints': ['exc. atoms', 'solv. rad.', 'solv. pol.',
                                          'shell atoms', 'solute', 'total'],
                           'SUM': ['potential', 'kinetic', 'total'],
                           'temperature': ['T (K)']}

        plot_yscroll = Scrollbar(self.plot_frame)
        plot_yscroll.grid(row=0, rowspan=10, column=3, sticky='nsw', padx=(0,10))
        self.plot_listbox = Listbox(self.plot_frame, yscrollcommand=plot_yscroll.set,
                                   width=20, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        plot_yscroll.config(command=self.plot_listbox.yview)
        self.plot_listbox.grid(row=0, rowspan=10, column = 0, columnspan=3, sticky='e')
        self.plot_listbox.config(font=tkinter.font.Font(family="Courier", size=12))
        self.plot_listbox.bind('<<ListboxSelect>>', self.list_plot_event)

        term_yscroll = Scrollbar(self.plot_frame)
        term_yscroll.grid(row=0, rowspan=10, column=7, sticky='nsw', padx=(0,10))
        self.term_listbox = Listbox(self.plot_frame, yscrollcommand=term_yscroll.set,
                                   width=15, height=10, highlightthickness=0, relief=GROOVE, selectmode=EXTENDED,
                                   exportselection=False)
        term_yscroll.config(command=self.term_listbox.yview)
        self.term_listbox.grid(row=0, rowspan=10, column = 4, columnspan=3, sticky='e')
        self.term_listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        for part in sorted(self.plot_parts.keys()):
            self.plot_listbox.insert(END, part)
        for term in self.plot_parts[self.plot_listbox.get(0)]:
            self.term_listbox.insert(END, term)

        self.plot_listbox.selection_set(0)

        plotit = Button(self.plot_frame, text='Plot it!', highlightbackground=self.main_color, command=self.plot_it)
        plotit.grid(row=0, rowspan=10, column=8)

        #Bottom frame Quit/save
        quit_button = Button(bottomframe, text='Close', highlightbackground=self.main_color, command=self.destroy)
        quit_button.grid(row=0, column=0)
