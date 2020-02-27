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

from tkinter import Label, TOP, Button, Listbox, Scrollbar, BROWSE, Spinbox, Entry, Text, Frame, \
    Toplevel, DISABLED, END, GROOVE, NORMAL, BOTH, IntVar, StringVar, Checkbutton
import tkinter.font
import numpy as np
import os
import shutil
import matplotlib
matplotlib.use('TkAgg')
#Implement default mpl key bindings
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from tkinter.filedialog import askopenfilenames, asksaveasfilename, askdirectory
matplotlib.rcParams['text.usetex'] = True



class EvbCalibration(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root

        #List to store all energy files:
        self.ene_files = []

        self.autoupdate = IntVar()
        self.show_target_lines = IntVar()

        #Trace changes in parameters to update table/plot
        self.alpha_var = StringVar()
        self.alpha_var.trace('w', self.parameters_changed)
        self.hij_var = StringVar()
        self.hij_var.trace('w', self.parameters_changed)
        self.target_dG_var = StringVar()
        self.target_dG_var.trace('w', self.target_changed)
        self.target_dG0_var = StringVar()
        self.target_dG0_var.trace('w', self.target_changed)


        self.dg_plot = None
        self.dG = []
        self.enegaps = []

        self.dialog_window()
        self.canvas.mpl_connect('key_press_event', self.on_key_event)

        #Insert default EVB parameter values:
        self.alpha_entry.delete(0, END)
        self.alpha_entry.insert(0, '0')
        self.hij_entry.delete(0,END)
        self.hij_entry.insert(0, '0')
        self.target_dG.delete(0, END)
        self.target_dG.insert(0, '0.0')
        self.target_dG0.delete(0, END)
        self.target_dG0.insert(0, '0.0')

        #Insert default FEP settings:
        self.temperature.insert(0, 300)
        self.bins.insert(0, 50)
        self.binpoints_min.insert(0, 30)
        self.skip_points.insert(0, 100)


    def on_key_event(self, event):
        print(('you pressed %s' % event.key))
        key_press_handler(event, self.canvas, self.toolbar)

    def import_all_from_dir(self):
        """
        Ask user for directory and finds all .en files.
        """
        dirname = askdirectory(parent=self, mustexist=False,
                               title='Select EVB reference directory', initialdir= self.app.workdir)

        if dirname:
            found_files = False
            ene_files = self.find_files(['.en'], dirname)
            if len(ene_files) > 0:
                found_files = True
                for ene_file in ene_files:
                    self.ene_files.append(ene_file)

            else:
                #Look for subfolders
                print('Seraching subdirectories')
                nr_dir = 1
                while True:
                    subdir = '%s/%d' % (dirname, nr_dir)

                    if os.path.isdir(subdir):
                        ene_files = self.find_files(['.en'], subdir)
                        if len(ene_files) > 0:
                            found_files = True
                            for ene_file in ene_files:
                                self.ene_files.append(ene_file)

                    else:
                        break

                    nr_dir += 1

            if not found_files:
                self.app.errorBox('Error', 'Could not find any energy files in directory,')
            else:
                self.update_table()
        else:
            return

        #self.update_plot()

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
                        if os.path.getsize('%s/%s' % (qdir, files)) > 0:
                            list_of_files.append('%s/%s' % (qdir, files))

        return list_of_files

    def make_qfep_input(self, enefiles=[], out_dir='workdir'):
        if len(enefiles) < 1:
            print('No energy files specified')
            return
        if out_dir == 'workdir':
            out_dir = self.app.workdir

        nr_files = len(enefiles)
        kT = 0.001987209 * float(self.temperature.get())
        skip = int(self.skip_points.get())
        bins = int(self.bins.get())
        #Too many bins will cause Qfep to fail
        if bins > 100:
            self.app.errorBox('Warning', 'Too many bins (%d). Setting bins to 100.' % bins)
            bins = 100
            self.bins.delete(0, END)
            self.bins.insert(0, bins)
        min_points = int(self.binpoints_min.get())
        alpha = float(self.alpha_entry.get())
        hij = float(self.hij_entry.get())
        linear_comb = self.linear_comb.get()

        qfepinp = open('%s/qfep.inp' % out_dir, 'w')

        qfepinp.write('%d\n' % nr_files)
        qfepinp.write('2  0\n')
        qfepinp.write('%.4f  %d\n' % (kT, skip))
        qfepinp.write('%d\n' % bins)
        qfepinp.write('%d\n' % min_points)
        qfepinp.write('%.2f\n' % alpha)
        qfepinp.write('1\n')
        qfepinp.write('1  2  %.2f  0.0 0 0.000\n' % hij)
        qfepinp.write('%s\n' % linear_comb)

        for enefile in enefiles:
            qfepinp.write('%s\n' % enefile)

        qfepinp.write('stop\n')
        qfepinp.close()

    def compute_evb(self, return_stats = False):
        if len(self.ene_files) < 1:
            print('No energy files specified')
            return
        #Energy files can not be in different directories...
        #Make tempdir:
        tempdir = '%s/%s' % (self.app.workdir, '.qfep_tempdir')
        if os.path.isdir(tempdir):
            shutil.rmtree(tempdir)
        os.mkdir(tempdir)
        enefiles = []

        #Copy all enenrgy files to temp. directory and rename files.
        for i in range(len(self.ene_files)):
            shutil.copy2(self.ene_files[i], '%s/%04d.en' % (tempdir, i))
            enefiles.append('%04d.en' % i)
        try:
            self.make_qfep_input(enefiles, tempdir)
            self.make_qfep_input(enefiles)
        except:
            return

        #Get qfep command from settings
        qfep = self.app.q_settings[ 'executables' ][2]

        #Change directory to tempdir:
        os.chdir(tempdir)

        #POPEN Qfep
        #tmpfile = open(tempdir + '/.tmpfile', 'wb')
        #self.session = Popen([qfep, '<','qfep.inp','>','qfep.out'], stdout=tmpfile, stdin=PIPE, preexec_fn=os.setsid)
        os.system('%s <qfep.inp>qfep.out' % qfep)

        #Copy qfep.out back to workdir:
        if os.path.isfile('%s/qfep.out' % self.app.workdir):
            os.remove('%s/qfep.out' % self.app.workdir)
        shutil.copy2('%s/qfep.out' % tempdir, '%s/qfep.out' % self.app.workdir)

        os.chdir(self.app.workdir)

        #Remove tempdir when done:
        shutil.rmtree(tempdir)

        self.dG = []
        self.enegaps = []

        #Get enegap and dG
        part3 = False
        count = 0
        with open('%s/qfep.out' % self.app.workdir, 'r') as qfepout:
            for line in qfepout:
                if '# Part 3: Bin-averaged summary' in line:
                    self.app.log(' ', '\n')
                    part3 = True
                if part3:
                    if line == '':
                        part3 = False
                        break
                    if count > 0:
                        self.app.log(' ',line)
                    if count > 1:
                        try:
                            self.enegaps.append(float(line.split()[1]))
                            self.dG.append(float(line.split()[3]))
                        except:
                            continue
                    else:
                        count += 1

        rsE, tsE, rE = self.find_local_min_max(self.dG)
        self.dG = [x - float(rsE) for x in self.dG]
        self.update_plot()

        self.update_ddG(tsE, rE)

        if return_stats:
            return tsE, rE

    def update_ddG(self, tsE, rE):
        """
        updates the computed values
        """

        #Get ddG values:
        try:
            target_dG = float(self.target_dG.get())
            target_dG0 = float(self.target_dG0.get())
            ddG = float(tsE) - target_dG
            ddG0 = float(rE) - target_dG0
        except:
            ddG = 'na'
            ddG0 = 'na'

        self.dG_real_txt.config(state=NORMAL)
        self.dG0_real_txt.config(state=NORMAL)
        self.ddG_real_txt.config(state=NORMAL)
        self.ddG0_real_txt.config(state=NORMAL)

        self.dG_real_txt.delete(0.0, END)
        self.dG0_real_txt.delete(0.0, END)
        self.ddG_real_txt.delete(0.0, END)
        self.ddG0_real_txt.delete(0.0, END)

        self.dG_real_txt.insert(0.0, tsE)
        self.dG0_real_txt.insert(0.0, rE)
        self.ddG_real_txt.insert(0.0, ddG)
        self.ddG0_real_txt.insert(0.0, ddG0)

        self.dG_real_txt.config(state=DISABLED)
        self.dG0_real_txt.config(state=DISABLED)
        self.ddG_real_txt.config(state=DISABLED)
        self.ddG0_real_txt.config(state=DISABLED)

    def find_local_min_max(self, dG):
        """
        Returns RS, TS and reaction energy
        """
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

    def linear_reg(self, x, y):
        """
        Performs linear regression for x and y (y = ax +b)
        returns a and b.
        """
        x = list(map(float, x))
        y = list(map(float, y))

        A = np.array([x, np.ones(len(x))])
        w = np.linalg.lstsq(A.T, y)[0]

        return w[0], w[1]

    def scan_evb(self):
        """
        Find optimal alpha and Hij to target dG values
        Iterates between dG_rxn(alpha) and dG_act(Hij)
        """
        self.autoupdate.set(0)
        if len(self.ene_files) < 1:
            return

        #Convergence treshold:
        treshold = 0.20

        #get target values:
        dG_target = float(self.target_dG.get())
        dG0_target = float(self.target_dG0.get())

        if int(dG_target) == 0:
            if int(dG0_target) == 0:
                self.app.errorBox('Error', 'Please specify target values before starting scan.')
                return

        #Find TS and RS energy:
        alpha = float(self.alpha_entry.get())
        hij = float(self.hij_entry.get())
        current_tsE, current_rE = self.compute_evb(True)

        a_step = 10
        h_step = 10
        counter = 0
        done = False
        while not done:
            self.update()
            counter += 1
            print(('RRC iteration %d' % counter))

            #dG_rxn(alpha):
            #Decide step sign:
            alpha += a_step
            self.insert_entry(self.alpha_entry, alpha)
            try:
                new_rE = float(self.compute_evb(True)[1])
            except:
                try:
                    a_step = -(a_step)
                    alpha += (2 * a_step)
                    self.insert_entry(self.alpha_entry, alpha)
                    new_rE = float(self.compute_evb(True)[1])
                except:
                    self.app.errorBox('Error', 'Could not fit alpha. Try different starting point.')
                    return

            if abs(new_rE - dG0_target) >= abs(current_rE - dG0_target):
                a_step = -a_step
                alpha += (2 * a_step)
                self.insert_entry(self.alpha_entry, alpha)
                new_rE = float(self.compute_evb(True)[1])

            current_rE = new_rE

            x = list()
            x.append(float(alpha))
            y = list()
            y.append(float(current_rE))
            for i in range(1):
                alpha += a_step
                self.insert_entry(self.alpha_entry, alpha)
                current_rE = self.compute_evb(True)[1]
                x.append(float(alpha))
                y.append(float(current_rE))
            a, b = self.linear_reg(x, y)

            alpha = (dG0_target - b) / a

            if b < 0:
                b = '- %.2f' % b
            else:
                b = '+ %.2f' % b
            self.app.log(' ', '\n' + 70 * '#' + '\n')
            self.app.log(' ', "\N{GREEK CAPITAL LETTER DELTA}G_rxn(\N{GREEK SMALL LETTER ALPHA}, Hij=%.2f) = %.2f"
                              "\N{GREEK SMALL LETTER ALPHA} %s" % (hij, a, b) + '\n')
            self.app.log(' ', '\n' + 70 * '#' + '\n')
            self.insert_entry(self.alpha_entry, alpha)

            #dG_act(Hij)
            #decide step sign:
            hij += h_step
            hij = abs(hij)
            self.insert_entry(self.hij_entry, hij)
            try:
                new_tsE, new_rE = float(self.compute_evb(True)[0:2])
            except:
                try:
                    h_step = -h_step
                    hij += (2 * h_step)
                    self.insert_entry(self.hij_entry, hij)
                    new_tsE = self.compute_evb(True)[0]
                except:
                    self.app.errorBox('Error', 'Could not fit Hij. Try different starting point.')
                    return

            if abs(new_tsE - dG_target) > abs(current_tsE - dG_target):
                h_step = -h_step
                hij += (2 * h_step)
                self.insert_entry(self.hij_entry, hij)
                new_tsE = self.compute_evb(True)[0]
            current_tsE = new_tsE

            x = list()
            y = list()
            x.append(hij)
            y.append(current_tsE)

            for i in range(1):
                hij += h_step
                self.insert_entry(self.hij_entry, hij)
                current_tsE = self.compute_evb(True)[0]
                x.append(hij)
                y.append(current_tsE)

            a, b = self.linear_reg(x, y)
            hij = (dG_target - b)/a
            if b < 0:
                b = '- %.2f' % b
            else:
                b = '+ %.2f' % b
            self.app.log(' ', '\n' + 70 * '#' + '\n')
            self.app.log(' ', "\N{GREEK CAPITAL LETTER DELTA}G_act(Hij, \N{GREEK SMALL LETTER ALPHA}=%.2f) = %.2f"
                              "Hij %s" % (alpha,a, b) + '\n')
            self.app.log(' ', '\n' + 70 * '#' + '\n')
            self.insert_entry(self.hij_entry, hij)

            #Check for convergence:
            current_tsE, current_rE = self.compute_evb(True)
            if abs(current_tsE - dG_target) < treshold:
                if abs(current_rE - dG0_target) < treshold:
                    done = True
                    self.app.log(' ', '\n\n' + 70 * '#' + '\n')
                    self.app.log('info', 'EVB calibration converged (treshold %.2f)' % treshold)
                    self.app.log(' ', 'Alpha = %.2f    Hij = %.2f' % (alpha, hij))
                    self.app.log(' ', '\n' + 70 * '#' + '\n')
            if counter == 2:
                a_step = (a_step / 2)
                h_step = (h_step / 2)

            if counter == 10:
                self.app.log(' ','\n\n'+ 70*'#'+'\nRRC did not converge.\n' + 70*'#'+ '\n')
                done = True

            #If iteration > 3: Do old search method in between:
            if counter in [2, 4, 6, 8] and not done:
                print('Running fine tuned search...')
                done, alpha, hij, treshold = self.classic_serch(dG_target, current_tsE, dG0_target, current_rE)
                self.insert_entry(self.alpha_entry, alpha)
                self.insert_entry(self.hij_entry, hij)
                current_tsE, current_rE = self.compute_evb(True)
                if done:
                    self.app.log(' ', '\n\n' + 70 * '#' + '\n')
                    self.app.log('info', 'EVB calibration converged (treshold %.2f)' % treshold)
                    self.app.log(' ', 'Alpha = %.2f    Hij = %.2f' % (alpha, hij))
                    self.app.log(' ', '\n' + 70 * '#' + '\n')

            if done:
                self.update_plot()
                break

            if counter > 5:
                treshold = (treshold + 0.2)
                print(('New convergence treshold = %f' % treshold))

    def classic_serch(self, dG_target, dG_current, dG0_target, dG0_current):
        """
        Performs a fine tuned stepping approach towards targe values.
        """
        done = False
        iteration = 1
        treshold = 0.2
        alpha_step = 1.0
        hij_step = 1.0

        ddG_target = (dG_target - dG0_target)
        current_ddG = (dG_current-dG0_current)

        alpha = float(self.alpha_entry.get())
        hij = float(self.hij_entry.get())

        while not done:
            self.update()
            #Use alpha to tune activation energy
            alpha_step = (alpha_step/abs(alpha_step)) * self.set_step_size(dG_target, dG_current)
            alpha += alpha_step
            self.insert_entry(self.alpha_entry, alpha)

            reverse_sign = False

            try:

                new_dG = float(self.compute_evb(True)[0])
                if abs(new_dG - dG_target) <= abs(dG_current - dG_target):
                    dG_current = new_dG
                else:
                    reverse_sign = True
            except:
                reverse_sign = True

            if reverse_sign:
                try:
                    alpha_step = -(alpha_step)
                    alpha += (2 * alpha_step)
                    self.insert_entry(self.alpha_entry, alpha)
                    dG_current = float(self.compute_evb(True)[0])
                except:
                    self.app.log('info', 'Could not fit alpha with classic approach. Try different starting point.')
                    break

            #Check for convergence:
            if abs(dG_current - dG_target) < treshold:
                print('Activation energy now within treshold')
                dG0_current = float(self.compute_evb(True)[1])
                if abs(dG0_current - dG0_target) < treshold:
                    print('Reaction energy now within target.')
                    done = True
                    break
            #Use hij to tune dG(rxn)
            hij_step = (hij_step/abs(hij_step)) * self.set_step_size(dG0_target, dG0_current)
            hij += hij_step
            self.insert_entry(self.hij_entry, hij)

            reverse_sign = False
            try:
                new_dG, new_dG0 = self.compute_evb(True)[0:2]
                new_ddG = (new_dG - new_dG0)
                current_ddG = (dG_current - dG0_current)
                if abs(new_ddG - ddG_target) <= (abs(current_ddG - ddG_target + treshold)):
                    dG_current = new_dG
                    dG0_current = new_dG0
                    current_ddG = new_ddG
                else:
                    reverse_sign = True

            except:
                reverse_sign = True

            if reverse_sign:
                try:
                    hij_step = -hij_step
                    hij += (2 * hij_step)
                    self.insert_entry(self.hij_entry, hij)
                    dG_current, dG0_current = self.compute_evb(True)[0:2]
                    current_ddG = (dG_current - dG0_current)
                except:
                    self.app.log('info', 'Could not fit Hij with classic approach. Try different starting point.')
                    break

            #Check for convergence:
            if abs(current_ddG - ddG_target) <= treshold:
                print('Reaction energy now within treshold')
                if abs(dG_current - dG_target) < treshold:
                    print('Activation energy now within treshold.')
                    done = True
                    break

            iteration += 1

            if iteration in [10, 20, 20]:
                treshold += 0.1
                print(('New treshold: %f' % treshold))

            if iteration > 40:
                break

            print((alpha, hij, dG_current, dG0_current))

        print((alpha, hij, dG_current, dG0_current))

        return done, alpha, hij, treshold

    def set_step_size(self, target, current):
        """
        Decide how large steps to take in classic search
        """
        ddE = abs(target - current)
        step = 1
        if ddE >= 30:
            step = 20
        if ddE >= 20 and ddE < 30:
            step = 10
        elif ddE >= 10 and ddE < 20:
            step = 5
        elif ddE >= 5 and ddE < 10:
            step = 3
        elif ddE >= 2 and ddE < 5:
            step = 2
        elif ddE >= 1 and ddE < 2:
            step = 1
        elif ddE < 1:
            step = 0.5

        return step

    def update_table(self):
        """
        Updates table
        """
        self.listbox.delete(0, END)

        if len(self.ene_files) > 0:
            for enefile in self.ene_files:
                self.listbox.insert(END, '../%s/%s' % (enefile.split('/')[-2], enefile.split('/')[-1]))

    def update_plot(self):
        """
        Update plot
        """
        if len(self.dG) < 1:
            return

        if self.dg_plot:
            self.dg_plot.clear()


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

        self.dg_plot.plot(self.enegaps, self.dG, 'k.', label='Reference')
        self.dg_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})

        self.dg_plot.plot(self.enegaps, self.dG, 'r-')
        self.dg_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}, numpoints=1)

        if self.show_target_lines.get() == 1:
            try:
                dG_target = float(self.target_dG.get())
                dG0_target = float(self.target_dG0.get())
                self.dg_plot.axhline(y=dG_target, ls='--', alpha=0.5)
                self.dg_plot.axhline(y=dG0_target, ls='--', alpha=0.5)
            except:
                pass

        self.canvas.draw()

    def parameters_changed(self, *args):
        """
        Update table/plot whenever LIE parameters are changed
        """
        if self.autoupdate.get() != 1:
            return

        self.compute_evb()

    def target_changed(self, *args):
        """
        Update computed values when target values is toggled
        """
        if len(self.dG) < 1:
            return

        tsE, rE = self.find_local_min_max(self.dG)[1:]
        self.update_ddG(tsE, rE)
        self.update_plot()

    def insert_entry(self, entry, value):
        value = '%.3f' % value
        entry.delete(0, END)
        entry.insert(0, value)

    def add_enefile(self):
        """
        Add EVB energy files
        """
        enefiles = askopenfilenames(parent = self, initialdir = self.app.workdir,
                                   filetypes=(("Energy file", "*.en"),("All files","*.*")))

        if enefiles:
            for enefile in enefiles:
                self.ene_files.append(enefile)
            self.update_table()

    def delete_selected(self):
        """
        Removes selected item from list
        """
        try:
            sel_index = int(self.listbox.curselection()[0])
        except:
            return

        del self.ene_files[sel_index]
        self.listbox.delete(sel_index)

    def sort_files(self):
        """
        Sort or reverse list by names.
        """

        sorted_list = sorted(self.ene_files, key=lambda x: x.split('/')[-1])

        #If file already is sorted, reverse it:
        if sorted_list == self.ene_files:
            self.ene_files.reverse()
        else:
            self.ene_files.sort(key=lambda x: x.split('/')[-1])

        self.update_table()

    def clear_table(self):
        """
        Clears all data in table (resets everything)
        """
        self.listbox.delete(0, END)
        self.ene_files = []
        if self.dg_plot:
            self.dg_plot.clear()
            self.canvas.draw()
            self.dg_plot = None

    def reset_parameters(self):
        self.hij_entry.delete(0,END)
        self.alpha_entry.delete(0, END)
        self.hij_entry.insert(0, '0')
        self.alpha_entry.insert(0, '0')
        #self.update_plot()

    def save_evb_parameters(self):
        """
        Save LIE fit results to file
        """
        savename = None
        savename = asksaveasfilename(parent=self, title='Save EVB parameters', initialdir=self.app.workdir,
                                      filetypes=(("Q EVB", "*.qevb"), ("All files","*.*")),
                                      initialfile = 'evb_rrc.qfep')
        if savename:
            if not savename.endswith('.qfep'):
                savename = savename.split('.')[0] + '.qfep'

        else:
            return

        qfepinp = open(savename, 'w')

        kT = 0.001987209 * float(self.temperature.get())
        skip = int(self.skip_points.get())
        bins = int(self.bins.get())
        #Too many bins will cause Qfep to fail
        if bins > 100:
            self.app.errorBox('Warning', 'Too many bins (%d). Setting bins to 100.' % bins)
            bins = 100
            self.bins.delete(0, END)
            self.bins.insert(0, bins)
        min_points = int(self.binpoints_min.get())
        alpha = float(self.alpha_entry.get())
        hij = float(self.hij_entry.get())
        linear_comb = self.linear_comb.get()

        qfepinp.write('#ene_files\n')
        qfepinp.write('2  0\n')
        qfepinp.write('%.4f  %d\n' % (kT, skip))
        qfepinp.write('%d\n' % bins)
        qfepinp.write('%d\n' % min_points)
        qfepinp.write('%.2f\n' % alpha)
        qfepinp.write('1\n')
        qfepinp.write('1  2  %.2f  0.0 0 0.000\n' % hij)
        qfepinp.write('%s\n\n' % linear_comb)

        qfepinp.close()
        self.app.log('info', 'EVB parameters saved to %s' % savename.split('/')[-1])
        self.app.log(' ', '(%s can be imported in the analyze EVB reaction energies)\n' % savename.split('/')[-1])

    def dialog_window(self):
        self.title('EVB Reference Reaction Calibration')
        self.mainframe = Frame(self, bg=self.main_color)
        self.mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame with import button etc
        frame1 = Frame(self.mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0, pady=(0,10))

        #Frame with FEP settings
        frame5 = Frame(self.mainframe, bg=self.main_color)
        frame5.grid(row=0, column=1)

        #Frame with EVB energy files to fit
        frame2 = Frame(self.mainframe, bg=self.main_color)
        frame2.grid(row=2, column=0)

        #Frame with plot:
        frame3 = Frame(self.mainframe, bg=self.main_color)
        frame3.grid(row=1, rowspan=2, column=1)
        frame_plot= Frame(frame3, bg=self.main_color)
        frame_plot.pack(side=TOP, padx=(10, 10), pady=(10, 0), fill=BOTH)

        #Frame in plot frame with dG values etc.
        frame6 = Frame(self.mainframe, bg=self.main_color)
        frame6.grid(row=3, column=1)

        #Frame with target dG values:
        frame7 = Frame(self.mainframe, bg=self.main_color)
        frame7.grid(row=3, column=0)

        #Frame quit/save
        frame4 = Frame(self.mainframe, bg=self.main_color)
        frame4.grid(row=4, column=0, columnspan=2, pady=(20,0))



        #Frame1:
        import_button = Button(frame1, text='Import from dir', highlightbackground=self.main_color,
                               command=self.import_all_from_dir)
        import_button.grid(row=0, column=0, columnspan=4)

        self.add_button = Button(frame1, text='Add', highlightbackground=self.main_color,
                                 command=self.add_enefile)
        self.add_button.grid(row=1, column=0)

        self.delete_selected = Button(frame1, text='Delete', highlightbackground=self.main_color,
                                      command=self.delete_selected)
        self.delete_selected.grid(row=1, column=1)

        self.edit_selected = Button(frame1, text='Sort', highlightbackground=self.main_color,
                                    command=self.sort_files)
        self.edit_selected.grid(row=1, column=2)

        self.clear_button = Button(frame1, text='Clear all', highlightbackground=self.main_color, command=self.clear_table)
        self.clear_button.grid(row=1, column=3)

        #Frame5 FEP settings:
        temperature = Label(frame5, text='     T     ', bg=self.main_color)
        temperature.grid(row=0, column=0)

        bins = Label(frame5, text='   Bins    ', bg=self.main_color)
        bins.grid(row=0, column=1)

        binpoints = Label(frame5, text='Min. points', bg=self.main_color)
        binpoints.grid(row=0, column=2)

        skip_points = Label(frame5, text='   Skip    ', bg=self.main_color)
        skip_points.grid(row=0, column=3)

        linear_combination = Label(frame5, text='Linear combination:', bg=self.main_color)
        linear_combination.grid(row=0, rowspan=2, column=4, padx=(10,5))

        state1 = Label(frame5, text="\N{GREEK CAPITAL LETTER PHI}(1)", bg=self.main_color)
        state1.grid(row=0, column=5)

        state2 = Label(frame5, text="\N{GREEK CAPITAL LETTER PHI}(2)", bg=self.main_color)
        state2.grid(row=0, column=6)

        self.temperature = Entry(frame5, width=5, highlightthickness=0)
        self.temperature.grid(row=1, column=0)

        self.bins = Entry(frame5, width=5, highlightthickness=0)
        self.bins.grid(row=1, column=1)

        self.binpoints_min = Entry(frame5, width=5, highlightthickness=0)
        self.binpoints_min.grid(row=1, column=2)

        self.skip_points = Entry(frame5, width=5, highlightthickness=0)
        self.skip_points.grid(row=1, column=3)

        self.linear_comb = Spinbox(frame5, width=6, highlightthickness=0, relief=GROOVE,
                                   values=(' 1    -1 ', '  -1    1'))
        self.linear_comb.grid(row=1, column=5, columnspan=2)

        #Frame2 (Table with data)
        table_heading = Text(frame2, width=40, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        table_heading.grid(row=0, column=0, columnspan=10)
        table_heading.config(font=tkinter.font.Font(family="Courier", size=12))
        table_heading.insert(0.0, 'Energy files:')

        listbox_scroll = Scrollbar(frame2)
        listbox_scroll.grid(row = 1, column = 10, sticky = 'nsw', padx=(0,10))
        self.listbox = Listbox(frame2, yscrollcommand = listbox_scroll.set, width=40, height=21,
                               highlightthickness=0, relief=GROOVE, selectmode=BROWSE)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row=1, column = 0, columnspan=9, sticky = 'e')
        self.listbox.config(font=tkinter.font.Font(family="Courier", size=12))

        alpha = Text(frame2, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        alpha.tag_configure("subscript", offset=-1)
        alpha.insert("insert","\N{GREEK SMALL LETTER ALPHA}","",'ij','subscript')
        alpha.grid(row=2, column=0)
        alpha.config(state=DISABLED)

        self.alpha_entry = Spinbox(frame2, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-500.00, to=500.00, increment=1, textvariable=self.alpha_var)
        self.alpha_entry.grid(row=3, column=0)

        hij = Text(frame2, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        hij.tag_configure("subscript", offset=-1)
        hij.insert("insert",'H',"",'ij','subscript')
        hij.grid(row=2, column=1)
        hij.config(state=DISABLED)

        autoupdate = Label(frame2, text='Auto update:', bg=self.main_color)
        autoupdate.grid(row=2, column=2, columnspan=2, sticky='w')

        self.autoupdate_check = Checkbutton(frame2, variable=self.autoupdate, bg=self.main_color)
        self.autoupdate_check.grid(row=2, column=3, sticky='e')

        self.hij_entry = Spinbox(frame2, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=0, to=200, increment=1.0, textvariable=self.hij_var)
        self.hij_entry.grid(row=3, column=1)


        resetbutton = Button(frame2, text='Reset', highlightbackground=self.main_color, command=self.reset_parameters)
        resetbutton.grid(row=3, column=2)

        compute_button = Button(frame2, text='Compute', highlightbackground=self.main_color, command=self.compute_evb)
        compute_button.grid(row=3, column=3)


        #Frame3 (plot)
        self.plot_window = Figure(figsize=(5.5,3), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.plot_window, master=frame_plot)
        self.plot_window.patch.set_facecolor('white')

        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame_plot)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        #Frame 6 (real dG values) u"\u2021" (doble dagger)
        computed_label = Label(frame6, text='Computed: ', bg=self.main_color)
        computed_label.grid(row=0, rowspan=2, column=0)

        real_dG = Text(frame6, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        real_dG.tag_configure("superscript", offset=3)
        real_dG.insert("insert","\N{GREEK CAPITAL LETTER DELTA}G","","\\u2021",'superscript')
        real_dG.grid(row=0, column=1)
        real_dG.config(state=DISABLED)

        self.dG_real_txt = Text(frame6, width=7, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.dG_real_txt.grid(row=1, column=1, pady=(5,0), padx=(10, 10))
        self.dG_real_txt.config(state=DISABLED)

        real_dG0 = Text(frame6, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        real_dG0.tag_configure("superscript", offset=3)
        real_dG0.insert("insert","\N{GREEK CAPITAL LETTER DELTA}G","", 'o', 'superscript')
        real_dG0.grid(row=0, column=2)
        real_dG0.config(state=DISABLED)

        self.dG0_real_txt = Text(frame6, width=7, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.dG0_real_txt.grid(row=1, column=2, pady=(5,0), padx=(10, 10))
        self.dG0_real_txt.config(state=DISABLED)

        real_ddG = Text(frame6, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        real_ddG.tag_configure("superscript", offset=3)
        real_ddG.insert("insert","\N{GREEK CAPITAL LETTER DELTA}\N{GREEK CAPITAL LETTER DELTA}G","","\\u2021",
                        'superscript')
        real_ddG.grid(row=0, column=3)
        real_ddG.config(state=DISABLED)

        self.ddG_real_txt = Text(frame6, width=7, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.ddG_real_txt.grid(row=1, column=3, pady=(5,0), padx=(10, 10))
        self.ddG_real_txt.config(state=DISABLED)

        real_ddG0 = Text(frame6, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        real_ddG0.tag_configure("superscript", offset=3)
        real_ddG0.insert("insert","\N{GREEK CAPITAL LETTER DELTA}\N{GREEK CAPITAL LETTER DELTA}G","",'o','superscript')
        real_ddG0.grid(row=0, column=4)
        real_ddG0.config(state=DISABLED)

        self.ddG0_real_txt = Text(frame6, width=7, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        self.ddG0_real_txt.grid(row=1, column=4, pady=(5,0), padx=(10, 10))
        self.ddG0_real_txt.config(state=DISABLED)

        #Frame7 target dG values:
        target_label = Label(frame7, text='Target: ', bg=self.main_color)
        target_label.grid(row=0, rowspan=2, column=0)

        target_dG = Text(frame7, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        target_dG.tag_configure("superscript", offset=3)
        target_dG.insert("insert","\N{GREEK CAPITAL LETTER DELTA}G","","\\u2021",'superscript')
        target_dG.grid(row=0, column=1)
        target_dG.config(state=DISABLED)

        self.target_dG = Spinbox(frame7, width=5, highlightthickness=0, relief=GROOVE,
                                  from_=-100.00, to=100.00, increment=1, textvariable=self.target_dG_var)
        self.target_dG.grid(row=1, column=1, pady=(5, 0), padx=(5, 5))

        target_dG0 = Text(frame7, width=5, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        target_dG0.tag_configure("superscript", offset=3)
        target_dG0.insert("insert","\N{GREEK CAPITAL LETTER DELTA}G","", 'o', 'superscript')
        target_dG0.grid(row=0, column=2)
        target_dG0.config(state=DISABLED)

        self.target_dG0 = Spinbox(frame7, width=5, highlightthickness=0, relief=GROOVE,
                                  from_=-100.00, to=100.00, increment=1, textvariable=self.target_dG0_var)
        self.target_dG0.grid(row=1, column=2, pady=(5, 0), padx=(5, 5))

        scan_button = Button(frame7, text='Scan', highlightbackground=self.main_color, command=self.scan_evb)
        scan_button.grid(row=0, rowspan=2, column=3, pady=(5,0), padx=(5, 5))

        display_line = Label(frame7, text='Show target lines:', bg=self.main_color)
        display_line.grid(row=2, column=1, columnspan=2)

        target_line_check = Checkbutton(frame7, variable=self.show_target_lines, bg=self.main_color,
                                        command=self.update_plot)
        target_line_check.grid(row=2, column=3)




        #Frame 4 (close/save)
        savebutton = Button(frame4, text='Save', highlightbackground=self.main_color, command=self.save_evb_parameters)
        savebutton.grid(row=0, column=0)

        closebutton = Button(frame4, text='Close', highlightbackground=self.main_color, command=self.destroy)
        closebutton.grid(row=0, column=1)
