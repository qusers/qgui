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

from Tkinter import Label, TOP, Button, Listbox, Scrollbar, BROWSE, Spinbox, Entry, LabelFrame, Text, Frame, \
    Toplevel, DISABLED, END, GROOVE, NORMAL, BOTH, IntVar, StringVar, Checkbutton

import tkFont
import numpy as np
import os

import matplotlib
matplotlib.use('TkAgg')
#Implement default mpl key bindings
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import  key_press_handler
from matplotlib.figure import Figure
from tkFileDialog import askopenfilename, asksaveasfilename
from edit_lie import EditLIE


class FitQlie(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root

        #List with LIE data to fill table:
        self.lie_data = []

        #Control what parameters to change:
        self.fit_alpha_var = IntVar()
        self.fit_alpha_var.set(1)
        self.fit_beta_var = IntVar()
        self.fit_beta_var.set(0)
        self.fit_gamma_var = IntVar()
        self.fit_gamma_var.set(1)


        #Trace changes in parameters to update table/plot
        self.alpha_var = StringVar()
        self.alpha_var.trace('w', self.parameters_changed)

        self.beta_var = StringVar()
        self.beta_var.trace('w', self.parameters_changed)

        self.gamma_var = StringVar()
        self.gamma_var.trace('w', self.parameters_changed)

        self.dg_plot = None

        self.dialog_window()
        self.canvas.mpl_connect('key_press_event', self.on_key_event)

        #Insert default LIE parameter values:
        self.alpha_entry.delete(0, END)
        self.alpha_entry.insert(0, '0.18')
        self.beta_entry.delete(0,END)
        self.beta_entry.insert(0, '0.50')
        self.gamma_entry.delete(0,END)
        self.gamma_entry.insert(0, '0.00')

    def on_key_event(self, event):
        print('you pressed %s' % event.key)
        key_press_handler(event, self.canvas, self.toolbar)

    def import_lie_results(self):
        """
        Ask user for *.qlie and fill table with content from file
        """
        filename = askopenfilename(parent = self, initialdir = self.app.app.workdir,
                                   filetypes=(("Q-LIE", "*.qlie"),("All files","*.*")))

        if filename:
            if not os.path.isfile(filename):
                return
            oldfile = open(filename, 'r').readlines()

            for i in range(len(oldfile)):
                if '#' in oldfile[i][0]:
                    if not '#NAME' in oldfile[i].split()[0]:
                        ligname = oldfile[i].split()[0]
                        d_el = (float(oldfile[i+3].split()[3]) + float(oldfile[i+3].split()[5])) - (
                            float(oldfile[i+1].split()[3]) + float(oldfile[i+1].split()[5]))
                        d_vdw = (float(oldfile[i+3].split()[4]) + float(oldfile[i+3].split()[6])) - (
                            float(oldfile[i+1].split()[4]) + float(oldfile[i+1].split()[6]))
                        alpha = float(oldfile[i + 1].split()[9])
                        beta = float(oldfile[i + 1].split()[10])
                        gamma = float(oldfile[i + 1].split()[11])
                        dG = float(oldfile[i+1].split()[12])
                        dg_exp = float(oldfile[i+1].split()[13])

                        try:
                            corr = float(oldfile[i+3].split()[13]) - float(oldfile[i+1].split()[13])
                            dG += corr
                        except:
                            corr = 0
                        tmplist = [ligname, beta, d_el, alpha, d_vdw, gamma, dG, dg_exp, corr]
                        self.lie_data.append(tmplist)
                    elif '#NAME' in oldfile[i].split()[0]:
                        ligname = oldfile[i+1].split()[0]
                        beta = float(oldfile[i+1].split()[1])
                        d_el = float(oldfile[i+1].split()[2])
                        alpha = float(oldfile[i+1].split()[3])
                        d_vdw = float(oldfile[i+1].split()[4])
                        gamma = float(oldfile[i+1].split()[5])
                        dG = float(oldfile[i+1].split()[6])
                        dg_exp = float(oldfile[i+1].split()[7])

                        try:
                            corr = float(oldfile[i+1].split()[8])
                            dG += corr
                        except:
                            corr = 0.00
                        tmplist = [ligname, beta, d_el, alpha, d_vdw, gamma, dG, dg_exp, corr]
                        self.lie_data.append(tmplist)

        else:
            return

        self.update_table()
        self.update_plot()

    def update_table(self):
        """
        Updates table
        """
        self.listbox.delete(0, END)
        for i in range(len(self.lie_data)):
            self.listbox.insert(END,'%11s %5.2f %8.2f %5.2f %8.2f %7.2f %7.2f %7.2f' % (
                self.lie_data[i][0].ljust(11), self.lie_data[i][1],self.lie_data[i][2], self.lie_data[i][3], self.lie_data[i][4],
                self.lie_data[i][5], self.lie_data[i][6], self.lie_data[i][7]))

    def linear_regression(self, x, y):
        """
        Performs linear regression f(x) = ax + b.
        returns a and b
        """
        x = map(float, x)
        y = map(float, y)

        A = np.array([x, np.ones(len(x))])
        w = np.linalg.lstsq(A.T, y)[0]

        return w[0], w[1]

    def sum_errors_squared(self, predicted, real):
        """
        Return the residual sum of errors (SSE or RSS) squared between fitted function
        """
        predicted = np.array(predicted)
        real = np.array(real)
        r_i = (predicted - real)

        sse = np.sum(map(lambda x: x ** 2, r_i))

        return sse

    def explained_sum_squares(self, predicted, real):
        """
        Returns the correlation coefficient
        The explained sum of squares (ESS) is the sum of the squares of the deviations
        of the predicted values from the mean value of a response variable.
        """
        ave_real = np.average(real)

        ess = np.sum(map(lambda x: (x - ave_real) **2, predicted))

        return ess

    def coefficient_determination(self, predicted, real):
        """
        Returns R^2
        TSS = RSS + ESS
        """
        ess = self.explained_sum_squares(predicted, real)
        rss = self.sum_errors_squared(predicted, real)

        #r_sq = 1.0 - (rss / (rss + ess))
        r_sq = ess / (rss+ess)
        return r_sq

    def update_plot(self, return_stat=False):
        """
        Update plot
        """
        if self.dg_plot:
            self.dg_plot.clear()
        titles = []
        dG = []
        dG_exp = []

        for i in range(len(self.lie_data)):
            titles.append(self.lie_data[i][0])
            dG.append(float(self.lie_data[i][-3]))
            dG_exp.append(float(self.lie_data[i][-2]))

        #Get regrssion line for points:
        if len(self.lie_data) > 0:
            a,b = self.linear_regression(dG_exp, dG)
            dG_reg = map(lambda x: (x * a) + b, dG_exp)

        #Set color cycle for plots:
        matplotlib.rcParams['axes.color_cycle'] = ['k', 'b', 'g', 'r', 'm', 'y', 'c', 'brown',
                                                   'burlyWood', 'cadetBlue', 'DarkGreen', 'DarkBlue',
                                                   'DarkMagenta', 'DarkSalmon', 'DimGray', 'Gold']

        #Create subplot
        self.dg_plot = self.plot_window.add_subplot(111, axisbg='white')
        self.plot_window.subplots_adjust(hspace=0.5)

        #X/Y labels
        self.dg_plot.set_xlabel(r'$\Delta G_{bind}^{exp}$ (kcal/mol)')
        self.dg_plot.set_ylabel(r'$\Delta G_{bind}^{calc}$ (kcal/mol)')

        #Move label box outside plot region
        box = self.dg_plot.get_position()
        self.dg_plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        #Fit subplot to figure/canvas
        #rect=(left,bottom,top,right)
        self.plot_window.tight_layout(rect=(0.005, 0, 0.8, 1))


        for point in range(len(dG)):
            self.dg_plot.plot(dG_exp[point], dG[point], 'o', label='%s' % titles[point].split('#')[-1])
            self.dg_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})



        if len(self.lie_data) > 0:
            #Plot regression line:
            self.dg_plot.plot(dG_exp, dG_reg, 'b--')
            #Plot diagonal line:
            self.dg_plot.plot(dG_exp, dG_exp, 'k-')
            #Insert SSE for dG vs dG_reg:
            sse = self.sum_errors_squared(dG_reg, dG)
            self.sum_se.delete(0.0, END)
            self.sum_se.insert(0.0, '%8.4f' % sse)
            sse_0 = self.sum_errors_squared(dG, dG_exp)
            self.sum_se0.delete(0.0, END)
            self.sum_se0.insert(0.0, '%8.4f' % sse_0)
            #Insert coefficients of determination
            r = self.coefficient_determination(dG, dG_reg)
            self.r_coef.delete(0.0, END)
            self.r_coef.insert(0.0, '%8.4f' % r)
            r_0 = self.coefficient_determination(dG, dG_exp)
            self.r0_coef.delete(0.0, END)
            self.r0_coef.insert(0.0, '%8.4f' % r_0)

        self.canvas.show()
        if return_stat:
            if len(self.lie_data) > 0:
                return sse_0

    def parameters_changed(self, *args):
        """
        Update table/plot whenever LIE parameters are changed
        """
        if len(self.lie_data) < 1:
            return

        lie_data_index = [3, 1, 5]
        new_parameters = ['alpha','beta','gamma']
        try:
            if self.fit_alpha_var.get() == 1:
                new_parameters[0] = float(self.alpha_entry.get())
            if self.fit_beta_var.get() == 1:
                new_parameters[1] = float(self.beta_entry.get())
            if self.fit_gamma_var.get() == 1:
                new_parameters[2] = float(self.gamma_entry.get())
        except:
            return

        for i in range(len(self.lie_data)):
            for j in range(len(new_parameters)):
                if str(new_parameters[j])[-1].isdigit():
                    self.lie_data[i][lie_data_index[j]] = new_parameters[j]
            self.lie_data[i][-3] = ( (self.lie_data[i][1] * self.lie_data[i][2]) + (
                self.lie_data[i][3] * self.lie_data[i][4] ) + (
                self.lie_data[i][5] + self.lie_data[i][-1]))


        self.update_table()
        self.update_plot()

    def insert_entry(self, entry, value):
        value = '%5.3f' % value
        entry.delete(0, END)
        entry.insert(0, value)

    def autofit(self):
        """
        Auto fits alpha/beta/gamma...
        Minimize sum of squared errors
        """
        if len(self.lie_data) < 1:
            return
        parameter_var = [self.fit_alpha_var.get(), self.fit_beta_var.get(), self.fit_gamma_var.get()]
        entries = [self.alpha_entry, self.beta_entry, self.gamma_entry]
        param_names = ['alpha', 'beta', 'gamma']
        #Get parameters to fit
        fit_entries = []
        fit_names = []
        for i in range(len(parameter_var)):
            if parameter_var[i] == 1:
                fit_entries.append(entries[i])
                fit_names.append(param_names[i])

        fit_names.reverse()
        fit_entries.reverse()


        for lie_param in range(len(fit_names)):
            self.app.app.log('info', 'Auto fit %s started' % fit_names[lie_param])
            param_value = float(fit_entries[lie_param].get())
            sse = self.update_plot(True)
            step = 0.02
            if fit_names[lie_param] == 'gamma':
                step = 0.4

            x = []
            y = []

            #Decide sign on step:
            param_value += step
            self.insert_entry(fit_entries[lie_param], param_value)
            new_sse = self.update_plot(True)
            if new_sse > sse:
                step = -step
                param_value += step
            else:
                sse = new_sse
                x.append(float(param_value))
                y.append(float(sse))
            counter = 0
            done = False
            while not done:
                param_value += step
                self.insert_entry(fit_entries[lie_param], param_value)
                sse_new = self.update_plot(True)

                if sse_new > sse:
                    counter += 1
                    if counter > 2:
                        done = True
                elif sse_new < sse:
                    sse = sse_new
                x.append(float(param_value))
                y.append(float(sse_new))

            #Make second order polynomial from points
            z = np.polyfit(x, y, 2)

            #Get parameter minima from dz/dparameter = 0
            param_value = -z[1] / (2 * z[0])

            self.app.app.log(' ', 'Sum of Errors squared minimum with %s = %5.3f\n' %
                                  (fit_names[lie_param], float(param_value)))

            #Update final data/plot
            self.insert_entry(fit_entries[lie_param], param_value)
            self.update_plot()

    def edit_or_add(self, edit=True):
        """
        Edit or add LIE data to table
        """
        name = 'NEW'
        beta = 0.50
        d_el = 0.00
        alpha = 0.18
        d_vdw = 0.00
        gamma = 0.00
        dG = 0.00
        dG_exp = 0.00
        sel_index = 'END'
        if edit:
            #Edit existing entry in list
            try:
                sel_index = int(self.listbox.curselection()[0])
            except:
                return
            name = self.lie_data[sel_index][0]
            beta = self.lie_data[sel_index][1]
            d_el = self.lie_data[sel_index][2]
            alpha = self.lie_data[sel_index][3]
            d_vdw = self.lie_data[sel_index][4]
            gamma = self.lie_data[sel_index][5]
            dG = self.lie_data[sel_index][6]
            dG_exp = self.lie_data[sel_index][7]

        data_to_insert = [name, beta, d_el, alpha, d_vdw, gamma, dG, dG_exp]
        self.lie_edit = EditLIE(self, self.root, edit, data_to_insert, sel_index)

    def delete_selected(self):
        """
        Removes selected item from list
        """
        try:
            sel_index = int(self.listbox.curselection()[0])
        except:
            return

        del self.lie_data[sel_index]
        self.listbox.delete(sel_index)
        self.update_plot()

    def clear_table(self):
        """
        Clears all data in table (resets everything)
        """
        self.listbox.delete(0, END)
        self.lie_data = []
        if self.dg_plot:
            self.dg_plot.clear()
            self.canvas.show()
            self.dg_plot = None

        self.sum_se.delete(0.0, END)
        self.sum_se0.delete(0.0, END)
        self.r_coef.delete(0.0, END)
        self.r0_coef.delete(0.0, END)

        self.alpha_entry.delete(0,END)
        self.beta_entry.delete(0,END)
        self.gamma_entry.delete(0,END)

        self.alpha_entry.insert(0, 0.18)
        self.beta_entry.insert(0, 0.50)
        self.gamma_entry.insert(0, 0.00)

    def reset_parameters(self):
        self.alpha_entry.delete(0, END)
        self.beta_entry.delete(0, END)
        self.gamma_entry.delete(0, END)

        self.alpha_entry.insert(0, 0.18)
        self.beta_entry.insert(0, 0.50)
        self.gamma_entry.insert(0, 0.00)

        self.update_plot()

    def save_lie(self):
        """
        Save LIE fit results to file
        """
        savename = asksaveasfilename(parent=self, title='Save LIE results', initialdir=self.app.app.workdir,
                                      filetypes=(("Q LIE", "*.qlie"), ("All files","*.*")),
                                      initialfile = 'results.qlie')
        if savename:
            if not savename.endswith('.qlie'):
                savename = savename.split('.')[0] + '.qlie'
            if os.path.isfile(savename):
                oldfile = open(savename, 'r').readlines()
            else:
                oldfile = []

            savefile = open(savename, 'w')
            #If file exist, write old stuff to file first
            if len(oldfile)>0:
                for line in oldfile:
                    savefile.write(line)


            heading = '\n#%10s %8s %8s %8s %8s %8s %8s %8s\n' % ('NAME'.ljust(10),
                                                'beta', 'd_el', 'alpha','d_vdW', 'gamma', 'dG', 'dG(exp)')
            for i in range(len(self.lie_data)):
                name = self.lie_data[i][0].split('#')[-1]
                beta = self.lie_data[i][1]
                d_el = self.lie_data[i][2]
                alpha = self.lie_data[i][3]
                d_vdw = self.lie_data[i][4]
                gamma = self.lie_data[i][5]
                dG = self.lie_data[i][6]
                dG_exp = self.lie_data[i][7]
                savefile.write(heading)
                savefile.write('%11s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' %
                               (name.ljust(11), beta, d_el, alpha, d_vdw, gamma, dG, dG_exp))

            savefile.close()
        else:
            return

    def dialog_window(self):
        self.title('Fit LIE parameters')
        self.mainframe = Frame(self, bg=self.main_color)
        self.mainframe.pack(fill='both', padx=(10, 10), pady=(10, 10))

        #Frame with import button etc
        frame1 = Frame(self.mainframe, bg=self.main_color)
        frame1.grid(row=0, column=0, pady=(0,10))

        #Frame with LIE results to fit
        frame2 = Frame(self.mainframe, bg=self.main_color)
        frame2.grid(row=2, column=0)

        #Frame with plot:
        frame3 = Frame(self.mainframe, bg=self.main_color)
        frame3.grid(row=2, column=1)
        frame_plot= Frame(frame3, bg=self.main_color)
        frame_plot.pack(side=TOP, padx=(10, 10), pady=(10, 10), fill=BOTH)

        #Frame with control panel for fitting:
        frame4 = Frame(self.mainframe, bg=self.main_color)
        frame4.grid(row=3, column=0, columnspan=2)



        #Frame1:
        import_button = Button(frame1, text='Import LIE results', highlightbackground=self.main_color,
                               command=self.import_lie_results)
        import_button.grid(row=0, column=0, columnspan=4)

        self.add_button = Button(frame1, text='Add', highlightbackground=self.main_color,
                                 command=lambda: self.edit_or_add(False))
        self.add_button.grid(row=1, column=0)

        self.edit_selected = Button(frame1, text='Edit selected', highlightbackground=self.main_color,
                                    command=lambda: self.edit_or_add(True))
        self.edit_selected.grid(row=1, column=1)

        self.delete_selected = Button(frame1, text='Delete selected', highlightbackground=self.main_color,
                                      command=self.delete_selected)
        self.delete_selected.grid(row=1, column=2)

        self.clear_button = Button(frame1, text='Clear table', highlightbackground=self.main_color, command=self.clear_table)
        self.clear_button.grid(row=1, column=3)



        #Frame2 (Table with data)
        table_heading = Text(frame2, width=70, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        table_heading.grid(row=0, column=0, columnspan=10)
        table_heading.config(font=tkFont.Font(family="Courier", size=12))
        table_heading.insert(0.0, u"#NAME         \N{GREEK SMALL LETTER BETA}     \N{GREEK CAPITAL LETTER DELTA}E(el)"
                                  u"   \N{GREEK SMALL LETTER ALPHA}    \N{GREEK CAPITAL LETTER DELTA}E(vdW)     "
                                  u"\N{GREEK SMALL LETTER GAMMA}    "
                                  u"  \N{GREEK CAPITAL LETTER DELTA}G     \N{GREEK CAPITAL LETTER DELTA}G(exp)")

        listbox_scroll = Scrollbar(frame2)
        listbox_scroll.grid(row = 1, column = 10, sticky = 'nsw', padx=(0,10))
        self.listbox = Listbox(frame2, yscrollcommand = listbox_scroll.set, width=70, height=15,
                               highlightthickness=0, relief=GROOVE, selectmode=BROWSE)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row=1, column = 0, columnspan=9, sticky = 'w')
        self.listbox.config(font=tkFont.Font(family="Courier", size=12))


        alpha_label = Label(frame2, text= u"\N{GREEK SMALL LETTER ALPHA}", bg=self.main_color)
        alpha_label.grid(row=2, column=0)
        self.alpha_entry = Spinbox(frame2, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-1.00, to=1.00, increment=0.01, textvariable=self.alpha_var)
        self.alpha_entry.grid(row=3, column=0)

        beta_label = Label(frame2, text= u"\N{GREEK SMALL LETTER BETA}", bg=self.main_color)
        beta_label.grid(row=2, column=1)
        self.beta_entry = Spinbox(frame2, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-1.00, to=1.00, increment=0.01, textvariable=self.beta_var)
        self.beta_entry.grid(row=3, column=1)

        gamma_label = Label(frame2, text= u"\N{GREEK SMALL LETTER GAMMA}", bg=self.main_color)
        gamma_label.grid(row=2, column=2)
        self.gamma_entry = Spinbox(frame2, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=-99.99, to=99.99, increment=0.05, textvariable=self.gamma_var)
        self.gamma_entry.grid(row=3, column=2)

        self.auto_fit_button = Button(frame2, text='Auto fit', highlightbackground=self.main_color,
                                      command=self.autofit)
        self.auto_fit_button.grid(row=3, column=3)

        resetbutton = Button(frame2, text='Reset', highlightbackground=self.main_color, command=self.reset_parameters)
        resetbutton.grid(row=3, column=4)

        alpha_check = Checkbutton(frame2, bg=self.main_color, variable=self.fit_alpha_var)
        alpha_check.grid(row=4, column=0)

        beta_check = Checkbutton(frame2, bg=self.main_color, variable=self.fit_beta_var)
        beta_check.grid(row=4, column=1)

        gamma_check = Checkbutton(frame2, bg=self.main_color, variable=self.fit_gamma_var)
        gamma_check.grid(row=4, column=2)

        sum_squarer_error = Text(frame2, width=6, height=1, bg=self.main_color,borderwidth=0, highlightthickness=0)
        sum_squarer_error.insert(END, 'SSE = ')
        sum_squarer_error.grid(row=5, column=0, sticky='e')
        sum_squarer_error.config(state=DISABLED)

        self.sum_se = Text(frame2, width=7, height=1, bg=self.main_color, highlightthickness=0, borderwidth=0)
        self.sum_se.grid(row=5, column=1)

        r_label = Text(frame2, width=7, height=1, bg=self.main_color,borderwidth=0, highlightthickness=0)
        r_label.insert(END, 'COD = ')
        r_label.grid(row=5, column=2,sticky='e')
        r_label.config(state=DISABLED)

        self.r_coef = Text(frame2, width=7, height=1, borderwidth=0, highlightthickness=0, bg=self.main_color)
        self.r_coef.grid(row=5, column=3)

        sum_se = Text(frame2, width=7, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        sum_se.tag_configure("subscript", offset=-3)
        sum_se.insert("insert",'SSE',"",'0','subscript')
        sum_se.insert(END,' = ')
        sum_se.grid(row=6, column=0, sticky='e')
        sum_se.config(state=DISABLED)

        self.sum_se0 = Text(frame2, width=7, height=1, bg=self.main_color, highlightthickness=0, borderwidth=0)
        self.sum_se0.grid(row=6, column=1)

        r0_label = Text(frame2, width=8, height=1, bg=self.main_color, borderwidth=0, highlightthickness=0)
        r0_label.tag_configure("subscript", offset=-3)
        r0_label.insert("insert",'COD',"",'0','subscript')
        r0_label.insert(END,' = ')
        r0_label.grid(row=6, column=2, sticky='e')
        r0_label.config(state=DISABLED)

        self.r0_coef = Text(frame2, width=7, height=1, borderwidth=0, highlightthickness=0, bg=self.main_color)
        self.r0_coef.grid(row=6, column=3)


        #Frame3 (plot)
        self.plot_window = Figure(figsize=(5.5,3), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.plot_window, master=frame_plot)
        self.plot_window.patch.set_facecolor('white')

        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, frame_plot)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        #Frame 4 (close/save)
        savebutton = Button(frame4, text='Save', highlightbackground=self.main_color, command=self.save_lie)
        savebutton.grid(row=0, column=0)

        closebutton = Button(frame4, text='Close', highlightbackground=self.main_color, command=self.destroy)
        closebutton.grid(row=0, column=1)
