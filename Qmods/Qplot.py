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

# BB: matplotlib changes on lines 64, 85, 103, 164
# NavigationToolbar2TkAgg -> NavigationToolbar2Tk

from Tkinter import X, BOTH, TOP, Button, Frame, Toplevel, Label, Spinbox, GROOVE


import matplotlib
matplotlib.use('TkAgg')
#Implement default mpl key bindings
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from cycler import cycler


class Qplot(Toplevel):
    def __init__(self, app, root, ylist=[], xlist=[], plot_titles=[], line_titles=[], xlabel='X', ylabel='Y'):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.ylist = ylist
        self.xlist = xlist
        self.plot_titles = plot_titles
        self.line_titles = line_titles
        self.xlabel = xlabel
        self.ylabel = ylabel

        self.plotlist = None

        self.dialog_window()
        self.canvas.mpl_connect('key_press_event', self.on_key_event)

        self.make_plot(self.xlist, self.ylist, self.line_titles)

    def make_plot(self, xlist, ylist, titles):

        if self.plotlist:
            for i in range(len(self.plotlist)):
                self.plotlist[i].clear()

        plots_to_make = 1
        try:
            if len(ylist[1][0]) > 0:
                plots_to_make = len(ylist)
        except:
            plots_to_make = 1


        #Default color cycle
        matplotlib.rcParams['axes.prop_cycle'] = cycler('color', ['k', 'b', 'g', 'r', 'm', 'y', 'c', 'brown',
                                                   'burlyWood', 'cadetBlue', 'DarkGreen', 'DarkBlue',
                                                   'DarkMagenta', 'DarkSalmon', 'DimGray', 'Gold'])

        #Find min X and max X to scale xrange:
        min_x = 100000.00
        max_x = -100000.00
        print len(xlist)
        print len(xlist[0])
        print len(xlist[0][0])
        for plot in range(len(xlist)):
            for types in range(len(xlist[plot])):
                if xlist[plot][types][0] < min_x:
                    min_x = xlist[plot][types][0]
                if xlist[plot][types][-1] > max_x:
                    max_x = xlist[plot][types][-1]

        self.plotlist = []

        for plot in range(1, plots_to_make + 1):
            placement = '%d1%d' % (plots_to_make, plot)
            self.plotlist.append(self.plot_window.add_subplot(int(placement), facecolor='white'))
            self.plot_window.subplots_adjust(hspace=0.5)

        for plots in range(len(self.plotlist)):
            self.plotlist[plots].set_xlabel(self.xlabel)
            self.plotlist[plots].set_ylabel(self.ylabel)
            self.plotlist[plots].set_title(self.plot_titles[plots])
            #Shring axis by 20% to fit legend outside plotting area
            box = self.plotlist[plots].get_position()
            self.plotlist[plots].set_position([box.x0, box.y0, box.width * 0.8, box.height])

            for types in range(len(ylist[plots])):
                self.plotlist[plots].plot(xlist[plots][types], ylist[plots][types],
                                     label='%s' % titles[plots][types])
                self.plotlist[plots].legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
                #plotlist[plots].set_autoscale_on(False)
                self.plotlist[plots].set_xlim([min_x, max_x])

        self.canvas.draw()

    def on_key_event(self, event):
        print('you pressed %s' % event.key)
        key_press_handler(event, self.canvas, self.toolbar)

    def update_plot(self):


        #self.plot_window.clear()

        plot_int = 1
        try:
            plot_int = int(self.interval.get())
        except:
            plot_int = 1

        titles = list()
        xlist = list()
        ylist = list()

        for plot in range(len(self.xlist)):
            tmp_x = list()
            tmp_y = list()
            titles.append(self.line_titles[plot])
            for term in range(len(self.xlist[plot])):
                tmp_tmp_x = list()
                tmp_tmp_y = list()
                for i in range(0, len(self.xlist[plot][term]), plot_int):
                    tmp_tmp_x.append(self.xlist[plot][term][i])
                    tmp_tmp_y.append(self.ylist[plot][term][i])
                tmp_x.append(tmp_tmp_x)
                tmp_y.append(tmp_tmp_y)
            xlist.append(tmp_x)
            ylist.append(tmp_y)

        self.make_plot(xlist, ylist, titles)

    def _quit(self):
        #self.quit()
        self.destroy()

    def dialog_window(self):
        """
        Plot window
        """
        self.title('Q-plot')
        #self.geometry("900x700+100+200")
        #self.maxsize(300,150)
        self.minsize(850, 650)

        frame = Frame(self, bg=self.main_color)
        frame.pack(side=TOP, padx=(10, 10), pady=(10, 10), fill=X)

        frame2 = Frame(self, bg=self.main_color)
        frame2.pack(side=TOP, padx=(10, 10), pady=(10, 10), fill=X)

        # a tk.DrawingArea
        self.plot_window = Figure(figsize=(8,5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.plot_window, master=frame)
        self.plot_window.patch.set_facecolor('white')
        self.canvas.draw() #self.canvas.show()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        int_label = Label(frame2, text='Interval:', bg=self.main_color)
        int_label.grid(row=1, column=0)

        self.interval = Spinbox(frame2, width=7, highlightthickness=0, relief=GROOVE,
                                  from_=1, to=10000, increment=1)
        self.interval.grid(row=1, column=1)

        update_plot = Button(frame2, text='Update plot', highlightbackground=self.main_color, command=self.update_plot)
        update_plot.grid(row=1, column=10)

        quit_button = Button(frame2, text='Quit', highlightbackground=self.main_color, command=self._quit)
        quit_button.grid(row=2, column=0, columnspan=10)








