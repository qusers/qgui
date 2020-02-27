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

from tkinter import PhotoImage, Toplevel, Canvas, Scrollbar, Listbox, Checkbutton, DISABLED, NORMAL, END, GROOVE, IntVar
import os
import sys

class SplashScreen(Toplevel):
    def __init__(self, app, root, image="/Qmods/qgui_box.gif", timeout=2000):
        """
        create a splash screen from a specified image file
        keep splash screen up for timeout milliseconds
        """
        Toplevel.__init__(self, root)
        self.root = root
        self.app = app

        if 'darwin' in sys.platform:
            try:
                os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python2" to true' ''')
            except:
                os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

        # don't show main window
        #self.main.withdraw()
        self.overrideredirect(1)

        # use Tkinter's PhotoImage for .gif files
        self.image = PhotoImage(file=image)
        self.after_idle(self.centerOnScreen)

        self.update()
        self.after(timeout, self.destroySplash)

    def centerOnScreen(self):
        self.update_idletasks()
        self.width, self.height = self.image.width(), self.image.height()

        xmax = self.winfo_screenwidth()
        ymax = self.winfo_screenheight()

        x0 = self.x0 = xmax/2 - self.width/2
        y0 = self.y0 = ymax/2 - self.height/2
        self.geometry("%dx%d+%d+%d" % (self.width, self.height, x0, y0))
        self.createSplash()

    def createSplash(self):
        # show the splash image
        self.canvas = Canvas(self, height=self.height, width=self.width)
        self.canvas.create_image(0,0, anchor='nw', image=self.image)
        self.canvas.pack()
        self.canvas.tag_raise('firstRect')

    def destroySplash(self):
        # bring back main window and destroy splash screen
        self.withdraw()
        self.destroy()
