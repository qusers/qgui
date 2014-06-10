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

from Tkinter import Text, Label, Button, Frame, X, CENTER, Scrollbar, YES, Menu, DISABLED, NORMAL, GROOVE, END, TOP, PhotoImage, Toplevel

class AboutQ(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.root = root

        self.dialog_window()

    def dialog_window(self):
        """
        Window
        """
        self.title('About Qgui')

        frame = Frame(self, bg='black')
        frame.pack(padx=(10, 10), pady=(10, 10))

        frame2 = Frame(self, bg='black')
        frame2.pack(padx=(10, 10), pady=(10, 10))

        photo = PhotoImage(file=self.app.settings_path + "/Qmods/qgui_box.gif")
        w = Label(frame, image=photo, bg='black', highlightthickness=2)
        w.photo = photo
        w.grid(row=0, column=5, sticky = 'nswe')

        txt = Text(frame2, width = 80, height=20, bg='WhiteSmoke')
        txt.grid(row=1, column=0, columnspan=10)
        txt.insert(0.0, """Qgui source file is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Qgui is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Qgui.  If not, see <http://www.gnu.org/licenses/>.

QGui is still under development and testing!
Please contact geir.isaksen@uit.no to report bugs or errors.""")
        txt.config(highlightthickness=2, relief=GROOVE, highlightbackground='lightblue', state=DISABLED)

        close = Button(frame2, text='OK', command=self.destroy)
        close.config(highlightbackground='black')
        close.grid(row=3,column = 9)



