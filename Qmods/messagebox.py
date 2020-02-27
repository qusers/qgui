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

# -*- coding: utf-8 -*-
from tkinter.ttk import Frame

import tkinter.messagebox

from tkinter import Button


class MessageBox(Frame):
    def __init__(self, root, msg): #Receives app and root from Qgui class.
        Frame.__init__(self, root) # app = self from qgui and root = self.root from qgui
        self.root = root
        self.message = msg

        self.message_box(self.message)

    def message_box(self, message):

        msg_box = tkinter.messagebox.showinfo('Error', self.message)


        button = Button(self, text = 'OK')

        self.pack()
