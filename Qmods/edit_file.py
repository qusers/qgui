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

from tkinter import Button, Text, Frame, Toplevel, Scrollbar, END, LEFT


class FileEdit(Toplevel):
    """Opens a file to that can be edited"""

    def __init__(self, root, file_to_read, main_color='Gray'):
        Toplevel.__init__(self, root)
        self.root = root

        self.main_color = main_color

        self.filename = file_to_read

        self.dialog_box()
        self.update_txt()

    def update_txt(self):
        """
        Reads text from file to window
        """
        self.txt.insert(END, open(self.filename).read())

    def save_file(self):
        """
        Saves changes to file
        """
        savefile = open(self.filename,'w')
        #Avoid the additional blank line at the end
        #that you get with savefile.write(self.txt.get(1.0, END))
        output = self.txt.get(1.0, END)
        for line in range(len(output) - 1):
            savefile.write(output[line])
        savefile.close()
        self.destroy()

    def dialog_box(self):
        """
        Layout
        """
        self.title('Edit file')
        self.config(background=self.main_color)

        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=10, padx=(10,10))

        self.txt = Text(left_frame, width =80, height=30)
        self.txt.grid(row=0, column=0, columnspan=10)

        scrollb = Scrollbar(left_frame, command = self.txt.yview)
        scrollb.grid(row = 0, column = 10, sticky = 'nsew')
        self.txt['yscrollcommand'] = scrollb.set

        cancel_button = Button(left_frame, text = 'Cancel', command = self.destroy)
        cancel_button.grid(row=1, column=9)
        cancel_button.config(highlightbackground=self.main_color)

        save_button = Button(left_frame, text = 'Save', command=self.save_file)
        save_button.grid(row=1, column=8)
        save_button.config(highlightbackground=self.main_color)
