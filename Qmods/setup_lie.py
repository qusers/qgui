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

from Tkinter import X, Label, Button, Spinbox, BOTTOM, Entry, LabelFrame, Frame, Toplevel, DISABLED, \
    NORMAL, END, GROOVE

from select_atoms import AtomSelectRange
import os
from tkFileDialog import askopenfilename
from setup_md import SetupMd
from topoprep import TopologyPrepare
from ttk import Progressbar
import random
from subprocess import call


class SetupLie(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.progress = 0

        self.workdir = self.app.workdir
        self.complex_pdb = None
        self.complex_top = None
        self.complex_md = False
        self.ligand_pdb = None
        self.ligand_top = None
        self.ligand_md = False

        self.md_ligand = False
        self.md_complex = False

        self.qgui_path = self.app.qgui_path

        self.dialog_window()
        self.app.log('info', 'Setup LIE session started.')

    def update_progress(self):
        """
        Updates progressbar, buttons and entries
        """
        self.progress = 0
        if self.complex_pdb:
            self.complex_entry.delete(0, END)
            self.complex_entry.insert(0, self.complex_pdb.split('/')[-1])
            self.progress += 12.5
            self.prep_complex_top.config(state=NORMAL)
            self.select_ligand.config(state=NORMAL)
        if self.complex_top:
            self.complex_top_entry.delete(0, END)
            self.complex_top_entry.insert(0, self.complex_top.split('/')[-1])
            self.progress += 12.5
        if self.complex_pdb and self.complex_top:
            self.setup_complex_md.config(state=NORMAL)
        if self.complex_md:
            self.progress += 25
        if self.ligand_pdb:
            self.ligand_entry.delete(0, END)
            self.ligand_entry.insert(0, self.ligand_pdb.split('/')[-1])
            self.progress += 12.5
            self.prep_ligand_top.config(state=NORMAL)
        if self.ligand_top:
            self.ligand_top_entry.delete(0, END)
            self.ligand_top_entry.insert(0, self.ligand_top.split('/')[-1])
            self.progress += 12.5
        if self.ligand_pdb and self.ligand_top:
            self.setup_ligand_md.config(state=NORMAL)
        if self.ligand_md:
            self.progress += 25

        if self.md_complex:
            self.progress += 25
        if self.md_ligand:
            self.progress += 25

        if self.progress == 100:
            self.write_button.config(state=NORMAL)
            self.run_lie_button.config(state=NORMAL)

        self.progressbar.config(value=self.progress)


    def open_pdb(self, entry=None):
        filename = askopenfilename(parent=self, initialdir = self.app.workdir, filetypes=(("pdb", "*.pdb"),("All files","*.*")))
        if filename != '':
            entry.delete(0,END)
            entry.insert(0,filename.split('/')[-1])
            if entry == self.complex_entry:
                self.complex_pdb = filename
            else:
                self.ligand_pdb = filename

        self.update_progress()

    def open_top(self, entry=None):
        filename = askopenfilename(parent=self, initialdir=self.app.workdir, filetypes=(("TOP", "*.top"), ("All files", '*.*')))
        if filename != '':
            entry.delete(0,END)
            entry.insert(0, filename.split('/')[-1])
            if entry == self.complex_top_entry:
                self.complex_top = filename
            else:
                self.ligand_top = filename

        self.update_progress()

    def setup_md(self, lie_type='complex'):
        """
        Opens the Setup MD for Qdyn
        """
        pdb = None
        self.md_written = False

        if lie_type == 'complex':
            top = self.complex_top
            self.add_title = 'complex'
            dirname = self.complex_dir.get().strip()
            if self.complex_pdb:
                pdb = self.complex_pdb

        else:
            top = self.ligand_top
            self.add_title = 'ligand'
            dirname = self.ligand_dir.get().strip()
            if self.ligand_top:
                pdb = self.ligand_pdb

        self.workdir = self.app.workdir + '/' + dirname + '/inputfiles'
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        os.system('cp %s %s' % (top, self.workdir))
        os.system('cp %s %s' % (pdb, self.workdir))

        self.setup_md_ = SetupMd(self,self.root, pdb, top, False, fep=False)
        self.setup_md_.configure(background = self.main_color)
        self.setup_md_.title('Configure MD for %s (LIE)' % self.add_title)
        self.setup_md_.resizable()

    def prepare_topology(self, lie_type='complex'):
        """
        Opens the Prepare topology window
        """
        pdb = None
        if lie_type == 'complex':
            self.add_title = 'complex'
            pdb = self.complex_pdb
        else:
            self.add_title = 'ligand'
            pdb = self.ligand_pdb

        prms = self.app.q_settings[1]
        libs = self.app.q_settings[2]

        self.topo_prepare = TopologyPrepare(self, self.root, pdb, prms, libs, False)
        self.topo_prepare.configure(background = self.main_color)
        self.topo_prepare.title('Prepare %s topology' % self.add_title)
        self.topo_prepare.resizable()

    def select_atoms(self):
        """
        opens dialog to select atoms and inserts first into entry 1 and last into entry 2
        """
        entry1 = self.ligand_entry
        entry2 = entry1

        self.select_atomrange = AtomSelectRange(self,self.root, self.complex_pdb, entry1, entry2, 'EXTENDED', True)
        self.select_atomrange.configure(bg = self.main_color)
        self.select_atomrange.title('Select atoms')
        self.select_atomrange.resizable()

    def write_lie(self, message=True):
        """
        Creates n directories (1-n; n = runs) in ligand and complex
        Writes run script to submit ligand and complex md!
        """

        submitfile = open(self.app.workdir + '/runLIE.sh', 'w')
        submitfile.write('#! /bin/bash\n\n')

        ligand_basename = self.ligand_entry.get().split('.')[0]
        complex_basename = self.complex_entry.get().split('.')[0]
        runs = int(self.runs.get())
        complex_dir = self.app.workdir + '/' + self.complex_dir.get().strip()
        ligand_dir = self.app.workdir + '/' + self.ligand_dir.get().strip()
        print complex_dir
        print ligand_dir
        for dir_ in [complex_dir, ligand_dir]:
            for i in range(1, runs + 1):
                if not os.path.exists('%s/%d' % (dir_, i)):
                    os.makedirs('%s/%d' % (dir_, i))
                os.system('cp %s/inputfiles/* %s/%d/' % (dir_, dir_, i))
                if dir_ == complex_dir:
                    basename = complex_basename
                elif dir_ == ligand_dir:
                    basename = ligand_basename
                try:
                    eq1 = open('%s/%d/%s_eq1.inp' % (dir_, i, basename), 'r').readlines()
                    eq1_new = open('%s/%d/%s_eq1.inp' % (dir_, i, basename), 'w')
                    for line in eq1:
                        if 'random_seed' in line:
                            newseed = random.randrange(1000, 9999)
                            eq1_new.write('%25s %d\n' % ('random_seed'.ljust(25), newseed))
                        else:
                            eq1_new.write(line)
                    eq1_new.close()
                except:
                    print 'Equilibration file with random seed not found!'
                submitfile.write('cd %s/%d\n' % (dir_, i))
                submitfile.write('%s %srun.sh\n' % (self.app.q_settings[4][1], basename))
                submitfile.write('cd ../../\n')
        submitfile.close()

        if message:
            self.app.errorBox('Info','Input files written. Use runLie.sh to submit jobs.')

    def run_lie(self):
        """
        Writes inputfiles, runscript and submission script
        and submits them
        """
        self.write_lie(False)
        os.chdir(self.app.workdir)
        tmpfile = open('.tmpfile', 'w')
        #os.system('bash runLIE.sh')
        call('bash runLIE.sh', shell=True, stdout=tmpfile, stderr=tmpfile)
        job_id = open(self.app.workdir + '/.tmpfile','r').readlines()
        self.app.log('info','Submitting LIE jobs ...')
        for line in job_id:
            self.app.main_window.update_txt(line)

        self.app.errorBox('Info','Jobs submitted!')



    def dialog_window(self):
        """
        Window
        """
        self.title('Setup LIE')

        #Complex
        frame = Frame(self, bg=self.main_color)
        frame.pack(padx=(10, 10), pady=(10, 10), fill=X)
        #Ligand
        frame2 = Frame(self, bg=self.main_color)
        frame2.pack(padx=(10, 10), pady=(10, 10), fill=X)
        #write/run
        frame3 = Frame(self, bg=self.main_color)
        frame3.pack(padx=(10, 10), pady=(10, 10), fill=X)
        #progressbar
        frame4 = Frame(self, bg=self.main_color)
        frame4.pack(fill=X)



        #COMPLEX ENTRY
        complex_label = LabelFrame(frame, text='COMPLEX', padx=10, bg=self.main_color)
        complex_label.grid()

        complex_structure = Label(frame, text='Structure: ', bg=self.main_color)
        complex_structure.grid(in_=complex_label, row=0, column=0)

        self.complex_entry = Entry(frame, width=20, highlightthickness=0)
        self.complex_entry.grid(in_=complex_label, row=0, column=1)

        load_complex = Button(frame, text='Load', highlightbackground=self.main_color, command=lambda: self.open_pdb
            (self.complex_entry))
        load_complex.grid(in_=complex_label, row=0, column=2)

        complex_top = Label(frame, text='Topology: ', bg=self.main_color)
        complex_top.grid(in_=complex_label, row=1,column=0)

        self.complex_top_entry = Entry(frame, width=20, highlightthickness=0)
        self.complex_top_entry.grid(in_=complex_label, row=1, column=1)

        load_complex_top = Button(frame, text='Load', highlightbackground=self.main_color, command=lambda: self.open_top(
            self.complex_top_entry))
        load_complex_top.grid(in_=complex_label, row=1, column=2)

        self.prep_complex_top = Button(frame, text='Prepare', highlightbackground=self.main_color,
                                       command=lambda: self.prepare_topology('complex'))
        self.prep_complex_top.grid(in_=complex_label, row=1, column=3)
        self.prep_complex_top.config(state=DISABLED)

        complex_dir = Label(frame, text='Dir. name:', bg=self.main_color)
        complex_dir.grid(in_=complex_label, row=2, column=0)

        self.complex_dir = Entry(frame, width=20, highlightthickness=0)
        self.complex_dir.grid(in_=complex_label, row=2, column=1)
        self.complex_dir.insert(0, 'complex')

        self.setup_complex_md = Button(frame, text='Configure MD', highlightbackground=self.main_color,
                                       command=lambda: self.setup_md('complex'))
        self.setup_complex_md.grid(in_=complex_label, row=3, column=0, columnspan=4)
        self.setup_complex_md.config(state=DISABLED)

        #LIGAND ENTRY
        ligand_label = LabelFrame(frame2, text='LIGAND', padx=10, bg=self.main_color)
        ligand_label.grid()

        ligand_structure = Label(frame2, text='Structure: ', bg=self.main_color)
        ligand_structure.grid(in_=ligand_label, row=0, column=0)

        self.ligand_entry = Entry(frame2, width=20, highlightthickness=0)
        self.ligand_entry.grid(in_=ligand_label, row=0, column=1)

        load_ligand = Button(frame2, text='Load', highlightbackground=self.main_color, command=lambda: self.open_pdb(self.ligand_entry))
        load_ligand.grid(in_=ligand_label, row=0, column=2)

        self.select_ligand = Button(frame2, text='Select', highlightbackground=self.main_color, width=6, command=self.select_atoms)
        self.select_ligand.grid(in_=ligand_label, row=0, column=3)
        self.select_ligand.config(state=DISABLED)

        ligand_top = Label(frame2, text='Topology: ', bg=self.main_color)
        ligand_top.grid(in_=ligand_label, row=1,column=0)

        self.ligand_top_entry = Entry(frame2, width=20, highlightthickness=0)
        self.ligand_top_entry.grid(in_=ligand_label, row=1, column=1)

        load_ligand_top = Button(frame2, text='Load', highlightbackground=self.main_color, command=lambda: self.open_top(
            self.ligand_top_entry))
        load_ligand_top.grid(in_=ligand_label, row=1, column=2)

        self.prep_ligand_top = Button(frame2, text='Prepare', highlightbackground=self.main_color,
                                      command=lambda: self.prepare_topology('ligand'))
        self.prep_ligand_top.grid(in_=ligand_label, row=1, column=3)
        self.prep_ligand_top.config(state=DISABLED)

        ligand_dir = Label(frame2, text='Dir. name:', bg=self.main_color)
        ligand_dir.grid(in_=ligand_label, row=2, column=0)

        self.ligand_dir = Entry(frame2, width=20, highlightthickness=0)
        self.ligand_dir.grid(in_=ligand_label, row=2, column=1)
        self.ligand_dir.insert(0, 'ligand')

        self.setup_ligand_md = Button(frame2, text='Configure MD', highlightbackground=self.main_color,
                                      command=lambda: self.setup_md('ligand'))
        self.setup_ligand_md.grid(in_=ligand_label, row=3, column=0, columnspan=4)
        self.setup_ligand_md.config(state=DISABLED)

        #Write/run
        runs = Label(frame3, text='Runs:', bg=self.main_color)
        runs.grid(row=0, column=0)

        self.runs = Spinbox(frame3, width=2, from_=1, to=99, highlightthickness = 0, relief = GROOVE)
        self.runs.grid(row=0, column=1, padx=(5,0))

        self.write_button = Button(frame3, text='Write', highlightbackground=self.main_color, command=self.write_lie)
        self.write_button.grid(row=0, column=2, sticky='e', padx=(60, 0))
        self.write_button.config(state=DISABLED)

        self.run_lie_button = Button(frame3, text='Run', highlightbackground=self.main_color, command=self.run_lie)
        self.run_lie_button.grid(row=0, column=3, sticky='e')
        self.run_lie_button.config(state=DISABLED)

        self.close_button = Button(frame3, text='Close', highlightbackground=self.main_color, command=self.destroy)
        self.close_button.grid(row=0, column=4, padx=(85, 0))

        #progressbar
        self.progressbar = Progressbar(frame4, mode='determinate')
        self.progressbar.pack(side=BOTTOM, expand=True, fill=X)
        self.progressbar.config(value=self.progress)
