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

from Tkinter import X, Label, Button, Spinbox, Entry, Radiobutton, LabelFrame, Text, Frame, Toplevel, DISABLED, \
    NORMAL, END, GROOVE, IntVar, StringVar, Checkbutton
import cPickle
from select_atoms import AtomSelectRange
from subprocess import call
from tkFileDialog import askopenfilename
import os
import time
import numpy as np


class CreatePrmLib(Toplevel):
    def __init__(self, app, root):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.q_settings = cPickle.load(open(self.app.settings_path + '/' + 'Qsettings','rb'))
        self.ffld_path = self.app.q_settings[ 'schrodinger path' ]
        if not self.ffld_path:
            self.app.errorBox('Error', 'Schrodinger path not in settings.')
            self.destroy()
        else:
            print self.ffld_path

        self.atomi_changed = StringVar()
        self.atomj_changed = StringVar()
        self.prm_id_var = StringVar()

        #Replace original pdb entry
        self.check_replace = IntVar()
        self.check_replace.set(1)

        self.dialog_window()
        self.prm_id_var.set('T')

        self.atomi_changed.trace('w', self.update_atomnr)
        self.atomj_changed.trace('w', self.update_atomnr)



        self.resname.insert(0,'LIG')
        self.resnr.insert(0, 1)
        self.atomnr.insert(0, 1)

        #Mass dictionary used in functions:
        self.massDict = {'H': 1.008, 'C': 12.01, 'N': 14.01,'O': 16.00,'F': 19.00,'Mg': 24.305, 'P': 30.97, 'S': 32.07,
                         'Cl': 35.45, 'Ca': 40.078, 'Fe':55.847, 'Zn': 65.38, 'Br':79.90,'I':126.90}

        #look in mae file for atom weight and return atom:
        self.atoms = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9:'F',12: 'Mg', 15:'P',16:'S',17:'Cl', 20: 'Ca', 26: 'Fe', 30: 'Zn', 35:'Br', 53:'I'}

    def select_atoms(self, entry1, entry2):
        """
        opens dialog to select atoms and inserts first into entry 1 and last into entry 2
        """

        self.select_atomrange = AtomSelectRange(self,self.root, self.app.pdb_id, entry1, entry2)
        self.select_atomrange.configure(bg = self.main_color)
        self.select_atomrange.title('Select atoms')
        self.select_atomrange.resizable()

    def update_atomnr(self,*args):
        """
        If pdb atoms are selected, automatically suggest
        resname, resnr and atomnumber
        """
        try:
            first_atom = int(self.atom_i.get())
            last_atom = int(self.atom_j.get())
        except:
            return
        if first_atom > last_atom:
            return

        resname = 'LIG'
        resnr = 1
        pdbfile = open(self.app.pdb_id).readlines()

        done = False
        line = 0
        while not done:
            if 'ATOM' in pdbfile[line] or 'HETATM' in pdbfile[line]:
                if first_atom == int(pdbfile[line].split()[1]):
                    resname = pdbfile[line].split()[3]
                    resnr = pdbfile[line].split()[4]
                    done = True
                else:
                    line += 1
            else:
                line += 1
            if line >= 50000:
                self.app.errorBox('Error','Could not find atoms in pdb file!')
                return

        self.resname.delete(0, END)
        self.resname.insert(0, resname)
        self.resnr.delete(0, END)
        self.resnr.insert(0, resnr)
        self.atomnr.delete(0, END)
        self.atomnr.insert(0, first_atom)

    def rename_atoms(self, pdb_dict):
        """
        Creates a single residue from multiple residues with new atomnames
        """
        atom_count = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0, 'Cl': 0, 'F': 0, 'Br': 0}

        resnr = self.resnr.get().strip()
        atomnr = int(self.atomnr.get())

        for i in sorted(pdb_dict.keys()):
            j = 13
            found_atom = False
            while not found_atom:
                atom = pdb_dict[i][j]
                if atom not in atom_count.keys():
                    if j > 16:
                        self.app.errorBox('Error', 'Could not translate all atomtypes. Pleas verify pdb file.')
                        pdb_dict = None
                        break
                    j += 1
                else:
                    if pdb_dict[i][j + 1].lower() == 'l':
                        if atom == 'C':
                            atom = 'Cl'
                    elif atom == 'B':
                        if pdb_dict[i][j + 1].lower() == 'r':
                            atom = 'Br'
                    found_atom = True
            if not pdb_dict:
                break

            atom_count[atom] += 1
            pdb_dict[i] = '%5s%6d  %1s%3s%4s%5s%s' % (pdb_dict[i][0:5], atomnr, atom, str(atom_count[atom]).ljust(3),
                                                  self.ligname.ljust(4), resnr, pdb_dict[i][26:])
            atomnr += 1
        return pdb_dict

    def create_pdb(self):
        """
        Creat a new pdb file from selection
        """

        try:
            first_atom = int(self.atom_i.get())
            last_atom = int(self.atom_j.get())
        except:
            return
        if first_atom > last_atom:
            return


        self.lig_atoms = dict()
        self.lig_atoms_renamed = dict()

        resnr = self.resnr.get()
        atomnr = int(self.atomnr.get())


        #create list of selected atom numbers:
        atom_numbers = []
        for i in range(first_atom, last_atom + 1):
            atom_numbers.append(i)

        pdbfile = open(self.app.pdb_id, 'r').readlines()
        self.ligname = self.resname.get().strip()


        #If more than 1 residue, the entire molecule must be redefined with single residue name and new atomnames!
        #If pdb atomname occures more than one time in singular residue, rename as well!
        residues = []
        pdbnames = []
        pdb_duplicate = False
        #Get selection from pdb file:
        count = 0
        for line in pdbfile:
            if 'ATOM' in line or 'HETATM' in line:
                if int(line.split()[1]) in atom_numbers:
                    print line
                    count += 1
                    atom = line[13]
                    if line[14].lower() == 'l':
                        if atom == 'C':
                            atom = 'Cl'
                    elif atom == 'B':
                        if line[14].lower() == 'r':
                            atom = 'Br'
                    pdbname = line[13:17].strip()
                    self.lig_atoms[count] = line.split('\n')[0].rstrip() + '  0.00  0.00          %2s\n' % atom
                    if line[22:26].strip() not in residues:
                        residues.append(line[22:26].strip())
                    if pdbname not in pdbnames:
                        pdbnames.append(pdbname)
                    elif pdbname in pdbnames:
                        pdbnames.append(pdbname)
                        #Check how many times name occures and compare to number of residues
                        for name in pdbnames:
                            if name == pdbname:
                                count += 1
                        if len(residues) < count:
                            pdb_duplicate = True


        if len(residues) > 1:
            self.app.errorBox('Info','More than 1 residue in selection. A new single residue will be generated')
            self.lig_atoms_renamed = self.rename_atoms(self.lig_atoms)

        else:
            if not pdb_duplicate:
                for i in sorted(self.lig_atoms.keys()):
                    self.lig_atoms_renamed[i] = '%5s%6d  %s%4s%5s%s' % (self.lig_atoms[i][0:5], atomnr, self.lig_atoms[i][13:17],
                                                          self.ligname.ljust(4), resnr, self.lig_atoms[i][26:])
                    atomnr += 1
            elif pdb_duplicate:
                self.app.errorBox('Info','Duplicate atomnames in residue! Renaming all atoms.')
                self.lig_atoms_renamed = self.rename_atoms(self.lig_atoms)
        #If replace original entry in pdbfile:
        if self.check_replace.get() == 1:
            res, current_res = 1, 1
            found_lig = False

            org_file = open(self.app.pdb_id, 'r').readlines()
            self.app.pdb_id = self.app.pdb_id.split('.pdb')[0] + '_new.pdb'
            self.app.main_window.set_entryfield(self.app.pdb_id.split('/')[-1])
            new_file = open(self.app.pdb_id, 'w')
            for line in org_file:

                if 'GAP' in line:
                        new_file.write(line)
                else:
                    if int(line.split()[1]) in atom_numbers:
                        #Increas res nr by 1 the first time new res is inserted:
                        if not found_lig:
                            res += 1
                            found_lig = True

                        #Find key to replace:
                        for i in sorted(self.lig_atoms.keys()):
                            print '%s %s' % (self.lig_atoms[i].split()[1], line.split()[1])
                            if int(self.lig_atoms[i].split()[1]) == int(line.split()[1]):
                                newline = '%s\n' % (self.lig_atoms_renamed[i][0:55])
                                new_file.write(newline)
                                break
                    else:
                        #Make sure that res nr are in order
                        if current_res != int(line[21:26]):
                            res += 1
                            current_res = int(line[21:26])
                        newline = '%s%5d%s' % (line[0:21], res, line[26:])
                        new_file.write(newline)

            new_file.close()

        #Write pdb for selected pdb entries
        ligfile = open(self.app.workdir + '/' + '%s.pdb' % self.ligname, 'w')
        for i in sorted(self.lig_atoms_renamed.keys()):
            ligfile.write(self.lig_atoms_renamed[i])

        ligfile.close()

    def getParam(self):
        """
        Uses $SCHRODINGER/utilities/ffld_server -ipdb -print_parameters -version 11 or version 14
        """
        self.impact_failed = False

        self.opls = self.forcefield.get().strip()
        if self.opls == 'OPLS2001':
            version = '11'
        #OPLS-2011 and up requires special license of course, so...
        elif self.opls == 'OPLS2011':
            version = '16'
        else:
            version = '14'

        os.chdir(self.app.workdir)

        tmpfile = open('.tmpfile','wb')

        input_type = 'pdb'
        #Check if maestro file exist
        maefile = self.ligname.split('.pdb')[0]+'.mae'
        if os.path.isfile(self.app.workdir+'/'+maefile):
            print 'Found maestro file, using this by default'
            input_type = 'mae'

        call('%s/utilities/ffld_server -i%s %s/%s.%s -print_parameters -version %s' %
             (self.ffld_path, input_type, self.app.workdir, self.ligname, input_type, version),
             shell=True, stdout=tmpfile, stderr=tmpfile)

        self.update()

        done = False

        self.app.log('info','Generating parameters ...')
        job_id = open(self.app.workdir + '/.tmpfile','r').readlines()
        for line in job_id:
            self.app.main_window.update_txt(line)
            if 'OPLSAA FORCE FIELD TYPE ASSIGNED' in line:
                done = True
            if 'exception' in line:
                done = True
                self.impact_failed = True
            if 'FATAL' in line:
                self.impact_failed = True
                done = True

        if not done:
            #Wait for ffld_server to finish if not done:
            self.app.main_window.txt.config(state=NORMAL)
            self.app.main_window.txt.insert(END, '\nWaiting for parameter assignment to finish\n')
            logline = len(self.app.main_window.txt.get(0.0, END).split('\n')) - 2
            logline = str(logline) + '.0'

        dotcount = 1
        count = 0
        while not done:
            job_id = open(self.app.workdir + '/.tmpfile','r').readlines()
            for line in job_id:
                if 'OPLSAA FORCE FIELD TYPE ASSIGNED' in line:
                    done = True
            #Update monitor so that user can see that something is happening:
            self.app.main_window.txt.delete(logline, END)
            self.app.main_window.txt.insert(END, '\nWaiting for parameter assignment to finish' + dotcount * '.' + '\n')
            self.update()
            time.sleep(0.2)
            dotcount += 1
            count += 1
            if dotcount >= 40:
                dotcount = 1
            if count > 100:
                self.impact_failed = True
                break

    def impQ(self,a1,a2,a3,a4):
        """
        This function will make sure that impropers from
        Impact are reordered in the correct manner for Q.
        """
        atoms = [a1, a2, a3, a4]
        bonds = []
        bondList = open(self.app.workdir + '/.tmpfile','r').readlines()
        for atom in range(len(atoms)):
            count = 0
            for line in range(self.bndStart,self.bndEnd + 1):
                if atoms[atom] in bondList[line]:
                    count += 1
            bonds.append(count)

        if 1 in bonds:
            firstAtom = atoms[bonds.index(1)]
            first = False
            while first == False:
                if firstAtom == atoms[0]:
                    first = True
                else:
                    atoms.insert(0,atoms.pop())
            q1 = atoms[0]
            q2 = atoms[1]
            q3 = atoms[2]
            q4 = atoms[3]
        elif 1 not in bonds:
            q1 = atoms[1]
            q2 = atoms[2]
            q3 = atoms[3]
            q4 = atoms[0]

        atoms = [q1,q2,q3,q4]

        for i in (0,1):
            for j in range(self.bndStart,self.bndEnd + 1):
                if atoms[i] in bondList[j]:
                    done = False
                    popped = 0
                    while done == False :
                        if (atoms[i] in bondList[line] and atoms[(i + 1)] in bondList[line]):
                            done = True
                        else:
                            atoms.insert((i + 1),atoms.pop())
                            popped += 1
                        if popped > 2:
                            line += 1
                            popped = 0

        return atoms[0], atoms[1], atoms[2], atoms[3]

    def ligprmDict(self,impactLog='%s/impact_param.log', pdb = 'lig.pdb'):
        """
        Creates a dictionary for keeping track of atom names in impact log file
        and newly created atom names in lig.pdb. Need this for tracking parameters
        to correct atoms in lig.lib. Uses impact_out.mae for tracking!

        Returns ligDict, prmDict, atomnrDict
        """

        log = open(impactLog,'r').readlines()   #Defining list start and end points/ track atoms
        lig = open(pdb,'r').readlines()         #tracking: lig <>outMae/LOG

        self.ligDict = dict()        #Dictionary to keep track of new and original atom names.
        self.prmDict = dict()

        self.nbndStart = None
        self.angStart = None
        self.angEnd = None
        self.torStart = None
        self.torEnd = None
        self.impStart = None

        count = 0
        lig_line = 0
        found_atoms = False
        for i in range(len(log)):
            if found_atoms:
                if len(log[i].split()) > 6:
                    orgName = lig[lig_line][13:17].strip()
                    newName = log[i].split()[0]
                    prmName = log[i].split()[3]           #Create unique parameter names for .prm
                    self.ligDict[newName] = orgName
                    self.prmDict[newName] = prmName
                    lig_line += 1

                else:
                    found_atoms = False

            if '-----------------' in log[i]:
                count += 1

            if count == 2:
                found_atoms = True
            if 'OPLSAA FORCE FIELD TYPE ASSIGNED' in log[i]:
                self.nbndStart = i + 3        #Non-bonded terms start here
            if count == 3:
                self.nbndEnd = i - 1
                count += 1
            if 'Stretch' in log[i]:
                self.bndStart = i+1         #Bonded terms start here
            if 'Bending' in log[i]:
                self.angStart = i+1         #Angle terms start here
                self.bndEnd = i - 2
            try:
                if 'proper' == log[i].split()[0]:
                    self.torStart = i+1         #Torsions/imps start here ('==>' in tors)
                    self.angEnd = i - 2
                if 'improper' == log[i].split()[0]:
                    self.impStart = i+1           #Torsions/Imps end here
                    self.torEnd = i - 2
            except:
                continue
            if self.torStart:
                if len(log[i].split()) < 5:
                    self.torEnd = i
        if not self.torEnd:
            self.torEnd = len(log) - 1

        return self.ligDict, self.prmDict

    def makeLib(self, impactLog = 'impact_param.log'):
        """
        Returns lig.lib
        """

        ligDict, prmDict = self.ligprmDict(self.impact_log, self.lig_pdb)    #Get references for tracing new to orig. atom names.

        self.log = open(impactLog,'r').readlines()

        res = self.resname.get().strip()

        libname = self.app.workdir + '/' + self.ligname +'.lib'
        self.liglib = open(libname,'w')

        atoms = 0
        chargesum = 0
        atomnames = []
        prmnames = []
        charges = []

        print self.nbndStart
        for i in range(self.nbndStart + 1, self.nbndEnd + 1):
            atoms += 1
            atomname = ligDict[self.log[i].split()[0]]
            atomnames.append(atomname)
            prmname = self.prm_id.get() + prmDict[self.log[i].split()[0]]
            prmnames.append(prmname)
            charge = self.log[i].split()[4]
            charges.append(charge)
            chargesum += float(charge)

        self.liglib.write('*----------------------------------------------------------------------\n')
        self.liglib.write('{%3s}                 ! %s: atoms %d charge %.4f\n\n' % (res, self.ligname, atoms, chargesum))
        self.liglib.write('[info]\n')
        self.liglib.write('SYBYLtype RESIDUE\n\n')
        self.liglib.write('[atoms]\n')

        for i in range(len(atomnames)):
            self.liglib.write('%3d      %4s    %4s     %7.4f\n' %\
                              (i+1,atomnames[i].ljust(4),prmnames[i].ljust(4),float(charges[i])))

        self.liglib.write('\n[bonds]\n')

        #Make atom bond list based on stretches from impact:
        for i in range(self.bndStart,self.bndEnd + 1):
            atom1 = ligDict[self.log[i].split()[0]]
            atom2 = ligDict[self.log[i].split()[1]]
            self.liglib.write('%4s  %4s\n' %(atom1.ljust(4), atom2.ljust(4)))

        if self.impStart:
            self.liglib.write('\n[impropers]\n')

            #Make improper list from impact log file:
            i = self.impStart
            done = False
            while not done:
                try:
                    if len(self.log[i].split()) > 5:
                        atom1 = self.log[i].split()[0]
                        atom2 = self.log[i].split()[1]
                        atom3 = self.log[i].split()[2]
                        atom4 = self.log[i].split()[3]
                        qimp1,qimp2,qimp3,qimp4 = self.impQ(atom1,atom2,atom3,atom4)
                        qimp1 = ligDict[qimp1]
                        qimp2 = ligDict[qimp2]
                        qimp3 = ligDict[qimp3]
                        qimp4 = ligDict[qimp4]

                        self.liglib.write('%4s  %4s  %4s  %4s\n'\
                                      % (qimp1.ljust(4), qimp2.ljust(4),qimp3.ljust(4), qimp4.ljust(4)))
                        i += 1
                    else:
                        done = True
                except:
                    done = True

        self.liglib.write('\n')
        self.liglib.close()

        print '%s written' % libname
        self.app.log('info','Q library %s successfully generated.' % libname.split('/')[-1])

    def calcNBpar(self, sigma = 0, epsilon = 0):
        """
        >> npbar(sigma, epsilon)

        Returns: avdw1,avdw2,avdw3,bvdw1,bvdw2
        """
        eps = float(epsilon)
        sig = float(sigma)

        avdw1 = np.sqrt(4*eps*sig**(12))
        avdw2 = avdw1
        avdw3 = np.sqrt(0.5) * avdw1

        bvdw1 = np.sqrt(4 * eps * sig**(6) )
        bvdw23 = np.sqrt(0.5) * bvdw1

        return (avdw1,avdw2,avdw3,bvdw1,bvdw23)

    def makeVdw(self, impactLog = 'impact_param.log'):
        """
        Writes  lig_vdw.prm
        """
        self.log = open(impactLog,'r').readlines()
        vdwFile = open('%s/%s_vdw.prm' % (self.app.workdir, self.ligname), 'w')

        vdwFile.write('! Ligand vdW parameters for %s\n' % self.ligname)

        #Remove redundancies
        prm_added = []

        prmDict = self.ligprmDict(self.impact_log, self.lig_pdb)[1]
        for i in range(self.nbndStart + 1, self.nbndEnd + 1):
            prmname = self.prm_id.get() + prmDict[self.log[i].split()[0]]
            if prmname not in prm_added:
                prm_added.append(prmname)
                sigma = float(self.log[i].split()[5])
                epsilon = float(self.log[i].split()[6])
                avdw1,avdw2,avdw3,bvdw1,bvdw23 = self.calcNBpar(sigma,epsilon)
                mass = self.massDict[self.log[i].split()[0][0]]
                comment = str(self.log[i][52:])
                vdwFile.write('%4s  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f !%s'\
                          % (prmname.ljust(4), avdw1, avdw2, bvdw1, avdw3, bvdw23, mass, comment))

        vdwFile.write('\n')
        vdwFile.close()

        self.app.log('info', '%s_vdw.prm successfully generated' % self.ligname)

    def makeBnd(self,impactLog = 'impact_param.log'):
        """
        Writes  lig_bnd.prm
        """
        self.log = open(impactLog,'r').readlines()
        bndFile = open('%s/%s_bnd.prm' % (self.app.workdir, self.ligname),'w')

        bndFile.write('! Ligand bond parameters for %s\n' % self.ligname)

        prmDict = self.ligprmDict(self.impact_log, self.lig_pdb)[1]
        bnd_added = []

        for i in range(self.bndStart, self.bndEnd + 1):
            atom1 = self.prm_id.get() + prmDict[self.log[i].split()[0]]
            atom2 = self.prm_id.get() + prmDict[self.log[i].split()[1]]
            if [atom1, atom2] not in bnd_added and [atom2, atom1] not in bnd_added:
                bnd_added.append([atom1, atom2])
                bnd_added.append([atom2, atom1])
                K = 2 * float(self.log[i].split()[2])
                Req =  float(self.log[i].split()[3])
                comment = str(self.log[i][39:].strip())
                if '\n' not in comment:
                    comment += '\n'
                bndFile.write('%4s  %4s  %7.2f  %7.2f !%s'\
                          % (atom1, atom2, K, Req, comment))

        bndFile.write('\n')
        bndFile.close()
        print '%s_bnd.prm written' % self.ligname
        self.app.log('info', '%s_bnd.prm successfully generated.' % self.ligname)

    def makeAng(self,impactLog = 'impact_param.log'):
        """
        Writes  lig_ang.prm
        """
        self.log = open(impactLog,'r').readlines()
        angFile = open('%s/%s_ang.prm' % (self.app.workdir, self.ligname),'w')

        angFile.write('! Ligand angle parameters for %s\n' % self.ligname)

        prmDict = self.ligprmDict(self.impact_log, self.lig_pdb)[1]
        ang_added = []

        for i in range(self.angStart, self.angEnd + 1):
            atom1 = self.prm_id.get() + prmDict[self.log[i].split()[0]]
            atom2 = self.prm_id.get() +  prmDict[self.log[i].split()[1]]
            atom3 = self.prm_id.get() + prmDict[self.log[i].split()[2]]
            if [atom1, atom2, atom3] not in ang_added and [atom3, atom2, atom1] not in ang_added:
                ang_added.append([atom1, atom2, atom3])
                ang_added.append([atom3, atom2, atom1])
                K = 2 * float(self.log[i].split()[3])
                Theta = float(self.log[i].split()[4])
                comment = str(self.log[i][48:].strip())
                if '\n' not in comment:
                    comment += '\n'
                angFile.write('%4s  %4s  %4s  %7.2f  %7.2f  !%s'\
                          % (atom1,atom2,atom3,K,Theta,comment))

        angFile.write('\n')
        angFile.close()
        print '%s_ang.prm written' % self.ligname
        self.app.log('info', '%s_ang.prm successfully generated.' % self.ligname)

    def makeTor(self,impactLog = 'impact_param.log'):
        """
        Writes  lig_tor.prm
        """
        self.log = open(impactLog,'r').readlines()
        torFile = open('%s/%s_tor.prm' % (self.app.workdir, self.ligname),'w')

        torFile.write('! Ligand torsion parameters for %s\n' % self.ligname)

        prmDict = self.ligprmDict(self.impact_log, self.lig_pdb)[1]

        tor_added = []
        for i in range(self.torStart, self.torEnd + 1):
            try:
                atom1 = self.prm_id.get() + prmDict[self.log[i].split()[0]]
                atom2 = self.prm_id.get() + prmDict[self.log[i].split()[1]]
                atom3 = self.prm_id.get() + prmDict[self.log[i].split()[2]]
                atom4 = self.prm_id.get() + prmDict[self.log[i].split()[3]]
                if [atom1, atom2, atom3, atom4] not in tor_added and [atom4, atom3, atom2, atom1] not in tor_added:
                    tor_added.append([atom1, atom2, atom3, atom4])
                    tor_added.append([atom4, atom3, atom2, atom1])
                    V1 = 0.5 * float(self.log[i].split()[4])
                    V2 = 0.5 * float(self.log[i].split()[5])
                    V3 = 0.5 * float(self.log[i].split()[6])
                    comment = str(self.log[i][68:].strip())
                    if not '\n' in comment:
                        comment += '\n'

                    if (V1 == 0 and V2 == 0  and V3 == 0):
                        phace = 1
                        torFile.write('%4s   %4s   %4s   %4s   %7.3f   %3d    0.000   1 !%s'\
                            % (atom1.ljust(4),atom2.ljust(4), atom3.ljust(4), atom4.ljust(4), V1, phace,comment))

                    if V1 != 0:
                        if V2 != 0:
                            phace = -1
                        else:
                            phace = 1
                        torFile.write('%4s   %4s   %4s   %4s   %7.3f   %3d    0.000   1 !%s'\
                            % (atom1.ljust(4),atom2.ljust(4), atom3.ljust(4), atom4.ljust(4), V1, phace,comment))
                    if V2 != 0:
                        if V3 != 0:
                            phace = -2
                        else:
                            phace = 2
                        torFile.write('%4s   %4s   %4s   %4s   %7.3f   %3d  180.000   1 !%s'\
                            % (atom1.ljust(4),atom2.ljust(4), atom3.ljust(4), atom4.ljust(4), V2, phace,comment))
                    if V3 != 0:
                        phace = 3
                        torFile.write('%4s   %4s   %4s   %4s   %7.3f   %3d    0.000   1 !%s'\
                            % (atom1.ljust(4),atom2.ljust(4), atom3.ljust(4), atom4.ljust(4), V3, phace,comment))
            except:
                continue
        torFile.write('\n')
        torFile.close()

        print '%s_tor.prm written' % self.ligname
        self.app.log('info','%s_tor.prm successfully generated.' % self.ligname)

    def makeImp(self, impactLog = 'impact_param.log'):
        """
        Writes  lig_imp.prm
        """
        self.log = open(impactLog,'r').readlines()
        impFile = open('%s/%s_imp.prm' %(self.app.workdir, self.ligname), 'w')

        impFile.write('! Ligand improper parameters for %s\n' % self.ligname)

        try:
            if self.ligPrm:
                prmDict = self.ligPrm
        except:
            prmDict = self.ligprmDict(self.impact_log, self.lig_pdb)[1]

        i = self.impStart
        imp_added = []
        done = False
        while not done:
            if len(self.log[i].split()) > 5:
                atom1 = self.log[i].split()[0]
                atom2 = self.log[i].split()[1]
                atom3 = self.log[i].split()[2]
                atom4 = self.log[i].split()[3]
                V2 = 0.5 * float(self.log[i].split()[4])
                comment = str(self.log[i][43:].strip())
                if '\n' not in comment:
                    comment += '\n'

                qimp1, qimp2, qimp3, qimp4 = self.impQ(atom1, atom2, atom3, atom4)

                qimp1 = self.prm_id.get() + prmDict[qimp1]
                qimp2 = self.prm_id.get() + prmDict[qimp2]
                qimp3 = self.prm_id.get() + prmDict[qimp3]
                qimp4 = self.prm_id.get() + prmDict[qimp4]
                if [qimp1, qimp2, qimp3, qimp4] not in imp_added and [qimp1, qimp2, qimp4, qimp3] not in imp_added:
                    imp_added.append([qimp1, qimp2, qimp3, qimp4])
                    imp_added.append([qimp1, qimp2, qimp4, qimp3])

                    impFile.write('%4s   %4s   %4s   %4s   %7.3f  180.000 !%s'\
                              % (qimp1.ljust(4), qimp2.ljust(4), qimp3.ljust(4), qimp4.ljust(4), V2, comment))
                i += 1
            else:
                done = True

        impFile.write('\n')
        impFile.close()
        print '%s_imp.prm written' % self.ligname
        self.app.log('info', '%s_imp.prm successfully generated' % self.ligname)

    def writePrm(self):
        """
        Merges vdw, bond, angle, tors and imp parameters to final file for Q.
        """
        qname = 'Q%s_%s.prm' % (self.opls, self.ligname)
        finalFile = open('%s/%s' % (self.app.workdir, qname), 'w')
        liglist = []
        try:
            ligvdw = open('%s/%s_vdw.prm' % (self.app.workdir, self.ligname), 'r').readlines()
            liglist.append(ligvdw)
            if self.bndStart:
                ligbnd = open('%s/%s_bnd.prm' % (self.app.workdir, self.ligname), 'r').readlines()
                liglist.append(ligbnd)
            if self.angStart:
                ligang = open('%s/%s_ang.prm'% (self.app.workdir, self.ligname), 'r').readlines()
                liglist.append(ligang)
            if self.torStart:
                ligtor = open('%s/%s_tor.prm' % (self.app.workdir, self.ligname), 'r').readlines()
                liglist.append(ligtor)
            if self.impStart:
                ligimp = open('%s/%s_imp.prm' % (self.app.workdir, self.ligname), 'r').readlines()
                liglist.append(ligimp)
        except:
            print '\nWarning!'
            print 'Can not find ligand parameter files in specified folder ...\n'
            print 'Qoplsaa_lig.prm will be written without ligand parameters!\n'
            return

        finalFile.write('*-------------------------------------------------------------\n*\n')
        finalFile.write('*Q-FF parameters for %s (auto-generated with Q-gui 1.00)\n' % self.ligname)
        finalFile.write('*Translation of %s force field parameters assigned by ffld_server (schrodinger)\n' % self.opls)
        finalFile.write('*Always use auto-generated parameters with care!\n')
        finalFile.write('*-------------------------------------------------------------\n')
        finalFile.write('[options] force-field options\n')
        finalFile.write('name    %s\n' % qname.split('.')[0])
        finalFile.write('type    AMBER\n')
        finalFile.write('vdw_rule        geometric !vdW combination rule\n')
        finalFile.write('scale_14        0.5     !electrostatic 1-4 scaling factor\n')
        finalFile.write('switch_atoms    on\n')
        finalFile.write('improper_potential      periodic\n')
        finalFile.write('improper_definition     explicit\n\n')

        for i in range(len(liglist)):
            if liglist[i] == ligvdw:
                finalFile.write('[atom_types] atom type definitions\n')
                finalFile.write('*iaci  Avdw1    Avdw2      Bvdw1   Avdw3    Bvdw2&3   mass   comment\n')
                finalFile.write('*-----------------------------------------------------------------\n')
            elif liglist[i] == ligbnd:
                finalFile.write('[bonds] bond type definitions\n')
                finalFile.write('*iaci iacj    forceK      R0   comment\n')
                finalFile.write('*-----------------------------------------------------------------\n')
            elif liglist[i] == ligang:
                finalFile.write('[angles] angle type definitions\n')
                finalFile.write('*iaci iacj iack   forceK      angle0   comment\n')
                finalFile.write('*-----------------------------------------------------------------\n')
            elif liglist[i] == ligtor:
                finalFile.write('[torsions] torsions type definitions\n')
                finalFile.write('*iaci   iacj    iack    iacl    forceK  minima  phase   paths   comment\n')
                finalFile.write('*-----------------------------------------------------------------\n')
            elif liglist[i] == ligimp:
                finalFile.write('[impropers] improper type definitions\n')
                finalFile.write('*iaci   iacj    iack    iacl     forceK  phase   comment\n')
                finalFile.write('*-----------------------------------------------------------------\n')
            for line in liglist[i]:
                finalFile.write(line)

        #Insert headings even though parameters do not exist:
        if not self.torStart:
            finalFile.write('[torsions] torsions type definitions\n')
            finalFile.write('*iaci   iacj    iack    iacl    forceK  minima  phase   paths   comment\n')
            finalFile.write('*-----------------------------------------------------------------\n\n')

        if not self.impStart:
            finalFile.write('[impropers] improper type definitions\n')
            finalFile.write('*iaci   iacj    iack    iacl     forceK  phase   comment\n')
            finalFile.write('*-----------------------------------------------------------------\n\n')


        finalFile.write('\n ')
        finalFile.close()

        self.app.log('info','%s successfully generated' % qname)
        self.app.errorBox('Info','%s.lib with the corresponding %s.pdb and %s were successfully generated.' %
                                 (self.ligname, self.ligname, qname))

    def run_parameterisation(self):
        """
        Try to generate parameters...
        """
        self.impact_log = '%s/.tmpfile' % self.app.workdir
        self.ligname = self.resname.get().strip()
        self.lig_pdb = '%s.pdb' % self.ligname

        #Only need PDB file... use ffld_server -print_parameters -version 11 (for opls2001) -version 15 (opls2005)

        #Create pdb with only selected atoms
        self.create_pdb()

        self.getParam()
        if not self.impact_failed:
            self.makeLib(self.app.workdir + '/.tmpfile')
            self.makeVdw(self.app.workdir + '/.tmpfile')
            if self.bndStart:
                self.makeBnd(self.app.workdir + '/.tmpfile')
            if self.angStart:
                self.makeAng(self.app.workdir + '/.tmpfile')
            if self.torStart:
                self.makeTor(self.app.workdir + '/.tmpfile')
            if self.impStart:
                self.makeImp(self.app.workdir + '/.tmpfile')
            self.writePrm()
            self.app.log('info', 'Parameter assignment completed')

        else:
            self.app.log('info','Parameter assignment failed!')
            self.app.log(' ', 'Please verify that you have a valid structure/file format.')

    def dialog_window(self):
        """
        Window
        """
        self.title('Generate Q parameters')

        frame = Frame(self, bg=self.main_color)
        frame.pack(padx=(10, 10), pady=(10, 10))

        frame2 = Frame(self, bg=self.main_color)
        frame2.pack(padx=(10, 10), pady=(10, 10), fill=X)

        info_label = LabelFrame(frame, text='INFO', padx=10, bg=self.main_color)
        info_label.grid()

        info_txt = Text(frame, width=70, height=2, borderwidth=0, relief=GROOVE, highlightthickness=0, bg=self.main_color)
        info_txt.grid(in_=info_label, row=0, column=0, columnspan=10)
        info_txt.insert(0.0, '* Schrodinger path must be set in settings.\n'
                             '* Multiple residues will be redefined into one residue')
        info_txt.config(state=DISABLED)

        replace = Label(frame2, text='Replace original entry in current pdb:', bg=self.main_color)
        replace.grid(row=1, column =0, columnspan=4, sticky='w')

        self.replace_check = Checkbutton(frame2, bg=self.main_color, variable=self.check_replace)
        self.replace_check.grid(row=1, column=4)

        use_pdb = Label(frame2, text='Select atoms from current pdb:', bg=self.main_color)
        use_pdb.grid(row=2, column=0, columnspan=4,sticky='w')

        #self.pdb_check = Radiobutton(frame2, variable=self.check_mae, value=2, bg=self.main_color,
        #                             command=self.trace_mae_pdb)
        #self.pdb_check.grid(row=2, column=4, sticky='w')

        #self.mae_button = Button(frame2, text = 'Load .mae', highlightbackground=self.main_color, width=8,
        #                         command=self.get_local_file)
        #self.mae_button.grid(row=1, column=5, sticky='w')

        #self.maestructure = Text(frame2, width=18, height=1, bg=self.main_color, relief=GROOVE, highlightthickness=0)
        #self.maestructure.grid(row=1, column=6, columnspan=2)

        self.select_button = Button(frame2, text='Select', highlightbackground=self.main_color, width=8,
                                    command=lambda: self.select_atoms(self.atom_i, self.atom_j))
        self.select_button.grid(row=2, column=5)

        self.atom_i = Entry(frame2, width=6, highlightthickness=0, relief=GROOVE, textvariable=self.atomi_changed)
        self.atom_i.grid(row=2, column=6)

        self.atom_j = Entry(frame2, width=6, highlightthickness=0, relief=GROOVE, textvariable=self.atomj_changed)
        self.atom_j.grid(row=2, column=7)

        resname = Label(frame2, text='Residue name:', bg=self.main_color)
        resname.grid(row=3, column=0, columnspan=3, sticky = 'w', pady=(10,0))

        self.resname = Entry(frame2, width = 6, highlightthickness=0, relief=GROOVE)
        self.resname.grid(row=3, column=3, pady=(10,0))

        resnr = Label(frame2, text= 'Residue nr.:', bg=self.main_color)
        resnr.grid(row=4, column=0, columnspan=3, sticky='w')

        self.resnr = Entry(frame2, width=6, highlightthickness=0, relief=GROOVE)
        self.resnr.grid(row=4, column=3)

        atomnr = Label(frame2, text = 'Atom nr. start:', bg=self.main_color)
        atomnr.grid(row=5, column=0, columnspan=3, sticky='w')

        self.atomnr = Entry(frame2, width=6, highlightthickness=0, relief=GROOVE)
        self.atomnr.grid(row=5, column=3)

        prmname = Label(frame2, text = 'Parameter id.:', bg=self.main_color)
        prmname.grid(row=6, column=0, columnspan=3, sticky='w')

        self.prm_id = Spinbox(frame2, width = 4, highlightthickness=0, relief=GROOVE,
                              values=('A','B','C','D','E','F','G','H','I',
                                                         'J','K','L','M','N','O','P','Q','R',
                                                         'S','T','U','V','W','X','Y','Z'),
                              textvariable=self.prm_id_var)
        self.prm_id.grid(row=6, column=3)

        forcefield = Label(frame2, text='Force field:', bg=self.main_color)
        forcefield.grid(row=3, rowspan=3, column=5, columnspan=1)

        self.forcefield = Spinbox(frame2, width=9, highlightthickness=0, relief=GROOVE, values=('OPLS2005'))
        self.forcefield.grid(row=3, rowspan=3, column=6,columnspan=2)

        write_pdb_button = Button(frame2, width=7, text='Write pdb', highlightbackground=self.main_color, command=self.create_pdb)
        write_pdb_button.grid(row=6, column=5, pady=(10,0))

        self.run_button = Button(frame2, width=4, text='Run', highlightbackground=self.main_color,
                                 command=self.run_parameterisation)
        self.run_button.grid(row=6, column=6, pady=(10, 0))

        self.close_button = Button(frame2, width=4, text='Close', highlightbackground=self.main_color, command=self.destroy)
        self.close_button.grid(row=6, column=7, pady=(10, 0))




