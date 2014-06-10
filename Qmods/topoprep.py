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

from Tkinter import Radiobutton, Spinbox, StringVar, Entry, Text, Label, Frame, Button, Scrollbar, Toplevel, Checkbutton, Listbox, DISABLED, NORMAL, END,  GROOVE, LEFT, RIGHT, IntVar
#from qgui import PDB_COLOR

# -*- coding: utf-8 -*-

import os
import tkFont
import prepareTopology as pt
from select_xyz import AtomSelect
from edit_file import FileEdit
import time
import numpy as np



class TopologyPrepare(Toplevel):
    """Implements a dialog-box when Prepare -> Topology is chosen from menubar.
    Has got methods to modify pdf file/structure such that it can be used by Q.
    TopologyPrepare has got reference to the MainWindow-class."""

    def __init__(self, app, root, pdbfile, prm, lib, qgui_parent=True):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.pdbfile = pdbfile
        self.check_variable = IntVar()
        self.neutralize = IntVar()
        self.xvar = StringVar()
        self.yvar = StringVar()
        self.zvar = StringVar()
        self.qgui_parent = qgui_parent

        #Initialize dictionaries:
        self.charge_off = {'AR+': 'ARG', 'NAR+': 'NARG', 'CAR+': 'CARG',
                           'AS-': 'ASP', 'NAS-': 'NASP', 'CAS-': 'CASP',
                           'GL-': 'GLU', 'NGL-': 'NGLU', 'CGL-': 'CGLU',
                           'LY+': 'LYS', 'NLY+': 'NLYS', 'CLY+': 'CLYS',}

        self.charge_on = {'ARG': 'AR+', 'NARG': 'NAR+', 'CARG': 'CAR+',
                          'ASP': 'AS-', 'NASP': 'NAS-', 'CASP': 'CAS-',
                          'GLU': 'GL-', 'NGLU': 'NGL-', 'CGLU': 'CGL-',
                          'LYS': 'LY+', 'NLYS': 'NLY+', 'CLYS': 'CLY+'}

        self.charg_on_atomtypes = {'ARG': 'CZ', 'NARG': 'CZ', 'CARG': 'CZ',
                          'ASP': 'CG', 'NASP': 'CG', 'CASP': 'CG',
                          'GLU': 'CD', 'NGLU': 'CD', 'CGLU': 'CD',
                          'LYS': 'CD', 'NLYS': 'CD', 'CLYS': 'CD'}

        if not self.qgui_parent:
            self.app.log = self.app.app.log
            self.app.errorBox = self.app.app.errorBox
            self.app.workdir = self.app.app.workdir

        if not pt.checkPdb(self.pdbfile):
            pdb_id_q = pt.convertPdb(self.pdbfile, self.app.workdir, 0, 0)
            self.app.log('info','%s converted to Q format pdb (%s).' % (self.pdbfile.split('/')[-1], pdb_id_q.split('/')[-1]))
            self.pdbfile = self.app.workdir + '/' + pdb_id_q

        self.app.log('info', 'Topology Prepare session started.')
        self.dialog_box()
        self.set_structure_entry()


        #Auto insert name for topology and topology pdb
        self.topology_entry.insert(0.0,self.pdbfile.split('/')[-1].split('.')[0]+'.top')
        self.topo_pdb_entry.insert(0.0,self.pdbfile.split('/')[-1].split('.')[0]+'_top.pdb')

        #Generate SS?
        self.makeSS = False

        #XYZ sim center:
        self.set_sim_center()

        #Trace coordinate changes to charge radiuses:
        self.xvar.trace('w', self.simcenter_changed)
        self.yvar.trace('w', self.simcenter_changed)
        self.zvar.trace('w', self.simcenter_changed)

        #Sphere radius
        system_radius = pt.findRadius(self.pdbfile)
        self.sphere_entry.delete(0, END)
        self.sphere_entry.insert(0, '%d' % int(round(system_radius + 5)))

        #Set default solvation
        self.check_solvation.set(1)

        #Check if HIS in file and convert to HID:
        self.his = pt.getResnrs(self.pdbfile, 'HIS')
        self.hid = pt.getResnrs(self.pdbfile, 'HID')
        self.hie = pt.getResnrs(self.pdbfile, 'HIE')
        self.hip = pt.getResnrs(self.pdbfile, 'HIP')

        if len(self.his) > 0:
            self.app.log('info', 'Found residue HIS in file. Autoconverting to HID/HIE')
            self.app.log(' ','    ---------------------\n')
            for entry in self.his:
                pt.convertOldNew(self.pdbfile, entry, 'HID')
                self.app.log(' ','    HIS %4s --> HID %4s\n' % (entry, entry))
            self.app.log(' ','    ---------------------\n')


        #ParameterList (get this from .settings)
        self.prm = prm

        #Library list (get this from .settings)
        self.lib = lib

        self.res_charge = dict()
        self.res_atoms = dict()
        #Check if residues in pdb file exist in library:
        self.checkLib()

        self.set_total_charge()

        self.updateList()

    def set_sim_center(self):
        """
        finds the center of current loaded pdb file and
        write xyz into labels.
        """
        self.xc, self.yc, self.zc = pt.centerofmass(self.pdbfile)
        self.center_x_entry.delete(0, END)
        self.center_x_entry.insert(0, self.xc)
        self.center_y_entry.delete(0, END)
        self.center_y_entry.insert(0, self.yc)
        self.center_z_entry.delete(0, END)
        self.center_z_entry.insert(0, self.zc)

    def simcenter_changed(self, *args):
        """
        Updates list with chargable residues and radiuses whenever
        simulation center is changed
        """
        try:
            self.updateList()
        except:
            return

    def checkLib(self):
        """
        #Check if pdb residues exist in current library:
        """
        status = 'Lib status: OK'

        #Collect info from lib files
        for lib in self.lib:
            with open(lib, 'r') as lib:
                found_res = False
                found_atoms = False
                charges = []
                res = 'Hello'
                for line in lib:
                    if found_res:
                        if '[bonds]' in line or '*-------' in line:
                            self.res_charge[res] = round(np.sum(charges), 1)
                            print '%4s: %.1f' % (res, round(np.sum(charges), 1))
                            found_atoms = False
                            found_res = False

                        if found_atoms:
                            if len(line.split()) > 3:
                                if line.split()[0] != '!':
                                    self.res_atoms[res].append(line.split()[1])
                                    charges.append(float(line.split()[3]))
                        if '[atoms]' in line:
                            found_atoms = True

                    if '{' in line and '}' in line:
                        res = line.split('{')[1].split('}')[0]
                        self.res_atoms[res] = list()
                        self.res_charge[res] = 0.0
                        charges = list()
                        found_res = True

        #Get info from pdb file and compare to existing lib entries
        with open(self.pdbfile, 'r') as pdb:
            res_nr = 0
            for line in pdb:
                if 'ATOM' in line or 'HETATM' in line:
                    res = line[17:21].strip()
                    atom_name = line[13:17].strip()
                    if res_nr != int(line[21:26]):
                        res_nr = int(line[21:26])
                        if not res in self.res_atoms.keys():
                            self.app.log(' ','WARNING: %4s %5d not found in lib entries!\n' % (res, res_nr))
                            if status != 'Lib entries missing':
                                status = 'Lib entries missing'

                    if res in self.res_atoms.keys():
                        if not atom_name in self.res_atoms[res]:
                            self.app.log(' ','WARNING: Atomname %4s not found in lib for %4s %5d!\n' %
                                             (atom_name, res, res_nr))
                            if status != 'Lib entries missing':
                                status = 'Lib entries missing'

        self.lib_status.config(state=NORMAL)
        self.lib_status.delete(0.0, END)
        self.lib_status.insert(0.0, status.rjust(19))
        self.lib_status.config(state=DISABLED)

    def writeTopology(self):
        qprepinp_name = self.pdbfile.split('/')[-1].split('.')[0]+'_Qprep.inp'
        qprepinp = open(self.app.workdir + '/' + qprepinp_name,'w')
        self.topname = self.topology_entry.get(0.0,END).strip()
        self.top_pdbname = self.topo_pdb_entry.get(0.0,END).strip()
        for entry in self.lib:
            qprepinp.write('readlib %s\n' % entry)

        #Can only use 1 prm file. If more than 1 is in settings, merge them!
        if len(self.prm) > 1:
            self.app.log(' ','\n****\nMore than 1 parameter file is set in Settings.\nThese files will'
                                     ' be merged. \n--> Force fields must be of same type.'
                                     '\n--> If atom types with the same name exist,'
                                     'they will be redefined based on which comes last.\n****\n')
            atom_index = []
            bond_index = []
            angle_index = []
            torsion_index = []
            improper_index = []

            merged_name = 'merged.prm'
            merged_file = open('%s/%s' % (self.app.workdir, merged_name), 'w')

            #Go through files and find indices:
            prmfiles = []
            for entry in self.prm:
                prmfile = open(entry,'r').readlines()
                for i in range(len(prmfile)):
                    if '[atom_types]' in prmfile[i]:
                        atom_index.append(i)
                    elif '[bonds]' in prmfile[i]:
                        bond_index.append(i)
                    elif '[angles]' in prmfile[i]:
                        angle_index.append(i)
                    elif '[torsions]' in prmfile[i]:
                        torsion_index.append(i)
                    elif '[impropers]' in prmfile[i]:
                        improper_index.append(i)
                    else:
                        continue
                prmfiles.append(prmfile)

            #Merge atom types
            for i in range(0, bond_index[0]):
                merged_file.write(prmfiles[0][i])
            for entry in range(1, len(prmfiles)):
                for j in range(atom_index[entry] + 1, bond_index[entry]):
                    if '*' in prmfiles[entry][j][0:2]:
                        continue
                    else:
                        merged_file.write(prmfiles[entry][j])
            #Merge bonds
            try:
                for i in range(bond_index[0], angle_index[0]):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(bond_index[entry] + 1 ,angle_index[entry]):
                        if '*' in prmfiles[entry][j][0:2]:
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            #Merge angles
            try:
                for i in range(angle_index[0], torsion_index[0]):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(angle_index[entry] + 1 ,torsion_index[entry]):
                        if '*' in prmfiles[entry][j][0:2]:
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            #Merge torsions
            try:
                for i in range(torsion_index[0], improper_index[0]):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(torsion_index[entry] + 1 ,improper_index[entry]):
                        if '*' in prmfiles[entry][j][0:2]:
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            #Merge impropers
            try:
                for i in range(improper_index[0], len(prmfiles[0])):
                    merged_file.write(prmfiles[0][i])
                for entry in range(1, len(prmfiles)):
                    for j in range(improper_index[entry] + 1 , len(prmfiles[entry])):
                        if '*' in prmfiles[entry][j][0:3]:
                            continue
                        else:
                            merged_file.write(prmfiles[entry][j])
            except:
                pass

            merged_file.close()
            qprepinp.write('readprm %s/%s\n' % (self.app.workdir, merged_name))
        else:
            qprepinp.write('readprm %s\n' % self.prm[0])

        qprepinp.write('readpdb %s\n' % self.pdbfile)

        if self.makeSS:
            for atompair in range(len(self.cys_residues)):
                qprepinp.write('addbond %s:SG %s:SG\n' % (self.cys_residues[atompair][0], self.cys_residues[atompair][1]))

        radius = self.sphere_entry.get()

        xc = self.center_x_entry.get()
        yc = self.center_y_entry.get()
        zc = self.center_z_entry.get()

        qprepinp.write('boundary 1 %s %s %s %s\n' % (xc,yc,zc,radius))

        if self.check_solvation.get() != 0:
            if self.check_solvation.get() == 1:
                solvation = 'HOH'
            else:
                solvation = 'SPC'
            qprepinp.write('solvate %s %s %s %s 1 %s\n' % (xc,yc, zc,radius,solvation))

        qprepinp.write('maketop %s topology\n' % self.pdbfile.split('.')[0])
        qprepinp.write('writetop %s/%s' % (self.app.workdir, self.topname))
        qprepinp.write('\n')
        qprepinp.write('writepdb %s/%s' % (self.app.workdir, self.top_pdbname))
        qprepinp.write('\n')
        qprepinp.write('y\n')
        qprepinp.write('q\n')
        qprepinp.close()

        self.app.log('info','Qprep input file written: %s' % qprepinp_name)

    def run_qprep(self):
        qprepinp_name = self.pdbfile.split('/')[-1].split('.')[0]+'_Qprep.inp'
        qprepout_name = self.pdbfile.split('/')[-1].split('.')[0]+'_Qprep.log'

        self.writeTopology()
        os.system('Qprep5 <%s/%s>%s/%s' % (self.app.workdir, qprepinp_name, self.app.workdir, qprepout_name))

        qprep_done = False
        count = 0
        max_count = 10
        while qprep_done == False:
            count += 1
            if os.path.isfile(self.app.workdir + '/' + qprepout_name):
                with open(self.app.workdir + '/' + qprepout_name, 'r') as qpreplog:
                    for line in qpreplog:
                        self.app.log(' ', line)
                        if 'Topology successfully generated' in line:
                            #Qprep run OK: Update entry fields in parent (Qgui main or LIE)
                            if not self.qgui_parent:
                                if self.app.add_title == 'complex':
                                    self.app.complex_pdb = self.app.workdir + '/' + self.top_pdbname
                                    self.app.complex_top = self.app.workdir + '/' + self.topname
                                elif self.app.add_title == 'ligand':
                                    self.app.ligand_pdb = self.app.workdir + '/' + self.top_pdbname
                                    self.app.ligand_top = self.app.workdir + '/' + self.topname
                                self.app.update_progress()

                            else:
                                self.app.pdb_id = self.app.workdir + '/' + self.top_pdbname
                                self.app.top_id = self.app.workdir + '/' + self.topname
                                self.app.update_pdb_id_entryfield()
                                self.app.main_window.set_topology_entryfield(self.topname.split('/')[-1])

                            self.app.log('info','Topology successfully generated')
                            qprep_done = True
                            break
                        elif 'ERROR:' in line:
                            self.app.log('info','Topology generation failed.')
                            self.app.errorBox('Warning','Topology contains errors. Check log file and parameters.')
                            qprep_done = True
                            break
                        elif 'topology is incomplete' in line:
                            self.app.log('info', 'Topology contains errors! Please check log file!')
                            self.app.errorBox('Warning','Topology contains errors. Check log file and parameters.')
                            qprep_done = True
                            break
                        elif 'Correct the PDB file' in line:
                            self.app.log('info', 'Topology contains errors! Please check log file!')
                            self.app.errorBox('Warning','Topology contains errors. Check log file and parameters.')
                            qprep_done = True
                            break
            else:
                time.sleep(1)

            if count >= max_count:
                self.app.log('info','Qprep failed to generate topology!')
                self.app.errorBox('Error','Could not run Qprep! Please verify installation.')
                qprep_done = True

    def set_structure_entry(self):
        """Inserts the name of the current file to structure_entry field. """
        self.structure_entry.config(state = NORMAL)
        self.structure_entry.insert(0.0, self.pdbfile.split('/')[-1].rjust(19))
        self.structure_entry.config(state = DISABLED)

    def set_total_charge(self):

        """Gets the overall charge from countCharges() function from prepareTopology
        file. Inserts the charge into total_c_entry field. """
        self.total_charge = 0
        res_count = dict()

        res_nr = 0
        with open(self.pdbfile, 'r') as pdb:
            for line in pdb:
                if 'ATOM' in line or 'HETATM' in line:
                    if res_nr != int(line[21:26]):
                        res_nr = int(line[21:26])
                        res = line[17:21].strip()
                        if res in self.res_charge.keys():
                            self.total_charge += self.res_charge[res]
                            if abs(self.res_charge[res]) > 0:
                                if res not in res_count:
                                    res_count[res] = self.res_charge[res]
                                else:
                                    res_count[res] += self.res_charge[res]
                        else:
                            print 'Can not set charge for %4s: Residue not in lib!' % res

        self.app.log(' ', '\n____________________________________________\n')
        for res in res_count.keys():
            self.app.log(' ', '    %4s   %7.2f\n' % (res, res_count[res]))
        self.app.log(' ', '____________________________________________\n')
        self.app.log(' ', 'Sum charge %7.2f\n' % self.total_charge)
        self.app.log(' ', '============================================\n')

        self.total_c_entry.config(state=NORMAL)
        self.total_c_entry.delete(0.0, END)
        self.total_c_entry.insert(0.0, self.total_charge)
        self.total_c_entry.config(state=DISABLED)

    def turn_off_charges_button_pressed(self, off=True):
        """This method is called when the charge_off_button is pressed.
        Calls for function toggleAllCharges() from prepareTopology file.
        Will toggle on/off all charges within 5/6r where r is the sim. sphere r."""
        if off:
            charge = self.charge_off
            toggle = 'off'
        else:
            charge = self.charge_on
            toggle = 'on'

        xc = self.center_x_entry.get()
        yc = self.center_y_entry.get()
        zc = self.center_z_entry.get()
        simrad = self.sphere_entry.get()

        pt.toggleAllCharges(self.pdbfile, charge, self.charg_on_atomtypes, toggle, simrad, xc, yc, zc)

        self.set_total_charge()
        self.updateList()

    def toggle_charge_selected(self):
        charges = False
        nterm = False
        cterm = False


        neutral = ['ARG','LYS','ASP','GLU']
        charged = ['AR+','LY+','AS-','GL-']
        his = ['HID','HIE','HIP']

        try:
            list_index = int(self.listbox.curselection()[0])
        except:
            return

        oldRes, resnr = self.listbox.get(list_index).split()[0:2]
        radius = self.listbox.get(list_index).split('|')[-1]
        if '%5s' % resnr in self.nterm_nr:
            nterm = True
        if '%5s' % resnr in self.cterm_nr:
            cterm = True
        if oldRes in neutral:
            charges = True
        if oldRes in charged:
            charges = True
        if oldRes in his:
            charges = True

        if charges:
            resnr = '%4d' % int(resnr)
            if oldRes in neutral:
                old = neutral
                new = charged
            elif oldRes in charged:
                old = charged
                new = neutral
            elif oldRes in his:
                old = his
                new = his
            else:
                return
            for i in range(len(old)):
                if oldRes == old[i]:
                    if old != new:
                        newRes = new[i]
                    else:
                        try:
                            newRes = new[i+1]
                        except:
                            newRes = new[0]

            pt.convertOldNew(self.pdbfile, resnr, newRes)
            self.app.log('info','%s %s --> %s %s' % (oldRes,resnr,newRes,resnr))
            self.listbox.insert(list_index, '%4s   %4d' % (newRes, int(resnr)))
            self.listbox.delete(list_index + 1 )
            self.listbox.selection_set(list_index)

        if nterm:
            resnr = '%5s' % resnr
            newRes = pt.convertNterminal(self.pdbfile,resnr)

        if cterm:
            resnr = '%5s' % resnr
            newRes = pt.convertCterminal(self.pdbfile, resnr)
        if not charges:
            self.app.log('info','%s %s --> %s %s' % (oldRes,resnr,newRes,resnr))
        self.listbox.insert(list_index, '%4s %4d    |%5.1f' % (newRes, int(resnr), float(radius)))
        self.listbox.delete(list_index + 1 )
        self.listbox.selection_set(list_index)

        self.set_total_charge()

    def create_ss_bonds_checkbutton_pressed(self):
        """Find potential S-S bridges. Returns:
        1. Arrays with atomnumber pairs to be used in Qprep5: makebond atom_i atom_j
        2. Array with residue number pairs.
        3. List with distances between S-S below treshold. """
        state = self.check_variable.get()
        if state == 0:
            self.makeSS = False
        if state == 1:
            self.makeSS = True

        self.sg_atompairs, residues, sg_dist = pt.findSS(self.pdbfile)
        self.cys_residues = []
        if self.makeSS:
            print "Create S-S bonds selected:"
            if len(residues) > 0:
                pt.convertOldNew(self.pdbfile, residues,'CYX')
                sg = 0
                self.app.log('info','Creating S-S bonds for:')
                self.app.log(' ','    -------------------------\n')
                for cyx in range(len(residues)-1,0, -2):
                    self.app.log(' ','    CYX%s - CYX%s (%.2f A)\n'
                    % (residues[cyx-1], residues[cyx], sg_dist[sg]))
                    sg += 1
                    tmp = [residues[cyx-1], residues[cyx]]
                    self.cys_residues.append(tmp)
                self.app.log(' ','    -------------------------\n')
        else:
            pt.convertOldNew(self.pdbfile, residues,'CYS')

    def updateList(self):
        """
        Function to update listbox with charged/neutral residues
        for toggle state selection.
        """

        self.listbox.delete(0, END)

        #Get residue numbers for charged residues:
        arg = pt.getResnrs(self.pdbfile, 'ARG')
        arg_pos = pt.getResnrs(self.pdbfile, 'AR+')
        lys = pt.getResnrs(self.pdbfile, 'LYS')
        lys_pos = pt.getResnrs(self.pdbfile, 'LY+')
        asp = pt.getResnrs(self.pdbfile, 'ASP')
        asp_neg = pt.getResnrs(self.pdbfile, 'AS-')
        glu = pt.getResnrs(self.pdbfile, 'GLU')
        glu_neg = pt.getResnrs(self.pdbfile, 'GL-')
        hid = pt.getResnrs(self.pdbfile, 'HID')
        hie = pt.getResnrs(self.pdbfile, 'HIE')
        hip = pt.getResnrs(self.pdbfile, 'HIP')

        reslist = [arg,arg_pos,lys,lys_pos,asp,asp_neg,glu,glu_neg,hid,hie,hip]

        #Get current x,y,z simulation center:
        xc = self.center_x_entry.get()
        yc = self.center_y_entry.get()
        zc = self.center_z_entry.get()

        #Get radius from center for charged residues (atom NZ, CZ, CG, CD, N, C, CE1)
        arg_r = []
        arg_pos_r = []
        lys_r = []
        lys_pos_r = []
        asp_r = []
        asp_neg_r = []
        glu_r = []
        glu_neg_r = []
        hid_r = []
        hie_r = []
        hip_r = []

        atomtype = 'CA'

        for residue in reslist:
            if residue == arg:
                radius = arg_r
                atomtype = 'CZ'
            elif residue == arg_pos:
                radius = arg_pos_r
                atomtype = 'CZ'
            elif residue == lys:
                radius = lys_r
                atomtype = 'NZ'
            elif residue == lys_pos:
                radius = lys_pos_r
                atomtype = 'NZ'
            elif residue == asp:
                radius = asp_r
                atomtype = 'CG'
            elif residue == asp_neg:
                radius = asp_neg_r
                atomtype = 'CG'
            elif residue == glu:
                radius = glu_r
                atomtype = 'CD'
            elif residue == glu_neg:
                radius = glu_neg_r
                atomtype = 'CD'
            elif residue == hid:
                radius = hid_r
                atomtype = 'CE1'
            elif residue == hie:
                radius = hie_r
                atomtype = 'CE1'
            elif residue == hip:
                radius = hip_r
                atomtype = 'CE1'

            for resnr in residue:
                radius.append(pt.calcDist(self.pdbfile, resnr, atomtype, xc, yc, zc))

        radiuslist = [arg_r,arg_pos_r,lys_r,lys_pos_r,asp_r,asp_neg_r,glu_r,glu_neg_r,hid_r,hie_r,hip_r]

        resnames = ['ARG','AR+','LYS','LY+','ASP','AS-','GLU','GL-','HID','HIE','HIP']

        #Insert residues that can be charged/neutral
        for residue in range(len(resnames)):
            if len(reslist[residue]) > 0:
                for resnr in range(len(reslist[residue])):
                    self.listbox.insert(END,'%4s %4d    |%5.1f' %
                                            (resnames[residue], int(reslist[residue][resnr]), radiuslist[residue][resnr]))

        #Find N and C terminals:
        self.nterm_nr, self.nterm_res, self.cterm_nr, self.cterm_res = pt.findTerminals(self.pdbfile)

        #Find terminal radius to center:
        nter_r = []
        cter_r = []
        for terminal in [self.nterm_nr, self.cterm_nr]:
            if terminal == self.nterm_nr:
                radius = nter_r
                atomtype = 'N'
            else:
                radius = cter_r
                atomtype = 'C'
            for resnr in terminal:
                radius.append(pt.calcDist(self.pdbfile, resnr, atomtype, xc, yc, zc))

        #Insert Terminals:
        if len(self.nterm_nr) > 0:
            for i in range(len(self.nterm_nr)):
                self.listbox.insert(END, '%4s %4d (N)|%5.1f' % (self.nterm_res[i], int(self.nterm_nr[i]), nter_r[i]))
        if len(self.cterm_nr) > 0:
            for i in range(len(self.cterm_nr)):
                self.listbox.insert(END, '%4s %4d (C)|%5.1f' % (self.cterm_res[i], int(self.cterm_nr[i]), cter_r[i]))

    def open_atom_select(self):
        """
        Opens a new window to select atom in pdb file as simulation center
        """
        self.select_pdb = AtomSelect(self, self.root, self.pdbfile,
                                     self.center_x_entry, self.center_y_entry, self.center_z_entry)
        self.select_pdb.configure(bg = self.main_color)
        self.select_pdb.title('Select simulation center')
        self.select_pdb.resizable()

    def not_available(self):
        """
        Function not available yet pop-up
        """
        self.app.errorBox('Warning','Sorry, this function is not implemented yet!')
        self.neutralize.set(0)

    def editPdb(self):
        """
        Opens up pdb file for editing
        """
        self.fileEdit = FileEdit(self, self.pdbfile)
        self.fileEdit.config(bg=self.main_color)
        self.fileEdit.title('Edit file')
        self.fileEdit.resizable()


    def dialog_box(self):
        """Defines the outlook of Topology Prepare window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """

        self.title('Topology Prepare')
        self.config(background=self.main_color)

        # Define frames
        left_frame = Frame(self, bg = self.main_color)
        left_frame.pack(side = LEFT,pady=10, padx=(10,0))

        right_frame = Frame(self, bg = self.main_color)
        right_frame.pack(side = RIGHT)


        # Define elements in the left_frame
        sphere_label = Label(left_frame, text = 'Simulation sphere:')
        sphere_label.grid(row = 0, column = 0, sticky='e')
        sphere_label.config(background = self.main_color)

        self.sphere_entry = Spinbox(left_frame, width = 5, highlightthickness = 0, relief = GROOVE, from_=1, to=200)
        self.sphere_entry.grid(row = 0, column = 1)

        aangstrom_label = Label(left_frame, text='%s' % u'\xc5')
        aangstrom_label.grid(row = 0, column = 2, sticky='w')
        aangstrom_label.config(background = self.main_color)

        centre_label = Label(left_frame, text = 'Simulation centre:')
        centre_label.grid(row = 1, column = 0, sticky = 'e')
        centre_label.config(background = self.main_color)

        self.center_x_entry = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE,
                                    textvariable=self.xvar)
        self.center_x_entry.grid(row = 1, column = 1)


        self.center_y_entry = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE,
                                    textvariable=self.yvar)
        self.center_y_entry.grid(row = 1, column = 2)

        self.center_z_entry = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE,
                                    textvariable=self.zvar)
        self.center_z_entry.grid(row = 1, column = 3)

        center_x_label = Label(left_frame, text = 'x')
        center_x_label.grid(row = 2, column = 1)
        center_x_label.config(background = self.main_color)

        center_y_label = Label(left_frame, text = 'y')
        center_y_label.grid(row = 2, column = 2)
        center_y_label.config(background = self.main_color)

        center_z_label = Label(left_frame, text = 'z')
        center_z_label.grid(row = 2, column = 3)
        center_z_label.config(background = self.main_color)

        sim_change_button = Button(left_frame, text = 'Change', command = self.open_atom_select)
        sim_change_button.grid(row = 1, column = 4)
        sim_change_button.config(highlightbackground = self.main_color)

        sim_center_button = Button(left_frame, text = 'Center', command = self.set_sim_center)
        sim_center_button.grid(row = 2, column = 4)
        sim_center_button.config(highlightbackground = self.main_color)

        solvate_label = Label(left_frame, text = 'Solvate:')
        solvate_label.grid(row = 3, column = 0, sticky = 'e')
        solvate_label.config(background = self.main_color)

        tip3p_label = Label(left_frame, text = ' TIP3P')
        tip3p_label.grid(row = 4, column = 0, sticky = 'e')
        tip3p_label.config(background = self.main_color)

        self.check_solvation = IntVar()
        self.tip3p_radiobutton = Radiobutton(left_frame, variable = self.check_solvation, value=1)
        self.tip3p_radiobutton.grid(row = 4, column = 1, sticky = 'w')
        self.tip3p_radiobutton.config(background = self.main_color)

        spc_label = Label(left_frame, text = 'SPC')
        spc_label.grid(row = 5, column = 0, sticky = 'e')
        spc_label.config(background = self.main_color)

        self.scp_radiobutton = Radiobutton(left_frame, variable = self.check_solvation, value=2)
        self.scp_radiobutton.grid(row = 5, column = 1, sticky = 'w')
        self.scp_radiobutton.config(background = self.main_color)

        none_label = Label(left_frame, text = 'None')
        none_label.grid(row = 6, column = 0, sticky = 'e')
        none_label.config(background = self.main_color)

        self.none_radiobutton = Radiobutton(left_frame, variable = self.check_solvation, value = 0)
        self.none_radiobutton.grid(row = 6, column = 1, sticky = 'w')
        self.none_radiobutton.config(background = self.main_color)

        ssbond_label = Label(left_frame, text = 'Create S-S bonds:')
        ssbond_label.grid(row = 7, column = 0, sticky = 'e')
        ssbond_label.config(background = self.main_color)

        ssbond_checkbutton = Checkbutton(left_frame, variable = self.check_variable, command = self.create_ss_bonds_checkbutton_pressed)
        ssbond_checkbutton.grid(row = 7, column = 1, sticky = 'w')
        ssbond_checkbutton.config(background = self.main_color)

        total_c_label = Label(left_frame, text = 'Total charge:')
        total_c_label.grid(row = 8, column = 0, sticky = 'e')
        total_c_label.config(background = self.main_color)

        self.total_c_entry = Text(left_frame, width = 5, height = 1)
        self.total_c_entry.grid(row = 8, column = 1)
        self.total_c_entry.config(state = DISABLED, highlightthickness = 0)

        charge_all = Label(left_frame, text = 'Toggle all charges:', bg=self.main_color)
        charge_all.grid(row=9, column=0, sticky='e')

        charge_on_button = Button(left_frame, text = 'ON ',
                                  command = (lambda: self.turn_off_charges_button_pressed(False)))
        charge_on_button.grid(row = 9, column = 1, sticky = 'e')
        charge_on_button.config(highlightbackground = self.main_color)

        charge_off_button = Button(left_frame, text = 'OFF',
                                   command = (lambda: self.turn_off_charges_button_pressed(True)))
        charge_off_button.grid(row = 9, column = 2, sticky = 'w')
        charge_off_button.config(highlightbackground = self.main_color)


        toggle_sel = Label(left_frame, text = 'Toggle state:', bg=self.main_color)
        toggle_sel.grid(row=10, column=0, sticky='e')

        toggle_state_button = Button(left_frame, text = 'Selected', command = self.toggle_charge_selected)
        toggle_state_button.grid(row = 10, column = 1, columnspan=2)
        toggle_state_button.config(highlightbackground = self.main_color)

        neutralize_label = Label(left_frame, text = 'Neutralize system with NaCl:')
        neutralize_label.grid(row = 11, column = 0, sticky = 'e')
        neutralize_label.config(background = self.main_color)

        neutralize_checkbutton = Checkbutton(left_frame, variable=self.neutralize, command = self.not_available)
        neutralize_checkbutton.grid(row = 11, column = 1, sticky = 'w')
        neutralize_checkbutton.config(background = self.main_color)

        nacl_label = Label(left_frame, text  ='(Na+ Cl-)')
        nacl_label.grid(row = 12, column = 0, sticky = 'e')
        nacl_label.config(background = self.main_color)

        topology_label = Label(left_frame, text = 'Topology name:')
        topology_label.grid(row = 13, column = 0, sticky = 'e',pady=(10,0))
        topology_label.config(background = self.main_color)

        self.topology_entry = Text(left_frame, width = 40, height = 1)
        self.topology_entry.grid(row = 13, column = 1, columnspan = 4, pady=(10,0), sticky = 'w')
        self.topology_entry.config(state = NORMAL, highlightthickness = 0, relief = GROOVE)

        topo_pdb_label = Label(left_frame, text = 'Topology PDB:')
        topo_pdb_label.grid(row = 14, column = 0, sticky = 'e')
        topo_pdb_label.config(background = self.main_color)

        self.topo_pdb_entry = Text(left_frame, width = 40, height = 1)
        self.topo_pdb_entry.grid(row = 14, column = 1, columnspan = 4, sticky='w')
        self.topo_pdb_entry.config(state = NORMAL, highlightthickness = 0, relief = GROOVE)

        edit_pdb = Button(left_frame, text='Edit', width=5, command=self.editPdb)
        edit_pdb.grid(row=4, column=4)
        edit_pdb.config(highlightbackground=self.main_color)

        self.structure_entry = Text(left_frame, width = 20, height = 1)
        self.structure_entry.grid(row = 4, column = 2, columnspan = 2, sticky='e')
        self.structure_entry.config(state = DISABLED, bg = self.main_color, highlightthickness = 0, relief = GROOVE)

        self.lib_status = Text(left_frame, width=20, height=1)
        self.lib_status.grid(row=5, column=2, columnspan=2, sticky='e')
        self.lib_status.config(state = DISABLED, bg = self.main_color, highlightthickness = 0, relief = GROOVE)

        checklib = Button(left_frame, text='Check', command=self.checkLib)
        checklib.grid(row=5, column=4)
        checklib.config(highlightbackground=self.main_color)

        listbox_scroll = Scrollbar(left_frame)
        listbox_scroll.grid(row = 8, rowspan = 5, column = 7, sticky = 'nsw')
        self.listbox = Listbox(left_frame, yscrollcommand = listbox_scroll.set, highlightthickness = 0, relief = GROOVE)
        listbox_scroll.config(command=self.listbox.yview)
        self.listbox.grid(row = 8, rowspan = 5, column = 3, columnspan = 4, sticky = 'e')
        self.listbox.config(font=tkFont.Font(family="Courier", size=12))


        create_top_button = Button(left_frame, text = 'Write', command = self.writeTopology)
        create_top_button.grid(row = 15, column = 2, pady=(20,10), sticky = 'e')
        create_top_button.config(highlightbackground = self.main_color)

        run_top_button = Button(left_frame, text = ' Run ', command = self.run_qprep)
        run_top_button.grid(row = 15, column = 3,  pady=(20,10))
        run_top_button.config(highlightbackground = self.main_color)

        cancel_button = Button(left_frame, text = 'Close', command = self.destroy)
        cancel_button.grid(row = 15, column = 4, pady=(20,10))
        cancel_button.config(highlightbackground = self.main_color)












