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

from Tkinter import StringVar, SUNKEN,Y, Spinbox, Scrollbar,Listbox, Button, Entry, Text, Label, Frame, Toplevel, \
    Checkbutton, DISABLED, NORMAL, END, GROOVE, LEFT, IntVar

import prepareTopology as pt
from select_atoms import AtomSelectRange
from md_restraints import AddSequenceRestraints, AddAtomRestraints, AddDistanceRestraints, AddWallRestraints
import random
import os
import cPickle
from subprocess import call

class SetupMd(Toplevel):
    """Implements a dialog-box when Setup -> MD is chosen from menubar.
    Has got methods to generate input files for MD simuations with Qdyn"""

    def __init__(self, app, root, pdbfile, topology, run_md=True, fep=False): #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.pdbfile = pdbfile
        self.topology = topology
        self.fep = fep
        self.run_md = run_md
        self.restart_file = None

        self.simtimeVar = StringVar()
        self.simfilesVar = StringVar()
        self.stepsizeVar = StringVar()
        self.shake_solvent = IntVar()
        self.shake_solute = IntVar()
        self.shake_hydrogens = IntVar()
        self.lrf_check = IntVar()
        self.polarisation_check = IntVar()
        self.q_atoms = IntVar()
        self.use_eq = IntVar()

        self.qatoms_list = []
        self.sequence_restraints = []
        self.distance_restraints = []
        self.atom_restraints = []
        self.wall_restraints = []

        if not run_md:
            #self.app is Qgui main window
            self.app.errorBox = self.app.app.errorBox
            self.app.log = self.app.app.log
            self.app.settings_path = self.app.app.settings_path
            #self.parent is LIE/FEP/EVB window
            self.parent = self.app

        self.dialog_box()

        self.simtimeVar.trace('w', self.inputfiles_written)
        self.simfilesVar.trace('w', self.inputfiles_written)
        self.stepsizeVar.trace('w', self.inputfiles_written)

        if self.run_md:
            self.app.log('info', 'Setup MD session started.')

        atomnumbers = pt.pdb_info(self.pdbfile)[0]

        if not self.fep:
            self.qmd_settings = {'simulation time':0.5,
                        'stepsize': 1,
                        'inputfiles':1,
                        'temperature': 300,
                        'bath coupling': 10,
                        'shake solvent': 1,
                        'shake solute': 0,
                        'shake hydrogens': 0,
                        'lrf': 1,
                        'solute-solute cutoff': 10,
                        'solvent-solvent cutoff': 10,
                        'solute-solvent cutoff': 10,
                        'Q-atoms cutoff': 99,
                        'lrf cutoff': 99,
                        'shell force': 10.0,
                        'radial force': 60.0,
                        'polarization restraint': 1,
                        'polarization force': 20.0,
                        'update non-bonded list': 25,
                        'update energy summary': 5,
                        'write energy file': 10,
                        'write trajectory': 100,
                        'trajectory atoms': 'all'}
        else:
            self.qmd_settings = {'simulation time': self.app.md_settings['simtime'],
                        'stepsize': self.app.md_settings['stepsize'],
                        'inputfiles':self.app.md_settings['inputfiles'],
                        'temperature': ' ',
                        'bath coupling': self.app.md_settings['bath_coupling'],
                        'shake solvent': self.app.md_settings['shake_solvent'],
                        'shake solute': self.app.md_settings['shake_solute'],
                        'shake hydrogens': self.app.md_settings['shake_hydrogens'],
                        'lrf': self.app.md_settings['lrf'],
                        'solute-solute cutoff': self.app.md_settings['solute_solute_cut'],
                        'solvent-solvent cutoff': self.app.md_settings['solvent_solvent_cut'],
                        'solute-solvent cutoff': self.app.md_settings['solute_solvent_cut'],
                        'Q-atoms cutoff': self.app.md_settings['q_atoms_cut'],
                        'lrf cutoff': self.app.md_settings['lrf_cut'],
                        'shell force': self.app.md_settings['shell_force'],
                        'radial force': self.app.md_settings['radial_force'],
                        'polarization restraint': self.app.md_settings['polarisation'],
                        'polarization force': self.app.md_settings['pol_force'],
                        'update non-bonded list': self.app.md_settings['nonbond_list'],
                        'update energy summary': self.app.md_settings['ene_summary'],
                        'write energy file': self.app.md_settings['ene_file'],
                        'write trajectory': self.app.md_settings['trajectory'],
                        'trajectory atoms': self.app.md_settings['trajectory atoms']}
            self.sequence_restraints = self.app.md_settings['seq_rest']
            self.atom_restraints = self.app.md_settings['atom_rest']
            self.distance_restraints = self.app.md_settings['dist_rest']
            self.wall_restraints = self.app.md_settings['wall_rest']

            #Get FEP states (use this in restraints setup)
            self.fep_states = self.app.evb_states.get()


        self.inputfiles_entry.delete(0,END)
        #Check if default MD settings exist, and use values from there:
        if os.path.isfile(self.app.settings_path + '/qmd_settings') and not self.fep:
            print 'file found'
            with open(self.app.settings_path + '/qmd_settings', 'r') as mdsettings:
                for line in mdsettings:
                    if '#' in line:
                        self.qmd_settings[line.split('#')[1].strip('\n')] = line.split('#')[0].strip()

        self.simtime_entry.insert(0,self.qmd_settings['simulation time'])
        self.inputfiles_entry.insert(0, self.qmd_settings['inputfiles'])

        self.stepsize_variable.insert(0,self.qmd_settings['stepsize'])
        if not self.fep:
            self.temperature.insert(0, self.qmd_settings['temperature'])
        self.nonbond_interval.insert(0, self.qmd_settings['update non-bonded list'])
        self.output.insert(0, self.qmd_settings['update energy summary'])
        self.energyfile.insert(0, self.qmd_settings['write energy file'])
        self.trajectory.insert(0, self.qmd_settings['write trajectory'])

        if self.fep:
            self.inputfiles_entry.config(state=DISABLED)
            self.temperature.config(state=DISABLED)
            atom1 = False
            for q in sorted(self.app.q_atom_nr.keys()):
                if atom1:
                    if int(self.app.q_atom_nr[q]) > (prev_atom + 1):
                        atom1 = False
                        self.qatom_list.insert(END, 'Q-atoms: %d - %d' % (atom_start, prev_atom))
                    elif q == sorted(self.app.q_atom_nr.keys())[-1]:
                        self.qatom_list.insert(END, 'Q-atoms: %d - %d' % (atom_start, self.app.q_atom_nr[q]))
                    else:
                        prev_atom = self.app.q_atom_nr[q]
                if not atom1:
                    atom_start = self.app.q_atom_nr[q]
                    atom1 = True
                    prev_atom = atom_start
            self.qatom_add_button.config(state=DISABLED)
            self.qatom_select_button.config(state=DISABLED)


        self.bath_entry.insert(0, self.qmd_settings['bath coupling'])
        self.shake_solvent.set(int(self.qmd_settings['shake solvent']))
        self.shake_solute.set(int(self.qmd_settings['shake solute']))
        self.shake_hydrogens.set(int(self.qmd_settings['shake hydrogens']))
        self.lrf_check.set(int(self.qmd_settings['lrf']))
        self.polarisation_check.set(int(self.qmd_settings['polarization restraint']))
        self.solute_solute.insert(0, self.qmd_settings['solute-solute cutoff'])
        self.solvent_solvent.insert(0,self.qmd_settings['solvent-solvent cutoff'])
        self.solute_solvent.insert(0, self.qmd_settings['solute-solvent cutoff'])
        self.qatom_cutoff.insert(0, self.qmd_settings['Q-atoms cutoff'])
        self.shell_force.insert(0, self.qmd_settings['shell force'])

        if self.qmd_settings['trajectory atoms'] == 'all':
            self.trj_atom1.insert(0,'%d' % atomnumbers[0])
            self.trj_atom2.insert(0,'%d' % atomnumbers[-1])
        else:
            self.trj_atom1.insert(0, 'included')
            self.trj_atom2.insert(0, 'included')

        self.set_file_steps()
        self.set_total_simtime()
        self.set_total_steps()
        self.check_lrf()
        self.check_polarisation()

        self.q_atoms.set(0)

        self.use_eq.set(1)

        if not self.run_md:
            self.run_button.grid_forget()
            if not self.fep:
                self.q_atoms.set(1)
            self.q_atom_check.config(state=DISABLED)


        self.check_qatoms()
        self.check_topology()

    def check_topology(self):
        """
        Checks that pdb file and topology matches.
        Get solvent radius and insert as default
        """
        topology = open(self.topology, 'r').readlines()
        pdb = open(self.pdbfile, 'r').readlines()

        #Atoms in pdb file
        atoms_pdb = 0

        #atoms defined in topology
        atoms_top = 0

        #solvent radius from topology
        solv_r = 0

        #Find number of atoms in pdb file:
        for line in pdb:
            try:
                if 'ATOM' in line.split()[0] or 'HETATM' in line.split()[1]:
                    atoms_pdb += 1
            except:
                pass

        # Find number of atoms defined in topology and solvent radius
        for i in range(len(topology)):
            if 'END        of header' in topology[i]:
                atoms_top = topology[i + 1].split()[0]
            if 'solvent type' in topology[i]:
                try:
                    solv_r = topology[i + 1].split()[0]
                except:
                    pass

        #Check if atom numbers match:
        if int(atoms_top) != int(atoms_pdb):
            self.app.errorBox('Warning', 'Found %d atoms in %s, but %s specifies %d atoms.\n'
                                         'This is not important for the simulation, but you can not'
                                         'use %s for visualizing the final trajectory because it does not'
                                         'match the topology.' %
                                         (int(atoms_pdb), self.pdbfile.split('/')[-1], self.topology.split('/')[-1],
                                          int(atoms_top), self.pdbfile.split('/')[-1]))

        #If solvent radius in top, insert it. If not, use radius of system - 3.
        if solv_r != 0:
            self.shell_radius.insert(0, '%.1f' % (float(float(solv_r) * 0.85)))
            print solv_r
        else:
            system_radius = pt.findRadius(self.pdbfile)
            self.shell_radius.insert(0, '%d' % int(round(system_radius - 3)))

    def writeInputs(self, message = True):
        """
        Write eq and md input files
        self.q_settings = ['workdir',[prm],[lib],[equilibration],[sub (script[1/0, cmd)], [executables]]
        """
        #Control that Q-atoms is selected!
        if not self.run_md and not self.fep:
            if len(self.qatoms_list) == 0:
                self.app.errorBox('Warning', 'No Q-atoms specified! Please select Q-atoms.')
                return

        base_name = self.pdbfile.split('/')[-1].split('.')[0]

        # Create FEP file:
        if not self.fep:
            fepname = self.app.workdir+'/'+base_name+'.fep'
            fepfile = open(fepname,'w')
            fepfile.write('[fep]\nstates 1\n')
            fepfile.write('[atoms]\n')
            if self.q_atoms.get() == 1:
                if len(self.qatoms_list) > 0:
                    for i in range(len(self.qatoms_list)):
                        fepfile.write('%3d  %6s\n' % (i + 1, self.qatoms_list[i]))
            fepfile.close()

        #If EVB or FEP, return MD settings to app and close:
        if self.fep:
            self.app.md_settings['simtime'] = float(self.simtime_entry.get())
            self.app.md_settings['stepsize'] = float(self.stepsize_variable.get())
            self.app.md_settings['bath_coupling'] = float(self.bath_entry.get())
            self.app.md_settings['shake_solvent'] = self.shake_solvent.get()
            self.app.md_settings['shake_solute'] = self.shake_solute.get()
            self.app.md_settings['shake_hydrogens'] = self.shake_hydrogens.get()
            self.app.md_settings['lrf'] = self.lrf_check.get()
            self.app.md_settings['solute_solute_cut'] = float(self.solute_solute.get())
            self.app.md_settings['solute_solvent_cut'] = float(self.solute_solvent.get())
            self.app.md_settings['solvent_solvent_cut'] = float(self.solvent_solvent.get())
            self.app.md_settings['q_atoms_cut'] = float(self.qatom_cutoff.get())
            self.app.md_settings['lrf_cut'] = float(self.lrf_expansion.get())
            self.app.md_settings['shell_force'] = float(self.shell_force.get())
            self.app.md_settings['shell_rad'] = float(self.shell_radius.get())
            self.app.md_settings['radial_force'] = float(self.radial_force.get())
            self.app.md_settings['polarisation'] = self.polarisation_check.get()
            self.app.md_settings['pol_force'] = float(self.polarisation.get())
            self.app.md_settings['nonbond_list'] = self.nonbond_interval.get()
            self.app.md_settings['ene_summary'] = self.output.get()
            self.app.md_settings['ene_file'] = self.energyfile.get()
            self.app.md_settings['trajectory'] = self.trajectory.get()
            if self.trj_atom1.get() == 'included' or self.trj_atom1.get() == 'not excluded':
                self.app.md_settings['trajectory atoms'] = 'not excluded'
            else:
                self.app.md_settings['trajectory atoms'] = '%6s %6s' % (self.trj_atom1.get(), self.trj_atom2.get())
            self.app.md_settings['seq_rest'] = self.sequence_restraints
            self.app.md_settings['atom_rest'] = self.atom_restraints
            self.app.md_settings['dist_rest'] = self.distance_restraints
            self.app.md_settings['wall_rest'] = self.wall_restraints

            self.app.log('info', 'EVB MD settings configured.')
            self.app.update_status()
            self.app.update_temperatures()
            self.destroy()
            return

        #Get global settings:
        q_settings = cPickle.load(open(self.app.settings_path + '/Qsettings','rb'))
        output_int = self.output.get()
        trj_int = self.trajectory.get()
        ene_int = self.energyfile.get()
        non_bond_int = self.nonbond_interval.get()
        radial_force = self.radial_force.get()
        polarisation = self.polarisation.get()
        shell_force = self.shell_force.get()
        shell_radius = self.shell_radius.get()
        md_steps = self.file_steps.get(0.0, END).strip()
        md_stepsize = self.stepsize_variable.get().strip()
        md_temp = self.temperature.get().strip()
        md_bath = self.bath_entry.get()
        if self.shake_solvent.get() == 1:
            shake_solvent = 'on'
        else:
            shake_solvent = 'off'
        if self.lrf_check.get() == 1:
            lrf = 'on'
            lrf_cutoff = self.lrf_expansion.get()
        else:
            lrf = 'off'
        if self.polarisation_check.get() == 1:
            use_pol = 'on'
        else:
            use_pol = 'off'

        #Find atom nrs for solute and solvent:
        pdbfile = open(self.pdbfile,'r').readlines()
        solute = []
        solvent = []
        all_atoms = []
        for line in pdbfile:
            if 'ATOM' in line:
                if 'HOH' in line or 'SPC' in line:
                    solvent.append(line.split()[1])
                else:
                    solute.append(line.split()[1])
                all_atoms.append(line.split()[1])


        #Open submission script:
        self.submit = base_name + 'run.sh'
        self.submitcommand = q_settings[ 'subscript' ][1]
        submitfile = open(self.app.workdir + '/' + self.submit,'w')
        if int(q_settings[ 'subscript' ][0]) == 1:
            if os.path.isfile(self.app.settings_path + '/qsubmit'):
                submissionscipt = open(self.app.settings_path + '/qsubmit','r').readlines()
            elif os.path.isfile(self.app.workdir + '/' + 'qsubmit'):
                submissionscipt = open(self.app.workdir + '/' + 'qsubmit','r').readlines()
            else:
                submissionscipt = ['#!/bin/bash\n#Qdyn I/O\n']
                print 'submission script not found! Please edit this in settings'
            for line in submissionscipt:
                if '#Qdyn I/O' in line:
                    break
                else:
                    submitfile.write(line)

        #Check if Qdyn is MPI run or not:
        qdyn = q_settings[ 'executables' ][1]
        if qdyn[-1] == 'p':
            qdyn = 'mpirun %s' % qdyn



        #Write defeault equilibration procedure
        if self.use_eq.get() == 1:
            self.app.log('info','Writing MD equilibration input files ...')
            random_seed = random.randrange(1000, 9999)
            count = 0
            for i in range(len(q_settings[ 'equilibration' ])):
                count += 1
                inputfile = self.app.workdir + '/' + base_name + '_eq%d.inp' % count
                logfile = self.app.workdir + '/' + base_name + '_eq%d.log' % count
                eq_file = open(inputfile, 'w')

                if q_settings[ 'equilibration' ][i][0] == 'End':
                    temp = self.temperature.get()
                else:
                    temp = q_settings[ 'equilibration' ][i][0]
                if int(q_settings[ 'subscript' ][0]) == 1:
                    inputfile = base_name + '_eq%d.inp' % count
                    logfile = base_name + '_eq%d.log' % count

                submitfile.write('%s %s > %s\n' % (qdyn, inputfile, logfile))

                #self.app.main_window.update_txt('%s\n' % inputfile)
                self.app.log(' ','%s\n' % inputfile)
                eq_file.write('[MD]\n')
                eq_file.write('%25s %s\n' % ('steps'.ljust(25), q_settings[ 'equilibration' ][i][5]))
                eq_file.write('%25s %s\n' % ('stepsize'.ljust(25), q_settings[ 'equilibration' ][i][4]))
                eq_file.write('%25s %s\n' % ('temperature'.ljust(25), temp))
                eq_file.write('%25s %s\n' % ('bath_coupling'.ljust(25), q_settings[ 'equilibration' ][i][1]))
                if count == 1:
                    eq_file.write('%25s %d\n' % ('random_seed'.ljust(25), random_seed))
                    eq_file.write('%25s %s\n' % ('initial_temperature'.ljust(25), q_settings[ 'equilibration' ][i][0]))
                    eq_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), 'on'))
                if count > 1:
                    eq_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), shake_solvent))
                if self.shake_hydrogens.get() == 1:
                    eq_file.write('%25s %s\n' % ('shake_hydrogens'.ljust(25), 'on'))
                if self.shake_solute.get() == 1:
                    eq_file.write('%25s %s\n' % ('shake_solute'.ljust(25), 'on'))
                eq_file.write('%25s %s\n' % ('lrf'.ljust(25), lrf))

                eq_file.write('\n[cut-offs]\n')
                eq_file.write('%25s %s\n' % ('solute_solvent'.ljust(25), self.solute_solute.get()))
                eq_file.write('%25s %s\n' % ('solute_solute'.ljust(25), self.solute_solute.get()))
                eq_file.write('%25s %s\n' % ('solvent_solvent'.ljust(25), self.solvent_solvent.get()))
                eq_file.write('%25s %s\n' % ('q_atom'.ljust(25), self.qatom_cutoff.get()))
                if lrf == 'on':
                    eq_file.write('%25s %s\n' % ('lrf'.ljust(25).ljust(25), lrf_cutoff))

                eq_file.write('\n[sphere]\n')
                eq_file.write('%25s %s\n' % ('shell_force'.ljust(25), shell_force))
                eq_file.write('%25s %s\n' % ('shell_radius'.ljust(25), shell_radius))

                eq_file.write(('\n[solvent]\n'))
                eq_file.write('%25s %s\n' % ('radial_force'.ljust(25), radial_force))
                eq_file.write(('%25s %s\n') % ('polarisation'.ljust(25), use_pol))
                if use_pol == 'on':
                    eq_file.write('%25s %s\n' % ('polarisation_force'.ljust(25), polarisation))

                eq_file.write('\n[intervals]\n')
                eq_file.write('%25s %s\n' % ('output'.ljust(25), output_int))
                eq_file.write('%25s %s\n' % ('trajectory'.ljust(25), trj_int))
                eq_file.write('%25s %s\n' % ('non_bond'.ljust(25), non_bond_int))

                eq_file.write('\n[files]\n')
                if int(q_settings[ 'subscript' ][0]) == 1:
                    topology = self.topology.split('/')[-1]
                    fepname = fepname.split('/')[-1]
                else:
                    topology = self.topology
                eq_file.write('%25s %s\n' % ('topology'.ljust(25), topology))
                eq_file.write('%25s %s_eq%d.dcd\n' % ('trajectory'.ljust(25), base_name, count))
                if count != 1:
                    eq_file.write('%25s %s\n' % ('restart'.ljust(25), self.restart_file))
                eq_file.write('%25s %s_eq%d.re\n' % ('final'.ljust(25), base_name, count))
                eq_file.write('%25s %s\n' % ('fep'.ljust(25), fepname))
                self.restart_file = '%s_eq%d.re' % (base_name, count)

                eq_file.write('\n[trajectory_atoms]\n')
                trj_atoms = self.trj_atom2.get()
                if trj_atoms == 'included' or trj_atoms == 'not excluded':
                    eq_file.write('not excluded\n')
                else:
                    eq_file.write('%6s  %6s\n' % (self.trj_atom1.get(), self.trj_atom2.get()))

                if self.fep:
                    eq_file.write('\n[fep]\n')


                eq_file.write('\n[sequence_restraints]\n')
                force = float(q_settings[ 'equilibration' ][i][3])
                if q_settings[ 'equilibration' ][i][2] != 'None':
                    if q_settings[ 'equilibration' ][i][2] == 'All':
                        atomlist = all_atoms
                        #not excluded does not work for this section. TODO! (implement in Q)
                        #eq_file.write(' not excluded  %4.1f 0  0\n' % force)
                    elif q_settings[ 'equilibration' ][i][2] == 'Solute':
                        atomlist = solute
                    elif q_settings[ 'equilibration' ][i][2] == 'Solvent':
                        atomlist = solvent

                    try:
                        atom_i = atomlist[0]
                        atom_j = atomlist[-1]

                        eq_file.write('%6s %6s %4.1f 0  0\n' % (atom_i, atom_j, force))
                    except:
                        continue
                if len(self.sequence_restraints) > 0:
                    for restraint in self.sequence_restraints:
                        eq_file.write(restraint)

                if len(self.distance_restraints) > 0:
                    eq_file.write('\n[distance_restraints]\n')
                    for restraint in self.distance_restraints:
                        eq_file.write('%s\n' % restraint)

                if len(self.atom_restraints) > 0:
                    eq_file.write('\n[atom_restraints]\n')
                    for restraint in self.atom_restraints:
                        eq_file.write('%s\n' % restraint)

                if len(self.wall_restraints) > 0:
                    eq_file.write('\n[wall_restraints]\n')
                    for restraint in self.wall_restraints:
                        eq_file.write('%s\n' % restraint)

        mdfiles = int(self.inputfiles_entry.get())
        if mdfiles == 1:
            self.app.log('info','Writing MD input file ...')
        else:
            self.app.log('info','Writing MD input files ...')

        count = 0
        for i in range(mdfiles):
            count += 1
            inputfile = self.app.workdir + '/' + base_name + '_md%03d.inp' % count
            logfile = self.app.workdir + '/' + base_name + '_md%03d.log' % count
            md_file = open(inputfile, 'w')
            if int(q_settings[ 'subscript' ][0]) == 1:
                inputfile = base_name + '_md%03d.inp' % count
                logfile = base_name + '_md%03d.log' % count

            submitfile.write('%s %s > %s\n' % (qdyn, inputfile, logfile))

            self.app.log('','%s\n' % inputfile)
            md_file.write('[MD]\n')
            md_file.write('%25s %s\n' % ('steps'.ljust(25), md_steps))
            md_file.write('%25s %s\n' % ('stepsize'.ljust(25), md_stepsize))
            md_file.write('%25s %s\n' % ('temperature'.ljust(25), md_temp))
            md_file.write('%25s %s\n' % ('bath_coupling'.ljust(25), md_bath))
            md_file.write('%25s %s\n' % ('shake_solvent'.ljust(25), shake_solvent))
            if self.shake_hydrogens.get() == 1:
                md_file.write('%25s %s\n' % ('shake_hydrogens'.ljust(25), 'on'))
            if self.shake_solute.get() == 1:
                md_file.write('%25s %s\n' % ('shake_solute'.ljust(25), 'on'))
            md_file.write('%25s %s\n' % ('lrf'.ljust(25), lrf))

            md_file.write('\n[cut-offs]\n')
            md_file.write('%25s %s\n' % ('solute_solvent'.ljust(25), self.solute_solvent.get()))
            md_file.write('%25s %s\n' % ('solute_solute'.ljust(25), self.solute_solute.get()))
            md_file.write('%25s %s\n' % ('solvent_solvent'.ljust(25), self.solvent_solvent.get()))
            md_file.write('%25s %s\n' % ('q_atom'.ljust(25), self.qatom_cutoff.get()))
            if lrf == 'on':
                md_file.write('%25s %s\n' % ('lrf'.ljust(25).ljust(25), lrf_cutoff))

            md_file.write('\n[sphere]\n')
            md_file.write('%25s %s\n' % ('shell_force'.ljust(25), shell_force))
            md_file.write('%25s %s\n' % ('shell_radius'.ljust(25), shell_radius))

            md_file.write(('\n[solvent]\n'))
            md_file.write('%25s %s\n' % ('radial_force'.ljust(25), radial_force))
            md_file.write(('%25s %s\n') % ('polarisation'.ljust(25), use_pol))
            if use_pol == 'on':
                md_file.write('%25s %s\n' % ('polarisation_force'.ljust(25), polarisation))

            md_file.write('\n[intervals]\n')
            md_file.write('%25s %s\n' % ('output'.ljust(25), output_int))
            md_file.write('%25s %s\n' % ('energy'.ljust(25), ene_int))
            md_file.write('%25s %s\n' % ('trajectory'.ljust(25), trj_int))
            md_file.write('%25s %s\n' % ('non_bond'.ljust(25), non_bond_int))

            md_file.write('\n[files]\n')
            if int(q_settings[ 'subscript' ][0]) == 1:
                topology = self.topology.split('/')[-1]
                fepname = fepname.split('/')[-1]
            else:
                topology = self.topology
            md_file.write('%25s %s\n' % ('topology'.ljust(25), topology))
            md_file.write('%25s %s_md%03d.dcd\n' % ('trajectory'.ljust(25), base_name, count))
            md_file.write('%25s %s\n' % ('restart'.ljust(25), self.restart_file))
            md_file.write('%25s %s_md%03d.en\n' % ('energy'.ljust(25), base_name, count))
            md_file.write('%25s %s_md%03d.re\n' % ('final'.ljust(25), base_name, count))
            md_file.write('%25s %s\n' % ('fep'.ljust(25), fepname))
            self.restart_file = '%s_md%03d.re\n' % (base_name, count)

            md_file.write('\n[trajectory_atoms]\n')
            trj_atoms = self.trj_atom2.get()
            if trj_atoms == 'included' or trj_atoms == 'not excluded':
                md_file.write('not excluded\n')
            else:
                md_file.write('%6s  %6s\n' % (self.trj_atom1.get(), self.trj_atom2.get()))

            if self.fep:
                md_file.write('\n[fep]\n')

            if len(self.sequence_restraints) > 0:
                md_file.write('\n[sequence_restraints]\n')
                for restraint in self.sequence_restraints:
                    md_file.write(restraint)

            if len(self.distance_restraints) > 0:
                md_file.write('\n[distance_restraints]\n')
                for restraint in self.distance_restraints:
                    md_file.write('%s\n' % restraint)

            if len(self.atom_restraints) > 0:
                md_file.write('\n[atom_restraints]\n')
                for restraint in self.atom_restraints:
                    md_file.write('%s\n' % restraint)

            if len(self.wall_restraints) > 0:
                md_file.write('\n[wall_restraints]\n')
                for restraint in self.wall_restraints:
                    md_file.write('%s\n' % restraint)
            md_file.close()

        #If use submission script, check for end statements (comes after #Qdyn I/O):
        if int(q_settings[ 'subscript' ][0]) == 1:
            write_end = False
            submissionscipt = open(self.app.settings_path + '/qsubmit','r').readlines()
            for k in range(len(submissionscipt)):
                if '#Qdyn I/O' in submissionscipt[k]:
                    end_statements_start = k + 1
                    write_end = True
            if write_end:
                for line in range(end_statements_start, len(submissionscipt)):
                    submitfile.write(submissionscipt[line])
        submitfile.close()

        self.app.log('info','Submission script written to: %srun.sh' % base_name)

        if message:
            self.app.errorBox('Info','MD input files written!')

        if not self.run_md:
            if self.parent.add_title == 'complex':
                self.parent.complex_md = True
            elif self.parent.add_title == 'ligand':
                self.parent.ligand_md = True
            self.parent.update_progress()
            self.destroy()

    def submit_md(self):
        """
        Writes MD input files and submits the submission script *run.sh
        """
        self.writeInputs(False)
        os.chdir(self.app.workdir)
        tmpfile = open('.tmpfile', 'w')

        #os.system('%s %s' % (self.submitcommand, self.submit))
        call('%s %s' % (self.submitcommand, self.submit), shell=True, stdout=tmpfile, stderr=tmpfile)
        job_id = open(self.app.workdir + '/.tmpfile','r').readlines()
        self.app.log('info','Submitting MD jobs ...')
        for line in job_id:
            self.app.main_window.update_txt(line)

        self.app.log('info', 'Jobs sumbitted')
        self.app.errorBox('Info', 'Jobs submitted')

    def check_qatoms(self):
        """
        Check if Q-atoms are to be used or not (fep file)
        """
        check = self.q_atoms.get()
        if check == 0:
            self.qatom1_entry.config(state = NORMAL)
            self.qatom2_entry.config(state = NORMAL)
            self.qatom1_entry.delete(0,END)
            self.qatom2_entry.delete(0,END)
            self.qatom1_entry.insert(0,'NA')
            self.qatom2_entry.insert(0,'NA')
            self.qatom1_entry.config(state = DISABLED)
            self.qatom2_entry.config(state = DISABLED)
            self.qatom_add_button.config(state = DISABLED)
            self.qatom_remove_button.config(state=DISABLED)
            self.qatom_select_button.config(state = DISABLED)
            if not self.fep:
                self.qatom_list.delete(0,END)
                self.qatom_list.insert(0,'Q-atoms...')
            self.qatom_list.config(state=DISABLED)
            self.qatoms_list = []
        elif check != 0:
            self.qatom1_entry.config(state = NORMAL)
            self.qatom2_entry.config(state = NORMAL)
            self.qatom1_entry.delete(0,END)
            self.qatom2_entry.delete(0,END)
            self.qatom_add_button.config(state = NORMAL)
            self.qatom_select_button.config(state = NORMAL)
            self.qatom_remove_button.config(state=NORMAL)
            self.qatom_list.config(state = NORMAL)
            self.qatom_list.delete(0,END)

    def check_equilibration(self):
        check = self.use_eq.get()
        if check == 0:
            print 'Not using equilibration procedure'
            self.app.errorBox('info', 'Sorry, you have to use the equilibration procedure for now...')
            self.use_eq.set(1)
        elif check == 1:
            print 'Using default equilibration procedure'

    def check_lrf(self):
        check = self.lrf_check.get()
        if check == 0:
            self.lrf_expansion.config(state = NORMAL)
            self.lrf_expansion.delete(0, END)
            self.lrf_expansion.insert(0,'NA')
            self.lrf_expansion.config(state = DISABLED)
        elif check != 0:
            self.lrf_expansion.config(state = NORMAL)
            self.lrf_expansion.delete(0, END)
            self.lrf_expansion.insert(0, self.qmd_settings['lrf cutoff'])
    
    def check_polarisation(self):
        check = self.polarisation_check.get()
        if check == 0:
            self.polarisation.config(state = NORMAL)
            self.polarisation.delete(0, END)
            self.polarisation.insert(0,'NA')
            self.polarisation.config(state=DISABLED)
        elif check != 0:
            self.polarisation.config(state=NORMAL)
            self.polarisation.delete(0,END)
            self.polarisation.insert(0, self.qmd_settings['polarization force'])

    def inputfiles_written(self,*args):
        self.set_total_simtime()
        self.set_file_steps()
        self.set_total_steps()

    def set_total_simtime(self):
        try:
            simtime = float(self.simtime_entry.get())
            inputfiles = int(self.inputfiles_entry.get())
            tot_time = float(simtime * inputfiles)
        except:
            tot_time = 'NaN'
        self.total_time.config(state=NORMAL)
        self.total_time.delete(0.0, END)
        self.total_time.insert(0.0, tot_time)
        self.total_time.config(state = DISABLED)
        
    def set_total_steps(self):
        try:
            stepsize = float(self.stepsize_variable.get())
            simtime = float(self.simtime_entry.get())
            files = int(self.inputfiles_entry.get())
            tot_steps = int(round((simtime*1000000.00)/stepsize)*files)
        except:
            tot_steps = 'NaN'
        self.total_steps.config(state = NORMAL)
        self.total_steps.delete(0.0,END)
        self.total_steps.insert(0.0,tot_steps)
        self.total_steps.config(state = DISABLED)

    def set_file_steps(self):
        try:
            stepsize = float(self.stepsize_variable.get())
            simtime = float(self.simtime_entry.get())
            steps = int(round((simtime*1000000.00)/stepsize))
        except:
            steps = 'NaN'
        self.file_steps.config(state = NORMAL)
        self.file_steps.delete(0.0,END)
        self.file_steps.insert(0.0,steps)
        self.file_steps.config(state = DISABLED)

    def select_atoms(self, entry1, entry2):
        """
        opens dialog to select atoms and inserts first into entry 1 and last into entry 2
        """

        self.select_atomrange = AtomSelectRange(self,self.root, self.pdbfile, entry1, entry2)
        self.select_atomrange.configure(bg = self.main_color)
        self.select_atomrange.title('Select atoms')
        self.select_atomrange.resizable()
    
    def add_restraints(self, newtitle = '',restraintlist = []):
        """
        opens dialog to add restraints
        """
        if newtitle == 'sequence': 
            new_window = AddSequenceRestraints
        elif newtitle == 'atom':
            new_window = AddAtomRestraints
        elif newtitle == 'distance':
            new_window = AddDistanceRestraints
        elif newtitle == 'wall':
            new_window = AddWallRestraints
        self.restraints_add = new_window(self, self.root, self.pdbfile, newtitle, restraintlist)
        self.restraints_add.config(bg = self.main_color)
        self.restraints_add.title('Edit %s restraints' % newtitle)
        self.restraints_add.resizable()


    def add_qatoms(self):
        """
        Adds Q atoms to list (for use in fep file generation)
        """
        first_atom = self.qatom1_entry.get()
        last_atom = self.qatom2_entry.get()

        self.qatom1_entry.delete(0, END)
        self.qatom2_entry.delete(0, END)

        try:
            for atom in range(int(first_atom), int(last_atom) + 1):
                if atom not in self.qatoms_list:
                    self.qatoms_list.append(atom)
        except:
            return

        self.qatom_list.insert(END, 'Q-atoms: %5d - %5d' % (int(first_atom), int(last_atom)))
        print self.qatoms_list

    def remove_qatoms(self):
        """
        Removes selected Q-atoms from list
        """
        try:
            selected_index = int(self.qatom_list.curselection()[0])
            atom1,atom2 = self.qatom_list.get(selected_index).split()[1], self.qatom_list.get(selected_index).split()[3]
        except:
            return

        for atom in range(int(atom1),int(atom2)+1):
            if atom in self.qatoms_list:
                self.qatoms_list.remove(atom)

        self.qatom_list.delete(selected_index)
        print self.qatoms_list

    def dialog_box(self):
        """Defines the outlook of Setup MD window.
        Uses Frame-widget to define left and right side of the window
        and uses grid to organize widgets inside the Frames. """
       
        #self.title('Setup MD')
        self.config(background=self.main_color)

        # Define frames
        left_frame = Frame(self, bg = self.main_color,bd=1,relief=SUNKEN)
        left_frame.pack(side = LEFT, padx=10,pady=10, fill=Y)

        right_frame = Frame(self, bg = self.main_color,bd=1,relief=SUNKEN)
        right_frame.pack(side = LEFT, padx=10, pady=10, fill=Y)

        md_settings_label = Label(left_frame, text = 'General MD settings')
        md_settings_label.grid(row=0, column=0, columnspan=12, pady=(0,5))
        md_settings_label.config(bg=self.main_color)

        # Define elements in the left_frame
        simtime_label = Label(left_frame, text = 'Simulation time:')
        simtime_label.grid(row = 1, column = 0, sticky = 'w')
        simtime_label.config(background = self.main_color)

        #Sim. time 
        self.simtime_entry = Entry(left_frame,width = 7, highlightthickness = 0, relief = GROOVE, textvariable = self.simtimeVar)
        self.simtime_entry.grid(row = 1, column = 1, columnspan=2)
        
        ns_label = Label(left_frame, text = 'ns/file')
        ns_label.grid(row = 1, column = 3, columnspan = 2, sticky = 'w')
        ns_label.config(background = self.main_color)

        #Step size
        stepsize_label = Label(left_frame, text = 'Stepsize:')
        stepsize_label.grid(row = 2, rowspan = 2, column = 0,sticky = 'w')
        stepsize_label.config(background = self.main_color)

        self.stepsize_variable = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE, textvariable = self.stepsizeVar)
        self.stepsize_variable.grid(row = 2, rowspan = 2, column = 1, columnspan = 2)

        fs_label = Label(left_frame, text ='fs')
        fs_label.grid(row = 2, rowspan = 2, column = 3, sticky = 'w')
        fs_label.config(bg = self.main_color)

        steps_label = Label(left_frame, text = 'steps/file:')
        steps_label.grid(row=2, rowspan = 2, column = 3,columnspan=4, sticky = 'e')
        steps_label.config(bg = self.main_color)

        self.file_steps = Text(left_frame, width = 10, height = 1)
        self.file_steps.grid(row = 2, rowspan = 2, column = 8, columnspan = 3)
        self.file_steps.config(state = DISABLED, highlightthickness = 0, bg=self.main_color)


        inputfiles_label = Label(left_frame, text = 'Input files:')
        inputfiles_label.config(background = self.main_color)
        inputfiles_label.grid(row=4, rowspan = 2, column = 0,sticky = 'w')

        self.inputfiles_entry = Spinbox(left_frame, width = 4, from_=1, to=50,highlightthickness = 0, relief = GROOVE,
                                      textvariable=self.simfilesVar)
        self.inputfiles_entry.grid(row=4,rowspan = 2,column = 1, columnspan = 2)
        self.inputfiles_entry.config(highlightthickness = 0)

        
        total_time_label = Label(left_frame, text = 'Total sim. time:')
        total_time_label.grid(row = 4, column = 3, columnspan = 4, sticky = 'E')
        total_time_label.config(background = self.main_color)

        ns2_label = Label(left_frame, text = 'ns')
        ns2_label.grid(row = 4, column = 12)
        ns2_label.config(bg = self.main_color)

        total_steps_label = Label(left_frame, text = 'Total steps:')
        total_steps_label.grid(row = 5, column = 3, columnspan = 4, sticky = 'E')
        total_steps_label.config(background = self.main_color)

        self.total_time = Text(left_frame, width = 10, height = 1)
        self.total_time.grid(row = 4, column = 8, columnspan = 3)
        self.total_time.config(state=DISABLED, highlightthickness=0, bg=self.main_color)

        self.total_steps = Text(left_frame, width = 10, height = 1)
        self.total_steps.grid(row = 5, column = 8, columnspan =3)
        self.total_steps.config(state = DISABLED, highlightthickness = 0, bg=self.main_color)

        temp_label = Label(left_frame, text = 'Temperature:')
        temp_label.grid(row = 6, column = 0, pady = (10,0), sticky = 'w')
        temp_label.config(bg = self.main_color)

        self.temperature = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE)
        self.temperature.grid(row=6, column = 1, columnspan = 2, pady = (10,0))

        kelvin_label = Label(left_frame, text = 'K')
        kelvin_label.grid(row = 6, column = 3, pady = (10,0), sticky = 'w')
        kelvin_label.config(bg =self.main_color)

        bath_label = Label(left_frame, text = 'Bath coupling:')
        bath_label.grid(row = 6, column = 4, columnspan = 2, pady = (10,0), sticky ='E')
        bath_label.config(bg = self.main_color)

        self.bath_entry = Entry(left_frame, width = 7, highlightthickness = 0, relief = GROOVE)
        self.bath_entry.grid(row = 6, column = 8, columnspan = 3, pady = (10,0))

        shake_label = Label(left_frame, text = 'SHAKE:')
        shake_label.grid(row = 7, column = 0,sticky = 'w')
        shake_label.config(bg=self.main_color)

        solvent_label = Label(left_frame, text = 'Solvent')
        solvent_label.grid(row = 7, column = 2)
        solvent_label.config(bg = self.main_color)

        shake_solvent = Checkbutton(left_frame, variable = self.shake_solvent)
        shake_solvent.grid(row = 7, column = 3, sticky = 'w')
        shake_solvent.config(bg =self.main_color)

        solute_label = Label(left_frame, text = 'Solute')
        solute_label.grid(row = 7, column = 4, columnspan = 2, sticky = 'E')
        solute_label.config(bg = self.main_color)

        shake_solute = Checkbutton(left_frame, variable = self.shake_solute)
        shake_solute.grid(row=7, column = 7, sticky = 'w')
        shake_solute.config(bg=self.main_color)

        hydrogens_label = Label(left_frame, text = 'Hydrogens')
        hydrogens_label.grid(row=7, column = 8, columnspan = 3)
        hydrogens_label.config(bg = self.main_color)

        shake_hydrogens = Checkbutton(left_frame, variable = self.shake_hydrogens)
        shake_hydrogens.grid(row = 7, column = 12, sticky = 'w')
        shake_hydrogens.config(bg = self.main_color)

        lrf_label = Label(left_frame, text = 'LRF Taylor expansion:')
        lrf_label.grid(row = 8, column = 0, columnspan=3,sticky = 'w')
        lrf_label.config(bg = self.main_color)

        lrf_check = Checkbutton(left_frame, variable = self.lrf_check, command = self.check_lrf)
        lrf_check.grid(row=8, column = 3, sticky = 'w')
        lrf_check.config(bg = self.main_color)

        cutoff_label = Label(left_frame, text = 'Cut-offs for non-bonded interactions')
        cutoff_label.grid(row = 9, column = 0, columnspan = 12, pady = (10,5))
        cutoff_label.configure(bg = self.main_color)

        solute_solute_label = Label(left_frame, text = 'Solute-Solute:')
        solute_solute_label.grid(row=10, column = 0, sticky='w')
        solute_solute_label.config(bg=self.main_color)
	
        self.solute_solute = Entry(left_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.solute_solute.grid(row=10, column = 1, columnspan = 2,sticky = 'E')

        aa_label = Label(left_frame, text = '%s' % u'\xc5')
        aa_label.grid(row=10, column = 3, sticky = 'w')
        aa_label.config(bg = self.main_color)

        solvent_solvent_label = Label(left_frame, text = 'Solvent-Solvent:')
        solvent_solvent_label.grid(row=10, column = 5, columnspan = 3,sticky = 'E')
        solvent_solvent_label.config(bg =self.main_color)

        self.solvent_solvent = Entry(left_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.solvent_solvent.grid(row=10, column = 8, columnspan = 2,sticky = 'E')

        aa2_label = Label(left_frame, text = '%s' % u'\xc5')
        aa2_label.grid(row = 10, column = 10, sticky = 'w')
        aa2_label.config(bg = self.main_color)

        solute_solvent_label = Label(left_frame, text = 'Solute-Solvent:')
        solute_solvent_label.grid(row=11, column = 0, sticky='w')
        solute_solvent_label.config(bg=self.main_color)

        self.solute_solvent = Entry(left_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.solute_solvent.grid(row=11, column = 1, columnspan = 2,sticky = 'E')

        aa3_label = Label(left_frame, text = '%s' % u'\xc5')
        aa3_label.grid(row=11, column = 3, sticky = 'w')
        aa3_label.config(bg = self.main_color)

        lrf_expansion_label = Label(left_frame, text = 'Q-atom:')
        lrf_expansion_label.grid(row=11, column = 5, columnspan = 3,sticky = 'E')
        lrf_expansion_label.config(bg =self.main_color)

        self.qatom_cutoff = Entry(left_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.qatom_cutoff.grid(row=11, column = 8, columnspan = 2,sticky = 'E')

        aa5_label = Label(left_frame, text = '%s' % u'\xc5')
        aa5_label.grid(row = 11, column = 10, sticky = 'w')
        aa5_label.config(bg = self.main_color) 

        lrf_expansion_label = Label(left_frame, text = 'LRF expansion:')
        lrf_expansion_label.grid(row=12, column = 0, sticky = 'w')
        lrf_expansion_label.config(bg =self.main_color)

        self.lrf_expansion = Entry(left_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.lrf_expansion.grid(row=12, column = 1, columnspan = 2,sticky = 'E')
        self.lrf_expansion.config(state=DISABLED)        

        aa4_label = Label(left_frame, text = '%s' % u'\xc5')
        aa4_label.grid(row = 12, column = 3, sticky = 'w')
        aa4_label.config(bg = self.main_color)

        sphere_settings_label = Label(left_frame, text = 'Simulation sphere settings')
        sphere_settings_label.grid(row=13, column = 0, columnspan = 12, pady = (10,5))
        sphere_settings_label.config(bg=self.main_color)

        shell_force_label = Label(left_frame, text = 'Shell force:')
        shell_force_label.grid(row=14, column = 0, sticky='w')
        shell_force_label.config(bg=self.main_color)

        self.shell_force = Entry(left_frame, width = 5, highlightthickness = 0, relief = GROOVE)
        self.shell_force.grid(row=14, column = 1, columnspan = 2,sticky = 'E')


        shell_radius_label = Label(left_frame, text = 'Shell radius:')
        shell_radius_label.grid(row=14, column = 5, columnspan = 3,sticky = 'E')
        shell_radius_label.config(bg =self.main_color)

        self.shell_radius = Entry(left_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.shell_radius.grid(row=14, column = 8, columnspan = 2,sticky = 'E')

        aa5_label = Label(left_frame, text = '%s' % u'\xc5')
        aa5_label.grid(row = 14, column = 10, sticky = 'w')
        aa5_label.config(bg = self.main_color)

        solvent_sphere_label = Label(left_frame, text = 'Solvent sphere boundary settings')
        solvent_sphere_label.grid(row=15, column = 0, columnspan = 12, pady = (10,5))
        solvent_sphere_label.config(bg=self.main_color)

        radial_force_label = Label(left_frame, text = 'Radial force:')
        radial_force_label.grid(row=16, column = 0, sticky = 'w')
        radial_force_label.config(bg=self.main_color)

        self.radial_force = Entry(left_frame, width = 5, highlightthickness = 0, relief = GROOVE)
        self.radial_force.grid(row = 16, column = 1, columnspan = 2, sticky = 'E')
        self.radial_force.insert(0,'60.0')

        polarisation_label = Label(left_frame, text = 'Use polarisation restraint')
        polarisation_label.grid(row = 17, column = 0, columnspan = 2, sticky = 'E')
        polarisation_label.config(bg = self.main_color)


        polarisation_check = Checkbutton(left_frame, variable = self.polarisation_check, command = self.check_polarisation)
        polarisation_check.grid(row=17, column = 3, sticky = 'w')
        polarisation_check.config(bg = self.main_color)

        pol_restraint = Label(left_frame, text = 'Force:')
        pol_restraint.grid(row=17, column = 4, columnspan =2)
        pol_restraint.config(bg = self.main_color)

        self.polarisation = Entry(left_frame, width = 5, highlightthickness = 0, relief = GROOVE)
        self.polarisation.grid(row=17, column = 6, columnspan = 4, sticky = 'w')

        #RIGHT FRAMe
        recording_label = Label(right_frame, text = 'Output recording intervals')
        recording_label.grid(row = 0, column = 0, columnspan = 8, pady=(0,5))
        recording_label.config(bg=self.main_color)
  
        nonbond_label = Label(right_frame, text = 'Non-bonded list:')
        nonbond_label.grid(row = 1, column = 0, columnspan = 2,sticky = 'w')
        nonbond_label.config(bg = self.main_color)

        self.nonbond_interval = Entry(right_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.nonbond_interval.grid(row=1, column = 2, columnspan = 2)


        output_label = Label(right_frame, text = 'Energy summary:')
        output_label.grid(row=1, column = 5, columnspan = 2, sticky = 'w', padx=(10,0))
        output_label.config(bg = self.main_color)

        self.output = Entry(right_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.output.grid(row = 1, column = 7, columnspan = 2)

	
        energyfile_label = Label(right_frame, text = 'Energy file:')
        energyfile_label.grid(row = 2, column = 0, columnspan = 2, sticky = 'w')
        energyfile_label.config(bg=self.main_color)

        self.energyfile = Entry(right_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.energyfile.grid(row = 2, column = 2, columnspan = 2)


        trajectory_label = Label(right_frame, text = 'Trajectory:')
        trajectory_label.grid(row=2, column = 5, columnspan = 2, sticky = 'w', padx = (10,0))
        trajectory_label.config(bg=self.main_color)

        self.trajectory = Entry(right_frame, width = 4, highlightthickness = 0, relief = GROOVE)
        self.trajectory.grid(row=2, column = 7, columnspan = 2)


        trj_atoms_label = Label(right_frame, text = 'Trajectory atoms')
        trj_atoms_label.grid(row = 3, column = 0, columnspan = 8, pady=(10,5))
        trj_atoms_label.config(bg=self.main_color)

        atom1_label = Label(right_frame, text = 'First atom:')
        atom1_label.grid(row=4, column = 0, columnspan = 2,sticky = 'w')
        atom1_label.config(bg=self.main_color)

        self.trj_atom1 = Entry(right_frame, width = 8, highlightthickness = 0, relief = GROOVE)
        self.trj_atom1.grid(row=4, column=2, columnspan =2)
        	
        atom2_label = Label(right_frame, text = 'Last atom:')
        atom2_label.grid(row=4, column =4, columnspan = 2, sticky = 'w', padx =(10,0))
        atom2_label.config(bg=self.main_color)

        self.trj_atom2 = Entry(right_frame, width = 8, highlightthickness = 0, relief = GROOVE)
        self.trj_atom2.grid(row=4, column = 6, columnspan = 2)

        trj_select_button = Button(right_frame, text = 'Select', command=lambda: self.select_atoms(
            self.trj_atom1, self.trj_atom2))
        trj_select_button.grid(row = 4, column = 9)
        trj_select_button.config(highlightbackground = self.main_color)

        if self.run_md:
            q_atom_label = Label(right_frame, text = 'Q-atoms (optional)')
        else:
            q_atom_label = Label(right_frame, text = 'Q-atoms')
        q_atom_label.grid(row=5, column =0, columnspan = 8, sticky='nswe', pady=(10,5))
        q_atom_label.config(bg=self.main_color)

        self.q_atom_check = Checkbutton(right_frame, variable = self.q_atoms, command = self.check_qatoms)
        self.q_atom_check.grid(row=5, column = 9, sticky = 'w', pady=(10,5))
        self.q_atom_check.config(bg=self.main_color)

        qatom1_label = Label(right_frame, text = 'First atom:')
        qatom1_label.grid(row=6, column = 0, columnspan = 2, sticky = 'w')
        qatom1_label.config(bg=self.main_color)

        self.qatom1_entry = Entry(right_frame, width = 8, highlightthickness = 0, relief = GROOVE)
        self.qatom1_entry.grid(row = 6, column = 2, columnspan = 2)

        qatom2_label = Label(right_frame, text = 'Last atom:')
        qatom2_label.grid(row = 6, column = 4, columnspan = 2, sticky = 'w',padx =(10,0))
        qatom2_label.config(bg = self.main_color)

        self.qatom2_entry = Entry(right_frame, width = 8, highlightthickness = 0, relief = GROOVE)
        self.qatom2_entry.grid(row = 6, column = 6, columnspan = 2)

        self.qatom_select_button = Button(right_frame, text = 'Select', command = lambda: self.select_atoms(
            self.qatom1_entry, self.qatom2_entry))
        self.qatom_select_button.grid(row = 6, column = 9)
        self.qatom_select_button.config(highlightbackground = self.main_color)

        self.qatom_add_button = Button(right_frame, text = 'Add', command = self.add_qatoms)
        self.qatom_add_button.grid(row = 7, column =3, columnspan = 2, sticky = 'w')
        self.qatom_add_button.config(highlightbackground = self.main_color)

        self.qatom_remove_button = Button(right_frame, text = 'Delete', command = self.remove_qatoms)
        self.qatom_remove_button.grid(row=7, column = 5, columnspan = 2)
        self.qatom_remove_button.config(highlightbackgroun = self.main_color)

        self.qatom_list = Listbox(right_frame, height = 3, width = 30, relief = GROOVE)
        self.qatom_list.grid(row = 8, rowspan = 3, column = 0, columnspan = 7, sticky = 'e')
        self.qatom_list.config(highlightthickness = 0)

        scrollb = Scrollbar(right_frame, command = self.qatom_list.yview)
        scrollb.grid(row = 8, rowspan=3, column = 7, sticky = 'ns')
        self.qatom_list['yscrollcommand'] = scrollb.set

        restraints_label = Label(right_frame, text = 'Setup Restraints')
        restraints_label.grid(row = 11, column = 0, columnspan = 8, pady = (10,5))
        restraints_label.config(bg = self.main_color)

        sequence_restraints = Button(right_frame, text = 'Sequence', command=lambda: self.add_restraints(
                                    'sequence', self.sequence_restraints))
        sequence_restraints.grid(row = 12, column = 1, columnspan =2)
        sequence_restraints.config(highlightbackground = self.main_color)

        atom_restraints = Button(right_frame, text = 'Atom', command=lambda: self.add_restraints(
                                    'atom', self.atom_restraints))
        atom_restraints.grid(row = 12, column = 3, columnspan = 2)
        atom_restraints.config(highlightbackground = self.main_color)

        distance_restraints = Button(right_frame, text = 'Distance', command=lambda: self.add_restraints(
                                    'distance', self.distance_restraints))
        distance_restraints.grid(row = 12, column = 5, columnspan = 2, sticky = 'w')
        distance_restraints.config(highlightbackground = self.main_color)

        wall_restraints = Button(right_frame, text = 'Wall', command=lambda: self.add_restraints(
                                    'wall', self.wall_restraints))
        wall_restraints.grid(row = 12, column = 7, columnspan = 3, sticky = 'w')
        wall_restraints.config(highlightbackground = self.main_color)

        equilibration_label = Label(right_frame, text='Use default equilibration procedure')
        equilibration_label.grid(row=13, column = 0, columnspan = 8, pady=(10,0))
        equilibration_label.config(bg=self.main_color)

        eq_check = Checkbutton(right_frame, variable = self.use_eq, command = self.check_equilibration)
        eq_check.grid(row = 13, column = 9, sticky = 'w', pady=(10,5))
        eq_check.config(bg=self.main_color)

        write_button = Button(right_frame, text = 'Write', command=self.writeInputs)
        write_button.grid(row = 14, column = 3, columnspan =2, sticky = 'e', pady = (20,5))
        write_button.config(highlightbackground = self.main_color)

        self.run_button = Button(right_frame, text = 'Run', command=self.submit_md)
        self.run_button.grid(row = 14, column = 5, columnspan = 2, sticky = 'w', pady = (20,5))
        self.run_button.config(highlightbackground = self.main_color)

        cancel_button = Button(right_frame, text = 'Close ', command = self.destroy)
        cancel_button.grid(row = 14, column = 7, columnspan = 3, pady = (20,5))
        cancel_button.config(highlightbackground = self.main_color)










