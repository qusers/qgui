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

from tkinter import Entry, EXTENDED, Text, Label, Frame, Button, Scrollbar, Toplevel, Checkbutton, Listbox, END, INSERT, GROOVE, RIDGE, TOP, IntVar, BROWSE


class PdbFilePrepare(Toplevel):
    """Implements a dialog-box when Prepare -> PDB is chosen from menubar.
    Has got methods to modify pdf file/structure such that it can be used by Q.
    PdbFilePrepare has got reference to the MainWindow-class."""

    def __init__(self, app, root, pdbfile, workdir):         #Receives app and root from Qgui-class.
        Toplevel.__init__(self, root)
        self.app = app
        self.main_color = self.app.main_color
        self.root = root
        self.pdbfile = pdbfile
        self.workdir = workdir
        self.savename = None            # This is the eventual custom chosen name for the finished PDB prepare file
        self.check_variable = IntVar() # This is used to get the state of the Checkbutton for removing waters

        self.app.log('info', 'PDB file prepare session started.')

        self.atoms_deletion = []
        self.atoms_deletion_final = []
        self.water_deletion_final = []
        self.final_deletion = []

        self.dialog_box()
        self.pdb_info()
        self.set_chain_number()
        self.set_residues_number()
        self.set_atoms_number()
        self.fill_chains_list()
        self.listbox_selected()


    def dialog_box(self):
        """Defines a dialog-box."""

        self.title('Prepare PDB') #TODO


        upper_frame = Frame(self, bg = self.main_color)
        upper_frame.pack(side = TOP)

        middel_frame = Frame(self, bg = self.main_color)
        middel_frame.pack(side = TOP)

        bottom_frame = Frame(self, bg = self.main_color)
        bottom_frame.pack(side = TOP)

        residues_label = Label(upper_frame, text = 'Residues total:')
        residues_label.grid(row = 0, column = 0, sticky = 'e', pady = 5)
        residues_label.config(background = self.main_color)
        self.residues_entry = Entry(upper_frame)
        self.residues_entry.grid(row = 0, column = 1, pady = 5)
        self.residues_entry.config(highlightthickness = 0, relief = GROOVE)
        self.residues_entry.insert(0, 'number of atoms')

        residues_water_label = Label(upper_frame, text = 'Water molecules:')
        residues_water_label.grid(row = 0, column = 2, sticky = 'e')
        residues_water_label.config(background = self.main_color)
        self.residues_water_entry = Entry(upper_frame)
        self.residues_water_entry.grid(row = 0, column = 3)
        self.residues_water_entry.config(highlightthickness = 0, relief = GROOVE)
        self.residues_water_entry.insert(0, 'number of waters')

        atoms_label = Label(upper_frame, text = 'Atoms total:')
        atoms_label.grid(row = 1, column = 0, sticky = 'e')
        atoms_label.config(background = self.main_color)
        self.atoms_entry = Entry(upper_frame)
        self.atoms_entry.grid(row = 1, column = 1)
        self.atoms_entry.config(highlightthickness = 0, relief = GROOVE)
        self.atoms_entry.insert(0, 'numer of atoms')

        chains_label = Label(upper_frame, text = 'Chains:')
        chains_label.grid(row = 1, column = 2, sticky = 'e')
        chains_label.config(background = self.main_color)
        self.chains_entry = Entry(upper_frame)
        self.chains_entry.grid(row = 1, column = 3)
        self.chains_entry.config(highlightthickness = 0, relief = GROOVE)
        self.chains_entry.insert(0, 'number of chains')

        delete_label = Label(upper_frame, text = 'Delete HOH:')
        delete_label.grid(row = 2, column = 0, sticky = 'e')
        delete_label.config(background = self.main_color)

        self.remove_checkbutton = Checkbutton(upper_frame, variable = self.check_variable, command = self.remove_water_checked)
        self.remove_checkbutton.config(background = self.main_color, relief = GROOVE)
        self.remove_checkbutton.grid(row  =2, column = 1)

        self.remove_label = Label(upper_frame, text = 'Check the checkbutton to remove all HOH.')
        self.remove_label.grid(row = 2, column = 2, columnspan = 2, sticky = 'w')
        self.remove_label.config(background = self.main_color)

        chains_label_list = Label(middel_frame, text = 'Chains:')
        chains_label_list.grid(row = 0, column = 0)
        chains_label_list.config(background = self.main_color)

        residues_label_list = Label(middel_frame, text = 'Residues:')
        residues_label_list.grid(row = 0, column = 2)
        residues_label_list.config(background = self.main_color)

        atoms_label_list = Label(middel_frame, text = 'Atoms:')
        atoms_label_list.grid(row = 0, column = 4)
        atoms_label_list.config(background = self.main_color)

        chains_scroll = Scrollbar(middel_frame)
        chains_scroll.grid(row = 1, column = 1, sticky = 'nsw')
        chains_scroll.config(bg = self.main_color, relief = GROOVE)
        self.chains_listbox = Listbox(middel_frame, yscrollcommand = chains_scroll.set, bd = 1, exportselection=0, selectmode = EXTENDED)
        self.chains_listbox.grid(row = 1, column = 0, columnspan = 1, padx = 5)
        self.chains_listbox.config(highlightthickness = 0)
        self.grid_rowconfigure(0, weight=1)
        chains_scroll.config(command = self.chains_listbox.yview)


        self.chains_listbox.bind('<<ListboxSelect>>', self.fill_residues_list)

        residues_scroll = Scrollbar(middel_frame)
        residues_scroll.grid(row = 1, column = 3, sticky = 'nsw')
        self.residues_listbox = Listbox(middel_frame, yscrollcommand = residues_scroll.set, bd = 1, exportselection=0, selectmode = EXTENDED)
        self.residues_listbox.grid(row = 1, column = 2, columnspan = 1, padx = 5)
        residues_scroll.config(command = self.residues_listbox.yview)
        self.residues_listbox.bind('<<ListboxSelect>>', self.fill_atoms_list)

        atoms_scroll = Scrollbar(middel_frame)
        atoms_scroll.grid(row = 1, column = 5, sticky = 'nsw', padx = 5)
        self.atoms_listbox = Listbox(middel_frame, yscrollcommand = atoms_scroll.set, bd = 1, selectmode = EXTENDED)
        self.atoms_listbox.grid(row = 1, column = 4, columnspan = 1, padx = 5)
        atoms_scroll.config(command = self.atoms_listbox.yview)

        self.add_selection_button = Button(middel_frame, text = 'Add to be deleted', command = self.add_selection_button_clicked)
        self.add_selection_button.grid(row = 2, column = 4)
        self.add_selection_button.config(highlightbackground = self.main_color, relief = GROOVE)

        label_txt = Label(bottom_frame, text = 'Items to be deleted:')
        label_txt.grid(row = 0, column = 0, sticky = 'w')
        label_txt.config(background = self.main_color)

        label_txt2 = Label(bottom_frame, text = 'Choose an item from one of the lists and click the button Add to be deleted.')
        label_txt2.grid(row = 3, column = 0, columnspan = 2, sticky = 'w')
        label_txt2.config(background = self.main_color)

        self.txt = Text(bottom_frame, height = 5)
        self.txt.config(font = ("consolas", 12), undo= True, wrap = 'word', highlightbackground = self.main_color, relief = RIDGE)
        self.txt.grid(row = 1, column = 0, columnspan = 2)

        create_button = Button(bottom_frame, text = 'Create PDB file', command = self.create_pdb_file_button_pressed)
        create_button.grid(row = 4, column = 1, sticky = 'w')
        create_button.config(highlightbackground = self.main_color, relief = GROOVE)

        create_label = Label(bottom_frame, text = 'Create a new PDB file where the selection has been deleted:')
        create_label.grid(row = 4, column = 0, sticky = 'w')
        create_label.config(background = self.main_color)

        cancel_button = Button(bottom_frame, text = 'Cancel', command = self.destroy)
        cancel_button.grid(row = 5, column = 1, sticky = 'w')
        cancel_button.config(highlightbackground = self.main_color, relief = GROOVE)

        self.clear_button = Button(bottom_frame, text = 'Clear selection', command = self.clear_selection_pressed)
        self.clear_button.grid(row = 5, column = 0, sticky = 'e')
        self.clear_button.config(highlightbackground = self.main_color, relief = GROOVE)


    def exit(self):
        pass

    def listbox_selected(self):
        atom_chosen = self.check_atom_chosen()
        residue_chosen = self.check_residue_chosen()
        chain_chosen = self.check_chain_chosen()

        if atom_chosen != ():
            self.residues_listbox.config(selectmode = BROWSE)
        else:
            self.residues_listbox.config(selectmode = EXTENDED)

        if residue_chosen != ():
            self.chains_listbox.config(selectmode = BROWSE)
        else:
            self.chains_listbox.config(selectmode = EXTENDED)

    def clear_selection_pressed(self):#
        """Clears the selection from the deletion txt-window.
        This is called when the self.clear_button is pressed."""
        # Clears up the lists for the final deletion
        self.atoms_deletion = []
        self.atoms_deletion_final = []
        self.water_deletion_final = []
        self.final_deletion = []

        # Delete text from the txt widget
        self.txt.delete(1.0, END)

        # Set the state=0 for HOH checkbutton
        self.check_variable.set(0)

    def set_chain_number(self):
        """This method sets the number of unique chains in self.chains list to
        self.chains_entry entry field. """

        old_chain = []
        temporary_chain = []

        for chain in self.chains:
            if chain not in old_chain:
                temporary_chain.append(chain)
                old_chain.append(chain)
            number_chains = len(temporary_chain)

        self.chains_entry.delete(0, 50) #Delete any entry from the entry field
        self.chains_entry.insert(0, number_chains) #Insert the length of the self.chains list

    def set_residues_number(self):
        """This method sets both the total number of residues and the number
        of water molecules in self.residues_entry field and
        self.residues_water_entry field.
         """
        whole_length = 0    #An empty variabel
        water_length = 0
        for chain in self.residues_chain: #Iterate over self.residues_chain
            length = len(chain)             #Check the length of every
            whole_length += length  #Add upp the length of all list in the list
            if 'HOH' in chain:
                length2 = len(chain)
                water_length += length2

        self.residues_entry.delete(0, 50)   #Delete any previous text from the entry
        self.residues_entry.insert(0, whole_length) #Insert the length in the entry

        self.residues_water_entry.delete(0, 50)
        self.residues_water_entry.insert(0, water_length)

    def set_atoms_number(self):
        """Sets the total number of atoms in the self.atoms_entry field. """
        self.atoms_entry.delete(0, 50)
        self.atoms_entry.insert(0, self.atom_count)

    def fill_chains_list(self):
        """Fills the self.chains_listbox with chains from the
        pdbfile."""

        self.chains_listbox.delete(0, 50)
        chain_count = 0
        for chain in self.chains:
            chain_count += 1
            self.chains_listbox.insert(END, 'Chain %d: %s' % (chain_count, chain))

    def fill_residues_list(self, index):
        """Gets the chosen chain in listbox from check_chain_chosen() method.
        Displays residues from this chain in residues listbox. """
        #self.chains_listbox.config(selectmode = BROWSE)
        chain = self.check_chain_chosen()
        index_list = []             # A list for indexes converted to integers
        for index in chain:
            index_chain = int(index)
            index_list.append(index_chain)

        self.atoms_listbox.delete(0, 100000)        #Deletes entry from atoms listbox when new chain is chosen

        chain_chosen_list = []
        chain_chosen_flat = []
        for index in index_list:
            chain_chosen = self.residues_chain[index]
            chain_chosen_list.append(chain_chosen)
            for chain in chain_chosen_list:
                result = isinstance(chain, list)
                if result == True:
                    chain_chosen_flat.extend(chain)

        residue_number_list = []
        residue_number_flat = []
        for index in index_list:
            residue_number = self.resnr_chain[index]
            residue_number_list.append(residue_number)
            for number_chain in residue_number_list:
                result = isinstance(number_chain, list)
                if result == True:
                    residue_number_flat.extend(number_chain)

        self.residues_listbox.delete(0, 100000)

        for residue, number in zip(chain_chosen_flat, residue_number_flat):
            if len(index_list) < 2:                                         # Fills up the residue list only if one chain is chosen, if several chains are chosen, no residues are displayed
                self.residues_listbox.insert(END, (residue, number))

        self.chain_chosen_flat = chain_chosen_flat


    def fill_atoms_list(self, index):
        """Gets the chosen residue from check_residue_chosen() method.
        Displays atoms from the chosen residue. """

        chain = self.check_chain_chosen()
        index_list_chain = []
        for index in chain:
            index_chain = int(index)
            index_list_chain.append(index_chain)

        residue = self.check_residue_chosen()
        index_list_residue = []
        for index in residue:
            residue_index = int(index)
            index_list_residue.append(residue_index)

        self.atoms_listbox.delete(0, 100000)

        atoms_flat = []
        for index_chain, index_residue in zip(index_list_chain, index_list_residue):
            atoms  = self.atomtypes_chain_residues[index_chain][index_residue]
            result = isinstance(atoms, list)
            if result == True:
                atoms_flat.extend(atoms)

        self.atoms_listbox.delete(0, 100000)

        for atom in atoms_flat:
            if len(index_list_residue) < 2:
                self.atoms_listbox.insert(END, atom)


    def check_chain_chosen(self):
        """Gets the current selection from the chains listbox. """
        return self.chains_listbox.curselection()


    def check_residue_chosen(self):
        """Gets the current selection from residues listbox. """
        return self.residues_listbox.curselection()


    def check_atom_chosen(self):
        """Gets the current selection from atoms listbox. """
        return self.atoms_listbox.curselection()


    def add_selection_button_clicked(self):
        """This method is called when the Add to be deleted-button is pressed.
        There are three possible outcomes:

        1. Only one chain is selected from a listbox and it is selected to be deleted from the pdbfile
        2. Both a chain and a residue is chosen from the
        listboxes which means that a residue is selected to be deleted
        3. Chain, residue and atom are chosen which means that an atom is selected to be deleted

        Makes a list of items to be deleted.

        The deletetion itself is invoked when Create a PDB file-button is pressed,
        look at create_pdb_file_button_pressed-method."""

        atom_chosen = self.check_atom_chosen()
        residue_chosen = self.check_residue_chosen()
        chain_chosen = self.check_chain_chosen()

        # Three possibilities:
        # 1.only chain chosen
        if (chain_chosen != ()) and (residue_chosen == ()) and (atom_chosen == ()):
            for index in chain_chosen:
                chain_index = int(index)
                chain_count = chain_index + 1
                self.txt.insert(INSERT, 'Chain %s:[%s]; ' % (chain_count, self.chains[chain_index]))
                self.atoms_deletion.append(self.atoms_chain_residues[chain_index])

        # 2.chain and residue chosen
        elif (chain_chosen != ()) and (residue_chosen != ()) and (atom_chosen == ()):
            for index_chain in chain_chosen:
                chain_index = int(index_chain)
                chain_count = chain_index + 1
                for index_residue in residue_chosen:
                    residue_index = int(index_residue)
                    residue_chosen = self.residues_chain[chain_index][residue_index]
                    self.txt.insert(INSERT, 'Residue %s: [%s] in Chain %s:[%s];' % (self.resnr_chain[chain_index][residue_index], residue_chosen, chain_count, self.chains[chain_index]))
                    self.atoms_deletion.append(self.atoms_chain_residues[chain_index][residue_index])

        # 3. all three chosen
        elif (chain_chosen != ()) and (residue_chosen != ()) and (atom_chosen != ()):
            for index in chain_chosen:
                chain_index = int(index)
                chain_count = chain_index + 1
                for index_residue in residue_chosen:
                    residue_index = int(index_residue)
                    residue_chosen = self.residues_chain[chain_index][residue_index]
                    for index_atom in atom_chosen:
                        atom_index = int(index_atom)
                        atom_chosen = self.atomtypes_chain_residues[chain_index][residue_index][atom_index]
                        self.txt.insert(INSERT, 'Atom: [%s] in Residue %s: [%s] in Chain %s:[%s]; ' % (atom_chosen, self.resnr_chain[chain_index][residue_index], residue_chosen, chain_count, self.chains[chain_index]))
                        self.atoms_deletion.append(self.atoms_chain_residues[chain_index][residue_index][atom_index])

        else:
            pass #TODO
            #error

        #The following creates a flatlist from the self.atoms_deletion array.
        #Creates a new list called self.atoms_deletion_final which is used
        # at delete_atoms() method.

        flatlist = []
        flatlist2 = []
        flatlist3 = []

        for entry in self.atoms_deletion:
            result = isinstance(entry, list)
            if result == True:
                flatlist.extend(entry)
                for entry2 in entry:
                    result = isinstance(entry2, list)
                    if result == True:
                        flatlist2.extend(entry2)
                        for entry3 in entry2:
                            result = isinstance(entry3, list)
                            if result == True:
                                flatlist3.extend(entry3)
                            elif result == False:
                                self.atoms_deletion_final.append(entry3)

                    elif result == False:
                        self.atoms_deletion_final.append(entry2)

            elif result == False:
                self.atoms_deletion_final.append(entry)


    def create_pdb_file_button_pressed(self):
        """This method is called when the Create PDB file-button is pressed.
        This method calls delete_atoms() method to make a new PDB file, excluding
        the atoms chosen to be deleted."""
        self.delete_atoms()

    def remove_water_checked(self):#TODO
        """This method is called when the Remove waters button is checked.
        It gets remove waters checkbuttons state which can be either 0(not checked)
        or 1(checked).

        When state==1 the index of HOH is retreived from self.residues_chain and
        these indexis are appended to delete_list.

        Entries from self.atoms_chain_residues which correspond to delete_list
        indexis are appended to yet an other list, self.water_deletion.
        This list is made flat and put in to variable self.water_deletion_final.

        self.water_deletion_final is used is delete_atoms() method."""
        state = self.check_variable.get()

        delete_list = []

        if state == 1:
            for index_chain, chain in enumerate(self.residues_chain):
                for index_residue, residue in enumerate(chain):
                    if residue == 'HOH':
                        delete_list.append((index_chain, index_residue))
            self.txt.insert(INSERT, 'Delete all HOH;')

        self.water_deletion = []

        for index_chain, index_residue in delete_list:
            self.water_deletion.append(self.atoms_chain_residues[index_chain][index_residue])

        for entry in self.water_deletion:
            result = isinstance(entry, list)
            if result == True:
                self.water_deletion_final.extend(entry)

        if state == 0:
            a = self.txt.get(0.0, END)
            self.txt.delete(0.0, END)

            self.txt.insert(0.0, (' ').join(a.split('Delete all HOH;')))




    def delete_atoms(self):
        """
        Takes a standard PDB file and deletes selected chains/residues/atoms.

        """
        finished = False
        # Two list are combined to one final deletion list. self.final_deletion is used
        # to tell which lines from the PDB file are not written to the new PDB file.
        self.final_deletion = (self.atoms_deletion_final + self.water_deletion_final)

        if len(self.final_deletion) == 0:
            self.app.log('info', 'No atoms specified. Please, choose chains, residues, atoms or HOH to be deleted.')
            self.app.pdb_delete_atoms_msg()
            return

        try:
            file = open(self.pdbfile,'r').readlines()
        except:
            self.app.log('info', 'Could not open PDB file')
            self.app.pdb_delete_atoms_msg2()
            return

        #Generate name for new pdb if not specified:
        if not self.savename:
            if not '_' in self.pdbfile.split('/')[-1]:
                self.savename = self.pdbfile.split('/')[-1].split('.')[0]+'_2_.pdb'
            else:
                try:
                    newnumber = int(self.pdbfile.split('/')[-1].split('_')[1]) + 1
                    self.savename = '%s_%d_%s' % (self.pdbfile.split('/')[-1].split('_')[0], newnumber,
                                                  self.pdbfile.split('_')[2])
                    if not self.savename.endswith('.pdb'):
                        self.savename += '.pdb'

                except:
                    self.savename = self.pdbfile.split('/')[-1].split('.')[0]+'_new.pdb'
        newpdb = []

        for line in file:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                current_atomnr = int(line[6:11])
                if current_atomnr not in self.final_deletion:
                    newpdb.append(line)
                    finished = True
        self.app.pdb_id = self.workdir + '/' + self.savename
        newfile = open(self.app.pdb_id,'w')
        for line in newpdb:
            newfile.write(line)

        if finished == True:
            self.app.log('info', '%s created from PDB file prepare' % self.savename)
            self.destroy()
            self.app.update_pdb_id_entryfield(self.savename)

    def pdb_info(self):
        """Returns list of chains, sequences, sequence numbers
        and total number of atoms.
        Takes self.pdbfile as input and returns chains, sequences, seqnr, atoms
        """

        try:
            file = open(self.pdbfile,'r').readlines()
        except:
            print('Please specify a valid PDB file')
        chains = []         #chains
        residues = []
        residues_chain = []         #sequences[chain] --> residues in each chain (Array)
        resnr = []
        resnr_chain = []        #seqnr[chain] --> residue nr in each chain (Array)
        atoms_chain_residues = []
        atoms_residues = []
        atoms = []          #atoms[chain][residue] --> Atom numbers
        atomtypes_chain_residues = []
        atomtypes_residues = []
        atomtypes = []      #atomtypes[chain][residue] --> Atom types

        chain = 'toFill'
        old_resnr = 0
        atom_count = 0

        for line in file:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                current_chain = line[21:22]
                current_residue = line[17:21].strip()
                current_resnr = int(line[22:26])
                current_atom = int(line[6:11])
                current_atomtype = line[11:17].strip()
                atom_count += 1
                if atom_count == 1:
                    chain = current_chain
                    old_resnr = current_resnr
                    chains.append(current_chain)
                    residues.append(current_residue)
                    resnr.append(current_resnr)
                if chain != current_chain:
                    chain = current_chain
                    chains.append(chain)
                    old_resnr = current_resnr
                    if len(residues) > 0:
                        residues_chain.append(residues)
                        residues = []
                        residues.append(current_residue)
                        resnr_chain.append(resnr)
                        resnr = []
                        resnr.append(current_resnr)
                        atoms_residues.append(atoms)
                        atoms = []
                        atoms_chain_residues.append(atoms_residues)
                        atoms_residues = []
                        atomtypes_residues.append(atomtypes)
                        atomtypes = []
                        atomtypes_chain_residues.append(atomtypes_residues)
                        atomtypes_residues = []

                        atoms.append(current_atom)
                        atomtypes.append(current_atomtype)
                elif current_resnr != old_resnr:
                    residues.append(current_residue)
                    resnr.append(current_resnr)
                    atoms_residues.append(atoms)
                    atomtypes_residues.append(atomtypes)
                    atoms = []
                    atomtypes = []
                    atoms.append(current_atom)
                    atomtypes.append(current_atomtype)
                    old_resnr = current_resnr
                elif current_resnr == old_resnr:
                    atoms.append(current_atom)
                    atomtypes.append(current_atomtype)
        if len(residues) > 0:
            residues_chain.append(residues)
            residues = []
            resnr_chain.append(resnr)
            resnr = []
            atoms_residues.append(atoms)
            atoms = []
            atoms_chain_residues.append(atoms_residues)
            atoms_residues = []
            atomtypes_residues.append(atomtypes)
            atomtypes = []
            atomtypes_chain_residues.append(atomtypes_residues)
            atomtypes_residues = []

        self.chains = chains
        self.residues_chain = residues_chain
        self.resnr_chain = resnr_chain
        self.atoms_chain_residues = atoms_chain_residues
        self.atomtypes_chain_residues = atomtypes_chain_residues
        self.atom_count = atom_count


        #return chains, residues_chain, resnr_chain, atoms_chain_residues, atomtypes_chain_residues, atom_count
