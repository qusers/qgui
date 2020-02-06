Qgui (depricated)
================================================================================

Gui for Q, The Molecular Dynamics Package
--------------------------------------------------------------------------------

**Qgui** is a python2 package written by Geir Isaksen at the University of Tromso
to ease the interaction and submission of jobs to the specialized Molecular
Dynamics program **Q**. We say specialized as **Q** is tailored specifically
to tackle the challenge of predicting accurate binding free energies,
solvation free energies, and activation free energies since these are quantities
of particular interest as they are the direct result of thermodynamic or
kinetic experiments.


###Installation

**Qgui** is easy to install. The only demand on the user is that of making sure
of having a working version of **Q**, python, and the matplotlib and numpy python
packages used in the analysis modules.

The program can be cloned from github issuing the following command in the user
terminal:

```git
git clone https://github.com/qusers/qgui.git
```

This will ask for your github username and password.

As **Qgui** is written in python and uses its standard tkinter graphic libraries
no more dependencies need to be fullfiled apart from installing the free
for academics version of *Maestro* to be able to generate force-field libraries
and parameters for ligand modeling. A working installation of pymol also comes
handy for molecular visualization but it's not required.

Once downloaded the program can be invoked directly by typing the following
command in the terminal:

```bash
bash-3.2$ python qgui.py
```

There is also an installation script to automatically set up Qgui on your MAC/LINUX.
Simply run this by typing the following command in the terminal:
```bash
python2 INSTALL.py
``` 
or
```bash
./INSTALL.py
``` 
and follow the instructions given therin.

Alternatively, the user can also create an alias in their .cshrc or .bashrc file. For example
for the case of .bashrc:

```bash
alias qgui="python /Users/username/software/qgui/qgui.py"
```

In the previous command you will have to replace *username* by your own and
also give the correct path to where you have cloned **Qgui**. 


###Examples
The following are simple step-by-step examples of usage of **Qgui** from the
ground-up.

####n-butane
Now that you've installed **Qgui** and made sure that all additional packages
are installed we can start with a simple example.

First you will have to create a PDB (Protein Data Bank) formatted coordinate
file of n-butane, you can use *pymol* or *avogadro*, or download a coordinate file
online and convert it to the correct format. A popular molecule repository is
pubchem (https://pubchem.ncbi.nlm.nih.gov/). You can just go to pubchem and
use the search box to query for n-butane scroll-down to the 3D-Conformer
area and download the file in .sdf format. Later you can open the downloaded file
in pymol and save the file as .pdb. You will still have to do some editing
and delete all rows starting with CONNECT and replace HETATM with ATOM making
sure to include two spaces to replace the T, and M, letters in HETATM.
Finally you will have to replace UNK with LIG. Now you finally have a clean
file you can work with.

Next step is to go ahead and fire-up **Qgui** from a terminal with your
brand new alias:

```bash
qgui
```

Now you will have to make sure that the settings point to the right places.
In the upper side of your screen go to File/Settings and make sure that
the name of the executables for **Q** correspond to the ones in your system
(e.g. qprep5, qdyn5). Also set-up the path to the schrodinger installation,
usually it is automatically installed at:

```
/opt/schrodinger/suites2014-1
```

Also add the path to the oplsaa parameter and library files and finally click
on save.

As you will also want to have all files in the same folder holding your
*n-butane.pdb* file go to File/Change workdir and make sure to select the
folder holding your n-butane pdb file.

Now you can proceed to load your PDB file by clicking on Load, and Browse to
the folder where your file is located.

Once loaded you will want to generate parameters for you ligand. Go to
Prepare/Parameters, select the atoms in your pdb file and select a force-field.
Then click on Run. If succesful a new window will pop-up telling you that a
library file with a corresponding pdb, and parameter file has been generated.

Now go back to File/Setttings and load the LIG.lib library file and the
QOPLS2001_LIG.prm parameter file which have just been generated.
**Qgui** takes care of assembling a new parameter file which merges
the oplsaa parameters and the newly generated ligand parameters into a new file
called merged.prm. You must always check that file to make sure that the parameter
assignment makes chemical sense.

Now you can generate a topology file. Go to Prepare/Topology, make sure all
options are as you want them, for example, TIP3 water model or SPC, etc. Click
on Run to launch Qprep, or Write if you just want to write the input files.

Now you are ready to run simulations.

To run an MD simulation just go to Setup/MD and look at the options with
care. At this step if you have an experienced **Q** user close by it would
be useful so that someone can aid on understanding the meaning of the
options provided in the Setup MD window, of course, you can always turn
to the **Q** users manual.

After making sure that you've set-up the option most appropiate for your
simulation you can click on Run, and, if no error messages pop-up, wait for 
approximately 10 minutes running the non-parallel version of **qdyn5** on a 
3.3GHz intel core i7 processor with 16Gb of RAM. Alternatively, you can press
Write to write the files to disk, and submit them from command line.

The analysis tools in **Qgui** are quite powerful and customizable as they
rely on the matplotlib and numpy python libraries. Usage examples of the
analysis tools will soon be explored here.




