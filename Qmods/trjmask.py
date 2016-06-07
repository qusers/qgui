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

from Tkinter import Tk, NORMAL, DISABLED, Entry, Button, Grid, Label, END, Toplevel, Frame
from tkFileDialog   import askopenfilename, asksaveasfilename
import os
import time
import subprocess
import shutil

class TrjMask(Toplevel):
	def __init__(self, app, root):
		Toplevel.__init__(self, root)
		self.app = app
		self.root = root
		self.workdir=self.app.workdir
		self.main_color = self.app.main_color

		self.topname = None
		if self.app.top_id: #Inherits topology file from Qgui
			self.topname = self.app.top_id
		self.trjname = None
		self.pdbname = None
		self.app.log('info', 'Trajectory mask session started ...')
		self.makewindow()
	
	# Merges file name and extension. If name is in use the method adds "_new" to the name,
	# and calls itself recursively until unique name is generated. Then it returns the name.
	def newname(self,fname,fext):	
		if os.path.isfile(fname+'.'+fext) == True:
			return self.newname(fname+'_new',fext)
		else:
			return fname+'.'+fext

	def makeinp(self,top,trj,inp,out):	#Writes input file for Qprep.
		f = open(inp, 'w')
		f.write('rt %s\n' % top)
		f.write('trj %s\n' % trj)
		f.write('y\n')
		f.write('wp %s\n' % out)
		f.write('y\n')
		f.write('quit\n')
		f.close()

	def testlog(self,fname):	#Tries to open log file and scans for "PDB-written" string.
		try:
			f = open(fname, 'r')
			for line in f:
				if line == 'PDB file successfully written.\n':
					f.close()
					return True
			f.close()
			return False
		except:
			return False

	def cleanup(self,fname):	#Removes file, if it exists.
		if os.path.isfile(fname):
			os.remove(fname)
		
	def runqprep(self,inp,log):	#Launches Qprep with given input and creates logfile.
		logfile=open(log,'w')                                         
		tmp=subprocess.Popen(['Qprep5',inp], stdout=logfile)
		tmp.communicate()
		logfile.close()



	### Main Method
	def makemask(self,topsource,trjsource,out):
		inp = self.newname('mask','inp')		#Sets name for input file
		log = self.newname('mask','log')		#Sets name for log file
		top = self.newname('mask','top')		#Sets name for top file
		trj = self.newname('mask','dcd')		#Sets name for trj file
		shutil.copyfile(topsource,top)			#Copies top file to workdir
		shutil.copyfile(trjsource,trj)			#Copies trj file to workdir
		self.makeinp(top,trj,inp,out)			#Creates input file
		try:
			self.runqprep(inp,log)				#Tries to run Qprep
		except:
			self.app.log('info', 'Qprep5 failed to run.')					
		self.cleanup(inp)						#Removes input file
		self.cleanup(top)						#Removes top file
		self.cleanup(trj)						#Removes trj file
		if self.testlog(log): #Checks log file to see if PDB was written, then returns True/False
			self.cleanup(log) #Removes log file
			return True
		else:
			self.cleanup(log)
			return False


	
	#Opens file dialog for selecting top/trj files, then sets global file name.
	def getfile(self, whereto, ftype=[], *args):
		#print args #debug
		ftype.append(('All files','.*'))
		#print ftype #debug
		opts = {}
		opts['filetypes'] = ftype
		opts['parent'] = self
		opts['title'] = 'Select file'
		opts['initialdir'] = self.workdir
		name = askopenfilename(**opts)
		if name:
			if whereto == self.E1:
				self.topname = name
			else:
				self.trjname = name
				
			self.insertFilename(name, whereto)
	
	#Inserts file name into entry box.
	def insertFilename(self, name, whereto):
		whereto.config(state = NORMAL)
		whereto.delete(0, END)
		whereto.insert(0, name.split('/')[-1])
		whereto.config(state = DISABLED)
	
	#Opens file dialog to choose pdb file name, and then creates trjmask with given input.
	def buttonRun(self):
		if self.E1.get() and self.E2.get():
			opts = {}
			opts['defaultextension'] = '.pdb'
			opts['initialfile'] = 'mask.pdb'
			opts['title'] = 'Choose name for pdb file'
			opts['parent'] = self
			opts['initialdir'] = self.workdir
			name = None
			name = asksaveasfilename(**opts)
			if name:
				self.pdbname = name
				#Creates trjmask pdb and returns file holder to Qgui if successful.
				if self.makemask(self.topname,self.trjname,self.pdbname):
					self.app.pdb_id = self.pdbname
					self.app.main_window.set_entryfield(self.pdbname.split('/')[-1])
					self.app.top_id = self.topname
					self.app.main_window.set_topology_entryfield(self.topname.split('/')[-1])
					self.app.log('info', 'Trajectory mask created sucessfully.')
					self.destroy()
				else:
					self.app.log('info', 'Trajectory mask failed.')
		else:
			self.app.log('info', 'Please select files before creating PDB.')
	
	def buttonQuit(self): #Quits the application.
		self.app.log('info', 'Trajectory mask session ended.')
		self.destroy()  
		
	#Main window.
	def makewindow(self): 
		self.title('Make trj mask PDB')
		self.mainframe = Frame(self, bg=self.main_color)
		self.mainframe.pack(fill='both', padx=(10,10), pady=(10,10))
		
		#Entry line 1, Select topology
		L1 = Label(self.mainframe, text="Topology file:", width = 12, bg = self.main_color)
		L1.grid(row=0, column=0)		
		self.E1 = Entry(self.mainframe, bd = 5, state = DISABLED)
		self.E1.config(disabledbackground='white', disabledforeground='gray50', 
						highlightbackground=self.main_color)
		self.E1.grid(row=0, column=1)
		if self.topname: #If topology is inherited from Qgui it is inserted into entry box here
			self.insertFilename(self.topname, self.E1)
		self.buttonTop = Button(self.mainframe, text='Select', 
							command = lambda: self.getfile(self.E1,
							[('Topology files','.top')]))
		self.buttonTop.config(highlightbackground=self.main_color)
		self.buttonTop.grid(row=0, column=2)
		
		#Entry line 2, Select trajectory
		L2 = Label(self.mainframe, text="Trajectory file:", width=12, bg = self.main_color)
		L2.grid(row=1, column=0)			
		self.E2 = Entry(self.mainframe, bd = 5, state = DISABLED)
		self.E2.config(disabledbackground='white', disabledforeground='gray50', 
						highlightbackground=self.main_color)
		self.E2.grid(row=1, column=1)
		self.buttonTrj = Button(self.mainframe, text='Select', 
							command = lambda: self.getfile(self.E2,
							[('Trajectory files','.dcd')]))
		self.buttonTrj.config(highlightbackground=self.main_color)
		self.buttonTrj.grid(row=1, column=2)

		# Create and Quit buttons
		self.button1 = Button(self.mainframe, command = self.buttonRun, text='Create')
		self.button1.config(highlightbackground=self.main_color)
		self.button1.grid(row=2, column=0, columnspan=2)	
		self.button2 = Button(self.mainframe, command = self.buttonQuit, text='Quit')
		self.button2.config(highlightbackground=self.main_color)
		self.button2.grid(row=2, column=1, columnspan=2)	


