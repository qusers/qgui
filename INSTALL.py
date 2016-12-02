#!/usr/bin/env python2


###############EDITABLE REGION################################
#Path where Qgui will be installed. 
#Default is to create QGUI in /home/$user
#If you modify 'default', use the entire adress ('/Users/../')

install_path='default' 

#Pointer to bash file. This is to include Qgui in environment.
#Default is .bash_profile

bashfile='.bash_profile'

###############INSTALL SCRIPT###################################
import os
import shutil

if install_path == 'default':
    install_path='%s/QGUI' % os.path.expanduser("~")
elif install_path.split('/')[-1] != 'QGUI':
    install_path = '%s/QGUI' % install_path

bashfile='%s/%s' % (os.path.expanduser("~"), bashfile)

if not os.path.isdir(install_path):
    os.makedirs(install_path)
    print('Created Qgui directory: %s' % install_path)
else:
    print('Install path %s exists' % install_path)
    print('Existing files will be overwritten!\n')
    print('Press Y to continue or N to cancel:')
    if raw_input().lower() != 'y':
        sys.exit()

if not os.path.isdir('%s/Qmods' % install_path):
    os.makedirs('%s/Qmods' % install_path)

org_path = os.path.dirname(os.path.realpath(__file__))
shutil.copy2('%s/qgui.py' % org_path, '%s/Qgui' % install_path)

for f in os.listdir('Qmods'):
    shutil.copy2('%s/Qmods/%s' % (org_path, f), '%s/Qmods/%s' % (install_path, f)) 

if not os.path.isdir('%s/FF' % install_path):
    os.makedirs('%s/FF' % install_path)

for f in os.listdir('FF'):
    if os.path.exists('%s/FF/%s' % (install_path, f)):
        shutil.rmtree('%s/FF/%s' % (install_path, f))
    shutil.copytree('%s/FF/%s' % (org_path, f), '%s/FF/%s' % (install_path, f))

#Write to bash file so that Qgui can be started from command line.
if os.path.isfile(bashfile):
    oldfile = open(bashfile, 'r').readlines()
else:
    oldfile = list()
    print('%s created' % bashfile)

newfile = open(bashfile, 'w')

export_exist = False
export_line = 'PATH=$PATH:%s' % install_path

for line in oldfile:
    if export_line in line:
        export_exist = True
    elif export_line.split('/')[-1] in line:
        line = '\nPATH=$PATH:%s\n' % install_path
        export_exist = True
    newfile.write(line)

if not export_exist:
    newfile.write('\nPATH=$PATH:%s\n' % install_path)
    newfile.write('export PATH\n')

newfile.close()
print('%s updated' % bashfile)

print('\nInstallation successful!')
print('Restart terminal and type "Qgui" to launch')

