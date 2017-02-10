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

"""
Build.py
Masoud Kazemi Sep 2016
kazemimsod@gmail.com
This module provides the functions for:
    1) building by atoms/group 
    2) adjusting the bond, angle, and torsions by atoms/group

core functions: 
    rotate(p,rp,rv,theta),translate(p,dv,d)

auxillary functions:
    getfragment(fraglib,fragname)

user interface functions:
    BuildByAtom(flag,p1,p2=None,p3=None,p4=None,d=None,thetaA=None,thetaT=None)
    BuildByGroup(flag,p1,p2=None,p3=None,p4=None,d=None,thetaA=None,thetaT=None,fragname=None,fraglib=None,p1group=None)
    measure(v1,v2,v3=None,v4=None)

BuildByAtom() BuildByGroup() and measure() interact with external user 
"""

import numpy as np
import os 

def rotate(p,rp,rv,theta):
    """ rotate a point based on a given rotation vector (rv)
        that passes trough a rotation point (rp).
        input:  'rv'      rotation vector as numpy array
                'rp'      rotation point as numpy array
                'p'       the point to be rotated
                'theta'   the rotation angle in degrees
        output: 'q'      the transformed point """ 

    #the output point
    q = np.array((0.,0.,0.))

    #Convert the angle
    theta = np.deg2rad(theta)

    # making sure the rotation vector is normalized
    rv = rv/np.linalg.norm(rv)

    # set the values
    a, b, c = rp
    u, v, w = rv
    x, y, z = p

    # transformation based on the operation matrix T(-p1)T(-xz)T(-z)R(Z(theta))T(z)T(xz)T(p1)
    q[0] =  ((a*(v*v+w*w))-u*(b*v+c*w-u*x-v*y-w*z))*(1-np.cos(theta))+(x*np.cos(theta))+((-c*v+b*w-w*y+v*z)*np.sin(theta))
    q[1] =  ((b*(u*u+w*w))-v*(a*u+c*w-u*x-v*y-w*z))*(1-np.cos(theta))+(y*np.cos(theta))+((c*u-a*w+w*x-u*z)*np.sin(theta))
    q[2] =  ((c*(u*u+v*v))-w*(a*u+b*v-u*x-v*y-w*z))*(1-np.cos(theta))+(z*np.cos(theta))+((-b*u+a*v-v*x+u*y)*np.sin(theta))
    return q

def translate(p,dv,d):
    """ This function translate a point (p) according
        to a direction vector (dv) by a distance (d)
        input:      'dv'      rotation vector as numpy array
                    'p'       the point to be translated
                    'd'       translation distance in angstrom
        output:     'q'      the transformed point """

    #the output point
    q = np.array((0.,0.,0.))

    # making sure the rotation vector is normalized
    dv = dv/np.linalg.norm(dv)

    # tranlation
    q = np.array(p) + d*dv
    return q

def measure(v1,v2,v3=None,v4=None):
    """ This function measures:
        1) bond lenght if v1 and v2 are given
        2) angle if v1,v2,v3 are given. the v2 is the center atom
        3) torsion if v1,v2,v3,v4 are given
        input:  v1,v2,v3,v4 are vectors
        output: scaler"""
    #measure bond lenght
    if (v3 is None) and (v4 is None):
        return np.linalg.norm(np.array(v2) - np.array(v1))

    #measure v1-v2-v3 angle
    elif (v4 is None):
        d21 = np.array(v1) - np.array(v2)
        d23 = np.array(v3) - np.array(v2)
        d21 = d21/np.linalg.norm(d21)
        d23 = d23/np.linalg.norm(d23)
        d21Dotd23   = np.dot(d21,d23)
        d21crossd23 = np.cross(d21,d23)
        return np.rad2deg(np.arctan2(np.linalg.norm(d21crossd23),d21Dotd23))

    #measure v1-v2-v3-v4 torsion (rotation around x-v2-v3-x) in IUPAC 
    #convension (cis=0 and trans=180) between -180 to 180
    else:
        d12 = np.array(v2) - np.array(v1)
        d23 = np.array(v3) - np.array(v2)
        d34 = np.array(v4) - np.array(v3)
        n123 = np.cross(d12,d23)
        n234 = np.cross(d23,d34)
        n123 = n123/np.linalg.norm(n123)
        n234 = n234/np.linalg.norm(n234)
        n123Dotn234   = np.dot(n123,n234)
        n123crossn234 = np.cross(n123,n234)
        torsion = np.rad2deg(np.arctan2(np.linalg.norm(n123crossn234),n123Dotn234))

        # fix clockwise anticlockwise
        if (np.dot(n123,d34) >= 0):
            return torsion
        else:
            return -1*torsion

def getfragment(fraglib,fragname):
    """ parse the fragment library.
        input:  fraglib      fragment library path
                fragname     fragment name
        output: p1group      the dictionary of fragment
                                                           """
    p1group = dict()
    found = False
    if not os.path.isfile(fraglib):
        print 'No fragment library file was found'
        return None
    else:
        with open(fraglib, 'r') as f:
            for line in f:
                if len(line) : #ignor empty lines
                    if fragname in line: 
                        found = True
                        continue
                    elif found and "END" in line:
                        break
                    elif found:
                        field = line.split()
                        p1group.update({int(field[0]):{'name':field[1],'xyz':list(map(float,field[2:5])),'symbol':field[5]}})
            if not found:
                print 'The fragment was not found in the library'
                return None
            else:
                return p1group

def BuildByAtom(flag,p1,p2=None,p3=None,p4=None,d=None,thetaA=None,thetaT=None):
    """ This function is the interface to 1) get the data
        2) decide which operation should be done 3) prepare 
        the required inputs and call the basic functions for 
        calculations 4) return the results and errors.
        inputs:     "flag"                  define operation mode
                     "p1","p2","p3","p4"    atom dictionaris
                     "d"                    scaler for bond distance in Angstrom
                     "thetaA"               scaler for angle in degree
                     "thetaT"               scaler for diheral angle in degree and in
                                            IUPAC standard -180 to 180. thetaT > 0 clockwise and thetaT < 0 anticlockwise
        output       "q"                    atom dictionary

    the bhavior of function will change depeinding on the flag and input combination
    Mode a1: (flag=a,p1=x1,d=r) :
              atom adition (q) to the p1 with bond lenght r in the Z axis direction
    Mode a2: (flag=a,p1=x1,p2=x2,d=r,thetaA=t) :
              atom adition (q) to the p1 with bond lenght r with the t for "q-p1-p2" angle. the plane of q-p1-p2 is choosen at random.
    Mode a3: (flag=a,p1=x1,p2=x2,p3=x3,d=r,thetaA=t,thetaT=phi) :
              atom adition (q) to the p1 with bond lenght r with the t for "q-p1-p2" angle and phi for q-p1-p2-p3 torsion angle.
    Mode b: (flag=b,p1=x1,p2=x2,d=r) :
              adjust bond lenght p1-p2 to r by moving p1
    Mode c: (flag=c,p1=x1,p2=x2,p3=x3,thetaA=t) :
              adjust angle p1-p2-p3 to t by moving p1
    Mode d: (flag=d,p1=x1,p2=x2,p3=x3,p4=x4,thetaT=phi) :
              adjust torsion p1-p2-p3-p4 to phi by moving p1
                                                                    """

    if (flag == 'a') :
        q = dict()
        # Mode a1: atom addition with one atom specified (p1). 
        # The new atom is added to p1 in the defult diretion is (0,0,1).
        if (p1 is not None) :

            #check the distance
            if d == None or d <= 0:
                print 'Bad bond lenght'
                return None
            else:
                # transformation
                dv = np.array((1.,0.,0.0))
                q['xyz'] = list(translate(p1['xyz'],dv,d))

            # Mode a2: atom addition with two atoms specified (p1 ,p2). 
            # The new atom is added to p1 with the angle thetaA q-p1-p2
            if (p2 is not None) :

                #check the distance and angle
                if thetaA is None:
                    print 'Bad angle'
                    return None

                else:
                    # choose a rotation vector 
                    dq1 = np.array(q['xyz']) - np.array(p1['xyz'])
                    d21 = np.array(p2['xyz']) - np.array(p1['xyz'])
                    rv = np.cross(d21,dq1)
                    #in case atoms are linear
                    if np.linalg.norm(rv) == 0 : 
                        rv = np.cross(np.random.rand(3), d21)
                    current_thetaA = measure(q['xyz'],p1['xyz'],p2['xyz'])
                    q['xyz']  = list(rotate(q['xyz'],p1['xyz'],rv,(thetaA - current_thetaA)))

                # Mode a3: atom addition with three atoms specified (p1 ,p2, p3). 
                # The new atom is added to p1 with the angle (q, p1, p2)
                # thetaA and torsion (q, p1, p2, p3) thetaT
                if (p3 is not None) :

                    #check the distance and angle
                    if thetaT is None:
                        print 'Bad torsion'
                        return None

                    else:
                        #rotate around q the trosron bond (q-p1-p2-p3).
                        rv = np.array(p1['xyz']) - np.array(p2['xyz'])
                        current_thetaT = measure(q['xyz'],p1['xyz'],p2['xyz'],p3['xyz'])
                        q['xyz'] = list(rotate(q['xyz'],p1['xyz'],rv,(thetaT - current_thetaT)))

            #End (if p1 not None)
            print 'New atom coordinates returned'
            return q

    # Mode b: adjust a bond the moving atom is p1
    elif (flag == 'b') :
        if (p2 is not None) and (p3 is None) and (p4 is None):

            #check the distance
            if d == None or d <= 0:
                print 'Bad bond lenght'
                return None

            else:
                dv = np.array(p1['xyz']) - np.array(p2['xyz'])
                current_d = measure(v1=p1['xyz'],v2=p2['xyz'])
                p1['xyz']  = list(translate(p1['xyz'],dv,d-current_d))
                print 'bond lenght adjusted'
                return p1
        else :
            print 'Not proper parameters for bond adjustement'
            return None

    # Mode c: adjust an angle (p1,p2,p3) the center is p2, the moving atom is p1
    elif (flag == 'c') :
        if (p2 is not None) and (p3 is not None) and (p4 is None):

            #check the angle
            if thetaA is None:
                print 'Bad angle'
                return None

            #input OK then add an atom in the direction of p1 p2 bond
            else:
                # the plane is defined by p1,p2,p3 right handed
                d21 = np.array(p1['xyz']) - np.array(p2['xyz'])
                d23 = np.array(p3['xyz']) - np.array(p2['xyz'])
                rv = np.cross(d23,d21)
                current_thetaA = measure(p1['xyz'],p2['xyz'],p3['xyz'])
                p1['xyz']  = list(rotate(p1['xyz'],p2['xyz'],rv,(thetaA - current_thetaA)))
                print 'Angle adjusted'
                return p1
        else :
            print 'Not proper parameters for Angle adjustement'
            return None 

    # Mode d: adjustin torsion (p1,p2,p3,p4) move p1
    if (flag == 'd'):
        if (p2 is not None) and (p3 is not None) and (p4 is not None):

            #check the angle
            if thetaT is None:
                print 'Torsion angle'
                return None 

            else:
                #rotate around the trosron bond (p2-p3). the p1 should move
                rv = np.array(p2['xyz']) - np.array(p3['xyz'])
                current_thetaT = measure(v1=p1['xyz'],v2=p2['xyz'],v3=p3['xyz'],v4=p4['xyz'])
                deg2rotate = thetaT - current_thetaT
                p1['xyz'] = list(rotate(p1['xyz'],p2['xyz'],rv,(deg2rotate)))
                print 'Torsion adjusted'
                return p1 
        else :
            print 'Not proper parameters for Torsion adjustement'
            return None

    else: #over flags
        print 'Operation is not defined'
        return None


def BuildByGroup(flag,p1,p2=None,p3=None,p4=None,d=None,thetaA=None,thetaT=None,fragname=None,fraglib=None,p1group=None):
    """ This function is the interface to 1) get the data
        2) decide which operation should be done 3) prepare 
        the required inputs and call the basic functions for 
        calculations 4) return the results and errors.
        inputs:     "flag"                  define operation mode
                     "p1","p2","p3","p4"    atom dictionaries
                     "d"                    scaler for bond distance in Angstrom
                     "thetaA"               scaler for angle in degree
                     "thetaT"               scaler for diheral angle in degree and in
                                            IUPAC standard -180 to 180. thetaT > 0 clockwise and thetaT < 0 anticlockwise
        output       "p1group"              dictionary {1:{name:'H1',xyz:[], symbol:'H'}}

    the bhavior of function will change depeinding on the flag and input combination
    in the addition mode the  atom 1 in the library is the avtive atom (q)
    Mode a1: (flag=a,p1=x1,d=r,fragname,fraglib) :
              group adition (p1group) to the p1 with bond lenght r releive to q in the library defined direction
    Mode a2: (flag=a,p1=x1,p2=x2,d=r,thetaA=t,fragname,fraglib) :
              group adition (p1group) to the p1 with bond lenght r releive to q in the library defined direction and t angle got q-p1-p2
    Mode a3: (flag=a,p1=x1,p2=x2,p3=x3,d=r,thetaA=t,thetaT=phi,fragname,fraglib) :
              group adition (p1group) to the p1 with bond lenght r releive to q in the library defined direction and t angle got q-p1-p2 and
              phi torsion for q-p1-p2-p3
    NOTE:
    in the b,c,d mode the active atom is the p1, but the operation is done on the p1group. if p1 should change it must be included in the 
    p1group. this gives the freedom to operate on a group of atoms by defining arbitaray bond, angle, and torsions. 
    
    Mode b: (flag=b,p1=x1,p2=x2,d=r,p1group) :
              adjust bond lenght along p1-p2 to r by moving p1group
    Mode c: (flag=c,p1=x1,p2=x2,p3=x3,thetaA=t,p1group) :
              adjust angle p1-p2-p3 to t by moving p1group
    Mode d: (flag=d,p1=x1,p2=x2,p3=x3,p4=x4,thetaT=phi,p1group) :
              adjust torsion p1-p2-p3-p4 to phi by moving p1group
                                                                    """

    #Addition mode 
    if (flag == 'a'):

        if fraglib  == None or fragname == None:
            print 'Need fragment-library path and fragment name'
            return None

        #load fragment library and fetch the coordinates
        p1group = getfragment(fraglib,fragname)
        if p1group == None: return None

        #Mode a1. add a fragmet to p1 atom from library by distance d > 0. 
        if (p1 is not None):
            if d < 0 :
                print 'Bad distance'
                return None

            #tanslation vector for bringing the group to place by the dummy atom(D0)
            dv = np.array(p1['xyz']) - np.array(p1group[0]['xyz'])
            for atom in p1group:
                p1group[atom]['xyz'] = dv + np.array(p1group[atom]['xyz'])

            #tranlation to adjust d p1-p1group[1]
            current_d = measure( np.array(p1['xyz']),np.array(p1group[1]['xyz']) )
            dv = np.array(p1group[1]['xyz']) - np.array(p1['xyz'])
            for atom in p1group:
                p1group[atom]['xyz'] = list(translate( np.array(p1group[atom]['xyz']), dv, d-current_d ))

            #Mode a2. add a fragmet to p1 atom from library by distance d > 0 and angle thetaA (q1-p1-p2).
            #q1 is the first atom in the p1group
            if (p2 is not None) :
                if thetaA == None:
                    print 'Bad Angle'
                    return None

                #rotate the p1group by angle thetaA p1group[1]-p1-p2
                dq1 = np.array(p1group[1]['xyz']) - np.array(p1['xyz'])
                d21 = np.array(p2['xyz']) - np.array(p1['xyz'])
                rv = np.cross(d21,dq1)

                #in case atoms are linear
                if np.linalg.norm(rv) == 0 : 
                    rv = np.cross(d21,np.random.rand(3))
                current_thetaA = measure(np.array(p1group[1]['xyz']),np.array(p1['xyz']),np.array(p2['xyz']))
                for atom in p1group:
                    p1group[atom]['xyz'] = list(rotate(np.array(p1group[atom]['xyz']),np.array(p1['xyz']),rv,(thetaA - current_thetaA)))

                #Mode a3. rotate by torsion thetaT(p1group[1]-p1-p2-p3)
                if (p3 is not None) :
                    if thetaT == None:
                        print 'Bad Torsion angle'
                        return None

                    # The q1 first atom of fragment is the hot atom
                    rv = np.array(p1['xyz']) - np.array(p2['xyz'])
                    current_thetaT = measure(np.array(p1group[1]['xyz']),np.array(p1['xyz']),np.array(p2['xyz']),np.array(p3['xyz']))
                    for atom in p1group:
                        p1group[atom]['xyz'] = list(rotate(np.array(p1group[atom]['xyz']),np.array(p1['xyz']),rv,(thetaT - current_thetaT)))

            #remove the dummy atom before return
            del p1group[0]
            print 'Group added'
            return p1group

    if (flag == 'b'):
        if (p2 is not None) and (p1group is not None) and (p3 is None) and (p4 is None):
            if d < 0 :
                print 'Bad distance'
                return None

            #tranlation to adjust d
            current_d = measure(np.array(p1['xyz']),np.array(p2['xyz']))
            dv = np.array(p1['xyz']) - np.array(p2['xyz'])
            for atom in p1group:
                p1group[atom]['xyz'] = list(translate(p1group[atom]['xyz'], dv, d-current_d))
            print 'Bond fixed'
            return p1group
        else:
            print 'Not proper parameters for bond adjustement'
            return None

    if (flag == 'c'):
        if (p2 is not None) and (p1group is not None) and (p3 is not None) and (p4 is None):
            if thetaA == None:
                print 'Bad Angle'
                return None

            #rotate the p1group by p1-p2-p3 atom
            d21 = np.array(p1['xyz']) - np.array(p2['xyz'])
            d23 = np.array(p3['xyz']) - np.array(p2['xyz'])
            rv = np.cross(d23, d21)

            #in case atoms are linear
            if np.linalg.norm(rv) == 0 : 
                rv = np.cross(np.random.rand(3),d23)

            current_thetaA = measure(np.array(p1['xyz']),np.array(p2['xyz']),np.array(p3['xyz']))
            for atom in p1group:
                p1group[atom]['xyz'] = list(rotate(p1group[atom]['xyz'],np.array(p2['xyz']),rv,(thetaA - current_thetaA)))
            print 'Angle fixed'
            return p1group
        else:
            print 'Not proper parameters for angle adjustement'
            return None

    if (flag == 'd'):
        if (p2 is not None) and (p1group is not None) and (p3 is not None) and (p4 is not None):
            if thetaT == None:
                return (None,'Bad torsion')

            #rotate around the trosron bond (p2-p3). the p1 should move
            rv = np.array(p2['xyz']) - np.array(p3['xyz'])
            current_thetaT = measure(np.array(p1['xyz']),np.array(p2['xyz']),np.array(p3['xyz']),np.array(p4['xyz']))
            for atom in p1group:
                p1group[atom]['xyz'] = list(rotate(p1group[atom]['xyz'],np.array(p2['xyz']),rv,(thetaT - current_thetaT)))
            print 'Torsion fixed'
            return p1group
        else :
            print 'Not proper parameters for Torsion adjustement'
            return None


if __name__ == '__main__':

    #TEST building tyrosin side chain
    #        h3        h2
    #          \      /
    #           c3--c2       hb1    ha1
    #          /      \     /      /
    #ho--o4--c4        c1--cb-- + ca--ha3
    #          \      /     \      \
    #           c5--c6       hb2    ha2
    #          /      \
    #         h5       h6             

    #define a random point
    mol = dict()
    c1 = np.random.rand(3)
    atom = {'atomnr':15, 'xyz':list(c1),'name':'C1','symbol':'C'}
    mol['1'] = atom

    # add c2 in d = 1.5 A
    mol['2'] = BuildByAtom('a',p1=mol['1'],d=1.5)
    mol['2']['symbol'] , mol['2']['name'] , mol['2']['atomnr']= 'c', 'C2', 16

    # add c3 to c2 with 120 angle fro c3-c2-c1
    mol['3'] = BuildByAtom('a',p1=mol['2'],p2=mol['1'],d=1.5,thetaA=120.0)
    mol['3']['symbol'] , mol['3']['name'] , mol['3']['atomnr']= 'c', 'C3', 17

    # add c4 to c3 with 120 angle for (c4-c3-c2) and 0 degre torsion for (c4-c3-c2-c1)
    mol['4'] = BuildByAtom("a",p1=mol['3'],p2=mol['2'],p3=mol['1'],d=1.5,thetaA=120.0,thetaT=0.0)
    mol['4']['symbol'] , mol['4']['name'] , mol['4']['atomnr']= 'c', 'C4', 18

    mol['5'] = BuildByAtom("a",p1=mol['4'],p2=mol['3'],p3=mol['2'],d=1.5,thetaA=120.0,thetaT=0.0)
    mol['5']['symbol'] , mol['5']['name'] , mol['5']['atomnr']= 'c', 'C5', 19

    mol['6'] = BuildByAtom("a",p1=mol['5'],p2=mol['4'],p3=mol['3'],d=1.5,thetaA=120.0,thetaT=0.0)
    mol['6']['symbol'] , mol['6']['name'] , mol['6']['atomnr']= 'c', 'C6', 20

    mol['7'] = BuildByAtom("a",p1=mol['4'],p2=mol['3'],p3=mol['2'],d=1.4,thetaA=120.0,thetaT=180.0)
    mol['7']['symbol'] , mol['7']['name'] , mol['7']['atomnr']= 'o', 'O4', 21

    mol['8'] = BuildByAtom("a",p1=mol['1'],p2=mol['2'],p3=mol['3'],d=1.5,thetaA=120.0,thetaT=180.0)
    mol['8']['symbol'] , mol['8']['name'] , mol['8']['atomnr']= 'c', 'Cb', 22

    mol['9'] = BuildByAtom("a",p1=mol['8'],p2=mol['1'],p3=mol['2'],d=1.0,thetaA=120.0,thetaT=-30.0)
    mol['9']['symbol'] , mol['9']['name'] , mol['9']['atomnr']= 'h', 'Hb1', 23

    mol['10'] = BuildByAtom("a",p1=mol['8'],p2=mol['1'],p3=mol['2'],d=1.0,thetaA=109.0,thetaT=-140.0)
    mol['10']['symbol'] , mol['10']['name'] , mol['10']['atomnr']= 'h', 'Hb1', 24

    mol['11'] = BuildByAtom("a",p1=mol['2'],p2=mol['3'],p3=mol['4'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['11']['symbol'] , mol['11']['name'] , mol['11']['atomnr']= 'h', 'H2', 25

    mol['12'] = BuildByAtom("a",p1=mol['3'],p2=mol['4'],p3=mol['5'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['12']['symbol'] , mol['12']['name'] , mol['12']['atomnr']= 'h', 'H3', 26

    mol['13'] = BuildByAtom("a",p1=mol['5'],p2=mol['6'],p3=mol['1'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['13']['symbol'] , mol['13']['name'] , mol['13']['atomnr']= 'h', 'H5', 27

    mol['14'] = BuildByAtom("a",p1=mol['6'],p2=mol['5'],p3=mol['4'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['14']['symbol'] , mol['14']['name'] , mol['14']['atomnr']= 'h', 'H6', 28

    mol['15'] = BuildByAtom("a",p1=mol['7'],p2=mol['4'],p3=mol['3'],d=1.0,thetaA=109.0,thetaT=90.0)
    mol['15']['symbol'] , mol['15']['name'] , mol['15']['atomnr']= 'h', 'HO', 29





    with open('01.xyz','w') as f:
        f.write('%d \n' %(len(mol)))
        for key in mol: 
            print mol[key]
            f.write('%s %.4f %.4f %.4f \n' %(mol[key]['name'],mol[key]['xyz'][0],mol[key]['xyz'][1],mol[key]['xyz'][2]))

    #TEST adjust bond,angle,torsion by atom
    mol['15'] = BuildByAtom("b",p1=mol['15'],p2=mol['7'],d=1.5)
    mol['15'] = BuildByAtom("c",p1=mol['15'],p2=mol['7'],p3=mol['4'],d=1.5,thetaA=120.)
    mol['15'] = BuildByAtom("d",p1=mol['15'],p2=mol['7'],p3=mol['4'],p4=mol['3'],d=1.5,thetaA=120.,thetaT=-90.)

    with open('02.xyz','w') as f:
        f.write('%d \n' %(len(mol)))
        for key in mol: 
            print mol[key]
            f.write('%s %.4f %.4f %.4f \n' %(mol[key]['name'],mol[key]['xyz'][0],mol[key]['xyz'][1],mol[key]['xyz'][2]))


    #TEST add by group
    fraglib = os.path.join(os.getcwd(),'fragments.dat')

    #adding methyl group to the cb
    mtl = BuildByGroup("a",p1=mol['8'],p2=mol['1'],p3=mol['2'],d=1.5,thetaA=109.0,thetaT=90.0,fraglib=fraglib,fragname='MTL')
    #adding the rgening group to the O4
    arg = BuildByGroup("a",p1=mol['7'],p2=mol['4'],p3=mol['3'],d=1.5,thetaA=109.0,thetaT=90.0,fraglib=fraglib,fragname='ARG')
    
    with open('03.xyz','w') as f:
        f.write('%d \n' %(len(mol)+len(mtl)+len(arg)))
        for key in sorted(mol): 
            print key, mol[key]
            f.write('%s %.4f %.4f %.4f \n' %(mol[key]['name'],mol[key]['xyz'][0],mol[key]['xyz'][1],mol[key]['xyz'][2]))
        for key in sorted(mtl): 
            print key, mtl[key]
            f.write('%s %.4f %.4f %.4f \n' %(mtl[key]['name'],mtl[key]['xyz'][0],mtl[key]['xyz'][1],mtl[key]['xyz'][2]))
        for key in sorted(arg): 
            print key, arg[key]
            f.write('%s %.4f %.4f %.4f \n' %(arg[key]['name'],arg[key]['xyz'][0],arg[key]['xyz'][1],arg[key]['xyz'][2]))


    #TEST adjustmet by group. p1 must be included in p1group to be operated on
    #Rotating the ARG side chain around C2-C1--O4-HO
    arg = BuildByGroup("d",p1=arg['4'],p2=arg['1'],p3=mol['7'],p4=mol['15'],thetaT=0.,p1group=arg)


    with open('04.xyz','w') as f:
        f.write('%d \n' %(len(mol)+len(mtl)+len(arg)))
        for key in sorted(mol): 
            print key, mol[key]
            f.write('%s %.4f %.4f %.4f \n' %(mol[key]['name'],mol[key]['xyz'][0],mol[key]['xyz'][1],mol[key]['xyz'][2]))
        for key in sorted(mtl): 
            print key, mtl[key]
            f.write('%s %.4f %.4f %.4f \n' %(mtl[key]['name'],mtl[key]['xyz'][0],mtl[key]['xyz'][1],mtl[key]['xyz'][2]))
        for key in sorted(arg): 
            print key, arg[key]
            f.write('%s %.4f %.4f %.4f \n' %(arg[key]['name'],arg[key]['xyz'][0],arg[key]['xyz'][1],arg[key]['xyz'][2]))












