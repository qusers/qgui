############################################################################################################################
# Build.py
# Masoud Kazemi Sep 2016 
# kazemimsod@gmail.com
# This module provides the functions for 1) building by atoms 
#                                        2) adjusting the bond, angle, and torsions by atoms
# core functions: 
#               rotate(p,rp,rv,theta),translate(p,dv,d)
# auxillary functions:
#               properpoint(p)
# user interface functions:
#               BuildByAtom(flag,p1,p2=None,p3=None,p4=None,d=None,thetaA=None,thetaT=None), measure(v1,v2,v3=None,v4=None)
# BuildByAtom() BuildByGroup() and measure() interact with external user since they accept only points, distance and angle as arguments
############################################################################################################################
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

def properpoint(p):
    """ This function get a point (p) and check the elements and lenght"""
    try:
        for i in range(3):
            float(p[i])
    except (ValueError, IndexError,TypeError): 
        return False
    else:
        return True


def BuildByAtom(flag,p1,p2=None,p3=None,p4=None,d=None,thetaA=None,thetaT=None):
    """ This function is the interface to 1) get the data
        2) decide which operation should be done 3) prepare 
        the required inputs and call the basic functions for 
        calculations 4) return the results and errors.
        inputs:     "flag"                  define operation mode
                     "p1","p2","p3","p4"    cartesian points
                     "d"                    scaler for bond distance in Angstrom
                     "thetaA"               scaler for angle in degree
                     "thetaT"               scaler for diheral angle in degree and in
                                            IUPAC standard -180 to 180. thetaT > 0 clockwise and thetaT < 0 anticlockwise
        output       "q"                    cartesian points

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
        # Mode a1: atom addition with one atom specified (p1). The new atom is added to p1 in the defult diretion is (0,0,1).
        if (p2 is None) and (p3 is None) and (p4 is None):

            #check the distance
            if d == None or d <= 0:
                print 'Bad bond lenght'
                return None
            else:
                # transformation
                dv = np.array((0.,0.,1.0))
                q['xyz'] = list(translate(p1['xyz'],dv,d))
                print 'One atom added in the direction of Z axis'
                return q 

        # Mode a2: atom addition with two atoms specified (p1 ,p2). The new atom is added to p1 with the angle thetaA reletive p2
        # The algle plane is choosen random to avoid any possible runtime error.
        if (p2 is not None) and (p3 is None) and (p4 is None):

            #check the distance and angle
            if d == None or d <= 0 or thetaA is None:
                print 'Bad bond lenght or angle'
                return None

            else:
                # step 1 add the new atom to p1 in line with p1 and p2
                dv = np.array(p1['xyz']) - np.array(p2['xyz'])
                q['xyz'] = np.array(translate(p1['xyz'],dv,d))

                #step 2 rotate the q reletive to p1 with angle thetaA reletive to p2
                # choose a rotation vector at random no plane is defined
                rv = np.random.rand(3)
                rv = np.cross(dv,rv)
                q['xyz'] = list(rotate(q['xyz'],p1['xyz'],rv,(180-thetaA)))
                print 'One atom added with the set angle'
                return q

        # Mode a3: atom addition with three atoms specified (p1 ,p2, p3). The new atom is added to p1 with the angle (q, p1, p2)
        # thetaA and torsion (q, p1, p2, p3) thetaT
        if (p2 is not None) and (p3 is not None) and (p4 is None):

            #check the distance and angle
            if d == None or d <= 0 or thetaA is None or thetaT is None:
                print 'Bad bond lenght or angle'
                return None

            else:
                # step 1 add the new atom to p1 in line with p1 and p2
                dv = np.array(p1['xyz']) - np.array(p2['xyz'])
                q['xyz']  = np.array(translate(p1['xyz'],dv,d))

                #step 2 rotate the q reletive to p1 with angle thetaA reletive to p2
                # choose a rotation vector at random no plane is defined
                rv = np.random.rand(3)
                rv = np.cross(dv,rv)
                q['xyz']  = rotate(q['xyz'],p1['xyz'],rv,(180-thetaA))

                #step 3 rotate the around the trosron bond. the q should move
                current_thetaT = measure(v1=q['xyz'],v2=p1['xyz'],v3=p2['xyz'],v4=p3['xyz'])
                deg2rotate = thetaT - current_thetaT
                q['xyz'] = list(rotate(q['xyz'],p1['xyz'],dv,(deg2rotate)))
                print 'One atom added with the set angle'
                return q
        else:
            print 'Not proper parameters for atom addition'
            return None 

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
                p1['xyz'] = rotate(p1['xyz'],p2['xyz'],rv,(deg2rotate))
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
                     "p1","p2","p3","p4"    cartesian points
                     "d"                    scaler for bond distance in Angstrom
                     "thetaA"               scaler for angle in degree
                     "thetaT"               scaler for diheral angle in degree and in
                                            IUPAC standard -180 to 180. thetaT > 0 clockwise and thetaT < 0 anticlockwise
        output       "p1group"              dictionary {1 ['H1', x, y, z, 'H']}

    the bhavior of function will change depeinding on the flag and input combination
    in the addition mode the  atom 1 in the library is the avtive atom (q)
    Mode a1: (flag=a,p1=x1,d=r,fragname,fraglib) :
              group adition (p1group) to the p1 with bond lenght r releive to q in the library defined direction
    Mode a2: (flag=a,p1=x1,p2=x2,d=r,thetaA=t,fragname,fraglib) :
              group adition (p1group) to the p1 with bond lenght r releive to q in the library defined direction and t angle got q-p1-p2
    Mode a3: (flag=a,p1=x1,p2=x2,p3=x3,d=r,thetaA=t,thetaT=phi,fragname,fraglib) :
              group adition (p1group) to the p1 with bond lenght r releive to q in the library defined direction and t angle got q-p1-p2 and
              phi torsion for q-p1-p2-p3
    in the b,c,d mode the active atom is the p1, but the operation is done on the p1group. if p1 should change it must also be included in the 
    p1 group. this gives the freedom to operate on a group of atoms by defining arbitaray bond, angle, and torsions. 
    TODO this can be chaged to the case where the first element of p1group is the active atom???
    Mode b: (flag=b,p1=x1,p2=x2,d=r,p1group) :
              adjust bond lenght along p1-p2 to r by moving p1group
    Mode c: (flag=c,p1=x1,p2=x2,p3=x3,thetaA=t,p1group) :
              adjust angle p1-p2-p3 to t by moving p1group
    Mode d: (flag=d,p1=x1,p2=x2,p3=x3,p4=x4,thetaT=phi,p1group) :
              adjust torsion p1-p2-p3-p4 to phi by moving p1group
                                                                    """

    #make sure inputs are float
    for parm in [p1,p2,p3,p4]:
        if parm is not None:
            if not properpoint(parm): return (None,'Bad cartesian coordinate')
            parm = np.array(parm) + 0.
    for parm in [d,thetaA,thetaT]:
            if parm  is not None :
                    try:
                        parm = parm + 0.
                    except (ValueError,TypeError): 
                        return (None,'Bad scaler parameters')

    #Addition mode 
    if (flag == 'a'):
        #load fragment library and fetch the coordinates
        p1group = dict()
        found = False
        if fraglib  == None: return (None, 'Need fragment library')
        if fragname == None: return (None, 'Need fragment Nmae')
        if not os.path.isfile(fraglib):
            return (None, 'No fragment library file was found')
        else:
            with open(fraglib, 'r') as f:
                for line in f:
                    field = line.split()
                    if len(field) : #ignor empty lines
                        if field[0] == fragname: 
                            found = True
                            continue
                        elif found and field[0] == "END":
                            break
                        elif found:
                            p1group.update({field[0]:field[1:]})
                if not found:
                    return (None, 'The fragment was not found in the library')

        #Mode a1. add a fragmet to p1 atom from library by distance d > 0. 
        if (p1 is not None):
            if d < 0 :
                return (None, 'Bad distance')

            #tanslation vector for bringing the group to place by the dummy atom(q0)
            q0 = np.array(map(float,p1group['0'][1:4])) 
            dv = p1 - q0

            for atom in p1group:
                p1group[atom][1:4] = dv + np.array(map(float,p1group[atom][1:4]))

            #tranlation to adjust d
            q1 = np.array(map(float,p1group['1'][1:4])) 
            current_d = measure(p1,q1)
            dv = np.array(map(float,p1group['1'][1:4])) - p1
            for atom in p1group:
                p1group[atom][1:4] = translate( np.array(map(float,p1group[atom][1:4])), dv, d-current_d )

            #Mode a2. add a fragmet to p1 atom from library by distance d > 0 and angle thetaA (q1-p1-p2).
            #q1 is the first atom in the p1group
            if (p2 is not None) :
                if thetaA == None:
                    return (None,'Bad Angle')
                #rotate the p1group by q1 atom
                # the plane is defined at random to avoid the case q1-p1-p2 are linear
                q1 = np.array(map(float,p1group['1'][1:4])) 
                dq1 = q1 - p1
                d21 = p2 - p1
                rv = np.cross(d21,dq1)
                #in case atoms are linear
                if np.linalg.norm(rv) == 0 : rv = np.random.rand(3)
                current_thetaA = measure(q1,p1,p2)
                for atom in p1group:
                    p1group[atom][1:4] = rotate( np.array(map(float,p1group[atom][1:4])) ,p1,rv,(thetaA - current_thetaA))

                #Mode a3. add a fragmet to p1 atom from library by distance d > 0 and angle thetaA (q1-p1-p2) and a torsion thetaT(q1-p1-p2-p3)
                if (p3 is not None) :
                    if thetaT == None:
                        return (None,'Bad Torsion angle')
                    # The q1 first atom of fragment is the hot atom
                    q1 = np.array(map(float,p1group['1'][1:4])) 
                    rv = p1 - p2
                    current_thetaT = measure(v1=q1,v2=p1,v3=p2,v4=p3)
                    for atom in p1group:
                        p1group[atom][1:4] = rotate( np.array(map(float,p1group[atom][1:4])) ,p1,rv,(thetaT - current_thetaT))

            #remove the dummy atom before return
            del p1group['0']
            return (p1group, 'Group added')

    if (flag == 'b'):
        if (p2 is not None) and (p1group is not None) and (p3 is None) and (p4 is None):
            if d < 0 :
                return (None, 'Bad distance')
            #tranlation to adjust d
            current_d = measure(p1,p2)
            dv = p1 - p2
            for atom in p1group:
                p1group[atom][1:4] = translate( np.array(map(float,p1group[atom][1:4])), dv, d-current_d )
            return (p1group, 'Bond fixed')
        else:
            return (None, 'Not proper parameters for bond adjustement')

    if (flag == 'c'):
        if (p2 is not None) and (p1group is not None) and (p3 is not None) and (p4 is None):
            if thetaA == None:
                return (None,'Bad Angle')
            #rotate the p1group by q2 atom
            #the plane is defined at random to avoid the case p1-p1-p2 are linear
            d21 = p1 - p2
            d23 = p3 - p2
            rv = np.cross(d23,d21)
            current_thetaA = measure(p1,p2,p3)
            #in case atoms are linear
            if np.linalg.norm(rv) == 0 : rv = np.random.rand(3)
            current_thetaA = measure(p1,p2,p3)
            for atom in p1group:
                p1group[atom][1:4] = rotate( np.array(map(float,p1group[atom][1:4])) ,p2,rv,(thetaA - current_thetaA))
            return (p1group, 'Angle fixed')
        else:
            return (None, 'Not proper parameters for angle adjustement')

    if (flag == 'd'):
        if (p2 is not None) and (p1group is not None) and (p3 is not None) and (p4 is not None):
            if thetaT == None:
                return (None,'Bad torsion')
            #rotate around the trosron bond (p2-p3). the p1 should move
            rv = p2 - p3
            current_thetaT = measure(v1=p1,v2=p2,v3=p3,v4=p4)
            print current_thetaT
            for atom in p1group:
                p1group[atom][1:4] = rotate( np.array(map(float,p1group[atom][1:4])) ,p2,rv,(thetaT - current_thetaT))
            return (p1group, 'Torsion fixed')
        else :
            return (None, 'Not proper parameters for Torsion adjustement')


if __name__ == '__main__':

    #TEST building tyrosin side chain
    #        h3        h2
    #          \      /
    #           c3--c2       hb1  ha1
    #          /      \     /    /
    #ho--o4--c4        c1--cb--ca--ha3
    #          \      /     \    \
    #           c5--c6       hb2  ha2
    #          /      \
    #         h5       h6             

    #define a random point
    mol = dict()
    c1 = np.random.rand(3)
    atom = {'atomnr':15, 'xyz':list(c1),'name':'c1','symbol':'c'}
    mol['1'] = atom

    # add c2 in d = 1.5 A
    mol['2'] = BuildByAtom('a',p1=mol['1'],d=1.5)
    mol['2']['symbol'] , mol['2']['name'] , mol['2']['atomnr']= 'c', 'c2', 16

    # add c3 to c2 with 120 angle fro c3-c2-c1
    mol['3'] = BuildByAtom('a',p1=mol['2'],p2=mol['1'],d=1.5,thetaA=120.0)
    mol['3']['symbol'] , mol['3']['name'] , mol['3']['atomnr']= 'c', 'c3', 17

    # add c4 to c3 with 120 angle for (c4-c3-c2) and 0 degre torsion for (c4-c3-c2-c1)
    mol['4'] = BuildByAtom("a",p1=mol['3'],p2=mol['2'],p3=mol['1'],d=1.5,thetaA=120.0,thetaT=0.0)
    mol['4']['symbol'] , mol['4']['name'] , mol['4']['atomnr']= 'c', 'c4', 18

    mol['5'] = BuildByAtom("a",p1=mol['4'],p2=mol['3'],p3=mol['2'],d=1.5,thetaA=120.0,thetaT=0.0)
    mol['5']['symbol'] , mol['5']['name'] , mol['5']['atomnr']= 'c', 'c5', 19

    mol['6'] = BuildByAtom("a",p1=mol['5'],p2=mol['4'],p3=mol['3'],d=1.5,thetaA=120.0,thetaT=0.0)
    mol['6']['symbol'] , mol['6']['name'] , mol['6']['atomnr']= 'c', 'c6', 20

    mol['7'] = BuildByAtom("a",p1=mol['4'],p2=mol['3'],p3=mol['2'],d=1.4,thetaA=120.0,thetaT=180.0)
    mol['7']['symbol'] , mol['7']['name'] , mol['7']['atomnr']= 'o', 'o4', 21

    mol['8'] = BuildByAtom("a",p1=mol['1'],p2=mol['2'],p3=mol['3'],d=1.5,thetaA=120.0,thetaT=180.0)
    mol['8']['symbol'] , mol['8']['name'] , mol['8']['atomnr']= 'c', 'cb', 22

    mol['9'] = BuildByAtom("a",p1=mol['8'],p2=mol['1'],p3=mol['2'],d=1.0,thetaA=120.0,thetaT=-30.0)
    mol['9']['symbol'] , mol['9']['name'] , mol['9']['atomnr']= 'h', 'hb1', 23

    mol['10'] = BuildByAtom("a",p1=mol['8'],p2=mol['1'],p3=mol['2'],d=1.0,thetaA=109.0,thetaT=-140.0)
    mol['10']['symbol'] , mol['10']['name'] , mol['10']['atomnr']= 'h', 'hb1', 24

    mol['11'] = BuildByAtom("a",p1=mol['2'],p2=mol['3'],p3=mol['4'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['11']['symbol'] , mol['11']['name'] , mol['11']['atomnr']= 'h', 'h2', 25

    mol['12'] = BuildByAtom("a",p1=mol['3'],p2=mol['4'],p3=mol['5'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['12']['symbol'] , mol['12']['name'] , mol['12']['atomnr']= 'h', 'h3', 26

    mol['13'] = BuildByAtom("a",p1=mol['5'],p2=mol['6'],p3=mol['1'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['13']['symbol'] , mol['13']['name'] , mol['13']['atomnr']= 'h', 'h5', 27

    mol['14'] = BuildByAtom("a",p1=mol['6'],p2=mol['5'],p3=mol['4'],d=1.0,thetaA=120.0,thetaT=180.0)
    mol['14']['symbol'] , mol['14']['name'] , mol['14']['atomnr']= 'h', 'h6', 28

    mol['15'] = BuildByAtom("a",p1=mol['7'],p2=mol['4'],p3=mol['3'],d=1.0,thetaA=109.0,thetaT=90.0)
    mol['15']['symbol'] , mol['15']['name'] , mol['15']['atomnr']= 'h', 'ho', 29





    with open('01.xyz','w') as f:
        f.write('%d \n' %(len(mol)))
        for key in mol: 
            print mol[key]
            f.write('%s %.4f %.4f %.4f \n' %(mol[key]['name'],mol[key]['xyz'][0],mol[key]['xyz'][1],mol[key]['xyz'][2]))

    #TEST adjust bond,angle,torsion by atom
    mol['15'] = BuildByAtom("b",p1=mol['15'],p2=mol['7'],d=1.5)
    mol['15'] = BuildByAtom("c",p1=mol['15'],p2=mol['7'],p3=mol['4'],d=1.5,thetaA=120.)
    mol['15'] = BuildByAtom("d",p1=mol['15'],p2=mol['7'],p3=mol['4'],p4=mol['3'],d=1.5,thetaA=120.,thetaT=0.)

    with open('02.xyz','w') as f:
        f.write('%d \n' %(len(mol)))
        for key in mol: 
            print mol[key]
            f.write('%s %.4f %.4f %.4f \n' %(mol[key]['name'],mol[key]['xyz'][0],mol[key]['xyz'][1],mol[key]['xyz'][2]))


