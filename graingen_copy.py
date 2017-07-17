#!/usr/bin/python

# Modified from voronoi-1.3.py Nov. 12th 2012

# Bruteforce approach to generate polycrystalline Silicon
# for a given GB distribution size. 

# TODO: catch overlapping atoms, delete
# TODO: If needed implement PBC, right now the box edges
# will result in GB atoms

from random import randint
from math import acos,pi
import numpy as np

def rotate(crys,seed,axis):
    '''
    Use Rodrigues rotation formula axis x vector = matrix
    then use Rodrigues to rotate vector about axis for theta

    Input: 
    crys - perfect crystal
    seed - seed for angle (unitless)
    axis - random axis to be rotated
    Output:
    rotcrys - rotated crystal with seed about random axis
    '''
    angle = seed*np.pi/2. #radians


    axisXvec = np.array([[0.0, -axis[2], axis[1]],
                         [axis[2], 0.0 , -axis[0]],
                         [-axis[1], axis[0], 0.0]])

    #Rodrigues Formula
    rotmat = (np.identity(3) + axisXvec * np.sin(angle) 
             + np.dot(axisXvec,axisXvec) * (1. - np.cos(angle)))


    #Rotate Crystal
    natoms = np.size(crys[:,0])
    rotcrys = np.ndarray((natoms,3),dtype=float)

    for i in xrange(natoms):
        rotcrys[i,:] = np.dot(rotmat,crys[i,:])

    return rotcrys


def lattice(a,size,basis):
    ''' 
        input:
        a - lattice parameter
        size - tuple of multiples to increase size
        basis - basis vectors of lattice
        output:
        crystal - self-explanitory
        '''

    #Diamond Lattice wycoff sites
    
    #basis = np.array([[0.0,0.0,0.0],
    #                  [0.5,0.5,0.0],
    #                  [0.5,0.0,0.5],
    #                  [0.0,0.5,0.5],
    #                  [0.25,0.25,0.25],
    #                  [0.75,0.75,0.25],
    #                  [0.25,0.75,0.75],
    #                  [0.75,0.25,0.75]],dtype=np.float64)
    
    natoms = basis.shape[0]
    ntimes = natoms*size[0]*size[1]*size[2]
    crystal = np.zeros((ntimes,3))
   
    #Just for working; make more efficient
    ind = 0 #index 
    for i in xrange(size[0]):
        for j in xrange(size[1]):
            for k in xrange(size[2]):
                for n in xrange(natoms):
                    crystal[ind,0] = a * (basis[n,0] + i) #/ size[0]
                    crystal[ind,1] = a * (basis[n,1] + j) #/ size[1]
                    crystal[ind,2] = a * (basis[n,2] + k) #/ size[2]
                    ind += 1
    return crystal

def gengrain(numberGrain,box):
    """
    Randomly generates the specified number of grain centers. To maintain
    periodicity, each point is mirrored to the neighboring 26 boxes.

    Args:
        numberGrain (int)   - the number of grains to be produced
        box         (list)  - box dimensions

    Returns:
        complete    (dict)  - a dictionary where each entry value is the set of
                              points corresponding to a single grain and its
                              images.
    """

    points = np.random.random_sample((numberGrain,3))
    #points = np.random.normal(gsize,0.1,(numberGrain,3))
    for k in range(3):
        points[:,k] *= box[k]

    complete = {}

    for i in range(numberGrain):
        complete[i] = np.empty((0,3))
    
    rng = [-1,0,1]
    for x in rng:
        for y in rng:
            for z in rng:
                toAdd = np.copy(points)
                toAdd[:,0] += x*box[0]
                toAdd[:,1] += y*box[1]
                toAdd[:,2] += z*box[2]

                for i in range(len(toAdd)):
                    complete[i] = np.concatenate((complete[i],toAdd[i].reshape(1,3)), axis=0)

    return complete

def in_grain(p,grain,pcenters):
    ''' if atom is closest to this grain center keep'''

    flag = True
    dist_square = np.sum(np.square(p-grain))
    for i,g in enumerate(pcenters):
        tmp_dist_square = np.sum(np.square(p-g))
        tmp1 = list(g)
        tmp2 = list(grain) # used for comparing list rather than np.arr
        if (tmp_dist_square < dist_square) and (tmp1 != tmp2):
             flag = False
    return flag


def fill(grainDict,numgrains,latp,sizec):

    # IMPORTANT: each basis set should be the coordinates of ONE type of atom
    basis1 = np.array([[0.0,0.0,0.0],
                      [0.5,0.5,0.0],
                      [0.5,0.0,0.5],
                      [0.0,0.5,0.5]],dtype=np.float64)

    basis2 = np.array([[0.25,0.25,0.25],
                      [0.75,0.75,0.25],
                      [0.25,0.75,0.75],
                      [0.75,0.25,0.75]],dtype=np.float64)

    allBases = [basis1, basis2]

    minbond = 1.10 #Angstroms H-H 
    box = latp*sizec
    
    numgrains = len(grainDict[0])

    newcrystal = None
    atomcrystalgids = None
    allCenters = [e1 for e2 in grainDict.values() for e1 in e2]

    atomType = 0

    # Iterate over all grains
    # Note: all images of a single grain should be rotated by the same amount
    for grainNum in range(len(grainDict)):
        atomType += 1
        points = grainDict[grainNum]
        seed = randint(1,9)*0.1
        axis = np.random.random_sample((3)) #Unit vector

        for rid,gcenter in enumerate(points):
            print "Grain id: ", grainNum
            print "Image num: ", rid

            for i in range(len(allBases)):
                atomsgrains = []
                atomgrainids = []
                basis = allBases[i]

                tmpcrystal = lattice(latp,3*sizec,basis)

                #Shift extended crystal by sizec
                tmpcrystal[:,0] -= (box[0])
                tmpcrystal[:,1] -= (box[1])
                tmpcrystal[:,2] -= (box[2])
                rotatecrys = rotate(tmpcrystal,seed,axis)
                natoms = len(rotatecrys[:,0])
                # name = 'data.tmp'+str(rid)
                # writedata(rotatecrys,box,name)

                keep = [] #Atom index 
                for at in xrange(natoms):
                    xyz = rotatecrys[at,:]
                    #check if atom is in box
                    if xyz[0] < 0.0 or xyz[0] > box[0]:
                       continue
                    elif xyz[1] < 0.0 or xyz[1] > box[1]:
                       continue
                    elif xyz[2] < 0.0 or xyz[2] > box[2]:
                       continue
                    
                    if in_grain(xyz,gcenter,allCenters) == True:
                        keep.append(at)

                graincrys = rotatecrys[keep,:]
                
                #graincrys = np.concatenate((graincrys, np.ones(len(graincrys),1)*atomType), axis=1)
                #name = 'data.tmp'+str(rid)
                #writedata(graincrys,box,name)

                atomsgrains.append(graincrys)
                tmparry = np.empty(graincrys[:,0].size,dtype=np.int32)
                #tmparry.fill(i+1) # TODO: should be atom type based on structure
                tmparry.fill(atomType)
                atomgrainids.append(tmparry)
            
                # Create newcrystal and atomcrystalgids if don't already exist
                if newcrystal is not None:
                    atomsgrains = np.concatenate(atomsgrains, axis=0)
                    atomgrainids = np.concatenate(atomgrainids, axis=0)
                    newcrystal = np.concatenate((newcrystal, atomsgrains),axis=0)
                    atomcrystalgids = np.concatenate((atomcrystalgids, atomgrainids),axis=0)
                else:
                    newcrystal = np.concatenate(atomsgrains,axis=0)
                    atomcrystalgids = np.concatenate(atomgrainids,axis=0)

            print "Center of Grain: ", points[rid]
            print "Number of Atoms in cell: ", graincrys[:,0].size
            print "Finished Grain: ",rid
            print "------------------------------------"

    #print("atomgrainids: ", atomgrainids)
    #print("atomcrystalgids: ", atomcrystalgids)
    return newcrystal,atomcrystalgids

def prune(a):
    ''' TODO: remove atoms overlapping assume PBC conditions'''
    x=a

def writedata(crys,grids,box,filname):

    f = open(filname,'w')

    natoms = len(crys[:,0])
    
    f.write('#Polycrystalline Si -Voronoi method \n') 
    f.write('\n')
    f.write('%i atoms \n' %(natoms))
    f.write('2 atom types \n\n')
        
    f.write('0.0 %f xlo xhi \n' %(box[0]))
    f.write('0.0 %f ylo yhi \n' %(box[1]))
    f.write('0.0 %f zlo zhi \n' %(box[2]))
    f.write('\n')
    f.write('\n')
    f.write('Masses\n\n')
    f.write('1 65.39\n')
    f.write('2 32.066\n\n')
    f.write('Atoms \n')
    f.write('\n')
    for i in xrange(natoms):    
        f.write('%i 0 %i 0.0 %f %f %f 0 0 0\n' %(i+1,grids[i],crys[i,0],crys[i,1],crys[i,2]))

    f.close()

if __name__ == "__main__":

    # Main Driver code
    a = 5.403

    p = 250
    x = int(np.ceil(p/a))
    y = int(np.ceil(p/a))
    #z = int(np.ceil(60/a))
    z = x
    
    # TODO: code breaks if x,y,z aren't ALL equal to each other

    print(x)
    print(y)
    print(z)

    size = np.array([x,y,z])
    ngrain = 4
    print "Box Size: ", a*size
    vbox = gengrain(ngrain,a*size)
    #Fill grains with atoms
    vcrys,vids = fill(vbox,ngrain,a,size)    
    writedata(vcrys,vids,a*size,'data.poly.big')
