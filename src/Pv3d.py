import numpy as np

from ase import Atoms

def inBox(atom, boxDims):

    ''' Returns 'true' if the point falls within the rectangular box

    Args:
        atom    (Atom)      -   atom to be checked
        boxDims (np.arr)    -   3x2 array of box dimensions

    Returns:
        true if point is in box '''

    x,y,z = atom.position

    if (x >= boxDims[0][0]) and (x <= boxDims[0][1]):
        if (y >= boxDims[1][0]) and (y <= boxDims[1][1]):
            if (z >= boxDims[2][0]) and (z <= boxDims[2][1]):
                return True
    return False
