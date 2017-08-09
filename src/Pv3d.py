import numpy as np
from ase import Atoms

def inBox(p, boxDims):
    '''Returns 'true' if the point falls within the rectangular box

    Args:
        p          (np.arr)    -   the point to check
        boxDims    (np.arr)    -   box dimensions; organized as {(xlo,xhi),
                                (ylo,yhi), (zlo,zhi)}

    Returns:
        True if point is in box'''

    if (p[1] >= boxDims[0][0]) and (p[1] <= boxDims[0][1]):
        if (p[2] >= boxDims[1][0]) and (p[2] <= boxDims[1][1]):
            if (p[3] >= boxDims[2][0]) and (p[3] <= boxDims[2][1])):
                return True
    return False

def inRegion(p, centers, regionId):
    '''Checks to see if point 'p' falls into the Voronoi tile specified by
    regionId.

    Args:
        p           (np.arr)    -   xyz coordinates of checked point (with atom type
                                    info)
        centers     (np.arr)    -   collection of all tile centers
        regionId    (int)       -   tile id to be checked

    Returns:
        True if point is in tile'''

    # TODO: duplicate centers in all directions
    # TODO: pull row of regionId
    # TODO: find first occurrence of minimum of row
