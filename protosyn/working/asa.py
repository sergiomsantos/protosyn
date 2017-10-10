import numpy as np
from collections import defaultdict


# https://github.com/boscoh/pdbremix/blob/master/pdbremix/asa.py

def make_boxes(a, d_max):
    '''
    Returns dictionary which keys are indecies of boxes (regions)
    with d_max length side and values
    are indicies of atoms belonging to these boxes
    '''
    b = defaultdict(list) # space divided into boxes
    for i in xrange(len(a)):
        atom = a[i]
        box_coor = tuple(int(np.floor(x / d_max)) for x in atom.xyz)
        b[box_coor].append(i)
    return b

def neighbor_atoms(b, box):
    '''
    Returns list of atoms from half of neighbouring boxes of the box
    another half is accounted when symmetric (opposite) boxes considered
    '''
    na = [] # list for neighboring atoms
    x, y, z = box # coordinates of the box
    # top layer consisting of 9 boxes
    if (x + 1, y + 1, z +1) in b: na.extend(b[(x + 1, y + 1, z +1)])
    if (x, y + 1, z +1) in b: na.extend(b[(x, y + 1, z +1)])
    if (x + 1, y, z +1) in b: na.extend(b[(x + 1, y, z +1)])
    if (x, y, z +1) in b: na.extend(b[(x, y, z +1)])
    if (x - 1, y + 1, z +1) in b: na.extend(b[(x - 1, y + 1, z +1)])
    if (x + 1, y - 1, z +1) in b: na.extend(b[(x + 1, y - 1, z +1)])
    if (x, y - 1, z +1) in b: na.extend(b[(x, y - 1, z +1)])
    if (x - 1, y, z +1) in b: na.extend(b[(x - 1, y, z +1)])
    if (x - 1, y - 1, z +1) in b: na.extend(b[(x - 1, y - 1, z +1)])
    # half of the middle layer excluding the box itself (4 boxes)
    if (x + 1, y + 1, z) in b: na.extend(b[(x + 1, y + 1, z)])
    if (x, y + 1, z) in b: na.extend(b[(x, y + 1, z)])
    if (x + 1, y, z) in b: na.extend(b[(x + 1, y, z)])
    if (x + 1, y - 1, z) in b: na.extend(b[(x + 1, y - 1, z)])
    return na

def adjacency_list(a, d_max):
    '''
    Returns adjacency list from coordinate file
    in O(len(a)) time
    '''
    b = make_boxes(a, d_max) # put atoms into the boxes with dmax length side
    # now go on boxes and check connections inside 3x3 superboxes
    conn = [[] for i in xrange(len(a))] # list of bond lengths each atom implicated
    for box in b:
        lb = len(b[box])
        for i in range(lb):
            a1 = b[box][i]
            # check possible connections inside the box
            for j in range(i+1, lb):
                a2 = b[box][j]
                add_bond(a, a1, a2, conn, d_max)
            # check connections with atoms from neighbouring boxes
            na = neighbor_atoms(b, box) # list of such atoms
            for a2 in na:
                add_bond(a, a1, a2, conn, d_max)
    return conn

def add_bond(a, a1, a2, conn, d_max):
    '''
    If distance between atoms a1 and a2 is less than d_max (neighboring atoms),
    add atoms a1 and a2 in adjacency list conn to each other
    '''
    atom1 = a[a1]
    atom2 = a[a2]
    if distance(atom1.xyz,atom2.xyz) <= d_max: # * d_max:  # connected
        conn[a1].append(a2)
        conn[a2].append(a1)



def generate_sphere_points(n):
    """
    Returns list of coordinates on a sphere using the Golden-
    Section Spiral algorithm.
    """
    points = []
    inc = np.pi * (3 - np.sqrt(5))
    offset = 2.0 / n
    for k in xrange(n):
        y = k * offset - 1 + (offset / 2)
        r = np.sqrt(1 - y*y)
        phi = k * inc
        points.append([np.cos(phi)*r, y, np.sin(phi)*r])
    return np.array(points)

def distance(a,b):
    ab = b-a
    return np.sqrt(np.dot(ab,ab))

def find_neighbor_indices(atoms, probe, k):
    pivot = atoms[k]
    radius = pivot.radius + 2*probe
    indices = range(k) + range(k+1, len(atoms))
    return [i for i in indices if distance(pivot.xyz,atoms[i].xyz) < (radius+atoms[i].radius)]


def find_neighbor_indices_modified(atoms, indices, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k. 
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k.radius + probe + probe
    for i in indices:
        if i == k: continue
        atom_i = atoms[i]
        # dist2 = v3.mag2(atom_k.pos - atom_i.pos) # ToAn
        # if dist2 < (radius + atom_i.radius) ** 2: # ToAn
        #     neighbor_indices.append(i)
        dist = distance(atom_k.xyz, atom_i.xyz) # ToAn
        if dist < (radius + atom_i.radius): # ToAn
            neighbor_indices.append(i)
    return neighbor_indices

def calculate_asa_optimized(atoms, probe, n_sphere_point=960):
    """
    Returns the accessible-surface areas of the atoms, by rolling a
    ball with probe radius over the atoms with their radius
    defined.
    """
    sphere_points = generate_sphere_points(n_sphere_point)
    points = []
    const = 4.0 * np.pi / len(sphere_points)
    areas = []
    # neighbor_list = adjacency_list(atoms, 2 * (probe + max(atoms, key=lambda p: p.radius).radius))
    neighbor_list = adjacency_list(atoms, 2 * (probe + 1.0))

    directions = []
    for i, atom_i in enumerate(atoms):
    
        neighbor_indices = [neig for neig in neighbor_list[i]]
        neighbor_indices = find_neighbor_indices_modified(atoms, neighbor_indices, probe, i) # even further narrow diapazon
        n_neighbor = len(neighbor_indices)
        j_closest_neighbor = 0
        radius = probe + atom_i.radius
        
        n_accessible_point = 0
        for point in sphere_points:
            is_accessible = True
            test_point = point*radius + atom_i.xyz#v3.scale(point, radius) + atom_i.pos
            cycled_indices = range(j_closest_neighbor, n_neighbor)
            cycled_indices.extend(range(j_closest_neighbor))
            
            for j in cycled_indices:
                atom_j = atoms[neighbor_indices[j]]
                r = atom_j.radius + probe
                #diff2 = v3.mag2(atom_j.pos - test_point)
                diff2 = distance(atom_j.xyz, test_point)
                
                if diff2 < r:
                    j_closest_neighbor = j
                    is_accessible = False
                    break
            if is_accessible:
                directions.append(point + test_point)
                points.append(test_point)
                n_accessible_point += 1
    
        area = const*n_accessible_point*radius*radius 
        areas.append(area)

    return points, directions #areas



                



if __name__ == '__main__':

    from molecule import Molecule
    from sys import argv
    import utils

    molecule = Molecule.LoadFromFile(argv[1])
    atoms = [at for at in molecule.iteratoms()]
    surf_pts, surf_dirs = calculate_asa_optimized(atoms, 1.4, n_sphere_point=150)
    # print surf_pts

    #print molecule
    #print utils.as_pdb(surf_pts)
    # with open('surf.xyz', 'w') as fout:
    #     print >> fout, '%5d'%len(surf_pts)
    #     print >> fout, 'surface points'
    #     print >> fout, '\n'.join(['X %12.3f %12.3f %12.3f'%(x,y,z) for x,y,z in surf_pts])
    
    # with open('surf_dir.xyz', 'w') as fout:
    #     print >> fout, '%5d'%len(surf_dirs)
    #     print >> fout, 'surface directions'
    #     print >> fout, '\n'.join(['X %12.3f %12.3f %12.3f'%(x,y,z) for x,y,z in surf_dirs])
    
    with open('surface.cgo', 'w') as fout:
        fmt = 'VERTEX, %.2f, %.2f, %.2f, VERTEX, %.2f, %.2f, %.2f,'
        print >> fout, 'from pymol.cgo import *'
        print >> fout, 'from pymol import cmd'
        print >> fout, 'obj = [BEGIN, LINES, COLOR, 1.0, 1.0, 1.0,'

        print >> fout, '\n'.join([fmt%(x1,y1,z1,x2,y2,z2) for (x1,y1,z1),(x2,y2,z2) in zip(surf_pts, surf_dirs)])
        
        print >> fout, 'END]'
        print >> fout, 'cmd.load_cgo(obj, "treta")'

