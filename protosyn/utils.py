import numpy as np

PDB_FMT = 'ATOM %6d %4s %-4sA %3d    %8.3f%8.3f%8.3f  1.00  0.00%12s'


def fit(target_residue, target_selection, mobile_residue, mobile_selection):
    
    k = min(len(target_selection), len(mobile_selection))
    target = target_residue.get_coordinates(*target_selection[:k])
    mobile = mobile_residue.get_coordinates(*mobile_selection[:k])

    centroid_target = np.mean(target, axis=0)
    centroid_mobile = np.mean(mobile, axis=0)

    target -= centroid_target
    mobile -= centroid_mobile
    
    H = np.dot(target.T, mobile)
    U,S,Vt = np.linalg.svd(H)
    
    R = np.dot(U, np.dot(np.diag([1, 1, np.linalg.det(np.dot(U, Vt))]), Vt))
    T = centroid_target - np.dot(R,centroid_mobile)

    xyz = mobile_residue.get_coordinates()
    xyz = np.dot(R, xyz.T).T + T
    mobile_residue.set_coordinates(xyz)

    return R,T


def calc_principal_axis(residue, atoms):
    xyz = residue.get_coordinates(*atoms)
    cov_mat = np.cov(xyz.T)
    eig_vals, eig_vecs = np.linalg.eig(cov_mat)

    eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]
    eig_pairs.sort(key=lambda x: x[0], reverse=True)
    return eig_pairs


def rotation_about_vector(matrix, angle, axis):
    theta = 0.5*angle
    q0 = np.cos(theta)
    q1,q2,q3 = np.sin(theta)*axis/np.linalg.norm(axis)
    Rmat = np.array([[1-2*q2*q2-2*q3*q3,   2*q1*q2-2*q0*q3,   2*q1*q3+2*q0*q2],
                    [   2*q2*q1+2*q0*q3, 1-2*q3*q3-2*q1*q1,   2*q2*q3-2*q0*q1],
                    [   2*q3*q1-2*q0*q2,   2*q3*q2+2*q0*q1, 1-2*q1*q1-2*q2*q2]], float)
    if len(np.shape(matrix)) == 1:
        return (np.dot(Rmat, np.array([matrix]).T).T)[0]
    return np.dot(Rmat, matrix.T).T

def closest_point_on_line(line_point, line_dir, point):
    return line_point + np.dot(line_point-point, line_dir)*line_dir


def as_pdb(residue):
    pdb = []
    for at in residue.iter_atoms():
        x,y,z = at.xyz
        name = '%-3s'%at.name
        pdb.append(PDB_FMT%(at.index+1, name, at.resname, at.resnum+1, x, y, z, at.symbol))
    
    if residue.next is None:
        pdb.append('TER')
    
    return '\n'.join(pdb)


def bonds_as_pdb(residue):
    pdb = []
    for at in residue.iter_atoms():
        others = at.get_connected_atoms()
        if others:
            atoms = [at] + others
            pdb.append('CONECT' + ''.join('%5d'%(atm.index+1) for atm in atoms))
            
    return '\n'.join(pdb)


def calculate_dihedral(x1,x2,x3,x4):
    v12 = x2 - x1
    v23 = x3 - x2
    v34 = x4 - x3
    v123 = np.cross(v12,v23)
    v234 = np.cross(v23,v34)

    angle = np.arctan2(np.dot(np.cross(v123,v234),v23)/np.sqrt(np.dot(v23,v23)),
                       np.dot(v123,v234))
    return angle


def dihedral_distance(theta1, theta2, k):
    '''
    theta1, theta2 - the angles in radians
    k - multiplicity
    '''
    symk = 2*np.pi/k
    rem12 = np.remainder(theta1-theta2, symk)
    if isinstance(theta1, np.ndarray):
        return np.min([rem12, symk-rem12], axis=0).sum()
    return min(rem12, symk - rem12)


# s = sequence { ALA ARG ASN ASP CYS GLN GLU GLY HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL }

def get_movable(pivot, exclude, movable_list):
    movable_list.append(pivot.rindex)
    for atom in pivot.get_connected_atoms(only_intra=True):
        if (atom != exclude) and (pivot.parent == atom.parent) and (atom.rindex not in movable_list):
            get_movable(atom, exclude, movable_list)
    return movable_list
