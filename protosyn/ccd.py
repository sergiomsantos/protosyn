from utils import rotation_about_vector, calculate_dihedral, closest_point_on_line
from dihedral import Dihedral, DihedralType
from numpy.linalg import norm
import numpy as np


# ver https://github.com/andre-wojtowicz/ccd-loop-closer/blob/master/CCDLoopCloser.py

class CCD:
    '''
    Cyclic Coordinate Descent
    '''
    def __init__(self, anchor_residue, target_residue, max_iter=100, threshold=0.5):
        
        self.anchor = anchor_residue
        self.target = target_residue
        
        self.max_iter  = max_iter
        self.threshold = threshold
    
    def _update_loop(self, loop, dihedrals):
        for residue,(phi,psi) in zip(loop,dihedrals.reshape(-1,2)):
            residue.dihedrals[DihedralType.PHI].radians = phi
            residue.dihedrals[DihedralType.PSI].radians = psi
    
    def _get_loop(self, polymer):
        
        if isinstance(self.anchor, int):
            self.anchor = polymer.residues[self.anchor]

        if isinstance(self.target, int):
            self.target = polymer.residues[self.target]

        loop = []
        cur = self.anchor
        while cur is not None:
            if cur is self.target:
                raise Exception('Target residue cannot belong to the chain containing the loop region')
            if cur is not self.anchor:
                loop.append(cur)
            cur = cur.next

        return loop


    def run(self, polymer, logger=None):

        loop = self._get_loop(polymer)

        bbone_coords = np.vstack([r.get_coordinates('N','CA','C') for r in loop])
        bbone_angles = np.array([[r.dihedrals[DihedralType.PHI].radians,
                                  r.dihedrals[DihedralType.PSI].radians] for r in loop[:-1]]).flatten()
        
        fixed_coords = self.target.get_coordinates('N','CA','C')
        
        iteration = 0
        loop_size = len(loop)
        while True:
            
            rmsd = np.sqrt(((fixed_coords - bbone_coords[-3:])**2).sum()/3)

            if logger:
                self._update_loop(loop, bbone_angles)
                logger(polymer)

            if (rmsd <= self.threshold) or (iteration > self.max_iter):
                if not logger:
                    self._update_loop(loop, bbone_angles)
                success = rmsd <= self.threshold
                return success, rmsd, iteration

            i = n = k = 0
            while n < loop_size-1:

                i = 3*n + k%2
                
                pivot = bbone_coords[i].copy()
                axis = bbone_coords[i+1] - pivot
                axis /= norm(axis)

                a = b = c = 0.0
                for j in range(-3,0):
                    Oj = closest_point_on_line(pivot, axis, bbone_coords[j])
                    Rj = bbone_coords[j] - Oj
                    Fj = fixed_coords[j] - Oj

                    rj = norm(Rj)
                    fj = norm(Fj)

                    vRj = Rj / rj
                    vSj = np.cross(vRj, axis)


                    #a += rj**2 + fj**2
                    b += 2 * rj * np.dot(Fj, vRj)
                    c += 2 * rj * np.dot(Fj, vSj)
                
                cos_alpha = b / np.sqrt(b**2 + c**2)
                sin_alpha = c / np.sqrt(b**2 + c**2)
                rot_angle = np.arctan2(sin_alpha, cos_alpha)
                bbone_angles[k] -= rot_angle

                tmp = bbone_coords[i:] - pivot
                bbone_coords[i:] = rotation_about_vector(tmp, -rot_angle, axis) + pivot

                k += 1
                if k%2 == 0:
                    n += 1
                
            iteration += 1

# class CCD:
#     '''
#     Cyclic Coordinate Descent
#     '''
#     def __init__(self, loop_start, loop_end, anchor, max_iter=100, threshold=0.5):
        
#         self.anchor = anchor
#         self.loop_end = loop_end
#         self.loop_start = loop_start
        
#         self.max_iter  = max_iter
#         self.threshold = threshold
    
#     def update(self, loop, dihedrals):
#         for residue,(phi,psi) in zip(loop,dihedrals.reshape(-1,2)):
#             residue.dihedrals[DihedralType.PHI].radians = phi
#             residue.dihedrals[DihedralType.PSI].radians = psi
    
#     def _setup(self, polymer):
#         cur = polymer[self.loop_start]
#         end = polymer[self.loop_end]
#         loop = []
#         while cur != None:
#             loop.append(cur)
#             if cur is end:
#                 break
#             if cur.next is None:
#                 raise Exception('loop_start is not connected to loop_end')
#             cur = cur.next
        
#         anchor = polymer[self.anchor]
#         if anchor in loop:
#             raise Exception('anchor cannot belong to loop')
#         return anchor, loop

#     def run(self, polymer, logger=None):
#         anchor, loop = self._setup(polymer)

#         # get movable atom coordinates
#         xloop = np.vstack([r.get_coordinates('N','CA','C') for r in loop])
#         angles = np.array([[r.dihedrals[DihedralType.PHI].radians,
#                             r.dihedrals[DihedralType.PSI].radians] for r in loop[:-1]]).flatten()
        
#         xanchor = anchor.get_coordinates('N','CA','C')

#         iteration = 0
#         loop_size = len(loop)
#         while True:
            
#             rmsd = np.sqrt(((xanchor - xloop[-3:])**2).sum()/3)


#             #print iteration, rmsd
#             if logger:
#                 self.update(loop, angles)
#                 logger(polymer)

#             if rmsd < self.threshold:
#                 if not logger:
#                     self.update(loop, angles)
#                 return True, rmsd, iteration

#             if iteration > self.max_iter:
#                 if not logger:
#                     self.update(loop, angles)
#                 return False, rmsd, iteration
            
#             i = n = k = 0
#             while n < loop_size-1:

#                 i = 3*n + k%2
                
#                 pivot = xloop[i].copy()
#                 axis = xloop[i+1] - pivot
#                 axis /= norm(axis)

#                 a = b = c = 0.0
#                 for j in range(-3,0):
#                     Oj = closest_point_on_line(pivot, axis, xloop[j])
#                     Rj = xloop[j] - Oj
#                     Fj = xanchor[j] - Oj

#                     rj = norm(Rj)
#                     fj = norm(Fj)

#                     vRj = Rj / rj
#                     vSj = np.cross(vRj, axis)


#                     #a += rj**2 + fj**2
#                     b += 2 * rj * np.dot(Fj, vRj)
#                     c += 2 * rj * np.dot(Fj, vSj)
                
#                 cos_alpha = b / np.sqrt(b**2 + c**2)
#                 sin_alpha = c / np.sqrt(b**2 + c**2)
#                 rot_angle = np.arctan2(sin_alpha, cos_alpha)
#                 angles[k] -= rot_angle

#                 tmp = xloop[i:] - pivot
#                 xloop[i:] = rotation_about_vector(tmp, -rot_angle, axis) + pivot

#                 k += 1
#                 if k%2 == 0:
#                     n += 1
                
#             iteration += 1

# def save(mol, mode='a'):
#     with open('loop.pdb', mode) as fout:
#         print >> fout, 'MODEL'
#         print >> fout, mol.as_pdb(include_bonds=mode=='w')
#         print >> fout, 'ENDMDL'

# def main():
#     # from dihedral import Dihedral
#     from molecule import Molecule
#     from sys import argv
    
#     fname = argv[1]

#     mol = Molecule.LoadFromFile(fname)
#     save(mol, 'w')
    
#     ccd = CCD(loop_start=24, loop_end=29, anchor=0)
#     ccd.threshold = 0.25
#     ccd.max_iter = 2000
#     success, rmsd, iteration = ccd.run(mol)

#     save(mol)

    
# if __name__ == '__main__':
    
#     main()
