from utils import rotation_about_vector, calculate_dihedral, closest_point_on_line
from dihedral import Dihedral, DihedralType
from numpy.linalg import norm
import numpy as np


# see https://github.com/andre-wojtowicz/ccd-loop-closer/blob/master/CCDLoopCloser.py

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
        for residue,(phi,psi,_) in zip(loop,dihedrals.reshape(-1,3)):
            residue.dihedrals[DihedralType.PHI].radians = phi
            residue.dihedrals[DihedralType.PSI].radians = psi
    
    def _get_loop(self, polymer):
        
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

        fixed_coords = self.target.get_coordinates('N','CA','C')
        bbone_coords = np.vstack([r.get_coordinates('N','CA','C') for r in loop])
        bbone_angles = np.array([[r.dihedrals[DihedralType.PHI].radians,
                                  r.dihedrals[DihedralType.PSI].radians, 0.0] for r in loop[:-1]]).flatten()
        
        iteration = 0
        loop_length = len(loop)
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
            
            indices = range(3*(loop_length-1))
            np.random.shuffle(indices)
            # for n in xrange(3*(loop_length-1)):
            # for n in xrange(3*(loop_length-1)-1,-1,-1):
            for n in indices:    
                # skip omega
                if (n+1)%3 == 0:
                    continue
                
                pivot = bbone_coords[n].copy()
                axis = bbone_coords[n+1] - pivot
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
                rot_angle = -np.arctan2(sin_alpha, cos_alpha)

                bbone_angles[n] += rot_angle
                tmp = bbone_coords[n:] - pivot
                bbone_coords[n:] = rotation_about_vector(tmp, +rot_angle, axis) + pivot

            iteration += 1