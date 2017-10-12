from enum import Enum
import numpy as np
import utils
import re

class DihedralType(Enum):
    OMEGA = -2
    PHI = -1
    PSI  = 0
    CHI1 = 1
    CHI2 = 2
    CHI3 = 3
    CHI4 = 4
    CHI5 = 5
    CHI6 = 6
    
    @staticmethod
    def get_sc():
        return [getattr(DihedralType, 'CHI%d'%n) for n in range(1,7)]
    
    @staticmethod
    def get_bb():
        return [DihedralType.OMEGA, DihedralType.PHI, DihedralType.PSI]



class Dihedral(object):
    
    def __init__(self, dtype, a1, a2, a3, a4, movable, parent):
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.type = dtype
        self.parent = parent
        self.movable = sorted(movable)
        
        self._value = 0.0
        #if parent:
        self.calculate()
        
    #==========================================================================
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return ' Dihedral(%s, %s, %s, %s, %s, %s, %s) = %8.3f'%(self.type.name,
                        self.a1, self.a2, self.a3, self.a4, str(self.movable), self.parent, self.degrees)
    
    def calculate(self):
        at1 = self.parent.get_atom_by_name(self.a1)
        at2 = self.parent.get_atom_by_name(self.a2)
        at3 = self.parent.get_atom_by_name(self.a3)
        at4 = self.parent.get_atom_by_name(self.a4)
        self._value = utils.calculate_dihedral(at1.xyz, at2.xyz, at3.xyz, at4.xyz)

    #==========================================================================
    # def _set_dihedral(self, value):
    #     residue = self.parent
    #     at2 = residue.get_atom_by_name(self.a2)    
    #     at3 = residue.get_atom_by_name(self.a3)
    #     origin = at2.xyz
    #     axis   = at3.xyz - origin
    #     angle  = value - self._value
    #     xyz    = residue.get_coordinates()
    #     xmov = xyz[self.movable] - origin
    #     xyz[self.movable] = utils.rotation_about_vector(xmov, angle, axis) + origin
    #     residue.set_coordinates(xyz)

    #     if self.type.value < DihedralType.CHI1.value:
    #         residue = residue.next
    #         while residue != None:
    #             xyz = residue.get_coordinates()
    #             residue.set_coordinates(utils.rotation_about_vector(xyz-origin, angle, axis) + origin)
    #             residue = residue.next
        
    #     self._value = value #np.remainder(np.pi+angle, 2*np.pi)-np.pi


    #==========================================================================
    def set_degrees(self, value):
        # self._set_dihedral(np.radians(value))
        self.set_radians(np.radians(value))
        
    def get_degrees(self):
        # return np.rad2deg(self._value)
        return np.rad2deg(self.get_radians())
    degrees = property(get_degrees, set_degrees)
    
    #==========================================================================
    def set_radians(self, value):
        #self._set_dihedral(value)
        current = self.get_radians()
        residue = self.parent
        at2 = residue.get_atom_by_name(self.a2)
        at3 = residue.get_atom_by_name(self.a3)
        origin = at2.xyz
        axis   = at3.xyz - origin
        angle  = value - current
        xyz    = residue.get_coordinates()
        xmov = xyz[self.movable] - origin
        xyz[self.movable] = utils.rotation_about_vector(xmov, angle, axis) + origin
        residue.set_coordinates(xyz)

        if self.type.value < DihedralType.CHI1.value:
            residue = residue.next
            while residue != None:
                xyz = residue.get_coordinates()
                residue.set_coordinates(utils.rotation_about_vector(xyz-origin, angle, axis) + origin)
                residue = residue.next

    def get_radians(self):
        #return self._value
        at1 = self.parent.get_atom_by_name(self.a1)
        at2 = self.parent.get_atom_by_name(self.a2)
        at3 = self.parent.get_atom_by_name(self.a3)
        at4 = self.parent.get_atom_by_name(self.a4)
        return utils.calculate_dihedral(at1.xyz, at2.xyz, at3.xyz, at4.xyz)
    radians = property(get_radians, set_radians)
    

    #==========================================================================
    @classmethod
    def get_dihedral(cls, dtype, a1, a2, a3, a4, residue):
        at2 = residue.get_atom_by_name(a2)
        at3 = residue.get_atom_by_name(a3)
        movable = utils.get_movable(at3, at2, [])
        return cls(dtype, a1, a2, a3, a4, movable, residue)

    # @classmethod
    # def setup_backbone_dihedrals(cls, residue):
    #     dihedrals = residue.dihedrals
    #     DType = DihedralType

    #     for dtype in DihedralType.get_bb():
    #         dihedrals.pop(dtype, None)
        
    #     if residue.prev is not None:
    #         dihedrals[DType.OMEGA] = cls.get_dihedral(DType.OMEGA, '-CA','-C', 'N','CA', residue)
    #         dihedrals[DType.PHI  ] = cls.get_dihedral(DType.PHI,    '-C', 'N','CA', 'C', residue)
        
    #     if residue.next is not None:
    #         dihedrals[DType.PSI] = cls.get_dihedral(DType.PSI, 'N','CA', 'C', '+N', residue)


    @classmethod
    def setup_sidechain_dihedrals(cls, residue):
        
        rname = residue.name.upper()

        # if this is proline, do nothing
        if rname == 'PRO':
            return
        
        # clear out all previous dihedrals (might no longer exist)
        dihedrals = residue.dihedrals
        for dtype in DihedralType.get_sc():
            dihedrals.pop(dtype, None)
        
        # determine the possible side chain sequences
        possible_queries = 'BGD' if rname in ('HIE','PHE','TRP','TYR') else 'BGDEZH'

        # format atom names to apply regex all at a time
        atnames = ':' + ':'.join(at.name for at in residue.iter_atoms()) + ':'
        
        path = ['C','CA']
        for chi,s in enumerate(possible_queries):
            match = re.match(r'(.*):(?P<name>[A-GI-Z]%s[1]?):(.*)'%s, atnames, re.I)
            if match:
                path.append(match.group('name'))
            if len(path) == 4:
                dtype = getattr(DihedralType, 'CHI%d'%chi)
                dihedrals[dtype] = cls.get_dihedral(dtype, path[0], path[1], path[2], path[3], residue)
                path.pop(0)


# if __name__ == '__main__':
#     d = Dihedral(DihedralType.OMEGA, 'a1','a2','a3','a4',[1,2,3],None)
#     print d
#     print DihedralType.get_bb()
#     print DihedralType.get_sc()
