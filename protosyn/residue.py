from dihedral import Dihedral, DihedralType
from atom import Atom
import dbase
import copy

import numpy as np
import utils
import maps

class Residue(object):
    
    #==========================================================================
    def __init__(self, name):
        self.name = name
        self.rindex = -1
        self.index = -1
        self.dihedrals = {}
        self._next = None
        self._prev = None
        self.atoms = []
        self.parent = None
        self.letter = maps.THREE_2_ONE.get(name, '?')
        
    def get_coordinates(self, *selection):
        if not selection:
            return np.array([at.xyz for at in self.atoms], float)
        return np.array([self[s].xyz for s in selection], float)
    
    def set_coordinates(self, values):
        for xyz,atom in zip(values, self.atoms):
            atom.xyz = xyz.copy()
    
    #==========================================================================
    def __str__(self):
        pn,pi = (self._prev.name,self._prev.index) if self._prev else (None,'?')
        nn,ni = (self._next.name,self._next.index) if self._next else (None,'?')
        return '<Residue:{0}:{1} {2}:{3}<-->{4}:{5}>'.format(self.name, self.index, pn, pi, nn, ni)
    
    def __repr__(self):
        return str(self)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.get_atom_by_name(key)
        return self.atoms[key]
    
    #==========================================================================
    def append_atom(self, atom):
        if atom not in self.atoms:
            atom.rindex = self.atoms[-1].rindex+1 if self.atoms else 0    
            self.atoms.append(atom)
            atom.parent = self

    # def remove_atom(self, atom):
    #     try:
    #         self.atoms.remove(atom)
    #     except ValueError:
    #         pass
    #     else:
    #         atom.parent = None
    #         for n,at in enumerate(self.atoms):
    #             at.rindex = n
            
        
    #==========================================================================
    def get_atom_by_name(self, name):
        if name.startswith('-'):
            res = self._prev
        elif name.startswith('+'):
            res = self._next
        else:
            res = self
        
        the_atom = None
        if res is not None:
            name = name.strip('+-')
            for atom in res.atoms:
                if atom.name == name:
                    the_atom = atom
                    break
        return the_atom
    
    #==========================================================================
    def compile(self, auto_setup=False, auto_bond=False, setup_sidechain_dihedrals=False, setup_backbone_dihedrals=False):
        if (auto_setup or auto_bond) and (self.name in dbase.BONDS):
            # get connection table
            table = dbase.BONDS.get(self.name)
            # clear out all previous connections
            for atom in self.iter_atoms():
                atom.clear_connected()
            # get new connections
            for atom in self.iter_atoms():
                for other in table.get(atom.name, []):
                    atom.connect_to(other if other[0] in '+-' else self.get_atom_by_name(other))

        if auto_setup or setup_sidechain_dihedrals:
            Dihedral.setup_sidechain_dihedrals(self)
            
        # if auto_setup or setup_backbone_dihedrals:
        #     Dihedral.setup_backbone_dihedrals(self)

    def has_next(self):
        return self._next is not None
    
    @property
    def next(self):
        return self._next
    
    @next.setter
    def next(self, value):
        self._next = value
        if value is None:
            self.dihedrals.pop(DihedralType.PSI, None)
        else:
            self.dihedrals[DihedralType.PSI] = Dihedral.get_dihedral(DihedralType.PSI, 'N','CA', 'C', '+N', self)
    
    def has_prev(self):
        return self._prev is not None
    
    @property
    def prev(self):
        return self._prev
    
    @prev.setter
    def prev(self, value):
        self._prev = value
        if value is None:
            self.dihedrals.pop(DihedralType.OMEGA, None)
            self.dihedrals.pop(DihedralType.PHI, None)
        else:
            self.dihedrals[DihedralType.OMEGA] = Dihedral.get_dihedral(DihedralType.OMEGA, '-CA','-C', 'N','CA', self)
            self.dihedrals[DihedralType.PHI] = Dihedral.get_dihedral(DihedralType.PHI, '-C', 'N','CA', 'C', self)
    
    #==========================================================================
    def iter_atoms(self):
        for atom in self.atoms:
            yield atom
    
    #==========================================================================
    def count_atoms(self):
        return len(self.atoms)
    
    #==========================================================================
    def as_pdb(self):
        return utils.as_pdb(self)
    
    def bonds_pdb(self):
        return utils.bonds_as_pdb(self)
        
    def copy(self):
        nxt = self._next
        prev = self._prev
        parent = self.parent
        self._next = self._prev = self.parent = None
        residue = copy.deepcopy(self)
        self.parent = parent
        self._prev = prev
        self._next = nxt
        
        return residue

