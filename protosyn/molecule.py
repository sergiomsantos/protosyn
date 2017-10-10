import numpy as np

from utils import fit

class Molecule(object):
    def __init__(self, name):
        self.residues = []
        self.name = name
        self.size = 0

    
    def get_coordinates(self):
        return np.array([at.xyz for at in self.iter_atoms()])
    
    def set_coordinates(self, values):
        if values.shape == (self.size,3):
            for n,at in enumerate(self.iter_atoms()):
                at.xyz = values[n].copy()
        else:
            raise Exception('Invalid shape for coordinate matrix: should be (%d,3)!'%self.size)
    
    def __str__(self):
        return '<Molecule:%s  Nresidues=%d>\n  '%(self.name, len(self.residues)) +\
                '\n  '.join((str(r) for r in self.residues))
    
    def __repr__(self):
        return str(self)
    
    def append_residue(self, residue, is_head=False):

        # remove the residue from its previous parent (if any)
        parent = residue.parent
        if parent:
            parent.remove_residue(residue)
        
        if self.residues:
            prev = self.residues[-1]

            if is_head:
                residue.next = residue.prev = None
            else:
                prev.next = residue
                residue.prev = prev
                residue.next = None
        else:
            residue.next = residue.prev = None
        
        self.residues.append(residue)
        residue.parent = self


    def remove_residue(self, residue):
        if residue.parent is self:
            if residue.prev:
                residue.prev.next = residue.next
            if residue.next:
                residue.next.prev = residue.prev
            residue.next = residue.prev = residue.parent = None
            try:
                self.residues.remove(residue)
            except ValueError:
                pass


    # def replace_residue(self, old_residue, new_residue):
    #     if old_residue.parent is self:
    #         assert new_residue.parent is None

    #         fit(old_residue, ('N','CA','C'), new_residue, ('N','CA','C'))

    #         new_residue.parent = self
    #         if old_residue.next:
    #             old_residue.next.prev = new_residue
    #             new_residue.next = old_residue.next
    #         if old_residue.prev:
    #             old_residue.prev.next = new_residue
    #             new_residue.prev = old_residue.prev

    #         i = self.residues.index(old_residue)
    #         self.remove_residue(old_residue)
    #         self.residues.insert(i, new_residue)
    #         new_residue.compile(setup_sidechain_dihedrals=True)


    def compile(self, auto_setup=False, auto_bond=False, setup_sidechain_dihedrals=False, setup_backbone_dihedrals=False):
        atm_counter = 0
        for i,residue in enumerate(self.iter_residues()):
            for j,atom in enumerate(residue.iter_atoms()):
                atom.index = atm_counter
                atom.rindex = j
                atm_counter += 1
            residue.index = i
            residue.compile(auto_setup, auto_bond, setup_sidechain_dihedrals, setup_backbone_dihedrals)
        self.size = atm_counter
    

    def as_pdb(self, include_bonds=False, *remarks, **kwargs):
        lines = ['TITLE  ' + self.name]
        if remarks:
            lines += ['REMARK '+remark for remark in remarks]
        lines += [r.as_pdb() for r in self.residues]
        if include_bonds:
            lines += [r.bonds_pdb() for r in self.residues]
        return '\n'.join(lines)
        
    def iter_atoms(self):
        for residue in self.residues:
            for atom in residue:
                yield atom
    
    def iter_residues(self):
        for residue in self.residues:
            yield residue
    
    @staticmethod
    def LoadFromFile(fname):
        from residue import Residue
        from atom import Atom

        def parse_pdb_line(s):
            # atname, resname, resid, xyz
            return s[12:16].strip(), \
                   s[17:21].strip(), \
                   int(s[22:26]), \
                   np.array([s[i:i+8] for i in (30,38,46)], dtype=float)
        rindex = 0
        pivot = None
        atoms = [None]
        is_head = False
        molecule = Molecule('mol')
        with open(fname, 'r') as fin:
            for line in fin:
                if line.startswith('ATOM'):
                    atname,resname,resid,xyz = parse_pdb_line(line)
                    if is_head or (resid != rindex):
                        residue = Residue(resname)
                        molecule.append_residue(residue, is_head)
                        is_head = False
                        rindex = resid
                    atom = Atom(atname, xyz)
                    atoms.append(atom)
                    residue.append_atom(atom)
                elif line.startswith('TER'):
                    is_head = True
                elif line.startswith('CONECT'):
                    idxs = map(int, line[6:].split())
                    pivot = atoms[idxs[0]]
                    for i in idxs[1:]:
                        pivot.connect_to(atoms[i])
        
        auto_setup = pivot is None
        molecule.compile(auto_setup=auto_setup)
        return molecule
        

    def get_segments(self):
        segments = []
        for residue in self.iter_residues():
            if residue.prev is None:
                segment = []
                segments.append(segment)
            segment.append(residue)
        return segments

    def count_segments(self):
        return sum((1 for r in self.iter_residues() if r.prev is None))

    def count_residues(self):
        return len(self.residues)
    
    def count_atoms(self):
        return self.size


if __name__ == '__main__':
    molecule = Molecule.LoadFromFile('assets/loop.pdb')
    print molecule.as_pdb()