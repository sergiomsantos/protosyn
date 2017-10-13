from utils import rotation_about_vector, fit
from dihedral import DihedralType
from molecule import Molecule
from atom import Atom

import numpy as np
import re

from maps import ONE_2_THREE

#import markdown

def grow(chain, *frag_names):
    '''
    grow(chain, *fragment_names)

    Grow a peptide chain by appending fragments to the
    C-terminus of the chain.
    '''
    if not frag_names:
        return chain
    
    if chain.count_residues() == 0:
        fragment = FragmentProvider.get_fragment(frag_names[0])
        chain.append_residue(fragment, is_head=True)
        return grow(chain, *frag_names[1:])
    
    tail = chain.residues[-1]
    for n,frag_name in enumerate(frag_names, start=1):
        
        fragment = FragmentProvider.get_fragment(frag_name)
        
        if n == 1:
            R,T = fit(fragment, ('N','CA','C'), tail, ('N','CA','C'))
            invR = np.linalg.inv(R)

        # calculate new fragment's coordinates (flipped and translated),
        xyz = fragment.get_coordinates()*[1,(-1)**n,1] + [n*3.638,0,0]
        xyz = np.dot(invR, np.transpose(xyz - T)).T
        fragment.set_coordinates(xyz)
        chain.append_residue(fragment)

    # take the tail residue back to its original position
    xyz = tail.get_coordinates()
    xyz = np.dot(invR, np.transpose(xyz - T)).T
    tail.set_coordinates(xyz)

    return chain

def cap_chain(chain, cap_c=True, cap_n=True):
    
    if isinstance(chain, Molecule):
        for segment in chain.get_segments():
            cap_chain(segment, cap_c=cap_c, cap_n=cap_n)
        return

    if len(chain) < 2:
        raise Exception('Cannot cap chain with less than 2 residues')
    
    # build N terminus
    res = chain[0]
    if cap_n and not re.match(r'^N[A-Z]{3}$', res.name, re.I):
        H  = res.get_atom_by_name('H')
        N  = res.get_atom_by_name('N' )
        CA = res.get_atom_by_name('CA')
        x1 = rotation_about_vector(H.xyz - N.xyz,  2*np.pi/3, N.xyz - CA.xyz) + N.xyz
        x2 = rotation_about_vector(H.xyz - N.xyz, -2*np.pi/3, N.xyz - CA.xyz) + N.xyz
        
        H.name = 'H1'
        res.name = 'N'+res.name
        res.append_atom(Atom('H2', x1))
        res.append_atom(Atom('H3', x2))
        res.compile(auto_setup=True)

    # build C terminus
    res = chain[-1]
    if cap_c and not re.match(r'^C[A-Z]{3}$', res.name, re.I):
        O  = res.get_atom_by_name('O')
        C  = res.get_atom_by_name('C')
        CA = res.get_atom_by_name('CA')
        x1 = rotation_about_vector(O.xyz - C.xyz,  np.pi, C.xyz - CA.xyz) + C.xyz
        
        O.name = 'OC1'
        res.name = 'C'+res.name
        res.append_atom(Atom('OC2', x1))
        res.compile(auto_setup=True)
    res.parent.renumber()


class FragmentProvider:
    FRAGMENTS = {}

    @staticmethod
    def load_fragments():
        from residue import Residue
        import os

        ref = Residue('REF')
        ref.append_atom(Atom( 'N', [0.000,  0.300, 0.000]))
        ref.append_atom(Atom('CA', [1.181, -0.540, 0.000]))
        ref.append_atom(Atom( 'C', [2.440,  0.300, 0.000]))
        
        # load all residues
        this_dir, this_filename = os.path.split(__file__)
        DATA_PATH = os.path.join(this_dir, "static", "residues.pdb")
        
        residue_container = Molecule.LoadFromFile(DATA_PATH)
        for fragment in residue_container.iter_residues():
            fit(ref, ('N','CA','C'), fragment, ('N','CA','C'))
            FragmentProvider.FRAGMENTS[fragment.name] = fragment
        
    @staticmethod
    def get_fragment(name):
        # lazy loading of fragments
        if not FragmentProvider.FRAGMENTS:
            FragmentProvider.load_fragments()
            
        if len(name) == 1:
            name = ONE_2_THREE.get(name, '?')

        frag = FragmentProvider.FRAGMENTS.get(name.upper(), None)
        return frag if frag is None else frag.copy()


def impose_secondary_structure(peptide, selection, ss):
    residues = peptide.residues[selection]
    for dtype, value in ss.iteritems():
        for residue in residues:
            if dtype in residue.dihedrals:
                residue.dihedrals[dtype].degrees = value


class ConfLib(object):
    alpha = {
        DihedralType.PHI:    -57.0,
        DihedralType.PSI:    -47.0,
        DihedralType.OMEGA:  180.0
    }
    abeta = {
        DihedralType.PHI:   -139.0,
        DihedralType.PSI:    135.0,
        DihedralType.OMEGA: -178.0
    }
    pbeta = {
        DihedralType.PHI:   -119.0,
        DihedralType.PSI:    113.0,
        DihedralType.OMEGA:  180.0
    }

