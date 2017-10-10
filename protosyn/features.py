from utils import rotation_about_vector, calc_principal_axis
from builder import FragmentProvider, fit

from itertools import permutations
from enum import Enum
import numpy as np


TWO_PI = 2* np.pi


class Feature(object):
    
    def get_generators(self):
        raise NotImplementedError()


class AcceptorFeature(Feature):
    COMPLEMENTS = {
        'gly': ['N','H'],
    }
    
    def __init__(self, acceptor, neighbor, dAX=1.0):
        self.a = acceptor
        self.x = neighbor
        self.dAX = dAX
    
    def get_generators(self):
        return [AcceptorGenerator(cname, cdonor, ch, self.a, self.x, self.dAX)
                    for cname,(cdonor,ch) in AcceptorFeature.COMPLEMENTS.iteritems()]


class AromaticFeature(Feature):
    COMPLEMENTS = {
        'phe' : ['CG','CD1','CE1','CZ','CE2','CD2'],
        'tyr' : ['CG','CD1','CE1','CZ','CE2','CD2'],
    }
    
    class Side(Enum):
        NEGATIVE = -1
        POSITIVE =  1
        BOTH = 2
    
    def __init__(self, side=Side.BOTH, *atoms):
        self.atoms = atoms
        self.side = side
    
    def get_generators(self):
        comps = self.COMPLEMENTS
        side = np.sign(self.side.value)

        generators = [AromaticGenerator(cname, catoms, self.atoms, side) for cname,catoms in comps.iteritems()]
        if self.side == self.Side.BOTH:
            generators += [AromaticGenerator(cname, catoms, self.atoms, -side) for cname,catoms in comps.iteritems()]
        return generators


class AcceptorGenerator(object):
    '''
      REC             LIG
    D-----H * * * * A-----X
    '''

    def __init__(self, rec_name, rec_d, rec_h, lig_a, lig_x, dAX=1.0):
        self.rec_name = rec_name
        self.rec_d = rec_d
        self.rec_h = rec_h
        self.lig_a = lig_a
        self.lig_x = lig_x
        self.dAX = dAX
    
    def __str__(self):
        return 'AcceptorGenerator:{rec_name} d[{rec_d}]--h[{rec_h}] ~~ a[{lig_a}]--x[{lig_x}]'.format(**self.__dict__)
        
    def __call__(self, ligand):
        vA = ligand.get_coordinates(self.lig_a)[0]
        vX = ligand.get_coordinates(self.lig_x)[0]
        vAX = vX - vA

        complement = FragmentProvider.get_fragment(self.rec_name)
        fit(ligand, [self.lig_x, self.lig_a], complement, [self.rec_h,self.rec_d])

        H = complement.get_coordinates(self.rec_h)[0]
        xyz = complement.get_coordinates() - H
        theta = np.random.uniform(TWO_PI)
        xyz = rotation_about_vector(xyz, theta, vAX) + vA - self.dAX*vAX/np.linalg.norm(vAX)
        complement.set_coordinates(xyz)
        
        return complement


class AromaticGenerator(object):
    '''
       2----------3
      /            \
     1              4
      \            /
       6----------5
    '''

    def __init__(self, rec_name, rec_atoms, lig_atoms, side):
        self.rec_name = rec_name
        self.rec_atoms = rec_atoms
        self.lig_atoms = lig_atoms
        self.side = side

    def __call__(self, ligand):
        
        # get ring normal
        _,ring_normal = calc_principal_axis(ligand, self.lig_atoms)[-1]

        # get ring center
        center = ligand.get_coordinates(*self.lig_atoms).mean(0)

        v12 = ligand[self.lig_atoms[2]].xyz - ligand[self.lig_atoms[1]].xyz
        v10 = ligand[self.lig_atoms[0]].xyz - ligand[self.lig_atoms[1]].xyz
        side = np.sign(np.dot(np.cross(v12, v10), ring_normal))
        side *= np.sign(self.side)
        offset = side * 3.0 * ring_normal

        complement = FragmentProvider.get_fragment(self.rec_name)
        fit(ligand, self.lig_atoms, complement, self.rec_atoms)

        xyz = complement.get_coordinates() - center
        xyz = rotation_about_vector(xyz, np.random.uniform(TWO_PI), ring_normal)

        complement.set_coordinates(xyz + center + offset)

        return complement


# class Feature(object):
#     # AROMATIC         = 1
#     # ACCEPTOR         = 2
#     # DONOR            = 3
#     # NEGIONIZABLE     = 4
#     # POSIONIZABLE     = 5
#     # HYDROPHOBE       = 6
#     # LUMPEDHYDROPHOBE = 7
#     # ZNBINDER         = 8



def get_generators(features, permute=True):
    all_generators = [[]]
    for feature in features:
        tmp = []
        for generator in feature.get_generators():
            tmp.extend([gen+[generator] for gen in all_generators])
        all_generators = tmp
    
    for generators in all_generators:
        if permute:
            for generator in permutations(generators):
                yield generator
        else:
            yield generators
