import numpy as np
import dbase

class Atom(object):
    
    #==========================================================================
    def __init__(self, name, xyz=[0,0,0]):
        if len(name) == 4 and name[3].isdigit() and not name[0].isdigit():
            name = name[3]+name[:3]
        self.name = name
        self.parent = None
        self.iconnect = set()
        self.connect = set()
        self.symbol = name.strip('1234567890')[0]
        self.rindex = -1
        self.index = -1

        self.xyz = np.array(xyz, float)

        # self.mass = dbase.MASS.get(self.symbol, 0)
        # self.radius = dbase.RADII.get(self.symbol, 1.9)

    #==========================================================================
    def __str__(self):
        return 'Atom({0}, ri={1}, i={2})'.format(self.name, self.rindex, self.index)
    
    def __repr__(self):
        return str(self)
    
    #==========================================================================
    def connect_to(self, other):
        if other is not None:
            if isinstance(other, str):
                self.iconnect.add(other)
            elif other not in self.connect:
                self.connect.add(other)
                other.connect_to(self)
    
    def get_connected_atoms(self, only_intra=False):
        connected = list(self.connect)
        if not only_intra:
            connected += filter(None, [self.parent.get_atom_by_name(s) for s in self.iconnect])
        return connected
    
    def clear_connected(self):
        self.connect.clear()
        self.iconnect.clear()
    
    #==========================================================================
    @property
    def resnum(self):
        return getattr(self.parent, 'index', -1)
    
    #==========================================================================
    @property
    def resname(self):
        return getattr(self.parent, 'name', 'UNK')
    