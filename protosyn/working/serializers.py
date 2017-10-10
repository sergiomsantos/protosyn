import numpy as np


class Serializer:
    def __init__(self, model):
        self.dtypes = [d for res in model for d in sorted(res.dihedrals.keys())]
    
    def serialize(self, dihedrals):
        return np.array([dset[dtype] for dset in dihedrals for dtype in sorted(dset.keys())])
    
    def deserialize(self, values):
        out = []
        prev = None
        for dtype,value in zip(self.dtypes, values):
            if (prev is None) or (dtype <= prev):
                out.append({})
                prev = dtype
            out[-1][dtype] = value
        return out
