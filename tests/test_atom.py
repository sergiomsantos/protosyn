import unittest
from context import protosyn

class TestAtomResidue(unittest.TestCase):
    
    def setUp(self):
        self.atom = protosyn.atom.Atom('A')
    
    def test_resname_is_unk(self):
        self.assertEquals(self.atom.resname, 'UNK')
    
    def test_resnum_is_neg1(self):
        self.assertEquals(self.atom.resnum, -1)
    
    def test_parent_is_none(self):
        self.assertEquals(self.atom.parent, None)
    
class TestAtomIdentifiers(unittest.TestCase):
    def setUp(self):
        self.atom = protosyn.atom.Atom('HG12')
    
    def test_name_is_2HG1(self):
        self.assertEquals(self.atom.name, '2HG1')
    
    def test_rindex_is_neg1(self):
        self.assertEquals(self.atom.rindex, -1)
    
    def test_index_is_neg1(self):
        self.assertEquals(self.atom.index, -1)
        
class TestAtomConnections(unittest.TestCase):
    def setUp(self):
        self.atom1 = protosyn.atom.Atom('A1')
        self.atom2 = protosyn.atom.Atom('A2')
        self.atom1.connect_to(self.atom2)

    def test_at1_connected_to_at2(self):
        self.assertIn(self.atom1, self.atom2.connect)
    
    def test_at2_connected_to_at1(self):
        self.assertIn(self.atom2, self.atom1.connect)
    
#   def test_0_pi_1(self):
#     self.assertAlmostEqual(dihedral_distance(0, np.pi, 1), np.pi)

#   def test_0_pi_2(self):
#       self.assertAlmostEqual(dihedral_distance(0, np.pi, 2), 0)
  
#   def test_0_pi_3(self):
#       self.assertAlmostEqual(dihedral_distance(0, np.pi, 3), np.pi/3)

if __name__ == '__main__':
    unittest.main()