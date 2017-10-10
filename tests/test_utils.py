from utils import dihedral_distance
import numpy as np
import unittest

class TestDihedralDistance(unittest.TestCase):
  def setUp(self):
    pass

  def test_0_pi_1(self):
    self.assertAlmostEqual(dihedral_distance(0, np.pi, 1), np.pi)

  def test_0_pi_2(self):
      self.assertAlmostEqual(dihedral_distance(0, np.pi, 2), 0)
  
  def test_0_pi_3(self):
      self.assertAlmostEqual(dihedral_distance(0, np.pi, 3), np.pi/3)


if __name__ == '__main__':
  unittest.main()