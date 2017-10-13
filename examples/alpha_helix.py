from context import protosyn
from protosyn.builder import grow, impose_secondary_structure, ConfLib
from protosyn.molecule import Molecule
from protosyn.ccd import CCD

import os

OUTPUT_FOLDER = './out'

# make sure the output folder exists
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)

# instantiate a new molecule 
peptide = Molecule('alpha')

# define the peptide sequence
seq = 10*['ala']

# build the peptide
grow(peptide, *seq)

# impose the secondary structure
impose_secondary_structure(peptide, slice(0, None), ConfLib.alpha)

# for r in peptide:
#     pass

# save the peptide
with open(OUTPUT_FOLDER + '/alpha.pdb', 'w') as fout:
    print >> fout, 'MODEL'
    print >> fout, peptide.as_pdb(include_bonds=True)
    print >> fout, 'ENDMDL'
    
    # xyz = peptide.get_coordinates()
    # for n in range(10):

    #     ccd = CCD(peptide.residues[0], peptide.residues[0], max_iter=500)
    #     ccd.run(peptide)
    #     print >> fout, 'MODEL'
    #     print >> fout, peptide.as_pdb()
    #     print >> fout, 'ENDMDL'
    #     peptide.set_coordinates(xyz)