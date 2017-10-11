from context import protosyn
from protosyn.builder import grow, impose_secondary_structure, ConfLib
from protosyn.molecule import Molecule


peptide = Molecule('alpha')
grow(peptide, *(20*['ala']))

peptide.compile(setup_backbone_dihedrals=True)
impose_secondary_structure(peptide, slice(0, None), ConfLib.alpha)

print peptide.as_pdb()