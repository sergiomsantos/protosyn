import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from protosyn.features import AcceptorFeature, AromaticFeature, AromaticGenerator, get_generators
from protosyn.builder import grow, cap_chain
from protosyn.molecule import Molecule
from protosyn.ccd import CCD

import random


# read the ligand from a PDB file and extract the first residue
# (it is expected that the ligand contains only one single residue)
ligand = Molecule.LoadFromFile('caffeine.pdb').residues[0]

# define the features for the required ligand. In this case,
# two H-bond acceptors and two aromatic rings that can be involved
# in pi-pi stacking on bonth sides
features = [
    AcceptorFeature(acceptor=0, neighbor=8),
    AcceptorFeature(acceptor=1, neighbor=9),
    AromaticFeature(AromaticFeature.Side.BOTH, 3, 6, 7, 5, 10),
    AromaticFeature(AromaticFeature.Side.BOTH, 2, 7, 6, 8, 4, 9),
]

# generate all possible generator combinations from the above
# features. If permute=True, then all possible permutations
# within a given generator set will be generated.
generators_iterator = get_generators(features, permute=True)

# define the residues that will makeup the loops connecting
# the feature complements produced by the generators 
linkers = 6*['ala']

# go over all generators
for k,generator_set in enumerate(generators_iterator, start=1):

    # skip sets that yield aromatic rings in the same side
    # reduce(lambda x,y: x*y, [getattr(g,'side',1) for g in generator_set], 1):
    s1,s2 = [g.side for g in generator_set if isinstance(g, AromaticGenerator)]
    if s1*s2 > 0: 
        continue

    # count how may generators exist (each one will yield a single
    # complement - aminoacid residue)
    size = len(generator_set)

    # -------------------------------------------------------------------------
    # instantiate a new molecule
    # -------------------------------------------------------------------------
    title = '-'.join(g.rec_name for g in generator_set)
    peptide = Molecule(title)
    
    # -------------------------------------------------------------------------
    # now go over each generator in this generator set
    # -------------------------------------------------------------------------
    for n,generator in enumerate(generator_set, start=1):

        # get the complement it generates and add it to the
        # molecule object
        complement = generator(ligand)
        peptide.append_residue(complement, is_head=True)
        
        # if this is not a terminal complement, then grow
        # this chain
        if n < size:
            linkers = random.randint(3,6)*['ala'] + ['gly']
            grow(peptide, *linkers)
    
    # setup the backbone dihedrals for this peptide
    peptide.compile(setup_backbone_dihedrals=True)

    # -------------------------------------------------------------------------
    # close all chains
    # -------------------------------------------------------------------------
    chains = peptide.get_segments()
    for chain1,chain2 in zip(chains[:-1],chains[1:]):
        anchor = chain1[0]
        target = chain2[0]
        # use CCD to close the loop, holding in position the first residue
        # of chain1 and targetting the first residue of chain2
        #ccd = CCD(anchor, target, max_iter=200, threshold=0.5)
        #success, rmsd, iteration = ccd.run(peptide)

        # the last residue of chain1 is a dummy one that
        # its sole purpose is to calculate the distance to
        # the target position. Hence, it must be removed
        dummy = chain1.pop()
        peptide.remove_residue(dummy)

        # connect chain1 and chain2
        tail = chain1[-1]
        head = chain2[0]
        tail.next = head
        head.prev = tail
    
    # -------------------------------------------------------------------------
    # save peptide to PDB file
    # -------------------------------------------------------------------------

    # add terminal protons to peptide chain
    cap_chain(peptide)
    # renumber peptide atoms
    peptide.compile()
    # save peptide to file
    with open('out/complement_%d.pdb'%k, 'w') as fout:
        print >> fout, peptide.as_pdb(include_bonds=True)
        print str(k), 'Done creating sequence', ''.join(r.letter for r in peptide.iter_residues())
