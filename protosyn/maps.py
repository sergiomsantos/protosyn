# import glob

# http://dock.compbio.ucsf.edu/DOCK_6/tutorials/amber_score/amber_score.htm


# lst='''Alanine  Ala     A
# Arginine        Arg     R
# Asparagine      Asn     N
# Aspartate       Asp     D
# Aspartate_Asparagine    Asx     B
# Cysteine        Cys     C
# Glutamate       Glu     E
# Glutamine       Gln     Q
# Glutamate_Glutamine     Glx     Z
# Glycine Gly     G
# Histidine       His     H
# Isoleucine      Ile     I
# Leucine Leu     L
# Lysine  Lys     K
# Methionine      Met     M
# Phenylalanine   Phe     F
# Proline Pro     P
# Serine  Ser     S
# Threonine       Thr     T
# Tryptophan      Trp     W
# Tyrosine        Tyr     Y
# Valine  Val     V'''

ONE_2_THREE = {
 'A': 'ALA',
# 'B': 'ASX',
 'C': 'CYS',
 'D': 'ASP',
 'E': 'GLU',
 'F': 'PHE',
 'G': 'GLY',
 'H': 'HIS',
 'I': 'ILE',
 'K': 'LYS',
 'L': 'LEU',
 'M': 'MET',
 'N': 'ASN',
 'P': 'PRO',
 'Q': 'GLN',
 'R': 'ARG',
 'S': 'SER',
 'T': 'THR',
 'V': 'VAL',
 'W': 'TRP',
 'Y': 'TYR',
 #'Z': 'GLX'
}

THREE_2_ONE = dict(((v,k) for k,v in ONE_2_THREE.items()))


FEATURE_COMPLEMENT = {
    'acceptor' : None,
    'aromatic' : {
        'phe' : ['CG','CD1','CD2','CE1','CE2','CZ'],
        'tyr' : ['CG','CD1','CD2','CE1','CE2','CZ'],
        'trp' : ['CG','CD1','CD2','NE1','CE2'],
    }
}