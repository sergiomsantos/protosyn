
MASS = {
    'H':  1,
    'C':  6,
    'N':  7,
    'O':  8,
    'S': 16
}


#==========================================================================
BONDS = {
    'ALA' : {
        'N': ('H','CA',),
        'C': ('O','CA',),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','HB3'),
    },
    'ARG' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('HG1','HG2','CD'),
        'CD': ('HD1','HD2','NE'),
        'NE': ('HE','CZ'),
        'CZ': ('NH1','NH2'),
        'NH1': ('1HH1','2HH1'),
        'NH2': ('1HH2','2HH2'),
    },
    'ASN' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('OD1','ND2'),
        'ND2': ('1HD2','2HD2'),
    },
    'ASP' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('OD1','OD2'),
    },
    'CYS' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','SG'),
        'SG': ('HG',),
    }, 
    'GLN' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('HG1','HG2','CD'),
        'CD': ('OE1','NE2'),
        'NE2': ('1HE2','2HE2'),
    },
    'GLU' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('HG1','HG2','CD'),
        'CD': ('OE1','OE2'),
    },
    'GLY' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA1','HA2'),
    },
    'ILE' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB','CG1','CG2'),
        'CG1': ('1HG1','2HG1','CD'),
        'CG2': ('1HG2','2HG2','3HG2'),
        'CD': ('HD1','HD2','HD3'),
    },
    'LEU' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('HG','CD1','CD2'),
        'CD1': ('1HD1','2HD1','3HD1'),
        'CD2': ('1HD2','2HD2','3HD2'),
    },
    'LYS' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('HG1','HG2','CD'),
        'CD': ('HD1','HD2','CE'),
        'CE': ('HE1','HE2','NZ'),
        'NZ': ('HZ1','HZ2','HZ3'),
    },
    'MET' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('HG1','HG2','SD'),
        'SD': ('CE',),
        'CE': ('HE1','HE2','HE3'),
    },
    'PHE' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('CD1','CD2'),
        'CD1': ('HD1','CE1'),
        'CD2': ('HD2','CE2'),
        'CE1': ('HE1','CZ'),
        'CE2': ('HE2','CZ'),
        'CZ': ('HZ',),
    },
    'PRO' : {
        'N': ('CA','CD'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('HG1','HG2'),
        'CD': ('HD1','HD2','CG'),
    },
    'SER' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','OG'),
        'OG': ('HG',),
    },
    'THR' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB','CG2','OG1'),
        'OG1': ('HG1',),
        'CG2': ('1HG2','2HG2','3HG2'),
    },
    'TRP' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('CD1','CD2'),
        'CD1': ('HD1','NE1'),
        'CD2': ('CE2','CE3'),
        'NE1': ('HE1','CE2'),
        'CE3': ('HE3','CZ3'),
        'CE2': ('CD2','CZ2'),
        'CZ3': ('HZ3','CH2'),
        'CZ2': ('HZ2','CH2'),
        'CH2': ('HH2',),
    },
    'TYR' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('CD1','CD2'),
        'CD1': ('HD1','CE1'),
        'CD2': ('HD2','CE2'),
        'CE1': ('HE1','CZ'),
        'CE2': ('HE2','CZ'),
        'CZ': ('OH',),
        'OH': ('HH',),
    },
    'VAL' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB','CG1','CG2'),
        'CG1': ('1HG1','2HG1','3HG1'),
        'CG2': ('1HG2','2HG2','3HG2'),
    },
    'HIE' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('ND1','CD2'),
        'ND1': ('CE1',),
        'CD2': ('HD2','NE2'),
        'CE1': ('HE1','NE2'),
        'NE2': ('HE2',),
    },
    'HID' : {
        'N': ('H','CA'),
        'C': ('O','CA'),
        'CA': ('HA','CB'),
        'CB': ('HB1','HB2','CG'),
        'CG': ('ND1','CD2'),
        'ND1': ('HD1','CE1'),
        'CD2': ('HD2','NE2'),
        'CE1': ('HE1','NE2'),
        #'NE2': ('HE2',),
    }
}

def complete():
    
    for resname in BONDS.keys():
        
        BONDS[resname]['N'] += ('-C',)
        BONDS[resname]['C'] += ('+N',)

        nres = {}
        if resname=='PRO':
            nres.update(BONDS[resname], N=('H1','H2','CA','CD'))
        else:
            nres.update(BONDS[resname], N=('H1','H2','H3','CA'))
            
        cres = {}
        cres.update(BONDS[resname], C=('CA','OC1','OC2'))
        
        BONDS['C'+resname] = cres
        BONDS['N'+resname] = nres

    # for table in BONDS.itervalues():
    #     table['expected'] = set(table.keys()+[s for v in table.values() for s in v])

complete()

#print BONDS['NPHE']
# for table in BONDS.itervalues():
#     table['N'] += ('H1','H2','H3')
#     table['C'] += ('OC1','OC2')

# for table in BONDS.itervalues():
#     table.update({
#         # 'N': ('H','CA'),
#         # 'C': ('O','CA'),
#         'VC1':('C',),
#         'VC2':('C',),
#         'VN1':('N',),
#         'VN2':('N',),
#     })

# #==========================================================================
# from dihedral import DihedralType

# MULTIPLICITY = {
#     'ALA': {
#     },
#     'ARG': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#         DihedralType.CHI3: 1,
#         DihedralType.CHI4: 1,
#         DihedralType.CHI5: 2,
#     },

#     'ASN': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#     },
#     'ASP': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 2,
#     },
#     'CYS': {
#         DihedralType.CHI1: 1,
#     },
#     'GLN': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#         DihedralType.CHI3: 1,
#     },
#     'GLU': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#         DihedralType.CHI3: 2,
#     },
#     'GLY': {
#     },
#     'HIE': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#     },
#     'ILE': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#     },
#     'LEU': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#     },
#     'LYS': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#         DihedralType.CHI3: 1,
#         DihedralType.CHI4: 1,
#     },
#     'MET': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#         DihedralType.CHI3: 1,
#     },
#     'PHE': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 2,
#     },
#     'PRO': {
#     },
#     'SER': {
#         DihedralType.CHI1: 1,
#     },
#     'THR': {
#         DihedralType.CHI1: 1,
#     },
#     'TRP': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 1,
#     },
#     'TYR': {
#         DihedralType.CHI1: 1,
#         DihedralType.CHI2: 2,
#     },
#     'VAL': {
#         DihedralType.CHI1: 1,
#     },

# }


# def update():
#     # for table in MULTIPLICITY.itervalues():
#     for aa,table in MULTIPLICITY.items():
#         table.update({
#             DihedralType.OMEGA: 1,
#             DihedralType.PSI: 1,
#             DihedralType.PHI: 1
#         })
#         MULTIPLICITY['C'+aa] = table
#         MULTIPLICITY['N'+aa] = table


# update()




RADII = {
    'H' : 1.20,
    'D' : 1.20,
    'C' : 1.90,
    'N' : 1.50,
    'O' : 1.40,
    'F' : 1.35,
    'P' : 1.90,
    'S' : 1.85,
    'I' : 2.15,
    'CL': 1.80,
    'Fe': 0.64,
    'Cu': 1.28,
    'Zn': 1.38,
    'Br': 1.95,
# 'default': 1.90
}


if __name__ == '__main__':
    import json
    
    print json.dumps(BONDS, indent=2)
    print json.dumps(MULTIPLICITY, indent=2)