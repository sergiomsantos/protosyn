import numpy as np
from pathlib2 import Path
import gzip
import re


ONE_2_THREE = {
    'A': 'ALA',
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
}

THREE_2_ONE = dict(((v,k) for k,v in ONE_2_THREE.items()))

# query = re.compile(r'^ATOM[0-9 ]{7}  (CA|C |N )  (ILE|ALA) [A-Z]')
# query = re.compile(r'^ATOM[0-9 ]{7}  (CA|C |N )  [A-Z]{3} [A-Z]')
query_ca = re.compile(r'^ATOM[0-9 ]{7}  CA  [A-Z]{3} [A-Z]')

def norm(v):
    return np.sqrt(np.dot(v,v))

def atom_from_string(s):
    return s[12:16].strip(), \
           s[17:21].strip(), \
           int(s[22:26]), \
           np.array([s[i:i+8] for i in (30,38,46)], dtype=float)



def calculate_dihedral(x1,x2,x3,x4):
    v12 = x2 - x1
    v23 = x3 - x2
    v34 = x4 - x3
    v123 = np.cross(v12,v23)
    v234 = np.cross(v23,v34)

    angle = np.arctan2(np.dot(np.cross(v123,v234),v23)/np.sqrt(np.dot(v23,v23)),
                np.dot(v123,v234))
    return angle



def calculate_angle(x1, x2, x3):
    v21 = x1 - x2
    v23 = x3 - x2
    return np.arccos(np.dot(v21,v23) / (norm(v21)*norm(v23)))



class AlphaCarbon(object):
    __slots__ = ('name','xyz','letter')
    def __init__(self, name, xyz):
        self.letter = THREE_2_ONE.get(name, '?')
        self.name = name
        self.xyz = xyz
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return self.name



class Aminoacid(object):
    __slots__ = ('name','N','CA','C','letter')
    def __init__(self, name):
        self.name = name[-3:]
        self.letter = THREE_2_ONE.get(self.name, '?')
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return self.name


def read_pdb(fname, chain_id):
    aminoacids = []
    the_file = Path(fname)
    query = re.compile(r'^ATOM[0-9 ]{7}  CA  [A-Z]{3} %s'%chain_id.upper())
    prev = None
    if the_file.is_file():
        with gzip.open(fname, 'rb') as fin:
            for line in fin:
                if query.match(line):
                    
                    atname,resname,resid,xyz = atom_from_string(line)
                    if prev != resid:
                        aminoacids.append(AlphaCarbon(resname, xyz))
                        prev = resid
                #elif line.startswith('TER'):
                #    aminoacids.append(None)
    return aminoacids

def parse_file(fname):
    aminoacids = []
    the_file = Path(fname)
    if the_file.is_file():
        with gzip.open(fname, 'rb') as fin:
            for line in fin:
                if query_ca.match(line):
                    atname,resname,resid,xyz = atom_from_string(line)
                    aminoacids.append(AlphaCarbon(resname, xyz))
                elif line.startswith('TER'):
                    aminoacids.append(None)
    return aminoacids



# def parse_file(fname):
#     prev = None
#     aminoacids = []
#     with gzip.open(fname, 'rb') as fin:
#         for line in fin:
#             if query.match(line):
#                 atname,resname,resid,xyz = atom_from_string(line)
#                 if prev != resid:
#                     aa = Aminoacid(resname)
#                     aminoacids.append(aa)
#                     prev = resid
#                 setattr(aa, atname, xyz)
#             elif line.startswith('TER'):
#                 aminoacids.append(None)
#                 prev = None
#     return aminoacids


def calc_theta_tau(carbons, histograms_theta, theta_binw, tau_binw):
    
    for i in xrange(1,len(carbons)-2):
        c1 = carbons[i-1]
        c2 = carbons[i]
        c3 = carbons[i+1]
        #c4 = carbons[i+2]
        if c1 and c2 and c3: # and c4:
            try:
                theta =  calculate_angle(c1.xyz, c2.xyz, c3.xyz)
                # tau = calculate_dihedral(c1.xyz, c2.xyz, c3.xyz, c4.xyz)
            except:
                pass
            else:
                ktheta = c1.letter + c2.letter + c3.letter
                if '?' not in ktheta:

                    j = int(theta/theta_binw)
                    histograms_theta[ktheta][j] += 1

                # ktau = ktheta + c4.letter
                # j = int(tau/tau_binw)


def calculate_dihedrals(aminoacids):
    Nbins = 12
    map3 = {}
    w = 2*np.pi/Nbins
    for i in xrange(1,len(aminoacids)-1):
        prv = aminoacids[i-1]
        cur = aminoacids[i]
        nxt = aminoacids[i+1]
        if prv and cur and nxt:
            phi   = calculate_dihedral(prv.C,  cur.N,  cur.CA, cur.C)
            psi   = calculate_dihedral(cur.N,  cur.CA, cur.C,  nxt.N)
            omega = calculate_dihedral(cur.CA, cur.C,  nxt.N,  nxt.CA)
            key = prv.letter + cur.letter + nxt.letter
            if key not in map3:
                map3[key] = np.zeros(Nbins)
            j = int(v/w)
            map3[key][j] += 1
            #x.append(phi)
            #y.append(psi)
    return map3

# x = []
# y = []
# map2 = {}
# for i in xrange(1,len(aminoacids)-1):
#     prv = aminoacids[i-1]
#     cur = aminoacids[i]
#     nxt = aminoacids[i+1]
#     if prv and cur and nxt:
#         phi = calculate_dihedral(prv.C, cur.N, cur.CA, cur.C)
#         psi = calculate_dihedral(cur.N, cur.CA, cur.C, nxt.N)
#         omega = calculate_dihedral(cur.CA, cur.C, nxt.N, nxt.CA)
#         # print cur.letter, cur, phi,psi,omega
#         k2 = cur.letter+nxt.letter
#         if k2 not in map2:
#             map2[k2] = []
#         map2[k2].append(omega)
#         x.append(phi)
#         y.append(psi)



# print map2.keys()
# print len(map2.keys())
# plt.figure()
# plt.hold(True)
# for k,values in map2.iteritems():
#     hist, edges = np.histogram(np.rad2deg(values), bins=20, range=(-180,180))
#     plt.plot(edges[0:-1], hist)
#     #plt.hist(np.rad2deg(map2[k]), bins=20, range=(-180,180))
# # plt.plot(np.rad2deg(x), np.rad2deg(y), 'o')
# plt.xlabel('phi')
# plt.ylabel('psi')
# plt.show()
# # awk '/^ATOM[0-9 ]{7}  (CA|C |N )  (ILE|ALA)/' 4z7e.pdb 
# # awk '/^ATOM[0-9 ]{7}  (CA|C |N )  (ILE|ALA)/ {print $3,$4,$5,substr($0,31,8),substr($0,39,8),substr($0,47,8)}' 4z7e.pdb | head

def score(sequence):
    keys = sorted(THREE_2_ONE.values())
    # print keys
    mass = dict(((k,float(n)) for n,k in enumerate(keys, start=1)))
    # print mass
    scores = []
    size = len(sequence)
    mat = np.zeros((size,size),float)
    for i,pivot in enumerate(sequence):
        # score = 11 + 8*mass[pivot]
        
        score = mass[pivot]
        mat[i,i] = score
        left = sequence[:i][::-1]
        right = sequence[i+1:]
        # print i,left,pivot,right
        for j,s in enumerate(left,start=2):
            # print '   ', j, s
            # score += j*11 + 8*mass[s]
            score += mass[s]/j
            mat[i,i-j+1] = mass[s]/j
        for j,s in enumerate(right,start=2):
            score += mass[s]/j
            mat[i,i+j-1] = mass[s]/j
        # print score
        scores.append(score)
    return scores, mat

def get_protein_id(fname):
    with open(fname) as fin:
        prots = [l.split()[0] for l in fin if  l.startswith('>') and 'mol:protein ' in l]
    return list(set([p.split('_')[0].replace('>','') for p in prots]))


def list_dir(path):
    import glob
    dirs = glob.glob(path)
    return [d.split('/')[1] for d in dirs]

def main2():
    import matplotlib.pyplot as plt
    MAX_SIZE = 300
    #>101m_A mol:protein length:154  MYOGLOBIN
    #MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG
    
    #>108l_A mol:protein length:164  T4 LYSOZYME
    #MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKIELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL
    myo_seq  = 2*'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG'
    lyso_seq = 'MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKIELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL'
    #myo_seq += myo_seq[::-1]
    myo,myo_mat = score(myo_seq)
    lyso,lyso_mat = score(lyso_seq)

    #xpt = 200
    #merge = score(myo_seq[:xpt]+lyso_seq[xpt:])
    # x = np.arange(len(scores))
    # y = np.array(scores)
    # plt.plot(np.arange(len(myo)), myo, np.arange(len(lyso)), lyso, np.arange(len(merge)), merge)
    x = np.arange(len(myo))
    plt.hold(True)
    As = [n for n,s in enumerate(myo_seq) if s=='A']
    for n in As:
        a = myo_mat[n]
        plt.plot(x, n*0.1 + a)
    # for n,a in enumerate(myo_mat[::5]):
    #     print n, n*0.1
    #     plt.plot(x, n*0.1 + 5*a)
    plt.plot(x, np.array(myo))
    
    # x = np.linspace(-10, MAX_SIZE, 10*MAX_SIZE)
    # y = np.zeros(10*MAX_SIZE, float)
    # plt.figure()
    # plt.hold(True)
    # for n,s in enumerate(scores):
    #     yy = s*np.exp(-10*(x-n)**2)
    #     plt.plot(x,yy,'k')
    #     y += s*np.exp(-10*(x-n)**2)
    # plt.stem(range(len(scores)), scores)
    # plt.plot(x,y)
    plt.show()

def main():
    from sys import stdout
    import itertools
    import json

    available_dirs = list_dir('pdb/*')
    prot_ids = get_protein_id('pdb_seqres.txt')
    files = ['pdb/%s/pdb%s.ent.gz'%(p[1:3],p) for p in prot_ids if p[1:3] in available_dirs]
    
    # parameters
    Nbins = 12
    theta_binw = np.pi/Nbins
    tau_binw = np.pi/Nbins

    # initialize histograms
    histograms_theta = {}
    letters = THREE_2_ONE.values()
    for a,b,c in itertools.product(letters,letters,letters):
        histograms_theta[a+b+c] = np.zeros(Nbins)
    
    Nfiles = len(files)
    print "Nfiles=",Nfiles

    for n,fname in enumerate(files):
        p = 100.0*n/Nfiles
        stdout.write('%-100s %.2f %% \r' % (int(p)*'=',p))
        # stdout.write('processing %.2f %% \r' % (100.0*n/Nfiles))
        stdout.flush()
        try:
            carbons = parse_file(fname)
        except Exception as e:
            print fname
            print e
        else:
            calc_theta_tau(carbons, histograms_theta, theta_binw, tau_binw)
    
    out = {}
    for k,v in histograms_theta.iteritems():
        #print k, v
        out[k] = v.tolist()
    print json.dumps(out, indent=2)

    
    import pickle
    with open('data.pkl', 'wb') as fout:
        pickle.dump(histograms_theta, fout)
    

    # Nbins = 12
    # bin_w = 2*np.pi/Nbins
    # maps = {}
    # for fname in files[:100]:
    #     aas = parse_file(fname)
    #     calculate_dihedrals(aas, maps, Nbins, bin_w)
    # print maps

if __name__ == '__main__':
    main2()