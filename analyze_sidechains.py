from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import *
import glob
import sys

MCATOMS = ['CA', 'C', 'N']
AMINOS = ['PRO', 'GLY', 'ALA', 'VAL', 'TRP',
          'TYR', 'ASP', 'ASN', 'GLN', 'GLU',
          'SER', 'THR', 'ILE', 'LEU', 'PHE',
          'HIS', 'LYS', 'ARG', 'MET', 'CYS']


class BBarrel:
    def __init__(s, filename, strandsfile):
        parser = PDBParser()
        s.stru = parser.get_structure('struct', filename)[0]
        s.sti = s.readstrands(strandsfile)
        s.sresis = s.get_strands_resis()
        s.c = s.get_center()
        s.sig = IOSignature(s)
        s.orients = s.getstrandsorient()
        # s.sig.normalize()
        print(s.sig)

    def readstrands(s, strandsfile):
        sti = []
        with open(strandsfile) as f:
            lines = f.readlines()
        for l in lines:
            ll = l.split()
            if len(ll) == 2:
                sti.append((int(ll[0]), int(ll[1])))
        return sti

    def get_strands_resis(s):
        subset = []
        for chain in s.stru:
            print(chain)
            for res in chain:
                for st in s.sti:
                    if res.id[1] >= st[0] and res.id[1] <= st[1]:
                        subset.append(res)
        return subset

    def get_center(s):
        xyz = (0, 0, 0)
        n = 0
        for res in s.sresis:
            for mca in MCATOMS:
                try:
                    coords = res[mca].get_coord()
                    xyz += coords
                    n += 1
                except KeyError as e:
                    print(e)
                    continue
        return (Atom.Atom('C', xyz / n, 0, 0, '', 'center', 0))

    def getstrandsorient(s):
        orients = ""
        lastid = None
        for res in s.sresis:
            if lastid is None:
                lastid = res.id[1]
            ca = res['CA']
            try:
                if res.get_resname() == 'GLY':
                    orient = 'G'
                else:
                    cb = res['CB']
                    da = abs(ca - s.c)
                    db = abs(cb - s.c)
                    if da > db:
                        orient = 'I'
                        s.sig.add_inner(res.get_resname())
                    else:
                        orient = 'O'
                        s.sig.add_outer(res.get_resname())
            except:
                orient = 'G'
            finally:
                orients += orient

            if abs(res.id[1] - lastid) > 1:
                orients += "\n"
            lastid = res.id[1]
        return orients


class IOSignature:
    def __init__(s, parent= None):
        s.p = parent
        s.freqs = dict()
        for a in AMINOS:
            s.freqs[a] = [0, 0, 0]

    def add_inner(s, name):
        s.freqs[name][0] += 1
        s.freqs[name][2] += 1

    def add_outer(s, name):
        s.freqs[name][1] += 1
        s.freqs[name][2] += 1

    def normalize(s):
        for k, v in s.freqs.items():
            if v[2] > 0:
                s.freqs[k] = [v[0] / float(v[2]), v[1] / float(v[2]), float(v[2])]

    def merge(self, other):
        for a in AMINOS:
            for i, f in enumerate (self.freqs[a]):
                self.freqs[a][i] += other.freqs[a][i]

    def __str__(s):
        rep = ""
        for k, v in s.freqs.items():
            rep += "key: {0}     value: {1}\n".format(k, v)
        return (rep)


def read_groups(groupfile):
    groups = dict()
    with open(groupfile) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        ll = l.split()
        ll = l.split()
        groups[i + 1] = ll

    return groups


def remove_shit_renumber(infile, outfile):
    with open(infile) as f:
        lines = f.readlines()
    newatomnum = 1
    resnum = None
    lastresnum = None
    newresnum = 1
    with open(outfile, "w") as f:
        for l in lines:
            rec = l[0:4]
            resname = l[17:20]
            if rec == 'ATOM' and resname in AMINOS:
                resnum = int(l[22:26])
                if lastresnum is None:
                    lastresnum = resnum
                if resnum != lastresnum:
                    newresnum += 1
                f.write(l[0:6] + '{:5d}'.format(newatomnum) + l[11:22] + '{:4d}'.format(newresnum) + l[26:])
                lastresnum = resnum
                newatomnum += 1

def remove_shit_norenumber(infile, outfile):
    with open(infile) as f:
        lines = f.readlines()
    newatomnum = 1
    resnum = None
    lastresnum = None
    newresnum = 1
    with open(outfile, "w") as f:
        for l in lines:
            rec = l[0:4]
            resname = l[17:20]
            if rec == 'ATOM' and resname in AMINOS:
                resnum = int(l[22:26])
                if lastresnum is None:
                    lastresnum = resnum
                if resnum != lastresnum:
                    newresnum += 1
                f.write(l)
                lastresnum = resnum
                newatomnum += 1
