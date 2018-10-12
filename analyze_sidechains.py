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
		s.sig.normalize()
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
				coords = res[mca].get_coord()
				xyz += coords
				n += 1
		return(Atom.Atom('C', xyz / n, 0, 0, '', 'center', 0))

	def getstrandsorient(s):
		orients = ""
		lastid = None
		for res in s.sresis:
			if lastid is None:
				lastid = res.id[1]
			ca = res['CA']
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
			orients += orient
			if abs(res.id[1]-lastid) > 1:
				orients += "\n"
			lastid = res.id[1]
		return orients

class IOSignature:
	def __init__(s, parent):
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
				s.freqs[k] = [v[0]/v[2], v[1]/v[2], v[2]]

	def __str__(s):
		rep = ""
		for k, v in s.freqs.items():
			rep += "key: {0}     value: {1}\n".format(k, v)
		return(rep)







				




# def compute_midas(ref_file, rad, traj_file, ion_id, ion_num):
# 	parser = PDBParser()
# 	sup = Superimposer()
# 	ref = parser.get_structure('ref', ref_file)
# 	ion = get_ion(ref, ion_id, ion_num)
# 	subset = get_CA_sphere(ref, ion, rad, ion_num)
# 	# print subset
# 	# return
# 	atm_list = ['CA', 'C', 'N', 'O']
# 	ref_atoms = get_atoms(ref[0], subset, atm_list)
# 	traj = parser.get_structure('traj', traj_file)
# 	count = 0
# 	for model in traj:
# 		model_atoms = get_atoms(model, subset, atm_list)
# 		sup.set_atoms(ref_atoms, model_atoms)
# 		print '{} {}'.format(count*10, sup.rms) # 10=step size in picosec
# 		count += 1

# def get_atoms(model, subset, atm_list):
# 	atoms = []
# 	for chain in model:
# 		for res in chain:
# 			if res.id[1] in subset:
# 				for atm in atm_list:
# 					atoms.append(res[atm])
# 	return atoms

# def get_ion(struct, ion_id, ion_num):
# 	model = struct[0]
# 	for chain in model:
# 		for res in chain:
# 			if res.id[1] == ion_num:
# 				return res[ion_id]
# 	print "did not find ion..."


# def compute_midas_pymol(ref_file, rad, ion_id, ion_num):
# 	parser = PDBParser()
# 	sup = Superimposer()
# 	ref = parser.get_structure('ref', ref_file)
# 	ion = get_ion(ref, ion_id, ion_num)
# 	subset = get_CA_sphere(ref, ion, rad, ion_num)
# 	print_pymol_cmd(subset)
# 	return

# def print_pymol_cmd(subset):
# 	cmd = "select mid_sp, "
# 	flag = True
# 	for i in subset:
# 		if flag:
# 			cmd = cmd + "i. {} ".format(i)
# 			flag = False
# 		else:
# 			cmd = cmd + "or i. {} ".format(i)
# 	print cmd
# 	return



if __name__ == "__main__":
	tester = BBarrel('bama.pdb', 'bama_trimmed.tmstrands')


	




























