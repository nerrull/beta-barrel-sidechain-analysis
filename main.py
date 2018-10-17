from analyze_sidechains import read_groups, BBarrel, AMINOS, IOSignature
from os.path import join
import matplotlib.pyplot as plt
if __name__ == "__main__":
    groups = read_groups("barrel_groups.txt")
    group_results = {}
    for key, group in groups.items():
        if key ==6:continue#rejected in paper
	if key ==3:continue#oligomeric
        group_results[key] = []
        for protein in group:
            directory = "data/{}/".format(protein)
            pdb =join(directory, '{}.pdb'.format(protein))
            strands = join(directory, '{}.tmstrands'.format(protein))
            tester = BBarrel( pdb, strands)
            group_results[key].append(tester)

    signals = {}
    metaSignal= IOSignature()
    for group, group_proteins in group_results.items():
        groupSignal = IOSignature()
        aminoDict = {}
        for protein in group_proteins:
            groupSignal.merge(protein.sig)

        metaSignal.merge(groupSignal)
        groupSignal.normalize()
        signals[group] = groupSignal
    metaSignal.normalize()

    print (group_results)
