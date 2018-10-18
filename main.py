from analyze_sidechains import read_groups, BBarrel, AMINOS, IOSignature
from os.path import join
import matplotlib.pyplot as plt
import random
from analyze_sidechains import AMINOS, AMINOGROUPS

def bootstrap(n_iter, prot_list):
    metaSignal = dict()
    for a in AMINOS:
        metaSignal[a] = [None, None, None]
    bootstrapped_sigs = []
    for i in range(n_iter):
        for j in range(len(prot_list)):
            sig = IOSignature()
            r = random.randint(0, len(prot_list)-1)
            sig.merge(prot_list[r].sig)
        #sig.normalize()
        bootstrapped_sigs.append(sig)
    for a in AMINOS:
        mean = 0
        total = 0
        for b in bootstrapped_sigs:
            n = b.freqs[a][2]
            n_inside = b.freqs[a][0]
            if n > 0:
                mean += n_inside/n
            total += 1
        mean = mean/total
        var = 0
        for b in bootstrapped_sigs:
            n = b.freqs[a][2]
            n_inside = b.freqs[a][0]
            if n > 0:
                diff = mean-(n_inside/n)**2
                var += diff
        var = var / total
        std = var**0.5
        metaSignal[a] = [mean, std, total]
    return metaSignal

def bootstrap2(n_iter, prot_list):
    metaSignal = dict()
    for g in AMINOGROUPS:
        metaSignal[g] = [None, None, None]
    bootstrapped_sigs = []
    for i in range(n_iter):
        for j in range(len(prot_list)):
            sig = IOSignature()
            r = random.randint(0, len(prot_list)-1)
            sig.merge(prot_list[r].sig)
        #sig.normalize()
        bootstrapped_sigs.append(sig)
    for g in AMINOGROUPS.keys():
        print(g)
        mean = 0
        total = 0
        for a in AMINOGROUPS[g]:
            print(a)
            for b in bootstrapped_sigs:
                n = b.freqs[a][2]
                n_inside = b.freqs[a][0]
                if n > 0:
                    mean += n_inside/n
                total += 1
        mean = mean/total
        var = 0
        for a in AMINOGROUPS[g]:
            for b in bootstrapped_sigs:
                n = b.freqs[a][2]
                n_inside = b.freqs[a][0]
                if n > 0:
                    diff = mean-(n_inside/n)**2
                    var += diff
        var = var / total
        std = var**0.5
        metaSignal[g] = [mean, std, total]
    return (metaSignal, bootstrapped_sigs)





if __name__ == "__main__":
    groups = read_groups("barrel_groups.txt")
    group_results = {}
    prot_list = []
    for key, group in groups.items():
        if key ==6:continue #rejected in paper
        if key ==3:continue #oligomeric
        group_results[key] = []
        for protein in group:
            directory = "data/{}/".format(protein)
            pdb =join(directory, '{}.pdb'.format(protein))
            strands = join(directory, '{}.tmstrands'.format(protein))
            tester = BBarrel( pdb, strands)
            tester2 = BBarrel(pdb, strands)
            group_results[key].append(tester)
            prot_list.append(tester2)
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
    metaSig2 = bootstrap(10000, prot_list)
    metaSig3, bssigs = bootstrap2(10000, prot_list)

    print (metaSignal)
    print (metaSig2)
    print (metaSig3)
    with open("iosig_bootstrap_full_dataframe.txt", "w") as f:
        f.write("boot_iter amino freq_in number_in number_total\n")
        for i, b in enumerate(bssigs):
            for a in AMINOS:
                n = b.freqs[a][2]
                n_inside = b.freqs[a][0]
                if n > 0:
                    f.write("{0} {1} {2} {3} {4}\n".format(i, a, n_inside/n, n_inside, n))
                else:
                    f.write("{0} {1} {2} {3} {4}\n".format(i, a, 'NA', n_inside, n))



















