import os
from analyze_sidechains import remove_shit_norenumber, AMINOS
os.chdir("data")

keep_extensions = [".pdb", ".strands", ".tmstrands", ".sidechain"]

files =  [os.path.join(dp, f) for dp, dn, fn in os.walk(os.curdir) for f in fn]
print (files)
for file in files:
    if ".pdb" in file:
        remove_shit_norenumber(file, file)
        continue

    elif any(ext in file for ext in keep_extensions):
        continue

    else:
        # print(file)
        os.remove(file)

