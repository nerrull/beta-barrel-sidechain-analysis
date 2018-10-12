import os

os.chdir("data")

keep_extensions = [".pdb", ".strands", ".tmstrands", ".sidechain"]

files =  [os.path.join(dp, f) for dp, dn, fn in os.walk(os.curdir) for f in fn]
print (files)
for file in files:
    if any(ext in file for ext in keep_extensions):
        continue
    else:
        # print(file)
        os.remove(file)

