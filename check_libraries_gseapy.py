# -*- coding: utf-8 -*-
"""RNAseq_analysis_postDE.ipynb

This part of code helps the user to see which library is present in the gsepy so they can add these to the parameters file
If the user already know which library to use, this is not required

"""


import gseapy as gp

"""View the library in gseapy"""
gp_version = gp.__version__

names = gp.get_library_name()
with open("library_names_gseapy.csv","w") as f:
  f.write(f"Library names GSEApy version {gp_version}\nThere are a total of {len(names)} libraries\n\n")
  for values in names:
    f.write(f'{values}\n')


"""Look at GO libraries only"""
print(f"\n   Displaying GO libraries only\n")
for x in names:
    if x.startswith("GO_"):
        print (x)

print(f"\n   Displaying KEGG libraries only\n")
for x in names:
    if x.startswith("KEGG_"):
        print (x)

print(f"\n   Displaying WikiPathway libraries only\n")
for x in names:
    if x.startswith("Wiki"):
        print (x)

print(f"\n\nAll libraries are saved in library_names_gseapy.csv file\n")
