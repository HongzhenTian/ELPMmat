#!/usr/bin/env python
#coding=utf-8
import spglib
import ReadStructure
import re
import numpy as np
import os
import store


###############get space group number######################
config, atomic_numbers = ReadStructure.import_from_pwmat("strain.config")
lattice, positions, numbers = spglib.standardize_cell(config,symprec=1e-1)
lattice_new = []
for i in lattice:
    lattice_new.append(tuple(i))
positions_new = []
for position in positions:
    positions_new.append(tuple(position))
config = (lattice_new, positions_new, numbers)
spacegroup = spglib.get_spacegroup(config).encode('utf-8')
pattern1 = re.compile(r'\([0-9]+\)')
spacegroupnum = pattern1.findall(spacegroup)[0]
pattern2 = re.compile(r'[0-9]+')
spgnum = int(pattern2.findall(spacegroupnum)[0])
#print int(spgnum[0])
###########################################################

csystem = 0
enumber = 0

# Oblique systemc
# 2D space group: p1, p2
#(6 independent elastic constants)
if spgnum >= 3 and spgnum <=15:
    enumber = 6
    csystem = 1
#Rectanglar system
#2D space group: pm, pg, p2mm, p2gg, p2mg, cm, c2mm
#(4 independent elastic constants)
elif spgnum <= 74:
    enumber = 4
    csystem = 2
#Square system
#2D space group: p4, p4mm, p4gm
#(3 independent elastic constants)
elif spgnum <= 142:
    enumber = 3
    csystem = 3
#Hexagon system
#2D space group: p3, p6, p3m1, p31m, p6m
#(2 independent elastic constants)
elif spgnum <=194:
    enumber = 2
    csystem = 4
else:
    print "space group error!"


strains = [-0.018, -0.015, -0.012, -0.009, -0.006, -0.003, 0.000,\
            0.003, 0.006, 0.009, 0.012, 0.015, 0.018]
#for strain in np.arange(-0.018,0.018,0.,003):
#    strains.append(strain)


#Oblique system (6 independent elastic constants)
# c11 c12  c16
# c12 c22  c26
# c16 c26  c66
def oblique():
    if not os.path.exists("./oblique"):
        os.system("mkdir ./oblique")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, 0.0, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [strain, 0.0, 0.0, 0.0, 0.0, strain],\
                         [0.0, strain, 0.0, 0.0, 0.0, strain]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config oblique/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))


#Rectanglar system (4 independent elastic constants)
# c11 c12  0
# c12 c22  0
# 0 0 0  c66
def rectanglar():
    if not os.path.exists("./rectanglar"):
        os.system("mkdir ./rectanglar")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, 0.0, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [strain, strain, 0.0, 0.0, 0.0, 0.0]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config rectanglar/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))


#square system (3 independent elastic constants)
# c11 c12 0
# c12 c11 0
# 0   0 c66
def square():
    if not os.path.exists("./square"):
        os.system("mkdir ./square")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [strain, 0.0, 0.0, 0.0, 0.0, 0.0]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config square/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))



#Hexagon system (2 independent elastic constants)
# c11 c12 0
# c12 c11 0
# 0 0 0 c66=(c11-c12)/2
def hexagon():
    if not os.path.exists("./hexagon"):
        os.system("mkdir ./hexagon")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config hexagon/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))






if csystem == 1:
    oblique()
elif csystem == 2:
    rectanglar()
elif csystem == 3:
    square()
elif csystem == 4:
    hexagon()
else:
    print "unrecongnize space group"

os.system("rm -rf CONTCAR")
os.system("rm -rf fort.3")
os.system("rm -rf OLDPOS")

