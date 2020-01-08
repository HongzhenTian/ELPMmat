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

# Triclinic systemc
if spgnum <= 2:
    enumber = 21
    csystem = 1
#Monoclnic system
elif spgnum <= 15:
    enumber = 13
    csystem = 2
#Orthorhombic system
elif spgnum <= 74:
    enumber = 9
    csystem = 3
#Tetragonal II system
elif spgnum <=88:
    enumber = 7
    csystem = 4
#Tetragonal I system
elif spgnum <=142:
    enumber = 6
    csystem = 5
#Trigonal II system
elif spgnum <= 148:
    enumber = 7
    csystem = 6
#Trigonal I system
elif spgnum <= 167:
    enumber = 6
    csystem = 7
#Hexagonal system
elif spgnum <= 194:
    enumber = 5
    csystem = 8
#Cubic system
else:
    enumber = 3
    csystem = 9



strains = [-0.018, -0.015, -0.012, -0.009, -0.006, -0.003, 0.000,\
            0.003, 0.006, 0.009, 0.012, 0.015, 0.018]
#for strain in np.arange(-0.018,0.018,0.,003):
#    strains.append(strain)


#Triclinic system (21 independent elastic constants)
# c11 c12 c13 c14 c15 c16
# c12 c22 c23 c24 c25 c26
# c13 c23 c33 c34 c35 c36
# c14 c24 c34 c44 c45 c46
# c15 c25 c35 c45 c55 c56
# c16 c26 c36 c46 c56 c66
def Triclinic():
    if not os.path.exists("./triclinic"):
        os.system("mkdir ./triclinic")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, 0.0, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, strain, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [strain, 0.0, strain, 0.0, 0.0, 0.0],\
                         [strain, 0.0, 0.0, strain, 0.0, 0.0],\
                         [strain, 0.0, 0.0, 0.0, strain, 0.0],\
                         [strain, 0.0, 0.0, 0.0, 0.0, strain],\
                         [0.0, strain, strain, 0.0, 0.0, 0.0],\
                         [0.0, strain, 0.0, strain, 0.0, 0.0],\
                         [0.0, strain, 0.0, 0.0, strain, 0.0],\
                         [0.0, strain, 0.0, 0.0, 0.0, strain],\
                         [0.0, 0.0, strain, strain, 0.0, 0.0],\
                         [0.0, 0.0, strain, 0.0, strain, 0.0],\
                         [0.0, 0.0, strain, 0.0, 0.0, strain],\
                         [0.0, 0.0, 0.0, strain, strain, 0.0],\
                         [0.0, 0.0, 0.0, strain, 0.0, strain],\
                         [0.0, 0.0, 0.0, 0.0, strain, strain]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config triclinic/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))


#Monoclinic system (13 independent elastic constants)
# c11 c12 c13 0 c15 0
# c12 c22 c23 0 c25 0
# c13 c23 c33 0 c35 0
# 0   0   0 c44 0 c46
# c15 c25 c35 0 c55 0
# 0 0 0 c46 0 c66
def Monoclinic():
    if not os.path.exists("./monoclinic"):
        os.system("mkdir ./monoclinic")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, 0.0, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, strain, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, strain, strain, 0.0, 0.0, 0.0],\
                         [strain, 0.0, strain, 0.0, 0.0, 0.0],\
                         [strain, 0.0, 0.0, 0.0, strain, 0.0],\
                         [0.0, strain, 0.0, 0.0, strain, 0.0],\
                         [0.0, 0.0, strain, 0.0, strain, 0.0],\
                         [0.0, 0.0, 0.0, strain, 0.0, strain]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config monoclinic/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))


#Orthorhombic system (9 independent elastic constants)
# c11 c12 c13 0 0 0
# c12 c22 c23 0 0 0
# c13 c23 c33 0 0 0
# 0   0   0 c44 0 0
# 0   0   0 0 c55 0
# 0   0   0 0 0 c66
def Orthorhombic():
    if not os.path.exists("./orthorhombic"):
        os.system("mkdir ./orthorhombic")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, 0.0, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, strain, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, strain, strain, 0.0, 0.0, 0.0],\
                         [strain, 0.0, strain, 0.0, 0.0, 0.0]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config orthorhombic/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))



#Tetragonal II system (7 independent elastic constants)
# c11 c12 c13 0 0 c16
# c12 c11 c13 0 0 -c16
# c13 c13 c33 0 0 0
# 0   0   0 c44 0 0
# 0   0   0 0 c44 0
# c16   -c16   0 0 0 c66
def TetragonalII():
    if not os.path.exists("./tetragonalII"):
        os.system("mkdir ./tetragonalII")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, strain, 0.0],\
                         [strain, strain, strain, 0.0, 0.0, 0.0],\
                         [0.0, strain, strain, 0.0, 0.0, 0.0],\
                         [strain, 0.0, 0.0, 0.0, 0.0, strain]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config tetragonalII/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))




#Tetragonal I system (6 independent elastic constants)
# c11 c12 c13 0 0 0
# c12 c11 c13 0 0 0
# c13 c13 c33 0 0 0
# 0   0   0 c44 0 0
# 0   0   0 0 c44 0
# 0   0   0 0 0 c66
def TetragonalI():
    if not os.path.exists("./tetragonalI"):
        os.system("mkdir ./tetragonalI")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, strain, 0.0],\
                         [strain, strain, strain, 0.0, 0.0, 0.0],\
                         [0.0, strain, strain, 0.0, 0.0, 0.0]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config tetragonalI/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))




#Rhombohedral II system (7 independent elastic constants)
# c11 c12 c13 c14 c15 0
# c12 c11 c13 -c14 -c15 0
# c13 c13 c33 0 0 0
# c14 -c14   0 c44 0 -c15
# c15 -c15   0 0 c44 c14
# 0   0   0 -c15 c14 c66=(c11-c12)/2
def RhombohedralII():
    if not os.path.exists("./rhombohedralII"):
        os.system("mkdir ./rhombohedralII")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, strain, 0.0],\
                         [strain, strain, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, strain, strain],\
                         [0.0, strain, 0.0, 0.0, 0.0, strain]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config rhombohedralII/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))




#Rhombohedral I system (6 independent elastic constants)
# c11 c12 c13 c14 0 0
# c12 c11 c13 -c14 0 0
# c13 c13 c33 0 0 0
# c14 -c14   0 c44 0 0
# 0 0   0 0 c44 c14
# 0   0   0 0 c14 c66=(c11-c12)/2
def RhombohedralI():
    if not os.path.exists("./rhombohedralI"):
        os.system("mkdir ./rhombohedralI")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, strain, 0.0],\
                         [strain, strain, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, strain, strain]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config rhombohedralI/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))





#Hexagonal system (5 independent elastic constants)
# c11 c12 c13 0 0 0
# c12 c11 c13 0 0 0
# c13 c13 c33 0 0 0
# 0   0   0 c44 0 0
# 0   0   0 0 c44 0
# 0   0   0 0 0 c66=(c11-c12)/2
def Hexagonal():
    if not os.path.exists("./hexagonal"):
        os.system("mkdir ./hexagonal")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, 0.0, 0.0, strain],\
                         [0.0, 0.0, strain, 0.0, 0.0, 0.0],\
                         [0.0, 0.0, 0.0, strain, strain, 0.0],\
                         [strain, strain, strain, 0.0, 0.0, 0.0]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config hexagonal/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))




#Cubic system (3 independent elastic constants)
# c11 c12 c12 0 0 0
# c12 c11 c12 0 0 0
# c12 c12 c11 0 0 0
# 0   0   0 c44 0 0
# 0   0   0 0 c44 0
# 0   0   0 0 0 c44
def Cubic():
    if not os.path.exists("./cubic"):
        os.system("mkdir ./cubic")
    atomic_symbol = store.AtomicSymbol()
    store.AtomToOLDPOS()
    for strain in strains:
        strain_matrix = [[0.0, 0.0, 0.0, strain, strain, strain],\
                         [strain, strain, 0.0, 0.0, 0.0, 0.0],\
                         [strain, strain, strain, 0.0, 0.0, 0.0]]
        for i in range(enumber):
            e = str(strain_matrix[i]).strip("[").strip("]")
            os.system("echo " + e + " | defvect.x >> /dev/null")
            os.system("cp fort.3 init" + "_" + str(strain) + "_" + str(i))
            os.system('''sed -i '6i ''' + atomic_symbol + '''' init''' + "_" + str(strain) + "_" + str(i))
            os.system("poscar2config.x < init" + "_" + str(strain) + "_" + str(i) + ">> /dev/null")
            os.system("mv atom.config cubic/atom.config" + "_" + str(strain) + "_" + str(i))
            os.system("rm init" + "_" + str(strain) + "_" + str(i))


if csystem == 1:
    Triclinic()
elif csystem == 2:
    Monoclinic()
elif csystem == 3:
    Orthorhombic()
elif csystem == 4:
    TetragonalII()
elif csystem == 5:
    TetragonalI()
elif csystem == 6:
    RhombohedralII()
elif csystem == 7:
    RhombohedralI()
elif csystem == 8:
    Hexagonal()
elif csystem == 9:
    Cubic()
else:
    print "unrecongnize space group"

os.system("rm -rf CONTCAR")
os.system("rm -rf fort.3")
os.system("rm -rf OLDPOS")
