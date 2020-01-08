#!/usr/bin/env python
#coding=utf-8
import numpy as np
import spglib
import re

def import_from_pwmat(filename):
    """
    Import structure data from a PWmat atom.config file.

    Args:
        filename (str): The file to import from.

    Returns:
        A Structure object.
    """
    with open(filename) as f:
        lines = f.readlines()

        number_of_atoms = int(lines[0].split()[0])

        lattice_vectors = []
        if lines[1].split()[0][0].lower() == 'l':
            for i in range(2, 5):
                lattice_vectors.append([float(x) for x in lines[i].split()[:3]])

        if lines[5].split()[0][0].lower() == 'p':
            line_index = 6

        positions = []
        atomic_numbers = []
        end_line_index = line_index + number_of_atoms
        for i in range(line_index, end_line_index):
            atomic_numbers.append(int(lines[i].split()[0]))
            positions.append([float(x) for x in lines[i].split()[1:4]])
        a = tuple(lattice_vectors[0])
        b = tuple(lattice_vectors[1])
        c = tuple(lattice_vectors[2])
        spg_lattice = [a,b,c]
        spg_positions = []
        for position in positions:
            spg_positions.append(tuple(position))
        spglib_type = (spg_lattice,spg_positions,atomic_numbers)
        return spglib_type, atomic_numbers
#spglib_type = import_from_pwmat("atom.config")
#spacegroup = spglib.get_spacegroup(spglib_type).encode('utf-8')
#print spacegroup
#pattern1 = re.compile(r'\([0-9]+\)')
#spacegroupnum = pattern1.findall(spacegroup)[0]
#pattern2 = re.compile(r'[0-9]+')
#spgnum = pattern2.findall(spacegroupnum)
#print int(spgnum[0])
if __name__=='__name__':
    import_from_pwmat()
