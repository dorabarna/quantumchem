# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 22:56:46 2017

@author: user
"""

import math
import numpy as np

## Reading in the atoms and their coordinates from a .xyz file
## The input file contains the number of atoms, their respective atomic coordinates 
## and their 3D Cartesian coordinates.

class Molecule(object):
    def __init__(self, n_atoms, atoms, coordinates):
        self.n_atoms = n_atoms
        self.atoms = atoms
        self.coordinates = coordinates
        
    def bond_length(self):
        """
        Calculates the bond length for each bond from the atomic Cartesian coordinates.
        """
        R = np.zeros((self.n_atoms, self.n_atoms))
        for i in range(self.n_atoms):
            for j in range(i):
                x_i = self.coordinates[i][0]
                y_i = self.coordinates[i][1]
                z_i = self.coordinates[i][2]
                x_j = self.coordinates[j][0]
                y_j = self.coordinates[j][1]
                z_j = self.coordinates[j][2]
                R[i, j] = math.sqrt((x_i - x_j)**2 + (y_i - y_j)**2 + (z_i - z_j)**2)
        return R
    
filename = input("Enter the filename: ")

file = open(filename, "r")

data = file.readlines()

n_atoms = int(data[0])

atoms = []
coordinates = []

for line in data[1:]:
    new_atom = line.split()
    atomic_number = int(new_atom[0])
    x = float(new_atom[1])
    y = float(new_atom[2])
    z = float(new_atom[3])
    atoms.append(atomic_number)
    coordinates.append((x, y, z))

new_molecule = Molecule(n_atoms, atoms, coordinates)
# R = new_molecule.bond_length()
#print(R[:,1])