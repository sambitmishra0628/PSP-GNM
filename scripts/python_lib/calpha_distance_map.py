#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 16:33:02 2018

@author: Sambit Mishra
"""

# Calculate the distance map between the C-alpha atoms in a protein. The input
# file is required to be a C_alpha coordinate file

import sys
import re
import numpy as np
import matplotlib.pyplot as plt


def get_ca_coordinates (filename):
    # parse the c-alpha coordinates from the PDB records
    # pdb_records is a list of lines, each line corresponding to a line entry
    # in a pdb file
    fh = open(filename, 'r')
    all_coords = []; # create a multi-dimensional array to store the coordinates
    for line_i in fh:
        if re.match('^\s*?$', line_i):
            pass
        elif re.match('^ATOM', line_i):
            line_i = line_i.rstrip()
            coords_i = line_i[26:54]
            coords_i = coords_i.split() # split by white space into individual elements
            # convert into integers
            coords_i = list(map(float,coords_i)) # convert from string to numeric
            all_coords.append(coords_i)
    fh.close()    
    # convert the multi-dimensional array into numpy array
    all_coords_ca = np.array(all_coords)
    return all_coords_ca
        

def calculate_ca_dist(ca_coords):
    # calculate c-alpha distances
    nres = len(ca_coords)
    dist_mat = np.zeros((nres,nres), dtype=float) # declare a 0 x 0 numpy matrix 
                                                  # to store the values
    for i in range(0,nres-1):
        for j in range(i+1,nres):
            diff_ij = ca_coords[i,:]-ca_coords[j,:];
            r_ij = np.linalg.norm(diff_ij)
            dist_mat[i,j] = r_ij
            dist_mat[j,i] = r_ij
    return dist_mat        



# The main script which will invoke the functions
filename = sys.argv[1]
all_coords_ca = get_ca_coordinates(filename)
dist_mat = calculate_ca_dist(all_coords_ca)
plt.figure()
plt.imshow(dist_mat, cmap='jet')
plt.show()     
