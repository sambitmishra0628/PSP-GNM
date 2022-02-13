#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 18:30:22 2018

@author: Sambit Mishra

"""
# Calculate the residue-residue contact based on an atomic structure
# Expects the structure to be processed i.e.,
# have no ligand atoms, have no water molecules and have addressed the ANISOU
# records and also the alternative locations A and B.
# In short, this script considers the provided PDB file to be pre-processed.
import numpy as np
import matplotlib.pyplot as plt
import re
import sys

# Python will only search the current directory for the list of user defined
# modules to import. We will need to explicitly include the path to the files
# we intend to import functions from.
sys.path.insert(0, '/home/8sm/Documents/Mitchell_Lab/Projects/PyLib')
import pdb_utils as ptils

### All subroutines are defined here ###

# validate whether a given pdb file is processed and has only coordinate
# information which is commenced by the keyword ATOM
def validatefile(filename):    
    rv = 0;
    fh = open(filename, 'r')
    for line in fh:
        if re.match('^ATOM',line):
            #pass
			continue
        else:
            rv = 1
            return rv
    fh.close()    
    return rv    


def calc_interresidue_shortest_dist(pdb_dict):
    all_coord = pdb_dict["COORD"]
    uniq_res_id = pdb_dict["UNIQ_RES_ID"]
    uniq_res = pdb_dict["UNIQ_RES"]
    all_res = pdb_dict["ALL_RES"]
    all_res_id = pdb_dict["ALL_RES_ID"]
    nres = len(uniq_res_id)
    atomic_dist_map = np.zeros((nres,nres)) # initialize the distance matrix
    
    # pick the shortest distance between any two heavy atoms of a residue pair
    for i in range(0,len(uniq_res_id)-1):
        for j in range(i+1,len(uniq_res_id)):
            res_i = uniq_res_id[i]
            res_j = uniq_res_id[j]
            #print("Res i = ",res_i,"Res j = ", res_j,"\n")
            shortest_dist = 100000 # arbitrarily assign a large valude
            
            # Find indices of all the atoms for res i and res j and find the
            # shortest distance between any two atom pairs from the two residues
            ind_i = np.where(all_res_id == res_i)
            ind_j = np.where(all_res_id == res_j)
            
            for k in np.nditer(ind_i):
                for l in np.nditer(ind_j):
                    #print("k=",k," l= ",l,"\n")
                    coord_ik = all_coord[k,:] # the x,y,z coordinate of the kth atom
                                              # of the ith residue
                    coord_jl = all_coord[l,:] # the x,y,z coordinate of the lth atom
                                              # of the jth residue                             
                    diff_ij = coord_ik-coord_jl;
                    r_ij = np.linalg.norm(diff_ij) # calculate the euclidean distance
                    #print ("i= ",i, " j= ",j, "k= ",k, "l= ",l, "rij = ",r_ij, "\n")
                    if r_ij < shortest_dist:
                        shortest_dist = r_ij                         
            atomic_dist_map[i,j] = shortest_dist
            atomic_dist_map[j,i] = shortest_dist
    return atomic_dist_map       


def calculate_contact_map(atomic_dist_map,dist_cutoff, offset):
    s = np.shape(atomic_dist_map)    
    nres = s[0]
    contact_map = np.zeros((nres,nres))
    for i in range(0,nres-1,1):
        for j in range(i+1,nres,1):
            if abs(i-j) > offset: # skip neighboring 4 residues
                rij = atomic_dist_map[i,j]
                if rij <= dist_cutoff:
                    contact_map[i,j] = 1
                    contact_map[j,i] = 1
                
    return contact_map

def all_atom_dist(pdb_dict):
    all_coord = pdb_dict["COORD"]
    uniq_res_id = pdb_dict["UNIQ_RES_ID"]
    uniq_res = pdb_dict["UNIQ_RES"]
    all_res = pdb_dict["ALL_RES"]
    all_res_id = pdb_dict["ALL_RES_ID"]
    natoms = len(all_res)    
    atomic_dist = np.zeros((natoms,natoms))
    for i in range(0,natoms-1,1):
        for j in range(i+1,natoms,1):
            coord_i = all_coord[i,:]
            coord_j = all_coord[j,:]
            diff_ij = coord_i-coord_j
            r_ij = np.linalg.norm(diff_ij)
            atomic_dist[i,j] = r_ij
            atomic_dist[j,i] = r_ij
    return atomic_dist



# The main function
filename = sys.argv[1]
rv = validatefile(filename)  # verify if the file fulfills the requirements
if rv: # if rv is true then the file is not processed
    sys.exit('Error! PDB file is not processed!')
else: 
    # parse the pdb file
    pdb_dict = ptils.parse_coordinates(filename)
#    atomic_dist = all_atom_dist(pdb_dict)
    
    # calculate the shortest distance between the heavy atoms of residues
    atomic_dist_map = calc_interresidue_shortest_dist(pdb_dict)
    dist_cutoff = 5 # Set the distance cutoff to decide on contacts
    offset = 4 # for a residue, we will not consider the contacts with 4 
                # neighboring residues
    contact_map = calculate_contact_map(atomic_dist_map,dist_cutoff, offset)
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(atomic_dist_map, cmap='jet')
    plt.gca().invert_yaxis() # ensure that both axes are numbered in the same way
    plt.show() 
    plt.subplot(1,2,2)
    plt.imshow(contact_map, cmap='gray_r')
    plt.gca().invert_yaxis() # ensure that both axes are numbered in the same way
    plt.show() 
# get the contact map
#[contact_map] = calc_contact_map()












