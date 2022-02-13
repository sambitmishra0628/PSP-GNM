# Ensemble of functions for parsing PDB files

import numpy as np
import re
import sys

def parse_coordinates(filename):
	# parse coordinate information from a PDB file and return information
	# as a dictionary
	coords = [] # multi dimensional list to hold the coordinate information for each residue
	all_residues = [] # store the residue names
	all_res_ids = [] # store the residue ids
	all_res_chain_ids = []
	unique_residues = [] # store the unique residue names
	unique_residue_ids = [] # store the unique residue ids
	unique_res_chain_ids = [] # store the chain-residue identifier (e.g., A12,A13, B13, B15 and so on where A and B are the chain IDs)
	parsed_info = {} # dictionary to store the parsed information
	fh = open (filename, 'r')
	for line in fh:
		if re.match('^ATOM',line):
			line = line.rstrip()
			chain = line[21]
			X_coord = line[30:38]
			Y_coord = line[38:46]
			Z_coord = line[46:54]
			
			# Strip white space
			X_coord = re.sub('\s+','',X_coord)
			Y_coord = re.sub('\s+','',Y_coord)
			Z_coord = re.sub('\s+','',Z_coord)
			
			# convert from string to numeric
			X_coord = float(X_coord)
			Y_coord = float(Y_coord) 
			Z_coord = float(Z_coord) 
			coords.append([X_coord, Y_coord, Z_coord])						
			res_name = line[17:20]
			res_id = line[22:26]
			
			res_name = re.sub('\s+','',res_name)
			res_id = re.sub('\s+','',res_id)
			res_chain_id = chain+res_id

			all_res_chain_ids.append(res_chain_id)
			res_id = int(res_id)
			
			# include all residue names and ids
			all_residues.append(res_name) 
			all_res_ids.append(res_id)
			
			# check for unique residue ids i.e., remove redundancy occuring from side chain
			if len(unique_residues) == 0: # if empty then include current residue name and id
				unique_residues.append(res_name)
				unique_residue_ids.append(res_id)
				unique_res_chain_ids.append(res_chain_id)
			else: # Include the residue id, name and chain-residue identfier if previous id in the list 
					# is different from the current one. 
				prev_uniq_res_id = unique_residue_ids[-1]
				if res_id != prev_uniq_res_id:
					unique_residues.append(res_name)
					unique_residue_ids.append(res_id)
					unique_res_chain_ids.append(res_chain_id)
		else:
			#pass
			continue
	coord_np = np.array(coords)
	
	# Store the parsed information in a dictionary and return to caller
	parsed_info["COORD"] = coord_np
	parsed_info["UNIQ_RES"] = np.array(unique_residues)
	parsed_info["UNIQ_RES_ID"] = np.array(unique_residue_ids)
	parsed_info["ALL_RES"] = np.array(all_residues)
	parsed_info["ALL_RES_ID"] = np.array(all_res_ids)
	parsed_info["UNIQ_RES_CHAIN_ID"] = np.array(unique_res_chain_ids)
	parsed_info["ALL_RES_CHAIN_ID"] = np.array(all_res_chain_ids)
	return parsed_info


# parse the PDB file as a dictionary 
def parse_pdb_as_dict(filename):
	# parse the contents of a PDB file as a dictionary. This is especially useful
	# if the PDB file has more than one chain.
	PDB_struct = {}
	fh = open(filename, 'r')
	for line in fh:
		if re.match('^ATOM', line):
			chain = line[21]
			X_coord = line[30:38]
			Y_coord = line[38:46]
			Z_coord = line[46:54]
			
			# Strip white space
			X_coord = re.sub('\s+','',X_coord)
			Y_coord = re.sub('\s+','',Y_coord)
			Z_coord = re.sub('\s+','',Z_coord)
			
			# convert from string to numeric
			X_coord = float(X_coord)
			Y_coord = float(Y_coord) 
			Z_coord = float(Z_coord)	
			
			res_name = line[17:20]
			res_id = line[22:26]
			
			res_name = re.sub('\s+','',res_name)
			res_id = re.sub('\s+','',res_id)
			res_id = int(res_id)
			atomname = line[13:15]
			if chain not in PDB_struct.keys():
				PDB_struct[chain] = dict()
				PDB_struct[chain]['COORD'] = []
				PDB_struct[chain]['RESNAME'] = []
				PDB_struct[chain]['RESID'] = []
				PDB_struct[chain]['ATOM_NAME'] = []
			PDB_struct[chain]['COORD'].append([X_coord, Y_coord, Z_coord])
			PDB_struct[chain]['RESNAME'].append(res_name)
			PDB_struct[chain]['RESID'].append(res_id)
			PDB_struct[chain]['ATOM_NAME'].append(atomname)
	fh.close()
	return PDB_struct

# For a PDB file having 2 chains, write each chain into a separate file.
def split_pdb_2_chains(inputfile, chain1, chain2, outfile1, outfile2):
	fh = open(inputfile, 'r')
	fh1 = open(outfile1, 'w')
	fh2 = open(outfile2, 'w')
	for line in fh:
		if line[21] == chain1:
			fh1.write(line)
		elif line[21] == chain2:
			fh2.write(line)
	fh.close()
	fh1.close()
	fh2.close()

# def parse_calpha_coordinates(pdbfile):
# 	"""
# 	Parse the X,Y and Z coordinates for the C-alpha
# 	atoms from the given pdb file
# 	"""
# 	coords = [] # multi dimensional list to hold the coordinate information for each residue
# 	res_names = [] # List to hold the residue names
# 	res_ids = [] # List to hold the residue ids

# 	fh = open (filename, 'r')
# 	for line in fh:
# 		if re.match('^ATOM',line):
# 			atomname = line[13:15]
# 			atomname = atomname.replace(' ', '')
# 			if atomname != 'CA':
# 				continue
# 			line = line.rstrip()
# 			X_coord = line[30:38]
# 			Y_coord = line[38:46]
# 			Z_coord = line[46:54]
			
# 			# Strip white space
# 			X_coord = re.sub('\s+','',X_coord)
# 			Y_coord = re.sub('\s+','',Y_coord)
# 			Z_coord = re.sub('\s+','',Z_coord)
			
# 			# convert from string to numeric
# 			X_coord = float(X_coord)
# 			Y_coord = float(Y_coord) 
# 			Z_coord = float(Z_coord) 
# 			coords.append([X_coord, Y_coord, Z_coord])						
# 			res_name = line[17:20]
# 			res_id = line[22:26]
			
# 			res_name = re.sub('\s+','',res_name)
# 			res_id = re.sub('\s+','',res_id)
# 			res_chain_id = chain+res_id

# 			all_res_chain_ids.append(res_chain_id)
# 			res_id = int(res_id)
			
# 			# include all residue names and ids
# 			all_residues.append(res_name) 
# 			all_res_ids.append(res_id)
			
# 			# check for unique residue ids i.e., remove redundancy occuring from side chain
# 			if len(unique_residues) == 0: # if empty then include current residue name and id
# 				unique_residues.append(res_name)
# 				unique_residue_ids.append(res_id)
# 				unique_res_chain_ids.append(res_chain_id)
# 			else: # Include the residue id, name and chain-residue identfier if previous id in the list 
# 					# is different from the current one. 
# 				prev_uniq_res_id = unique_residue_ids[-1]
# 				if res_id != prev_uniq_res_id:
# 					unique_residues.append(res_name)
# 					unique_residue_ids.append(res_id)
# 					unique_res_chain_ids.append(res_chain_id)
# 		else:
# 			#pass
# 			continue
# 	coord_np = np.array(coords)
	
# 	# Store the parsed information in a dictionary and return to caller
# 	parsed_info["COORD"] = coord_np
# 	parsed_info["UNIQ_RES"] = np.array(unique_residues)
# 	parsed_info["UNIQ_RES_ID"] = np.array(unique_residue_ids)
# 	parsed_info["ALL_RES"] = np.array(all_residues)
# 	parsed_info["ALL_RES_ID"] = np.array(all_res_ids)
# 	parsed_info["UNIQ_RES_CHAIN_ID"] = np.array(unique_res_chain_ids)
# 	parsed_info["ALL_RES_CHAIN_ID"] = np.array(all_res_chain_ids)
# 	return parsed_info
