"""
Protein stability change upon point mutations using a weighted Gaussian network model

Key steps:
    1. Coarse-grain a protein structure using only the alpha carbons
    2. Model the coarse-grained protein structure as a weighted 
Gaussian Network Model (GNM)
    3. The interactions between amino acid residues in such a model will
be weighted using statistical potentials
    4. Simulate unfolding of the wild-type structure by identifying residues
with most internal distance changes as the ones to break contact first
    5. Calculate the difference in energies obtained from the Miyazawa-Jernigan potential
    and entropies (given by the mean-squared fluctuations in residue distance) for the
    mutant and wild-type structures
"""

import numpy as np
from scipy import linalg
import os
import re
import sys
import pandas as pd
import click
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from sklearn.linear_model import LinearRegression

FILE_MAP = {'MJ': 'MIYS960101.txt'}
POT_MAT_DIR = 'potential_matrices' # Directory containing the potential matrix


def map_resname_to_id(res_code):
    """
    Convert the 3-lettered residue code to single letter
    
    Parameters:
        res_code: The three-lettered amino acid code
    Returns:
        The corresponding single-letter amino acid code
    """
    resname_2_id = {'ALA' : 'A', 'ARG' : 'R', 'ASN' : 'N', 'ASP' : 'D',
                    'CYS' : 'C', 'GLY' : 'G', 'GLN' : 'Q', 'GLU' : 'E',
                    'HIS' : 'H', 'ILE' : 'I', 'LEU' : 'L', 'LYS' : 'K',
                    'MET' : 'M', 'PRO' : 'P', 'PHE' : 'F', 'SER' : 'S',
                    'THR' : 'T', 'TRP' : 'W', 'TYR' : 'Y', 'VAL' : 'V',}
    return resname_2_id[res_code]

def process_wt_pdb (input_dir, output_dir, pdb_chain_list):
    """
    Process the atomic PDB structures
    
    Parameters:
        input_dir: Directory containing the raw PDB files

        output_dir: Directory to which the processed files will be written into
        
        pdb_chain_list: List of 4 letter PDBID and Chain ID as a single string 
    """

    if not input_dir.endswith('/'):
        input_dir = input_dir + '/'
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        
    not_found_list = [] # PDB_CHAIN entries that do not have PDB files in the input_dir    
    pdbfiles = glob.glob(input_dir + '*.pdb') # List of pdb files in the input directory
    
    if len(pdbfiles) == 0:
        sys.exit(f"Error! No *.pdb files found in {input_dir}")
        
    # Check if the names of the PDB files in the input directory include any chain information.
    pdb_id_lens = [len(pdb_i.split('/')[-1].split('.pdb')[0]) for pdb_i in pdbfiles]
    
    # If all pdb files have ids that are of length 4, then no chain
    # ID is included in the file name
    if len(set(pdb_id_lens)) == 1 and list(set(pdb_id_lens))[0] == 4:
        chain_flag = False
    elif len(set(pdb_id_lens)) == 1 and list(set(pdb_id_lens))[0] == 5:
        print ("PDB file names have length 5! Assuming the 5th alphabet is the chain ID!")
        chain_flag = True
    else:
        print (f"Error! Inconsitent PDB filenames! PDB file names (.pdb extension excluded) in {input_dir} "
               "should either have only the PDB ID (e.g., 1cei) or PDB ID with chain (e.g., 1ceiA). Acceptable "
               "formats include 1cei.pdb, 1CEI.pdb, 1ceiA.pdb or 1CEIA.pdb. All files should have the same naming format.")
        sys.exit()
        
    
    # We can expect more pdb files in the output directory compared 
    # to the input directory - separate pdb files for each chain
    if len(glob.glob(output_dir + '*.pdb')) >= len(pdbfiles): 
        print (f"Processed pdb files already present at {output_dir}! Nothing to do!")
        return
    total_files = len(pdb_chain_list)
    count = 0
    # print (f"Chain flag = {chain_flag}")
    # print (f"{pdbfiles}")
    for pdb_i in pdb_chain_list:
        # The pdb ids present in the input dir may differ from
        # those in the pdb_chain_list by their case. Consider such differences.
        if input_dir + pdb_i[0:4].upper() + '.pdb' in pdbfiles and chain_flag == False:
            pdbfile = input_dir +  pdb_i[0:4].upper() +  '.pdb'
        elif input_dir + pdb_i[0:4].lower() + '.pdb' in pdbfiles and chain_flag == False:
            pdbfile = input_dir +  pdb_i[0:4].lower() +  '.pdb'
        elif input_dir + pdb_i[0:5].upper() + '.pdb' in pdbfiles and chain_flag == True:
            pdbfile = input_dir +  pdb_i[0:5].upper() +  '.pdb'
        elif input_dir + pdb_i[0:5].lower() + '.pdb' in pdbfiles and chain_flag == True:
            pdbfile = input_dir +  pdb_i[0:5].lower() + '.pdb' 
        else:
            not_found_list.append(pdb_i[0:4] +  '.pdb')
            continue

        pdb_name = pdb_i
        chain_id = pdb_i[4]
        outfile = output_dir + pdb_name + '.pdb'
        # print (f"pdb_name = {pdb_name}, outfile = {outfile}")
        # return
        fh1 = open(pdbfile, 'r')
        fh2 = open(outfile, 'w')
        # print("Processing file ", pdbfile)
        is_nmr = False
        for line in fh1:
            if line.startswith('EXPDTA') and 'NMR' in line:
                is_nmr = True
                continue
            elif line.startswith('MODEL'):
                # By default we will use the first model
                is_nmr = True
                model_num = int(line.split()[1])
                fh2.write(line)
                continue
            elif re.match('^ATOM', line) and line[21] == chain_id: # Only get coordinates for the given chain ID
                alt_loc = line[16] # alternate location if any
                alt_loc = alt_loc.replace(' ','')
                if alt_loc == 'A': # we will consider by default the 'A' location
                    fh2.write(line)
                elif alt_loc == '':
                    fh2.write(line)
                else:
                    continue						
            elif re.match('^TER', line) and line[21] == chain_id:
                fh2.write(line)
            elif line.startswith('ENDMDL'):
                fh2.write(line)
                break
        fh1.close()
        fh2.close()
        # print("Done!\n")
        count += 1
        print (f"Processed {count}/{total_files} pdb files", end='\r', flush=True)
    print ("\nDone!")
    
    if len(not_found_list) > 0:
        not_found_str = ','.join(not_found_list)
        sys.exit(f"Error! Wildtype PDB files could not be found for {not_found_str} in {input_dir}! Cannot continue with prediction!")

def parse_calpha_pdb (pdbfile):
    """
    Parse the coordinates of the c-alpha atoms from the given pdb file
   
    Parameters:
        pdb_file: The PDB file with path
    Returns:
        PDB_struct: A dictionary containing the parsed information
    """ 
    PDB_struct = {}
    fh = open(pdbfile, 'r')
    for line in fh:
        if re.match('^ATOM', line):
            atomname = line[13:15]
            atomname = atomname.replace(' ', '')
            # Only consider the information for the c-alpha atoms
            if atomname != 'CA':
                continue
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

def parse_energy_pot(filename):    
    """
    Parse the Miyazawa-Jernigan potential matrix file into a dictionary
   
    Parameters:
        filename : The name of the file that includes the MJ potential as given
        in the AA Index database.
    Returns:
        aa_ind_dict: A dictionary containing the parsed MJ potential
    """   
    with open(POT_MAT_DIR + '/' + filename, 'r') as fh:
        fc = fh.readlines()
        aa_ind_str = ''.join(fc)
    aa_ind_dict = {}
    aa_ind_str = aa_ind_str.strip()
    aa_ind_list = aa_ind_str.split('\n')
    aa_ind_id = aa_ind_list[0]
    aa_ind_id = str(aa_ind_id[2:])

    aa_ind_dict['AA_INDEX_ID'] = aa_ind_id
    aa_order = ''; mat_row_ind = 0;
    for line_i in aa_ind_list:
        # Scan for the line starting with M and capture the order of amino acids used
        # in the contact matrix
        rx = re.compile('^M.*rows = (.*?),')
        if rx.match(line_i):
            aa_order = rx.match(line_i).group(1)
            aa_order = aa_order.replace(' ', '')
            #print (f" aa_order len = {len(aa_order)}")
            continue
        elif aa_order: # We have the order of aa. Now parse contact potential values.
            row_elems = line_i.split() # Get the elements in a single row of the matrix
            #if len(row_elems) != len(aa_order):
            #    return {}
            if mat_row_ind > len(aa_order)-1:
                #print ("Condition satisfied")
                break
            for j in range(0,mat_row_ind+1):
                #print(mat_row_ind, len(aa_order))
                res1 = aa_order[mat_row_ind]
                res2 = aa_order[j]
                val_12 = row_elems[j]
                # check if value can be converted into float or not
                try:
                    val_12 = float(val_12)
                except ValueError:
                    continue
                aa_ind_dict[res1 + ':' + res2] = val_12
                aa_ind_dict[res2 + ':' + res1] = val_12
            mat_row_ind += 1
    # print (f"aa ind dict = {aa_ind_dict}")        
    return aa_ind_dict

def energy_weighted_kirchoff(coord, res_codes, pot_type, cutoff, contact_matrix=None):
    """ 
    Get the Kirchoff's matrix in which the contact springs are assigned based on
    the energetic interaction between residues. If contact matrix is provided, then use it.
    Otherwise, calculate the contact matrix based on the residue coordinates.
    Parameters:
        coord : The C-alpha coordinates
        res_codes: The amino acid sequence of the PDB (as obtained from the PDB file)
        pot_type: By default it's the MJ potential
        cutoff: The C-alpha distance cutoff for interacting residues
        contact_matrix: The C-alpha residue-residue contact matrix.
    Returns:
        K: The weighted Kirchhoff matrix
        energy_matrix: The matrix of interaction energies for the interacting residues
    """
    pot_file = FILE_MAP[pot_type]
    pot_dict = parse_energy_pot(pot_file)
    del(pot_dict['AA_INDEX_ID'])

    # print (pot_dict)
    r,c = np.shape(coord)
    # Convert the potential dictionary into a matrix of size r x r 
    # having weights of interactions.
    P = np.zeros((r,r))
    # Have another matrix contain the interaction energies
    E_matrix = np.zeros((r,r))
    for i in range(0,r):
        for j in range(0,r):
            res_i, res_j = res_codes[i], res_codes[j]
            # res_i, res_j = map_resname_to_id(res_names[i]), map_resname_to_id(res_names[j])
            # For neighboring residues, the weight will be more than
            # the maximum weight based on the potential used.
            if abs(i-j) == 1:
                # Make the interaction energies of the adjacent c_alpha atoms stronger
                pot_ij = min(pot_dict.values())-1 
            else:    
                pot_ij = pot_dict[res_i + ':' + res_j]
            # The potential value is in energy terms. Convert to weights.
            pot_ij = round(np.exp(-pot_ij), 2)
            P[i,j] = pot_ij
    if contact_matrix.size != 0:
        # Typically useful when we want to provide a customized contact_matrix
        C = contact_matrix
    else:
        C, D = get_ca_contacts(coord, cutoff, r)
    energy_matrix = np.multiply(C,P)
    C = -C

    # We will weigh the residue-residue contacts by their energy potentials.
    K = np.multiply(C,P)
    K_cpy = K.copy()
    for i in range(0,r):
        # Diagonals in Kirchoff matrix are sum of the rows/columns
        # except the diagonal. 
        K[i,i] = -(np.sum(C[i,:] * P[i,:])-K_cpy[i,i])
    return K, energy_matrix 

def get_ca_contacts(coord, cutoff, num_res):
    """
    Get the contact matrix for the C-alpha atoms 
    
    Parameters:
        coord: The coordinates of the C-alpha atoms as N-by-3 numpy array.
        N is the total number of residues.
        
        cutoff: Distance cutoff of C-alpha atoms that defines an interacting pair
        
        num_res: Total number of residues in the protein
    Returns:
        C: C-alpha contact matrix
        D: C-alpha distance matrix
    """
    C = np.zeros((num_res,num_res)) # Contact matrix
    D = np.ones((num_res, num_res)) # Distance matrix
    for rn in range(0,num_res-1):
        #print (rn)
        coord_rn = coord[rn]
        coord_rem = coord[rn+1:]
        r2,c2 = np.shape(coord_rem)
        dist = np.linalg.norm(np.tile(coord_rn,(r2, 1))-coord_rem, axis=1)
        for i,j,d in zip([rn]*r2, list(range(rn+1,num_res)), dist):
            d = round(d,2)
            if i != j:
                D[i,j] = d
                D[j,i] = d
            if d <= cutoff:
                C[i,j] = 1
                C[j,i] = 1
    return C, D

def calc_residue_cross_corr(V,E,n_modes=20):
    """
    Calculate the cross correlations between residues 
    
    Parameters:
        V: The matrix including the modes (eigen values)
        E: The eigen vectors associated with each mode
        n_modes: Total number of modes
    Returns:
        C: Cross correlation matrix
        bfact: The theoretical residue mean-squared fluctuations
    """    

    n_res = len(V[:,0])
    nrow,ncol = np.shape(V)
    if n_modes > ncol:
        n_modes = ncol
    Hinv = np.zeros((n_res, n_res))
    for m in list(range(0,n_modes)):
        Hinv += (1/E[m]) *np.outer(V[:,m], V[:,m])
    C = np.zeros((n_res, n_res))
    for i in range(0,n_res):
        for j in range(0,n_res):
            C[i,j] = Hinv[i,j]/(Hinv[i,i]*Hinv[j,j])**0.5
    bfact = np.diagonal(Hinv)
    return C, bfact    
    
def calc_internal_dist_change(C):
    """ 
    Calculate the internal distance change (aka mean-squared fluctuation in distance 
    between a residue pair) given the cross-correlation matrix
    
    Parameters:
        C: The matrix including correlations and cross-correlations between the C-alpha atoms
    Returns:
        I: Internal distance change matrix
    """
    r,c = np.shape(C)
    I = np.zeros((r,r))
    for i in range(0,r):
        for j in range(0,r):
            # Round the internal distance change to 3 decimal points
            I[i,j] = round((C[i,i] + C[j,j] - 2*C[i,j]), 3)
    return I

def calc_gnm(coord, cutoff=9, num_modes=10, spring_type=None, res_codes=None, contact_matrix=None):
    """ 
    Run calculations for GNM in which the interactions are weighted
    between the residues.
    
    Parameters:
        coord: The C-alpha coordinates
        
        cutoff: Cutoff distance between residue C-alpha atoms that defines
        the interactions between residues
        
        num_modes: The total number of low frequency modes with non-zero
        eigen values to be returned
        
        spring_type: Type of weighting for the interactions between residues.
        The interaction strength is obtained from the Miyazawa-Jernigan contact
        potential. Default spring type in this implementation is 'MJ'
        
        res_codes: The single letter amino acid sequence of the protein that
        corresponds to the PDB sequence.
        
        contact_matrix: The C-alpha contact matrix
    Returns:
        V: Matrix of low-frequency modes (eigen vectors). The total number
        of modes is equal to the num_modes
        E: Eigen values for the num_modes low-frequency modes
        e_matrix: The potential matrix defining the interaction strengths of the
        interacting residues    
    """
    if spring_type == None:
        sys.exit("Error! spring_type cannot be None for weighted GNM")
    if type(res_codes) != list:
        sys.exit("Error! res_codes must be a list of 3-letter residue codes of all residues")

    stat_pot = ['MJ']
    if spring_type in stat_pot:
        K, e_matrix = energy_weighted_kirchoff(coord, res_codes, spring_type, cutoff, contact_matrix)
    else:
        sys.exit(f"Invalid spring type {spring_type}")
    
    E,V = linalg.eigh(K, turbo=True) # Returns eigen values and vectors    
    return V[:,1:num_modes+1], E[1:num_modes+1], e_matrix

def sanity_check(data_file, pdb_dir):
    """
    Check the data_file to see if the residues in the WILD_RES column (forward mutants) or
    MUTANT_RES (reverse_mutants) column of the data_file are indeed present in the PDB file
    at the positions specified by RES_NUM_PDB.

    Parameters:
        data_file : Name of the data file with information on PDB_CHAIN, wt_residue, mut_residue, residue_position
        pdb_dir : The directory including the atomic processed PDB files

    Returns: Total number of records that show correct mapping
    """
    if not pdb_dir.endswith('/'):
        pdb_dir += '/'
    df = pd.read_csv(data_file, encoding='utf8')
    total_recs = len(df)
    correct_recs = 0
    serial_resnum_list = []
    for pdb_id_i, res_i_wt, res_i_mut, category_i, pos_i in zip(df['PDB_CHAIN'], df['WILD_RES'], df['MUTANT_RES'], df['Category'], df['RES_NUM_PDB']):
        if (not np.isnan(pos_i)):
            # Get the amino acid sequence from the PDB file
            pdbfile_i = pdb_dir + pdb_id_i + '.pdb'
            ca_dict = parse_calpha_pdb(pdbfile_i)
            # Get the target chain from the pdb id
            chain = pdb_id_i[-1]
            # Get the single letter residue names 
            seq_i = [map_resname_to_id(res_i)  for res_i in ca_dict[chain]['RESNAME']]
            # Get the PDB residue numbers for all residues
            pdb_res_nums = ca_dict[chain]['RESID']
            serial_res_num = pdb_res_nums.index(int(pos_i)) # Map the pos_i (the pdb residue number) to serial index
            serial_resnum_list.append(serial_res_num)
            if category_i == 'Forward' and seq_i[int(serial_res_num)] ==  res_i_wt:
                correct_recs += 1
            elif category_i == 'Reverse' and seq_i[int(serial_res_num)] ==  res_i_mut:
                correct_recs += 1
            else:
                print (f"Residue mismatch for {category_i} mutant: {pdb_id_i}, {pos_i}, wt_res = {res_i_mut}, mut_res = {res_i_wt}") 
        else:
            sys.exit(f"RES_NUM_PDB is not defined in the data file for {pdb_id_i} WT:{res_i_wt}, MUT:{res_i_mut}!")        
    print (f"{correct_recs}/{total_recs} records show correct mutant position and amino acid match")
    return correct_recs, serial_resnum_list

def simulate_unfolding(ca_coord, res_codes, pdb_id, dist_cutoff, num_modes, mut_or_wt='wt',serial_res_num=None):
    """
    Simulate unfolding based on change in internal distance (mean-squared fluctuation in distance)
    between residues. We will simulate the unfolding until 50 percent contacts in the starting structure
    are broken.
    
    Parameters:
        ca_coord: The C-alpha coordinates
        
        res_codes: The single letter amino acid sequence of the protein that
        corresponds to the PDB sequence.
        pdb_id: The four-letterd PDB ID plus the chain ID
        
        dist_cutoff: Cutoff distance between residue C-alpha atoms that defines
        the interactions between residues
        
        num_modes: The total number of low frequency modes with non-zero
        eigen values to be returned
        mut_or_wt: Whether the simulation is being done for the mutant or wildtype protein
        
        serial_res_num: The serial index of the mutant position (not the PDB position).
    Returns:
        df_contact_breaks: A dataframe including information on the contacts broken
    """
    df_contact_breaks = pd.DataFrame()
    break_threshold = 0.5

    #dist_cutoff = 9 # Cutoff distance for interacting residues
    L_adj = 3 # Strength of interactions between adjacent residues.
              # This value will be reset based on the type of potential used.
              # If this value is less than the maximum value in the calculated
              # contact potential, then this value will be reset in the energy_weighted_kirchoff
              # subroutine.
    int_potential = 'MJ' # Statistical potential for residue-residue interactions
    #num_modes = 10
    num_res = len(res_codes)
    C,D = get_ca_contacts(ca_coord, dist_cutoff, num_res)
    tot_contacts_folded = C.sum()/2
    # print (f"Total contacts in folded protein = {tot_contacts_folded}")
    if mut_or_wt != 'wt':
        res_at_mut_pos = res_codes[serial_res_num] # Get the residue at the mutant position
   
    # When simulating for mutant, then get contacts for mutant position
    if mut_or_wt != 'wt':
        mut_pos_contacts = list(np.where(C[serial_res_num,:] == 1)[0])
        mut_pos_tot_contacts = len(mut_pos_contacts) * 2 # Consider both (i,j) and (j,i) pairs
    
    # Parse the potential matrix
    pot_file = FILE_MAP[int_potential]
    pot_dict = parse_energy_pot(pot_file)
    del(pot_dict['AA_INDEX_ID'])    # Remove the identifer key "AA_INDEX_ID"

    mut_pos_num_broken_contacts = 0 # Number of contacts made by the mutant position that are broken
    total_contacts_broken = 0
    contact_matrix = C
    iteration = 1
    # Run simulation until 50% of all the contacts are broken
    while (total_contacts_broken <= round(tot_contacts_folded*break_threshold)):
        # print (f"Simulating unfolding iteration {iteration}...", end='\r', flush=True)
        # Run GNM. For the first iteration, we will not calculate the contact_matrix
        # For subsequent iterations, we will calculate the contact_matrix based on
        # the internal distance changes.
        #print ("Running calc_gnm...", end='', flush=True)
        V,E, e_matrix = calc_gnm(ca_coord, dist_cutoff, num_modes, int_potential, res_codes, contact_matrix)
        #print ("Done!")
        
        # If the second eigen value (calc_gnm excludes the first eigen value and 
        # vector) is less than 0.00001 (close to 0) then stop the simulation
        # as the model is no longer stable
        if E[0] < 0.00001:
            break
        # Calculate atomic cross correlations
        #print ("Running calc_residue_cross_corr...", end='', flush=True)
        corr_matrix, bfact = calc_residue_cross_corr(V,E,num_modes)
        #print ("Done!")
        #print ("Running calculations for internal distance change...", end='', flush=True)
        int_dist_matrix = calc_internal_dist_change(corr_matrix)
        #print ("Done!")
        
        #print ("Getting coordinates for true contacts broken...", end='', flush=True)    
        # Get the product of internal dist matrix and the contact matrix
        P_mat = contact_matrix * int_dist_matrix
        
        # Get the coordinates of the maximum value in the P_mat
        max_int_dist_val = np.max(P_mat)
        max_ind_arr = np.where(P_mat == max_int_dist_val) 
        num_nonzero_contacts = len(max_ind_arr[0])
        #print ("Done!")
        # Break the contacts between the residue pairs
        # having the selected high internal distance change value.
        # The actual number of broken contacts is half since the contact matrix is symmetric
        #print ("Running break contacts loop...", flush=True, end='')
        for row_index, col_index in zip(max_ind_arr[0], max_ind_arr[1]):
            if contact_matrix[row_index, col_index] == 1:
                contact_matrix[row_index,col_index] = 0
                # If the simulation is being done for a mutant position
                if mut_or_wt != 'wt':
                    # If the broken contact is between the residue at the mutant position
                    # and another residue, then note it.
                    if (row_index == serial_res_num and col_index in mut_pos_contacts) or (col_index == serial_res_num and row_index in mut_pos_contacts):
                        if row_index == serial_res_num:
                            pair_res_num = col_index
                        else:
                            pair_res_num = row_index    
                        mut_pos_num_broken_contacts += 1
                        # contact_break_str = 'Iteration-' + str(iteration) + ',(' + str(row_index) + ',' + str(col_index) + ')'
                        res_code_i = res_codes[serial_res_num]
                        res_code_j = res_codes[pair_res_num]
                        
                        column_names = ['PDB_ID', 'WT_or_Mut', 'Mutation_position', 'Contact_Position', 'Res_at_Mut_Position', 'Res_at_Contact_Pos', 'Energy_MJ', 'Int_dist_change', 'Contact_break_rank']
                        df_tmp = pd.DataFrame([[pdb_id, mut_or_wt, serial_res_num, pair_res_num, res_code_i,
                            res_code_j, pot_dict[res_code_i + ':' + res_code_j], max_int_dist_val, iteration]], columns=column_names)
                        df_contact_breaks = df_contact_breaks.append(df_tmp)
                else:
                    # If simulation is being done for a wildtype structure
                    res_code_i = res_codes[row_index]
                    res_code_j = res_codes[col_index]
                    column_names = ['PDB_ID', 'WT_or_Mut', 'Mutation_position', 'Contact_Position', 'Res_at_Mut_Position', 'Res_at_Contact_Pos', 'Energy_MJ', 'Int_dist_change', 'Contact_break_rank']
                    df_tmp = pd.DataFrame([[pdb_id, mut_or_wt, row_index, col_index, res_code_i,
                            res_code_j, pot_dict[res_code_i + ':' + res_code_j], max_int_dist_val, iteration]], columns=column_names)
                    df_contact_breaks = df_contact_breaks.append(df_tmp)             
        total_contacts_broken += num_nonzero_contacts/2
        #print("Done!")
        # figfile = pdb_id + '_' + mut_or_wt + '_cont_mat_iter_' + str(iteration) + '.png'    
        # create_contact_map_fig(contact_matrix,'',figfile, res_codes)
        iteration += 1
        #print (f"Total contacts broken = {total_contacts_broken}/{tot_contacts_folded}")
    return df_contact_breaks

def calc_mut_energy_folded_unfolded(processed_pdb_dir, pdb_i, wt_res, mut_res, res_num, dist_cutoff, num_modes):
    """ 
    Calculate the interaction energies in the folded and unfolded mutant structure.
    
    Parameters:
        processed_pdb_dir: The directory having the processed PDB files
        
        pdb_i: The four-letterd PDB ID plus the chain ID
        wt_res: The amino acid (single letter) at the mutation position in the wildtype 
        structure
        mut_res: The amino acid (single letter) at the mutation position in the mutant
        structure
        res_num: The mutation postion number as in the PDB file
        dist_cutoff: Cutoff distance between residue C-alpha atoms that defines
        the interactions between residues
        
        num_modes: The total number of low frequency modes with non-zero
        eigen values to be used for calculations
    Returns:
        df_contact_breaks: A dataframe including information on the contacts broken
    """
    pdb_id = pdb_i
    pdb_id.replace('.pdb', '') # Just the pdb id for naming purpose
    if not pdb_i.endswith('.pdb'):
        pdb_i += '.pdb'
    pdbfile = processed_pdb_dir + pdb_i
    
    # Parse the c-alpha atom information
    ca_dict = parse_calpha_pdb(pdbfile)
    # print (f"c-alpha dictionary = {ca_dict}")
    
    # Get the target chain from the pdb id
    chain = pdb_i.split('.pdb')[0][-1]
    # print (ca_dict)
    # Get the single letter residue names 
    res_codes = [map_resname_to_id(res_i)  for res_i in ca_dict[chain]['RESNAME']]
    # print (res_codes)
    
    # Get the residue numbers from the pdb file
    pdb_res_nums = ca_dict[chain]['RESID']
    serial_res_num = pdb_res_nums.index(int(res_num)) # Convert the res_num (the pdb residue number) to serial number

    # Change the residue code at the serial_res_num to mut_res
    res_codes[serial_res_num] = mut_res

    # Simulate unfolding and calculate the energy associated with unfolding
    mut_tag = wt_res + str(res_num) + mut_res
    df_contact_breaks = simulate_unfolding(ca_dict[chain]['COORD'], res_codes, pdb_id, dist_cutoff, num_modes, mut_tag, serial_res_num)
    return df_contact_breaks

def calc_wt_energy_folded_unfolded(processed_pdb_dir, pdb_i, dist_cutoff, num_modes):
    """ 
    Calculate the interaction energies in the folded and unfolded wildtype structure.
    
    Parameters:
        processed_pdb_dir: The directory having the processed PDB files
        
        pdb_i: The four-letterd PDB ID plus the chain ID
        dist_cutoff: Cutoff distance between residue C-alpha atoms that defines
        the interactions between residues
        
        num_modes: The total number of low frequency modes with non-zero
        eigen values to be used for calculations
    Returns:
        df_contact_breaks: A dataframe including information on the contacts broken
    
    """
    pdb_id = pdb_i
    pdb_id = pdb_id.replace('.pdb', '') # Just the pdb id for naming purpose
    if not pdb_i.endswith('.pdb'):
        pdb_i += '.pdb'
    pdbfile = processed_pdb_dir + pdb_i
    
    # Parse the c-alpha atom information
    ca_dict = parse_calpha_pdb(pdbfile)
    # print (f"c-alpha dictionary = {ca_dict}")
    
    # Get the target chain from the pdb id
    chain = pdb_i.split('.pdb')[0][-1]
    
    # Get the single letter residue names 
    res_codes = [map_resname_to_id(res_i)  for res_i in ca_dict[chain]['RESNAME']]
    
    # Simulate unfolding and calculate the energy associated with unfolding
    df_contact_breaks = simulate_unfolding(ca_dict[chain]['COORD'], res_codes, pdb_id,  dist_cutoff, num_modes, 'wt', '')
    return df_contact_breaks

def run_ab_initio_stability_prediction_wildtype(pdb_i, outdir, processed_pdb_dir, dist_cutoff, num_modes):
    """ 
    Wrapper function for running calculations on wildtype proteins and write the 
    contact break information into a .csv file.
    
    Parameters:
        pdb_i: The four lettered PDB ID plus the chain
        outdir: The directory to which the files will be written to
        processed_pdb_dir: The directory having the processed PDB files
        
        dist_cutoff: Cutoff distance between residue C-alpha atoms that defines
        the interactions between residues
        
        num_modes: The total number of low frequency modes with non-zero
        eigen values to be used for calculations
    """
    outfile = pdb_i + '_wt_contact_breaks.csv'
    
    # Do not run calculations if outputfile already exists
    if os.path.isfile(outdir + outfile):
        print (f"Skipping unfolding simulation for wildtype {pdb_i} as contact-break file is already present!", flush=True)
        return
    print (f"Running calculations for {pdb_i} wild type, dist_cutoff = {dist_cutoff}, num_modes = {num_modes}...", flush=True)
    df_contact_breaks_wt = calc_wt_energy_folded_unfolded(processed_pdb_dir, pdb_i, dist_cutoff, num_modes)
    
    # Write the calculation output to a file
    df_contact_breaks_wt.to_csv(outdir + outfile, index=False)

def run_ab_initio_stability_prediction_mutant(row_i, outdir, processed_pdb_dir, dist_cutoff, num_modes):
    """ 
    Wrapper function for running calculations on mutant
    
    Parameters:
        row_i: The ith row as given in the benchmark .csv datafile
        outdir: The directory to which the files will be written to
        processed_pdb_dir: The directory having the processed PDB files
        
        dist_cutoff: Cutoff distance between residue C-alpha atoms that defines
        the interactions between residues
        
        num_modes: The total number of low frequency modes with non-zero
        eigen values to be used for calculations
    """
    pdb_i, wt_res, mut_res, res_num_pdb, mut_category = row_i['PDB_CHAIN'], row_i['WILD_RES'], row_i['MUTANT_RES'], row_i['RES_NUM_PDB'], row_i['Category']

    # Skip calculations if output file is already present
    outfile = pdb_i + '_' + row_i['WILD_RES'] + str(row_i['RES_NUM_PDB']) + row_i['MUTANT_RES'] + '_contact_breaks.csv'
    if mut_category == 'Forward' and os.path.isfile(outdir + outfile):
        print (f"Skipping unfolding simulation for mutant {pdb_i}: {row_i['WILD_RES']}{row_i['RES_NUM_PDB']}{row_i['MUTANT_RES']} as contact-break file is already present!")
        return
    if mut_category == 'Forward' and not os.path.isfile(outdir + outfile):
        print (f"Running calculations for {pdb_i}, mutation {wt_res}{res_num_pdb}{mut_res} dist_cutoff = {dist_cutoff}, num_modes = {num_modes}...",flush=True)
        df_contact_breaks_mut = calc_mut_energy_folded_unfolded(processed_pdb_dir, pdb_i, wt_res, mut_res, res_num_pdb, dist_cutoff, num_modes)
    
    # If this mutant is a reverse mutant, we also need to generate the contact breaks file for the forward mutant
    # if that is not already present
    outfile2 = pdb_i + '_' + row_i['MUTANT_RES'] + str(row_i['RES_NUM_PDB']) + row_i['WILD_RES'] + '_contact_breaks.csv'
    if mut_category == 'Reverse' and not os.path.isfile(outdir + outfile2):
        print (f"Running calculations for {pdb_i}, mutation {mut_res}{res_num_pdb}{wt_res} dist_cutoff = {dist_cutoff}, num_modes = {num_modes}...",flush=True)
        df_contact_breaks_mut_2 = calc_mut_energy_folded_unfolded(processed_pdb_dir, pdb_i, mut_res, wt_res, res_num_pdb, dist_cutoff, num_modes) # We switch the mutant and the wildtype
    
    # Add info on exp ddG if df_contact_breaks_mut variable exists
    if 'df_contact_breaks_mut' in locals():
        # Write the calculation output to a file
        df_contact_breaks_mut.to_csv(outdir + outfile, index=False)
    
    # If the mutation category is reverse, we will need the contact break info
    # on the forward mutation if not already present.
    if mut_category == 'Reverse' and 'df_contact_breaks_mut_2' in locals():
        df_contact_breaks_mut_2.to_csv(outdir + outfile2, index=False)
    return    


@click.command()
@click.option('--data_file', required=True, type=str, help='Name of the \
.csv file containing the information on ddG for the mutants')
@click.option('--outfile', required=True, type=str, help='Name of the file to \
    which the PSP-GNM-calculated energies and experimental energies will be written')
@click.option('--outdir', required=True, type=str, help='Name of the directory to \
    which the intermittent result files will be written to')    
@click.option('--wt_pdb_dir',required=True, type=str, help='Directory containing \
the wild type atomic pdb files')
@click.option('--num_jobs',required=True, type=str, help='Maximum number \
of jobs to be run in parallel')
@click.option('--dist_cutoff',required=True, default=9, type=str, help='Distance cutoff \
for interactions in GNM', show_default=True)
@click.option('--num_modes',required=True, default=10, type=str, help='Number \
of modes to be used', show_default=True)

def run_ab_initio_stability_prediction_wrapper(data_file, outfile, outdir, wt_pdb_dir, num_jobs, dist_cutoff, num_modes):
    # Wrapper function that parallely performs calculations for each 
    # First perform a sanity check on the mutant csv file- check how many records correctly
    # correspond to the residue position and whether the sequence included in the
    # mutant_csv_file has the specified residue at that particular position.

    # Input paramter definitions are same as described for the command-line arguments.

    # Writes the PSP-GNM-calculated ddG into outfile
    df_data = pd.read_csv(data_file, encoding='utf8')
    
    if not os.path.isdir(wt_pdb_dir):
        sys.exit(f"Error! {wt_pdb_dir} not found!")
    if not wt_pdb_dir.endswith('/'):
        wt_pdb_dir += '/'
    
    # Process the raw pdb files
    pdb_uniq = df_data['PDB_CHAIN'].unique().tolist()
    processed_pdb_dir = wt_pdb_dir[:-1] + '_processed/'
    process_wt_pdb(wt_pdb_dir, processed_pdb_dir, pdb_uniq)

    # First perform a sanity check on the mutant csv file- check how many records correctly
    # correspond to the residue position and whether the sequence included in the
    # mutant_csv_file has the specified residue at that particular position.
    print ("Running sanity check...", flush=True)    
    num_corr_map, serial_resnum_list = sanity_check(data_file, processed_pdb_dir)
    if num_corr_map < len(df_data):
        print (f"Only {num_corr_map}/{len(df_data)} records in {data_file} mapped correctly. Please fix the other records and then re-run.")
        sys.exit()
    else:
        print ("Done!", flush=True)

    # Include the serial residue numbers for mutant positions (as opposed to PDB Residue number)
    # obtained as a column in df_data
    df_data['RES_IND_SEQ'] = serial_resnum_list  

    if not outdir.endswith('/'):
        outdir += '/'
    # Create output directory if not already present
    if not os.path.isdir(outdir):
        os.mkdir(outdir)    
    
    # Convert num_jobs, dist_cutoff and num_modes to int
    num_jobs = int(num_jobs) 
    num_modes = int(num_modes)
    dist_cutoff = float(dist_cutoff)
    
    # Store the overall output here
    df_output_all = pd.DataFrame()

    # First run calculations for the wild type structures 
    Parallel(n_jobs=num_jobs)(delayed(run_ab_initio_stability_prediction_wildtype)(pdb_i, outdir, processed_pdb_dir, dist_cutoff, num_modes) for pdb_i in pdb_uniq)
    
    # Next run calculations for all the mutant rows
    Parallel(n_jobs=num_jobs)(delayed(run_ab_initio_stability_prediction_mutant)(row_i, outdir, processed_pdb_dir, dist_cutoff, num_modes) for idx, row_i in df_data.iterrows())
    
    # Define the columns in the output file
    col_list = list(df_data.columns) + ['Calc_ddG', 'Calc_ddI', 'Calc_ddG_mean', 'Calc_ddI_mean', 'No_contacts_broken']

    # Go through each row of protherm data and perform calculations using the 
    # intermediary contact breaks files.
    for idx_i, row_i in df_data.iterrows():
        # res_num is the residue number in the pdb file
        pdb_i, wt_res, mut_res, res_num_pdb, res_num_serial, mut_category = row_i['PDB_CHAIN'], row_i['WILD_RES'], row_i['MUTANT_RES'], row_i['RES_NUM_PDB'], row_i['RES_IND_SEQ'], row_i['Category']
        if mut_category == 'Reverse':
            # Switch the contact break files (i.e., wildtype contact break is now treated as the mutant
            # contact break file and the mutant contact break file as the wildtype)
            wt_cont_brk_file = pdb_i + '_' + row_i['MUTANT_RES'] + str(row_i['RES_NUM_PDB']) + row_i['WILD_RES'] + '_contact_breaks.csv'    
            mut_cont_brk_file = pdb_i + '_wt_contact_breaks.csv' 
        elif mut_category == 'Forward':    
            # If the mutation category is Forward, then read the original wild-type and mutant 
            # contact break files
            wt_cont_brk_file = pdb_i + '_wt_contact_breaks.csv'

            # Read the mutant contact breaks file
            mut_cont_brk_file = pdb_i + '_' + row_i['WILD_RES'] + str(row_i['RES_NUM_PDB']) + row_i['MUTANT_RES'] + '_contact_breaks.csv' 
        else:
            sys.exit(f"Invalid mutation category {mut_category} for {pdb_i} {wt_res}{res_num_pdb}{mut_res}. Mutation category must be either Forward or Reverse.")

        try:
            df_cont_brk_wt = pd.read_csv(outdir + wt_cont_brk_file)
        except pd.errors.EmptyDataError:
            print (f"{outdir + wt_cont_brk_file} is empty! No contacts involving mutation position broken during simulation! Will assign ddG value 0!")
            df_output_tmp = pd.DataFrame(data=[list(row_i) + [0, 0, 0, 0, True]],
                columns = col_list)
            df_output_all = df_output_all.append(df_output_tmp)
            continue
        
        try:        
            df_cont_brk_mut = pd.read_csv(outdir + mut_cont_brk_file)
        except pd.errors.EmptyDataError:
            print (f"{outdir + mut_cont_brk_file} is empty! No contacts involving mutation position broken during simulation! Will assign ddG value 0!")
            df_output_tmp = pd.DataFrame(data=[list(row_i) + [0, 0, 0, 0, True]],
                columns = col_list)
            df_output_all = df_output_all.append(df_output_tmp)
            continue    
        
        # Skip the record if no contacts are broken with the residue at the mutation position
        if len(df_cont_brk_mut) == 0:
            print (f"No contacts involving mutation position broken for {mut_category} mutant, {pdb_i} : {wt_res}{res_num_pdb}{mut_res}")
            continue
        elif len(df_cont_brk_wt) == 0:
            print (f"No contacts involving mutation position broken for wildtype of the {mut_category} mutant, {pdb_i} : {wt_res}{res_num_pdb}{mut_res}")
            continue

        # Calculate the theoretical ddG
        df_wt = df_cont_brk_wt.copy()
        df_mut = df_cont_brk_mut.copy()

        if mut_category == 'Forward':
            df_wt = df_wt.loc[(df_wt['PDB_ID'] == pdb_i) & (df_wt['WT_or_Mut'] == 'wt') & (df_wt['Mutation_position'] == res_num_serial)]
            df_mut = df_mut.loc[(df_mut['PDB_ID'] == pdb_i) & (df_mut['WT_or_Mut'] != 'wt') & (df_mut['Mutation_position'] == res_num_serial) & (df_mut['Res_at_Mut_Position'] == mut_res)]
        elif mut_category == 'Reverse':
            df_wt = df_wt.loc[(df_wt['PDB_ID'] == pdb_i) & (df_wt['Mutation_position'] == res_num_serial) ]
            df_mut = df_mut.loc[(df_mut['PDB_ID'] == pdb_i) & (df_mut['Mutation_position'] == res_num_serial) & (df_mut['Res_at_Mut_Position'] == mut_res)]

        # Drop duplicate rows
        df_wt.drop_duplicates(inplace=True)
        df_mut.drop_duplicates(inplace=True)
        # print (f"df_mut = {df_mut}")
        # print (f"df_wt = {df_wt}")
        if len(df_mut) == 0:
            print (f"No contacts involving mutant position broken while unfolding mutant structure of the {mut_category} mutant: {pdb_i} {wt_res}{res_num_pdb}{mut_res}. Will assign ddG value of 0!")
            df_output_tmp = pd.DataFrame(data=[list(row_i) + [0, 0, 0, 0, True]],
                columns = col_list)
            df_output_all = df_output_all.append(df_output_tmp)
            continue
        
        # If no contacts are broken in the wild type for the mutation position, then skip this position
        if len(df_wt) == 0:
            print (f"No contacts involving mutant position broken while unfolding wildtype structure of the {mut_category} mutant: {pdb_i} {wt_res}{res_num_pdb}{mut_res}. Will assign ddG value of 0!")
            df_output_tmp = pd.DataFrame(data=[list(row_i) + [0, 0, 0, 0, True]],
                columns = col_list)
            df_output_all = df_output_all.append(df_output_tmp)
            continue
        df_wt.reset_index(drop=True, inplace=True)
        df_mut.reset_index(drop=True, inplace=True)
        
        # Sort by contact break rank
        df_wt = df_wt.sort_values(by=['Contact_break_rank'])
        df_mut = df_mut.sort_values(by=['Contact_break_rank'])
        df_wt.reset_index(drop=True, inplace=True)
        df_mut.reset_index(drop=True, inplace=True)
        
        # Only consider the minimum number of contacts broken either
        # in the wild type or in the mutant
        min_len = len(df_wt)

        if len(df_mut) < min_len:
            min_len = len(df_mut)
        del_energy = df_mut['Energy_MJ'][0:min_len] - df_wt['Energy_MJ'][0:min_len]
        del_int_dist = df_mut['Int_dist_change'][0:min_len] - df_wt['Int_dist_change'][0:min_len]
        calc_ddG = sum(del_energy)
        calc_ddI = sum(del_int_dist)
        calc_ddG_mean = calc_ddG/len(del_energy)
        calc_ddI_mean = calc_ddI/len(del_int_dist)

        df_output_tmp = pd.DataFrame(data=[list(row_i) + [calc_ddG, calc_ddI, calc_ddG_mean, calc_ddI_mean, False]],
            columns = col_list)
        df_output_all = df_output_all.append(df_output_tmp)
    # Scale the calculated energy and ddI using the coefficients obtained by fitting
    # the calculated energy to the experimental energy
    
    # For the forward mutants we will scale using the coefficients obtained from the S350 data
    # Fit linear regression model for forward mutations
    df_output_all_fw = df_output_all.copy()
    df_output_all_fw = df_output_all_fw.loc[df_output_all_fw['Category'] == 'Forward']
    if len(df_output_all_fw) > 0 :
        #coeff = 0.11
        #intercept = -0.85
        if dist_cutoff == 9 and num_modes == 10:
            coeff = 0.08
            intercept = -1.01
        elif dist_cutoff == 9 and num_modes == 20:
            coeff = 0.08
            intercept = -0.99
        else:
            print (f"Scaling coefficients and intercepts for user-defined num_modes={num_modes} and dist_cutoff={dist_cutoff} unavailable.\n"
            "Using coefficients for num_modes = 20, dist_cutoff = 9 instead!")
            coeff = 0.08
            intercept = -0.99    
               
        calc_ddG_unscaled = -(df_output_all_fw['Calc_ddG']-df_output_all_fw['Calc_ddI'])
        ddG_PSP_GNM_fw = np.array(list(calc_ddG_unscaled))*coeff + intercept
        df_output_all_fw['ddG_PSP_GNM'] = ddG_PSP_GNM_fw

    # Scale the calculated ddG for the reverse mutations
    df_output_all_rev = df_output_all.copy()
    df_output_all_rev = df_output_all_rev.loc[df_output_all_rev['Category'] == 'Reverse']

    if len(df_output_all_rev) > 0:
        #coeff = 0.11
        #intercept = 0.85
        if dist_cutoff == 9 and num_modes == 10:
            coeff = 0.08
            intercept = -1.01
        elif dist_cutoff == 9 and num_modes == 20:
            coeff = 0.08
            intercept = -0.99
        else:
            print (f"Scaling coefficients and intercepts for user-defined num_modes={num_modes} and dist_cutoff={dist_cutoff} unavailable.\n"
            "Using coefficients for num_modes = 20, dist_cutoff = 9 instead!")
            coeff = 0.08
            intercept = -0.99    
        calc_ddG_unscaled = -(df_output_all_rev['Calc_ddG']-df_output_all_rev['Calc_ddI'])
        ddG_PSP_GNM_rev = np.array(list(calc_ddG_unscaled))*coeff + intercept
        df_output_all_rev['ddG_PSP_GNM'] = ddG_PSP_GNM_rev
        
    if len(df_output_all_fw) > 0 and len(df_output_all_rev) > 0:
        df_output_all_new = pd.concat([df_output_all_fw, df_output_all_rev])
    elif len(df_output_all_fw) > 0:
        df_output_all_new = df_output_all_fw
    else:
        df_output_all_new = df_output_all_rev                   

    df_output_all_new.to_csv(outfile, index=False)
    print (f"Wrote all calculations to {outfile}")


if __name__ == '__main__':
    run_ab_initio_stability_prediction_wrapper()       
