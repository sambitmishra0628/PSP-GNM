"""
Ab initio modeling of protein stability change upon point mutations using
a weighted Gaussian network model

- Coarse-grain a protein structure using only the alpha carbons

- Model the coarse-grained protein structure as a weighted 
Gaussian Network Model (GNM)

- The interactions between amino acid residues in such a model will
be weighted using statistical potentials

- Simulate unfolding of the wild-type structure by identifying residues
with most internal distance changes as the ones to break contact first

- Calculate the difference in energies (from the potential matrix) and
entropies (AACC matrix) for the first few contacts broken at the mutant position
for the mutant and wild-type structures (calc_ddG). This is where this code 
differs from the ab_initio_psp.py code. Assumption is that the unfolding energy
is dictated by the local contact energy change at the mutant site.

- Verify correlation between experimental ddG and calc_ddG

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

FILE_MAP = {'BT': 'BETM990101.txt', 'MJ': 'MIYS960101.txt', 'TS': 'TANS760101.txt', 'AACC': 'aacc.csv'}
POT_MAT_DIR = 'potential_matrices' # Directory containing the potential matrix
ENTROPY_MAT_DIR = 'entropy_matrices' # Directory containing the entropy matrix

def map_resname_to_id(res_code):
    """Convert the 3-lettered residue code to single letter"""
    resname_2_id = {'ALA' : 'A', 'ARG' : 'R', 'ASN' : 'N', 'ASP' : 'D',
                    'CYS' : 'C', 'GLY' : 'G', 'GLN' : 'Q', 'GLU' : 'E',
                    'HIS' : 'H', 'ILE' : 'I', 'LEU' : 'L', 'LYS' : 'K',
                    'MET' : 'M', 'PRO' : 'P', 'PHE' : 'F', 'SER' : 'S',
                    'THR' : 'T', 'TRP' : 'W', 'TYR' : 'Y', 'VAL' : 'V',}
    return resname_2_id[res_code]

def process_wt_pdb (input_dir, output_dir):
    """
    Process the pdb file of the wild types to filter alternative
    locations
    """
    if not input_dir.endswith('/'):
        input_dir = input_dir + '/'
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    pdbfiles = glob.glob(input_dir + '*.pdb')
    if len(pdbfiles) == len(glob.glob(output_dir + '*.pdb')):
        print (f"Processed pdb files already present at {output_dir}! Nothing to do!")
        return
    total_files = len(pdbfiles)
    count = 0    
    for pdbfile in pdbfiles:
        if not re.match('.*?\.pdb$', pdbfile):
            continue
        else:
            pdb_name = pdbfile.split('/')[-1]
            outfile = output_dir + pdb_name
            # print (f"pdb_name = {pdb_name}, outfile = {outfile}")
            # return
            fh1 = open(pdbfile, 'r')
            fh2 = open(outfile, 'w')
            # print("Processing file ", pdbfile)
            for line in fh1:
                if re.match('^ATOM', line):
                    alt_loc = line[16] # alternate location if any
                    alt_loc = alt_loc.replace(' ','')
                    if alt_loc == 'A': # we will consider by default the 'A' location
                        fh2.write(line)
                    elif alt_loc == '':
                        fh2.write(line)
                    else:
                        continue						
                elif re.match('^TER', line):
                    fh2.write(line)
                else:
                    continue
            fh1.close()
            fh2.close()
            # print("Done!\n")
            count += 1
            print (f"Processed {count}/{total_files} pdb files", end='\r', flush=True)

def parse_calpha_pdb (pdbfile):
    """
    Parse the coordinates of the c-alpha atoms from the given pdb file
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
    "Parse the potential file into a dictionary"
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

def parse_aacc(filename):
    "Parse the amino acid contact change file into a dictionary."
    df = pd.read_csv(ENTROPY_MAT_DIR + '/' + filename)
    aa_order = ['A','I','L','V','M','F','W','G','P','C','N','Q','S','T','Y','D','E','R','H','K']
    aacc_dict = {}
    for aa_x in aa_order:
        aa_vals = list(df[aa_x])
        for aa_y, cc_val in zip(aa_order, aa_vals):
            key1 = aa_x + ':' + aa_y
            key2 = aa_y + ':' + aa_x
            aacc_dict[key1] = cc_val
            aacc_dict[key2] = cc_val
    return aacc_dict

def energy_weighted_kirchoff(coord, res_codes, pot_type, cutoff, L, contact_matrix=None):
    """ Get the Kirchoff's matrix in which the contact springs are assigned based on
    the energetic interaction between residues.
    If contact matrix is provided, then use it. Otherwise, calculate the contact matrix based
    on the residue coordinates.
    """
    pot_file = FILE_MAP[pot_type]
    if pot_type != 'AACC':
        pot_dict = parse_energy_pot(pot_file)
        # Remove the identifer key "AA_INDEX_ID" from the
        # dictionary as it causes issues while sorting the 
        # dictionary.
        del(pot_dict['AA_INDEX_ID'])
    else:
        pot_dict = parse_aacc(pot_file)

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
                if pot_type == 'AACC':
                    pot_ij = (1/(min(pot_dict.values()))) + 10
                else:    
                    pot_ij = min(pot_dict.values())-1 
            else:    
                if pot_type == 'AACC':
                    pot_ij = (1/pot_dict[res_i + ':' + res_j])
                else:    
                    pot_ij = pot_dict[res_i + ':' + res_j]
            # The potential value is in energy terms. Convert to weights.
            if pot_type != 'AACC':
                pot_ij = round(np.exp(-pot_ij), 2)
            P[i,j] = pot_ij
    if contact_matrix.size != 0:
        # Typically useful when we want to provide a customized contact_matrix
        C = contact_matrix
    else:
        C, D = get_ca_contacts(coord, cutoff, r)
    energy_matrix = np.multiply(C,P)
    #D_inv = -(1/D)**2
    # energy_matrix = np.multiply()
    C = -C

    # We will weigh the residue-residue contacts by their energy potentials.
    K = np.multiply(C,P)
    # We will multiply the energy weights with inverse of distance
    # Take the Hadamard product of matrices D_inv and P.
    # This is basically the element-wise product of the 2 matrices
    #K = np.multiply(K,D_inv) 
    K_cpy = K.copy()
    # Distance dependent GNM
    #D_inv = -(1/D)
    #K = D_inv.copy()
    # print ("new")
    # print (np.shape(K))
    for i in range(0,r):
        # Diagonals in Kirchoff matrix are sum of the rows/columns
        # except the diagonal. 
        #K[i,i] = -(np.sum(C[i,:] * P[i,:] * D_inv[i,:])-K_cpy[i,i])
        K[i,i] = -(np.sum(C[i,:] * P[i,:])-K_cpy[i,i])
    return K, energy_matrix 

def get_ca_contacts(coord, cutoff, num_res):
    """
    Get the contact matrix between c-alpha atoms
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
    """Calculate the cross correlations between residues"""
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
    # Calculate the internal distance change given the cross-corr matrix
    r,c = np.shape(C)
    I = np.zeros((r,r))
    for i in range(0,r):
        for j in range(0,r):
            # Round the internal distance change to 3 decimal points
            I[i,j] = round((C[i,i] + C[j,j] - 2*C[i,j]), 3)
    return I

def calc_gnm(coord, cutoff=7.5, L=1, num_modes=20,
             weighted=False, spring_type=None, res_codes=None, 
             contact_matrix=None):
    """ Run calculations for GNM.
    The GNM can be either unweighted (default) or weighted (either by potentials or
                                                            entropy)
    """
    if weighted:
        """ Run calculations for GNM in which the interactions are weighted
        between the residues.
    
        spring_type = Type of weighting for the interactions between residues.
        The interaction spring can be either a potential-based (Betancourt-Thirumalai potential
                                                                or Miyazawa-Jernigan potential) or may be entropy-based) 
    
        spring_type: 'MJ' (default), 'BT'
                 'AACC'(entropy)
        """
        if spring_type == None:
            sys.exit("Error! spring_type cannot be None for weighted GNM")
        if type(res_codes) != list:
            # print(type(res_codes))
            sys.exit("Error! res_codes must be a list of 3-letter residue codes of all residues")

        stat_pot = ['BT', 'MJ', 'AACC']
        if spring_type in stat_pot:
            K, e_matrix = energy_weighted_kirchoff(coord, res_codes, spring_type, cutoff, L, contact_matrix)
        #elif spring_type == 'AACC':
        #    K, e_matrix = entropy_weighted_kirchoff(coord, res_codes, spring_type, cutoff)
        else:
            sys.exit(f"Invalid spring type {spring_type}")
    else:    
        K = get_kirchoff(coord, cutoff, L)
    E,V = linalg.eigh(K, turbo=True) # Returns eigen values and vectors
    
    #idx = E.argsort()[::1]
    #V = V[idx[1:]] # Skip the eigen vectors with 0 
    #V = V[idx]
    # print (E[0], E[1:10])
    #E = E[idx] # Skip the eigen values of 0
    return V[:,1:num_modes+1], E[1:num_modes+1], e_matrix

def get_kirchoff(coord, cutoff=3, L=1):
    """Calculates the unweighted GNM Kirchoff matrix.
    Parameters - 
        coord - N x 3 matrix of coordinates
        cutoff - distance cutoff 
        L - lambda, a fixed integer corresponding to the stiffness of the 
            spring
    """
    # Declare empty numpy array 
    r,c = np.shape(coord)
    C, D = get_contacts(coord, cutoff, r)
    K = -C*L # Off diagonal elements in contact are allocated -L
    print (np.shape(K))
    print (f"cutoff= {cutoff}, L = {L}")
    for i in range(0,r):
        K[i,i] = sum(C[i])*L # Diagonals in Kirchoff matrix are sum of the rows/columns
        print ("Changed...")
    return K     

def sanity_check(protherm_data_file):
    """
    Check how many records in the mutant csv file have correct correspondence
    between the mutated amino acid and its position in the wildtype sequence.

    Parameters
    ----------
    mutant_csv_file : str
        The name of the csv file containing information about the mutants

    Returns
    -------
    A list containing match status (MATCH or MISMATCH) for each record, based on
    whether a residue at a given position could be mapped on to the PDB sequence
    or not.

    """
    df = pd.read_csv(protherm_data_file)
    total_recs = len(df)
    correct_recs = 0
    match_status = []
    for res_i,pos_i,seq_i in zip(df['WILD_RES'], df['RES_IND_SEQ'], df['PDB_SEQUENCE']):
        if (not np.isnan(pos_i) and seq_i[int(pos_i)] ==  res_i):
            correct_recs += 1
            match_status.append('MATCH')
        else:
            match_status.append('MISMATCH')
    print (f"{correct_recs}/{total_recs} records show correct mutant position and amino acid match")
    print (len(match_status))
    return match_status

def create_contact_map_fig(C,D,figfile, res_codes):
    plot_rows, plot_cols = 1,2
    num_res = len(res_codes)
    plt.figure()
    if len(D) == 0:
        plot_rows, plot_cols = 1,1
    if len(C) > 0:
        plt.subplot(plot_rows, plot_cols, 1)
        ax1 = sns.heatmap(C,cmap='jet', square=True, cbar=False)
        ax1.invert_yaxis()
        plt.title('C-alpha Contact Map')
        x_ticks = list(np.arange(0, len(res_codes), 10))
        y_ticks = list(np.arange(0, len(res_codes), 10))
        ax1.set_xticks(range(0,num_res, 10))
        ax1.set_xticklabels(x_ticks)
        ax1.set_yticks(range(0,num_res, 10))
        ax1.set_yticklabels(y_ticks)
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)
    if len(D) > 0:
        plt.subplot(plot_rows, plot_cols, 2)
        ax2 = sns.heatmap(D,cmap='jet', square=True, cbar=False)
        ax2.invert_yaxis()
        plt.title('C-alpha Distance Map')
        x_ticks = list(np.arange(0, len(res_codes), 10))
        y_ticks = list(np.arange(0, len(res_codes), 10))
        ax2.set_xticks(range(0,num_res, 10))
        ax2.set_xticklabels(x_ticks)
        ax2.set_yticks(range(0,num_res, 10))
        ax2.set_yticklabels(y_ticks)
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)
    plt.tight_layout()
    plt.savefig(figfile, dpi=300)
    plt.close()

def simulate_unfolding(ca_coord, res_codes, pdb_id, dist_cutoff, num_modes, mut_or_wt='wt',serial_res_num=None):
    """
    Simulate unfolding based on change in internal distance between residues.
    We will simulate the unfolding until 50 percent contacts made by the mutant position
    are broken.
    """
    print (f"Num modes = {num_modes}, dist_cutoff = {dist_cutoff}")
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
    print (f"Total contacts in folded protein = {tot_contacts_folded}")
    if mut_or_wt != 'wt':
        res_at_mut_pos = res_codes[serial_res_num] # Get the residue at the mutant position
   
    # When simulating for mutant, then get contacts for mutant position
    if mut_or_wt != 'wt':
        mut_pos_contacts = list(np.where(C[serial_res_num,:] == 1)[0])
        mut_pos_tot_contacts = len(mut_pos_contacts) * 2 # Consider both (i,j) and (j,i) pairs
    
    # Parse the potential matrix
    if int_potential == 'AACC':
        pot_file = FILE_MAP['MJ'] # Contingency fix. Change later.
    else:    
        pot_file = FILE_MAP[int_potential]
    pot_dict = parse_energy_pot(pot_file)
    del(pot_dict['AA_INDEX_ID'])    # Remove the identifer key "AA_INDEX_ID"

    # Parse the entropy matrix
    entropy_mat = 'AACC'
    entropy_file = FILE_MAP[entropy_mat]
    entropy_dict = parse_aacc(entropy_file)

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
        V,E, e_matrix = calc_gnm(ca_coord, dist_cutoff, L_adj, num_modes, True, int_potential, res_codes, contact_matrix)
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
                        
                        column_names = ['PDB_ID', 'WT_or_Mut', 'Mutation_position', 'Contact_Position', 'Res_at_Mut_Position', 'Res_at_Contact_Pos', 'Energy_MJ', 'Entropy_AACC', 'Int_dist_change', 'Contact_break_rank']
                        df_tmp = pd.DataFrame([[pdb_id, mut_or_wt, serial_res_num, pair_res_num, res_code_i,
                            res_code_j, pot_dict[res_code_i + ':' + res_code_j], entropy_dict[res_code_i + ':' + res_code_j],  max_int_dist_val, iteration]], columns=column_names)
                        df_contact_breaks = df_contact_breaks.append(df_tmp)
                else:
                    # If simulation is being done for a wildtype structure
                    res_code_i = res_codes[row_index]
                    res_code_j = res_codes[col_index]
                    column_names = ['PDB_ID', 'WT_or_Mut', 'Mutation_position', 'Contact_Position', 'Res_at_Mut_Position', 'Res_at_Contact_Pos', 'Energy_MJ', 'Entropy_AACC', 'Int_dist_change', 'Contact_break_rank']
                    df_tmp = pd.DataFrame([[pdb_id, mut_or_wt, row_index, col_index, res_code_i,
                            res_code_j, pot_dict[res_code_i + ':' + res_code_j], entropy_dict[res_code_i + ':' + res_code_j], max_int_dist_val, iteration]], columns=column_names)
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
    Calculate the interaction energy of the folded and unfolded mutant structure.
    res_num is the residue number in the pdb file
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
    Calculate the energy of the wild-type structure in the folded and 
    partially unfolded states
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
    Run calculations first for all the wild-types
    """
    outfile = pdb_i + '_wt_contact_breaks.csv'
    
    # Do not run calculations if outputfile already exists
    if os.path.isfile(outdir + outfile):
        print (f"Skipping calculations for wildtype {pdb_i} as file is already present ")
        return
    print (f"Running calculations for {pdb_i} wild type, dist_cutoff = {dist_cutoff}, num_modes = {num_modes}...", end='', flush=True)
    df_contact_breaks_wt = calc_wt_energy_folded_unfolded(processed_pdb_dir, pdb_i,dist_cutoff, num_modes)
    
    # Add info on exp ddG
    df_contact_breaks_wt['EXP_ddG'] = ['']*len(df_contact_breaks_wt)
    print ("Done!")

    # Write the calculation output to a file
    df_contact_breaks_wt.to_csv(outdir + outfile, index=False)

def run_ab_initio_stability_prediction_mutant(row_i, outdir, processed_pdb_dir, dist_cutoff, num_modes):
    """
    Run calculations for the mutant position
    """
    pdb_i, wt_res, mut_res, res_num_pdb, res_num_serial = row_i['PDB_CHAIN'], row_i['WILD_RES'], row_i['MUTANT_RES'], row_i['RES_NUM_PDB'], row_i['RES_IND_SEQ']

    # Skip calculations if output file is already present
    outfile = pdb_i + '_' + row_i['WILD_RES'] + str(row_i['RES_NUM_PDB']) + row_i['MUTANT_RES'] + '_contact_breaks.csv'
    if os.path.isfile(outdir + outfile):
        print (f"Skipping calculations for mutant {pdb_i}: {row_i['WILD_RES']}{row_i['RES_NUM_PDB']}{row_i['MUTANT_RES']} as file is already present!")
        return
    print (f"Running calculations for {pdb_i}, mutation {wt_res}{res_num_pdb}{mut_res} dist_cutoff = {dist_cutoff}, num_modes = {num_modes}...", end='', flush=True)
    df_contact_breaks_mut = calc_mut_energy_folded_unfolded(processed_pdb_dir, pdb_i, wt_res, mut_res, res_num_pdb, dist_cutoff, num_modes)
    
    # Add info on exp ddG
    df_contact_breaks_mut['EXP_ddG'] = [row_i['EXP_DDG']]*len(df_contact_breaks_mut)
    print ("Done!")

    # Write the calculation output to a file
    df_contact_breaks_mut.to_csv(outdir + outfile, index=False)

@click.command()
@click.option('--protherm_data_file', required=True, type=str, help='Name of the \
.csv file containing the information on ddG for the mutants')
@click.option('--outfile', required=True, type=str, help='Name of the file to \
    which the theoretical energies and experimental energies will be written')
@click.option('--outdir', required=True, type=str, help='Name of the directory to \
    which the intermittent result files will be written to')    
@click.option('--wt_pdb_dir',required=True, type=str, help='Directory containing \
the wild type atomic pdb files')
@click.option('--num_jobs',required=True, type=str, help='Maximum number \
of jobs to be run in parallel')
@click.option('--dist_cutoff',required=True, type=str, help='Distance cutoff \
for interactions in GNM')
@click.option('--num_modes',required=True, type=str, help='Number \
of modes to be used')

def run_ab_initio_stability_prediction_wrapper(protherm_data_file, outfile, outdir, wt_pdb_dir, num_jobs, dist_cutoff, num_modes):
    # Wrapper function that parallely performs calculations for each 
    # First perform a sanity check on the mutant csv file- check how many records correctly
    # correspond to the residue position and whether the sequence included in the
    # mutant_csv_file has the specified residue at that particular position.
    print ("Running sanity check...")    
    match_status = sanity_check(protherm_data_file)
    print  ("Done!")
    if not os.path.isdir(wt_pdb_dir):
        sys.exit(f"Error! {wt_pdb_dir} not found!")
    if not wt_pdb_dir.endswith('/'):
        wt_pdb_dir += '/'
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

    # Process the raw pdb files
    processed_pdb_dir = wt_pdb_dir[:-1] + '_processed/'
    process_wt_pdb(wt_pdb_dir, processed_pdb_dir)

    # Make predictions for each mutant type
    df_protherm_data = pd.read_csv(protherm_data_file)
    column_names = ['PDB', 'WT_res', 'Mut_res', 'Mut_Position', 'PH', 
    'Temperature', 'Experimental_ddG', 'Calc_energy_folded_WT', 'Calc_energy_unfolded_WT', 'Calc_dG_WT', 
    'Calc_energy_folded_mutant', 'Calc_energy_unfolded_mutant', 'Calc_dG_mutant', 'Calc_ddG']
    count = 0
    pdb_uniq = df_protherm_data['PDB_CHAIN'].unique().tolist()
    print (pdb_uniq)
    # First run calculations for the wild type structures 
    Parallel(n_jobs=num_jobs)(delayed(run_ab_initio_stability_prediction_wildtype)(pdb_i, outdir, processed_pdb_dir, dist_cutoff, num_modes) for pdb_i in pdb_uniq)
    
    # Next run calculations for all the mutant rows
    Parallel(n_jobs=num_jobs)(delayed(run_ab_initio_stability_prediction_mutant)(row_i, outdir, processed_pdb_dir, dist_cutoff, num_modes) for idx, row_i in df_protherm_data.iterrows())

    # Go through each row of protherm data and perform calculations using the 
    # intermediary contact breaks files.
    for idx_i, row_i in df_protherm_data.iterrows():
        # res_num is the residue number in the pdb file
        pdb_i, wt_res, mut_res, res_num_pdb, res_num_serial = row_i['PDB_CHAIN'], row_i['WILD_RES'], row_i['MUTANT_RES'], row_i['RES_NUM_PDB'], row_i['RES_IND_SEQ']
        
        # Read the wild-type contact breaks file
        wt_cont_brk_file = pdb_i + '_wt_contact_breaks.csv' 
        df_cont_brk_wt = pd.read_csv(outdir + wt_cont_brk_file)

        # Read the mutant contact breaks file
        mut_cont_brk_file = pdb_i + '_' + row_i['WILD_RES'] + str(row_i['RES_NUM_PDB']) + row_i['MUTANT_RES'] + '_contact_breaks.csv'
        df_cont_brk_mut = pd.read_csv(outdir + mut_cont_brk_file)
        
        # Skip the record if no contacts are broken with the residue at the mutation position
        if len(df_cont_brk_mut) == 0:
            continue
        # Calculate the theoretical ddG
        df_wt = df_cont_brk_wt.copy()
        df_mut = df_cont_brk_mut.copy()
        df_wt = df_wt.loc[(df_wt['PDB_ID'] == pdb_i) & (df_wt['WT_or_Mut'] == 'wt') & (df_wt['Mutation_position'] == res_num_serial)]
        df_mut = df_mut.loc[(df_mut['PDB_ID'] == pdb_i) & (df_mut['WT_or_Mut'] != 'wt') & (df_mut['Mutation_position'] == res_num_serial) & (df_mut['Res_at_Mut_Position'] == mut_res)]
        
        # If no contacts are broken in the wild type for the mutation position, then skip this position
        if len(df_wt) == 0:
            continue
            
        # Drop duplicate rows
        df_wt.drop_duplicates(inplace=True)
        df_mut.drop_duplicates(inplace=True)
        print (f"df_mut = {df_mut}")
        print (f"df_wt = {df_wt}")
        if len(df_mut) == 0:
            print (f"No contacts involving mutant residue broken while unfolding mutant structure {pdb_i} {wt_res}{res_num_pdb}{mut_res}")
            continue

        df_wt.reset_index(drop=True, inplace=True)
        df_mut.reset_index(drop=True, inplace=True)
        exp_ddG = df_mut['EXP_ddG'][0]

        # Sort by contact break rank
        df_wt = df_wt.sort_values(by=['Contact_break_rank'])
        df_mut = df_mut.sort_values(by=['Contact_break_rank'])
        df_wt.reset_index(drop=True, inplace=True)
        df_mut.reset_index(drop=True, inplace=True)

        # print (df_wt)
        # print (df_mut)
        
        # Only consider the minimum number of contacts broken either
        # in the wild type or in the mutant
        min_len = len(df_wt)
        #if min_len > 6:
        #    min_len = 6
        if len(df_mut) < min_len:
            min_len = len(df_mut)
        del_energy = df_mut['Energy_MJ'][0:min_len] - df_wt['Energy_MJ'][0:min_len]
        del_entropy = df_mut['Entropy_AACC'][0:min_len] -df_wt['Entropy_AACC'][0:min_len]
        del_int_dist = df_mut['Int_dist_change'][0:min_len] - df_wt['Int_dist_change'][0:min_len]
        calc_ddG = sum(del_energy)
        calc_ddE = sum(del_entropy)
        calc_ddI = sum(del_int_dist)
        calc_ddG_mean = calc_ddG/len(del_energy)
        calc_ddE_mean = calc_ddE/len(del_entropy)
        calc_ddI_mean = calc_ddI/len(del_int_dist)
        
        df_output_tmp = pd.DataFrame(data=[[pdb_i, wt_res, mut_res, res_num_pdb, res_num_serial, row_i['PH'], row_i['TEMPERATURE'], exp_ddG, calc_ddG, calc_ddE, calc_ddI, calc_ddG_mean, calc_ddE_mean, calc_ddI_mean]],
            columns = ['PDB_ID', 'WT_Residue', 'Mutant_Residue', 'Res_Num_PDB', 'Res_Num_Serial', 'PH', 'Temperature','Exp_ddG', 'Calc_ddG', 'Calc_ddE', 'Calc_ddI', 'Calc_ddG_mean', 'Calc_ddE_mean', 'Calc_ddI_mean'])
        df_output_all = df_output_all.append(df_output_tmp)
    df_output_all.to_csv(outfile, index=False)
    print (f"Wrote all calculations to {outfile}")


if __name__ == '__main__':
    run_ab_initio_stability_prediction_wrapper()       
