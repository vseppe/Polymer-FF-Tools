#!/usr/bin/env python3
# title             : link_interactions.py
# description       : Writes the links.ff file for the polyply input. Either two residues or argument "same" if bonding to itself. Looks in topology files for extra bonded parameters. Supply 3- to 4-lettered resnames or same if same residue.
# use               : python link_interactions_ff.py $itp_res1_uncapped.ff $itp_res2_uncapped/"same" $resname1 $resname2/"same" $connection_atomnames.txt $parent_forcefield.itp/.dat
# input file        : $itp_res1_uncapped.ff, $itp_res2_uncapped.ff
# output file       : links.ff

import numpy as np
import sys
from pathlib import Path 
import datetime
import os
import copy
import re
from collections import Counter
from remove_caps_ff import section_block, read_itp_res, find_nonint_lines 

# Functions
# Inputs:
def load_input_command():
    """ Process the arguments given in the command line. """
    if len(sys.argv) != 7:
        print(f"Usage: {sys.argv[0]} <itp_res1_uncapped.ff> <itp_res2_uncapped.ff>/\"same\" <resname1> <resname2>/\"same\" <connection_atomnames.txt> <parent_forcefield.itp/.dat>")
        sys.exit(1)
    
    itp_res1_file = Path(sys.argv[1])
    itp_res2_file = Path(sys.argv[2])
    res1_name = sys.argv[3]
    res2_name = sys.argv[4]
    connection_atomnames_file = Path(sys.argv[5])
    parent_ff_file = Path(sys.argv[6])

    if not itp_res1_file.is_file():
        print(f"Error: itp_res1_uncapped.ff input file '{itp_res1_file}' does not exist.")
        sys.exit(1)
    
    if not itp_res2_file.is_file():
        if str(itp_res2_file) == "same":
            itp_res2_file = itp_res1_file
        else:
            print(f"Error: itp_res2_uncapped.ff input file '{itp_res2_file}' does not exist. Either give the name or write \"same\" if within the same residue.")
            sys.exit(1)

    if len(res1_name) > 4:
        print(f"Error: Residue names may not exceed a length of 4 characters, not \'{res1_name}\' with length {len(res1_name)}")
        sys.exit(1)                    

    if len(res2_name) > 4:
        print(f"Error: Residue names may not exceed a length of 4 characters, not \'{res2_name}\' with length {len(res2_name)}")
        sys.exit(1)                    

    if not connection_atomnames_file.is_file():
        print(f"Error: connection_atomnames_file.txt input file '{connection_atomnames_file}' does not exist.")
        sys.exit(1)

    if not parent_ff_file.is_file():
        print(f"Error: {parent_ff_file} is not an existing file. Make sure it is either a .dat file in AMBER style for the bonded (bonds, angles, dihedrals) parameters, or a GROMACS ffbonded.itp type of file (found in gromacs/share/gromacs/top directory)")
        sys.exit(1)

    return itp_res1_file, itp_res2_file, res1_name, res2_name, connection_atomnames_file, parent_ff_file

def ff_type(parent_ff_file: str) -> str:
    """ Determines the type of parent force field file. Either AMBER (.dat) or GROMACS (.itp)."""

    extension_parent_ff = str(parent_ff_file)[-4:]

    if extension_parent_ff == ".dat":
        ff_type = "amber"
    elif extension_parent_ff == ".itp":
        ff_type = "gromacs"
    else: # No correct force field
        print(f"Parent force field has an unrecognizable extension \'{extension_parent_ff}\'.")
        print(".dat should be used for AMBER-style force fields, and .itp for GROMACS-style ffbonded.itp force field files.")
        print("Make sure that the parameters for bonds, angles, and dihedrals are in there.")
        print("Possibly consider changing the extension.")
        print("Exiting.")
        sys.exit(1)

    return ff_type

def read_connection_file(connection_atomnames_file: str, atom_names_res1: str, atom_names_res2: str) -> list:
    """ Read the input file connection_atomnames.txt, structured as: atomname res 1 atomname res 2 (eg C1 O1).
        Returns list of tuples (expected only 1 connection). 
        Also checks whether atom names in this file match with the read in atom names of the itp files. """

    # Bring everything to a list of the file
    with open(connection_atomnames_file, "r") as f:
        connection_lines = f.readlines()

    # Add atom connections to the list
    connections = []
    for line_i, line in enumerate(connection_lines):
        if not line.lstrip():
            # empty line, do not do anything
            pass
        else:
            line_list = line.split()
            if line_list[0][0] == "#":
                # comment, do not read
                pass
            elif len(line_list) > 2 and line_list[2][0] != "#": # stop if more than two columns and no comments
                print(f"There are more than two arguments in some lines of the connection file {connection_atomnames_file}, namely")
                print(f"\'{line_list[2:]}\'. Exiting.")
                sys.exit(1)
            else:
                atom1 = str(line_list[0])
                atom2 = str(line_list[1])

                if atom1.isdigit() or atom2.isdigit(): # stop if they are indices instead of atom names
                    print(f"There are plain integers in the connection file {connection_atomnames_file}, such as in line \'{line}\'.")
                    print(f"This file should not contain indices, but atom names. Exiting.")
                    sys.exit(1)
                else:   
                    connections.append([atom1, atom2])

    # Veryify whether atom names appear in .ff
    for connection in connections:
        if not connection[0] in atom_names_res1:
            print(f"The atom name {connection[0]} does not appear in the .ff file of the first residue. Exiting.")
            sys.exit(1)
        if not connection[1] in atom_names_res2: 
            print(f"The atom name {connection[1]} does not appear in the .ff file of the second residue. Exiting.")
            sys.exit(1)

    return connections

def read_itp_atomsbondspairs(itp_res_file: str):
    """ Reads the atoms, bonds and pairs blocks in an itp file and returns the data (without titles or comments).
        This function is dependent on functions in remove_caps_ff.py. """ 

    # Read the sections from the itp file. 
    section_names = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals"]
    itp_lines, section_linenrs, sections = read_itp_res(itp_res_file, section_names)

    # Only take the lines with data (int from integer), not the comments etc
    block_int = []
    for interaction_type in ["atoms", "bonds", "pairs"]:
        full_block = section_block(itp_lines, sections, interaction_type)
        nonint_lines = find_nonint_lines(full_block)
        block_int.append([ x for i, x in enumerate(full_block) if i not in nonint_lines ]) # indices of list match atom indices-1 
    
    atoms_block_int = block_int[0]
    bonds_block_int = block_int[1]
    pairs_block_int = block_int[2]
    
    # Dictionary of atom names to atom index, atom types array
    atom_types = []
    atom_names = []
    for line in atoms_block_int:
        line_list = line.split()
        atom_types.append(line_list[1])
        atom_names.append(line_list[4])

    namestoidx = dict(zip(atom_names, range(1, len(atom_names)+1)))

    return  atoms_block_int, bonds_block_int, pairs_block_int, namestoidx, atom_types, atom_names

def write_connectivity(connectivity: np.array, bonds_apart: int, atom_idx1: int, atom_idx2: int):
    """ In the connectivity array, writes the newly given connectivity and updates it. 
        Writes for both lower and upper triangle.
        Writes the lowest connectivity between the given value and the value in the matrix.
        The atom index is the input: the coordinate is -1, since 1-indexed in Gromacs."""

    # From 1-basis to 0-basis
    atom_idx1 -= 1    
    atom_idx2 -= 1    

    # Write the minimum in the matrix, both upper and lower
    bonds_apart_value = min(bonds_apart, connectivity[atom_idx1, atom_idx2])
    connectivity[atom_idx1, atom_idx2] = bonds_apart_value
    connectivity[atom_idx2, atom_idx1] = bonds_apart_value
    return connectivity 

def determine_connectivity(itp_file: str, atoms_block: list, bonds_block: list) -> np.array:
    """ For every atom in the itp file, determine the connectivity, i.e. how many bonds are in between atom A and B.
        This will be stored as a triangular matrix with each element ij the number of bonds separating atom i and atom j. 0 on diagonals.
        It uses the gromacs itp bonds section to determine what is bonded and what is not."""

    # Make square connectivity matrix
    N_atoms = len(atoms_block)
    connectivity = np.empty((N_atoms, N_atoms)) # initialize empty
    connectivity[:] = np.nan

    for i in range(N_atoms):
        connectivity[i, i] = 0 # diagonal

    # Cycle through the elements to find atoms that are bonded directly: 1 in matrix
    for line in bonds_block: # the atom indices are 1-based in gromacs, here the array indices are 0-based
        line_list = line.split()
        atom_idx1 = int(line_list[0])
        atom_idx2 = int(line_list[1])
        connectivity = write_connectivity(connectivity, 1, atom_idx1, atom_idx2)
        
    # Cycle through the matrix per increment of bond_distance (eg 2 apart, 3 apart etc)
    bonds_apart = 1
    while not np.all(~np.isnan(connectivity)): # While there are still NaNs (i.e. not been assigned a value)
        for atom_1 in range(N_atoms): # Using python 0-indexing here! 
            specified_atoms = np.where(connectivity[atom_1, :] == bonds_apart)
            for atom_2 in specified_atoms[0]:
                # There are atoms at this number of bonds_apart of atom_idx_1
                atoms_bonded_atom2 = np.where(connectivity[atom_2, :] == 1)
                for atom_bonded_atom2 in atoms_bonded_atom2[0]:
                    # Atom 2 has atoms bonded to it directly (should do!)
                    connectivity = write_connectivity(connectivity, bonds_apart + 1, atom_1+1, atom_bonded_atom2+1) 
        bonds_apart += 1 

    return connectivity

def determine_pair_thresholds(pairs_block: list, connectivity: np.array):
    """ Determine for which bond count there are pairs in the pairs section. """
    bond_counts = []
    for line in pairs_block:
        line_list = line.split()
        atom_i = int(line_list[0]) - 1
        atom_j = int(line_list[1]) - 1
        bond_counts.append(int(connectivity[atom_i, atom_j]))
        
    pair_threshold_min = min(bond_counts)
    pair_threshold_max = max(bond_counts)

    return (pair_threshold_min, pair_threshold_max)

def link_atom_names(connections: list, namestoidx_res1: dict, namestoidx_res2: dict, connectivity_res1: np.array, connectivity_res2: np.array, atom_names_res1: list, atom_names_res2: list, bonds_block_int_res1: list, bonds_block_int_res2: list, pair_thresholds: tuple):
    """ Write the atom names involved in the inter-residue interactions, i.e. those needed to be written to the link.ff file. """

    # ATOMS
    atoms_names = []
    for connection in connections:
        atom_1 = connection[0]
        atom_2 = "+" + connection[1] # Prepend +
        atoms_names.append(atom_1)
        atoms_names.append(atom_2)
 
    # BONDS
    bonds_names = copy.deepcopy(connections) # This is the same, but a plus in front of the second atom name. Deepcopy since nested lists.
    for atom in bonds_names:
        atom[1] = "+" + atom[1] # Prepending + for link file

    # EDGES
    edges_names = []
    bonds_blocks = [bonds_block_int_res1, bonds_block_int_res2]
    atom_names = [atom_names_res1, atom_names_res2]
    plusornot = "" # for res1, no plus
    for i, bonds_block in enumerate(bonds_blocks):
        for line in bonds_block:
            line_list = line.split()
            atom_name_1 = plusornot + atom_names[i][int(line_list[0])-1] # the bonds block has 1-based indices
            atom_name_2 = plusornot + atom_names[i][int(line_list[1])-1]
            edges_names.append([atom_name_1, atom_name_2])

        if i == 0:
            edges_names.extend(bonds_names) # Append the inter-residue bonds in between res1 and res2
            plusornot = "+" # plus for res2 coming up 


    # Determining cross-residue connectivity. Should work if multiple connections, but not expected.
    # In 0-based indices
    all_res1_border1 = []
    all_res2_border1 = []
    all_res1_border2 = []
    all_res2_border2 = []
    for connection in connections:    
        atom1_idx = namestoidx_res1[connection[0]] - 1 # atomX_idx is 1-based, but resX_borderX_idx is 0-based since matrix coords. Making both 0-based here.
        atom2_idx = namestoidx_res2[connection[1]] - 1
    
        # Directly bound to edge atoms  Z--A-|-B--C.
        res1_border1_idx = np.where(connectivity_res1[atom1_idx, :] == 1)[0]
        res2_border1_idx =  np.where(connectivity_res2[atom2_idx, :] == 1)[0]

        # Atoms bound to those atoms, thus two away from the edges Y--Z--A-|-B--C--D
        res1_border2 = []
        res2_border2 = []
        for res1_border1 in res1_border1_idx:
            res1_border2_element = np.where(connectivity_res1[res1_border1, :] == 1)[0]
            if not atom1_idx in res1_border1_idx: # Atom A (see schematic) must not be included in a A--Z--A angle etc
                res1_border2.append(res1_border2_element)            

        for res2_border1 in res2_border1_idx:
            res2_border2_element = np.where(connectivity_res2[res2_border1, :] == 1)[0]
            if not atom2_idx in res2_border1_idx: # Atom B (see schematic) must not be included in a B--C--B angle etc
                res2_border2.append(res2_border2_element)            
 
        # so eg [[0, 1]] in all_res1_border1 and [[2, 3], [4, 5]] in all_res1_border2 if only 1 connection 
        all_res1_border1.append(res1_border1_idx)
        all_res2_border1.append(res2_border1_idx)
        all_res1_border2.append(res1_border2)
        all_res2_border2.append(res2_border2)      

   
    # ANGLES
    angles_names = []
    # per connection, retrieve Z or C in Z--A-|-B, A-|-B--C
    for i, connection in enumerate(connections):
        atom_name_center_res1 = connection[0]
        atom_name_center_res2 = connection[1]
        atom_idx_center_res1 = namestoidx_res1[atom_name_center_res1] - 1 # making 0-based
        atom_idx_center_res2 = namestoidx_res2[atom_name_center_res2] - 1 
        atom_name_center_res2 = "+" + atom_name_center_res2 # res2 gets a + in front
        for res1_border1_el in all_res1_border1[i]: # Z--A-|-B
            atom_name_border_res1 = atom_names_res1[res1_border1_el]
            angle_atom_names = [atom_name_border_res1, atom_name_center_res1, atom_name_center_res2]
            angles_names.append(angle_atom_names)

        for res2_border1_el in all_res2_border1[i]: # same for res2, A-|-B--C
            atom_name_border_res2 = "+" + atom_names_res2[res2_border1_el] # + in front
            angle_atom_names = [atom_name_center_res1, atom_name_center_res2, atom_name_border_res2]
            angles_names.append(angle_atom_names)


    # DIHEDRALS: PROPER
    dihedrals_names = []
    # per connection, retrieve: Y and Z in Y--Z--A-|-B, Z and C in Z--A-|-B--C, or C and D in A-|-B--C--D
    for i, connection in enumerate(connections):
        atom_name_center_res1 = connection[0]
        atom_name_center_res2 = connection[1]
        atom_idx_center_res1 = namestoidx_res1[atom_name_center_res1] - 1 # making 0-based
        atom_idx_center_res2 = namestoidx_res2[atom_name_center_res2] - 1
        atom_name_center_res2 = "+" + atom_name_center_res2  # res2 gets a + in front

        for border1_i, res1_border1_el in enumerate(all_res1_border1[i]): #Z--A-|-B
            atom_name_border1_res1 = atom_names_res1[res1_border1_el]
            for res1_border2_el in all_res1_border2[i][border1_i]: #Y--Z--A-|-B
                atom_name_border2_res1 = atom_names_res1[res1_border2_el]
                if not atom_name_border2_res1 in [atom_name_border1_res1, atom_name_center_res1, atom_name_center_res2]: # if Y is not A
                    dihedral_atom_names = [atom_name_border2_res1, atom_name_border1_res1, atom_name_center_res1, atom_name_center_res2]
                    dihedrals_names.append(dihedral_atom_names)

            for res2_border1_el in all_res2_border1[i]: #Z--A-|-B--C
                atom_name_border2_res1 = "+" + atom_names_res2[res2_border1_el]
                dihedral_atom_names = [atom_name_border1_res1, atom_name_center_res1, atom_name_center_res2, atom_name_border2_res1]
                dihedrals_names.append(dihedral_atom_names)
            
        for border1_i, res2_border1_el in enumerate(all_res2_border1[i]): # A-|-B--C--D
            atom_name_border1_res2 = "+" + atom_names_res2[res2_border1_el]
            for res2_border2_el in all_res2_border2[i][border1_i]: # element D 
                atom_name_border2_res2 = "+" + atom_names_res2[res2_border2_el]
                if not atom_name_border2_res2 in [atom_name_center_res1, atom_name_center_res2, atom_name_border1_res2]: # if D is not B
                    dihedral_atom_names = [atom_name_center_res1, atom_name_center_res2, atom_name_border1_res2, atom_name_border2_res2]
                    dihedrals_names.append(dihedral_atom_names)

    # DIHEDRALS: IMPROPER
    imp_dihedrals_names = []
    # (central atom is the THIRD here, Amber-style) per connection, retrieve: Y and Z in Y--^Z-A-|-B,  C and D in A-|-^C-B--D if resp. A or B has only 3 bonding partners
    for i, connection in enumerate(connections):
        atom_name_center_res1 = connection[0]
        atom_name_center_res2 = connection[1]
        atom_idx_center_res1 = namestoidx_res1[atom_name_center_res1] - 1 # making 0-based
        atom_idx_center_res2 = namestoidx_res2[atom_name_center_res2] - 1
        atom_name_center_res2 = "+" + atom_name_center_res2  # res2 gets a + in front

        # Determine bonding partners
        res1_border1_idx = np.where(connectivity_res1[atom_idx_center_res1, :] == 1)[0]
        res2_border1_idx =  np.where(connectivity_res2[atom_idx_center_res2, :] == 1)[0]

        if len(res1_border1_idx) == 2: # Y--^Z--A-|-B and B as third connection
            atom_name_border_1 = atom_names_res1[res1_border1_idx[0]] 
            atom_name_border_2 = atom_names_res1[res1_border1_idx[1]] 
            imp_dihedral_atom_names = [atom_name_border_1, atom_name_border_2, atom_name_center_res1, atom_name_center_res2] # The third atom is the center
            imp_dihedrals_names.append(imp_dihedral_atom_names)

        if len(res2_border1_idx) == 2: # A-|--^C-B--D and A as third connection
            atom_name_border_1 = "+" + atom_names_res2[res2_border1_idx[0]] # res2 atoms get a + in front 
            atom_name_border_2 = "+" + atom_names_res2[res2_border1_idx[1]] 
            imp_dihedral_atom_names = [atom_name_border_1, atom_name_border_2, atom_name_center_res2, atom_name_center_res1] # The third atom is the center
            imp_dihedrals_names.append(imp_dihedral_atom_names)

    # PAIRS
    pairs_names = []
    pairs_bond_distance = list(range(pair_thresholds[0], pair_thresholds[1]+1)) # In this force field, the number of bond lengths between two atoms in a pair

    for connection in connections:
        atom_idx_center_res1 = namestoidx_res1[connection[0]] - 1 # idx-1 from 1-based to array coords to 0-based for arrays
        atom_idx_center_res2 = namestoidx_res2[connection[1]] - 1 

        connect_center_res1 = connectivity_res1[atom_idx_center_res1, :]
        connect_center_res2 = connectivity_res2[atom_idx_center_res2, :]

        for atomidx_res1 in range(len(connect_center_res1)):
            for atomidx_res2 in range(len(connect_center_res2)):
                bonds_apart = connect_center_res1[atomidx_res1] + connect_center_res2[atomidx_res2] + 1 # bonds between atom1 and atom2
                if bonds_apart in pairs_bond_distance: # if the number of bonds between is in the range of allowed ones
                    atom_name_res1 = atom_names_res1[atomidx_res1]
                    atom_name_res2 = "+" + atom_names_res2[atomidx_res2]
                    pairs_names.append([atom_name_res1, atom_name_res2])

    return atoms_names, bonds_names, angles_names, dihedrals_names, imp_dihedrals_names, pairs_names, edges_names

def get_atom_types_from_names(bonds_names: list, angles_names: list, dihedrals_names: list, imp_dihedrals_names: list, atom_types_res1: list, atom_types_res2: list, namestoidx_res1: dict, namestoidx_res2: dict):
    """ Gets the atom types for all atoms involved in bonded interactions (bonds, angles, dihedrals) from their atom name."""
    
    # Deep copy the interaction_names and sequentially replace each name by its type
    bonds_types, angles_types, dihedrals_types, imp_dihedrals_types = copy.deepcopy(bonds_names), copy.deepcopy(angles_names), copy.deepcopy(dihedrals_names), copy.deepcopy(imp_dihedrals_names)
    inttypes = [bonds_types, angles_types, dihedrals_types, imp_dihedrals_types]

    for inttype in inttypes: # For bonds, angles, dihedrals
       for term in inttype: # For every "term", eg one bond, one angle, or one dihedral
            for j, name2type in enumerate(term): # For each atom in that term
                if name2type[0] == "+":  # res2 starts with a plus. Don't use this + in dict
                    name2type = atom_types_res2[namestoidx_res2[name2type[1:]] - 1] # Replace the name by the type. -1 since this index is 1-based.
                    term[j] = name2type
                else:
                    name2type = atom_types_res1[namestoidx_res1[name2type] - 1]
                    term[j] = name2type

    return bonds_types, angles_types, dihedrals_types, imp_dihedrals_types

def retrieve_bonded_parameters(bonds_types: list, angles_types: list, dihedrals_types: list, imp_dihedrals_types: list, parent_ff_file: str, ff_type: str):
    """ Looks into the parent force field and retrieves the parameters for bonds, angles, dihedrals for the corresponding atom types. """

    bonds_parameters, angles_parameters, dihedrals_parameters, imp_dihedrals_parameters = [], [], [], []
    extension_parent_ff = str(parent_ff_file)[-4:]

    inttypes = [bonds_types, angles_types, dihedrals_types, imp_dihedrals_types]    
    strnames_inttypes = [[], [], [], [], []] # The formatted name strings of each interaction type to be sought. Eg 'c3-c3-cc-cd'. For GROMACS this is only used for printing purposes, parameters are used. 
    inttypes_names = ["bond", "angle", "proper dihedral", "improper dihedral"]
    ff_lines = [[], [], [], []] # The relevant force field lines for each interaction. A list of lists since multiple dihedrals are possible, but also list of lists for bonds and angles, yet warning when more.

    with open(parent_ff_file, 'r') as f: # Store the file in a list
        parent_ff_lines = f.readlines()

    if ff_type == "amber": # AMBER style of force field

        for i, inttype in enumerate(inttypes): 
            # Reword in AMBER style: eg c3-c3-cc-cd
            for atoms in inttype:
                types_line = ""
                for atom in atoms:
                    if len(atom) == 1: # If an atom is one character, append a space
                        atom = atom + " "
                    types_line = types_line + "-" + atom
                types_line = types_line[1:] # Remove the leading " - "
                strnames_inttypes[i].append(types_line)

            # Find the line in the force field file and append it to a list
            for s, strname in enumerate(strnames_inttypes[i]):
                ff_lines_element = []

                # AMBER atom types seem to be left aligned. If that is not the case, uncomment this to also try other possibilities.
                #if " " in strname: # one-letter atom type in AMBER can have different alignment
                #    allspaces_names_set = {strname} # Set not to have any doubles even while generating
                #    whitespace_indices = [ i for i, character in enumerate(strname) if character == " "]

                #    for whitespace_idx in whitespace_indices:
                #        # Iterate over each whitespace
                #        for space_name in list(allspaces_names_set): # Iterate over a copy as to add onto for the second space etc (exponential 2 - 4 - 8 etc) 
                #            chars = list(space_name)

                #            if whitespace_idx == len(chars)-1:
                #                whitespace_idx_shift = -2 # last character of line, so check if two before is indeed -
                #            else:
                #                whitespace_idx_shift = 1 # check if next character is -
                #                
                #            if chars[whitespace_idx + whitespace_idx_shift]  == "-": #X -YZ so also try " X-YZ"
                #                chars[whitespace_idx] = chars[whitespace_idx-1] 
                #                chars[whitespace_idx - 1] = " "
                #            else: #" X-YZ" so also try "X -YZ"
                #                chars[whitespace_idx] = chars[whitespace_idx+1] 
                #                chars[whitespace_idx + 1] = " "

                #            allspaces_names_set.add("".join(chars))
 
                #else:
                #    allspaces_names_set = {strname}
                allspaces_names_set = {strname}
                
                if i != 3: # Bonds, angles, proper dihedrals
                    for order_name in list(allspaces_names_set): # Rearrange the order as well, so that eg cs-sy and sy-cs are looked for. cs-os-sy & sy-os-cs but cs-os-sy-h1 and h1-sy-os-cs.
                        nchars_order_name = len(order_name)
                        if nchars_order_name == 5: # bond eg c3-os
                            new_order_name = order_name[3:5] + "-" + order_name[0:2]
                        elif nchars_order_name == 8: # angle eg h1-c3-os
                            new_order_name = order_name[6:8] + "-" + order_name[3:5] + "-" + order_name[0:2]
                        elif nchars_order_name == 11: # dihedral eg h1-c3-os-c3
                            new_order_name = order_name[9:11] + "-" + order_name[6:8] + "-" + order_name[3:5] + "-" + order_name[0:2]
                        else:
                            print(f"Erorr: unusual character count in AMBER-type force field string, {nchars_order_name} for {order_name}. This is a fault in the generation of possibilities, not the reading out.")
                            sys.exit(1) 

                        allspaces_names_set.add(new_order_name)

                else: # Improper dihedrals: rearrange all but third atom
                    for order_name in list(allspaces_names_set): # Rearrange the order as well, so that eg cs-sy and sy-cs are looked for. cs-os-sy & sy-os-cs but cs-os-sy-h1 and h1-sy-os-cs.
                        nchars_order_name = len(order_name)
                        if nchars_order_name != 11: # improper dihedral must have 11 like a proper
                            print(f"Erorr: unusual character count in AMBER-type force field string, {nchars_order_name} for {order_name}, for an improper dihedral (should be 11). This is a fault in the generation of possibilities, not the reading out.")
                            sys.exit(1) 
                        else: # if 11 characters
                            new_order_names = []
                            # Manually define the other combinations (the first one is already in there). Third atom is maintained as center (Amber-convention)
                            new_order_name1 = order_name[0:2] + "-" + order_name[9:11] + "-" + order_name[6:8] + "-" + order_name[3:5]
                            new_order_name2 = order_name[3:5] + "-" + order_name[0:2] + "-" + order_name[6:8] + "-" + order_name[9:11]
                            new_order_name3 = order_name[3:5] + "-" + order_name[9:11] + "-" + order_name[6:8] + "-" + order_name[0:2]
                            new_order_name4 = order_name[9:11] + "-" + order_name[3:5] + "-" + order_name[6:8] + "-" + order_name[0:2]
                            new_order_name5 = order_name[9:11] + "-" + order_name[0:2] + "-" + order_name[6:8] + "-" + order_name[3:5]

                            for new_order_name in [new_order_name1, new_order_name2, new_order_name3, new_order_name4, new_order_name5]:
                                allspaces_names_set.add(new_order_name)

                formattedname_list = list(allspaces_names_set)

                for formattedname_option in formattedname_list: # Cycle through all posibilities, instead of predicting how AMBER aligns the types. 
                    for line in parent_ff_lines: # Read through the file by reading each list item
                        if formattedname_option in line:
                            if "-" not in line[len(formattedname_option):11]: # eg c3-sy-sy-c3 has len 11. If bond c3-sy matches, it should not be the dihedral, so no extra - in those 11 characters.
                                appropriate_params = True

                                # For proper dihedrals, check whether this is not the parameters from an improper dihedral.
                                if i == 2: 
                                    nowhitespace_subline = re.sub(r"\s*-\s*", "-", line[:11]) + line[11:] # Get out the spaces before hyphens, only the first 11 characters: atom types
                                    fields = nowhitespace_subline.split()
                                    if int(float(fields[1])) != float(fields[1]):
                                        # The first field must be idivf, integer
                                        appropriate_params = False
                                    try:
                                        int(float(fields[4])) # The fourth field must be pn, integer
                                    except ValueError:
                                        appropriate_params = False
                                        
                                # Only proper dihedrals are possibly not accepted now
                                if appropriate_params: 
                                    ff_lines_element.append(line)
                                    formattedname_definite = formattedname_option # Save the exact formatted name that matched

                if not ff_lines_element: # if empty (i.e. no parameter line found in force field file), report
                    noparam_error = True # Boolean for error at the end, since possibly interfering improper dihedrals make it ok
                    try_wildcards = False # Boolean to see if wildcards were tried

                    if i == 2: # Proper dihedrals: look if involved in an improper dihedral. If not, try wildcards on edge atoms, eg X -c3-c3-X

                        # Check if matching atoms with an improper dihedral
                        atom_types_matches = []
                        for imp_dihedrals_type in imp_dihedrals_types:
                            atom_types_matches.append(sum((Counter(imp_dihedrals_type) & Counter(dihedrals_types[s])).values())) # Count number of matches between dihedral types and the imp dihedral types

                        # From 3 matching atom types on, very probable that an improper dihedral is more adequate. Thus neglect this proper one.
                        if any(atom_types_match >= 3 for atom_types_match in atom_types_matches):
                            print(f"Note: proper dihedral {formattedname_option} was not found in {parent_ff_file}, but at least three atom types are involved in an improper dihedral. No action required.")
                            ff_lines_element.append("8888 8888 8888 8888 8888 8888 8888 8888 8888 8888") # Appending junk numbers not to mess up order, but not written in ff (8888 detected when writing).
                            noparam_error = False

                        # These atoms are not involved in an improper dihedral. Only checking for proper dihedrals for these general dihedrals X ... X
                        else: 
                            try_wildcards = True
                            transferrable_dihedral_names = set()
                            for formattedname_option in formattedname_list: # For all "options", make the names with wildcards
                                transferrable_dihedral_name = "X -" + formattedname_option[3:5] + "-" + formattedname_option[6:8] + "-X " 
                                transferrable_dihedral_names.add(transferrable_dihedral_name) # a set with all possibilities


                    elif i == 3: # Improper dihedrals: check for wildcards like eg X -X -c -X
                        try_wildcards = True
                        transferrable_dihedral_names = set()
                        for formattedname_option in formattedname_list: # For all "options", make the names with wildcards
                            # Very manual and not efficient
                            name1 = "X -" + formattedname_option[3:5] + "-" + formattedname_option[6:8] + "-" + formattedname_option[9:11]   # X 2 3 4
                            name2 = formattedname_option[0:2] + "-" + formattedname_option[3:5] + "-" + formattedname_option[6:8] + "-X "    # 1 2 3 X
                            name3 = formattedname_option[0:2] + "-X -" + formattedname_option[6:8] + formattedname_option[9:11]              # 1 X 3 4
                            #name4 = "X -" + formattedname_option[3:5] + "-" + formattedname_option[6:8] + "-X "                             # X 2 3 X # This is reserved for proper dihedrals, no confusion
                            name5 = "X -X -" + formattedname_option[6:8] + "-" + formattedname_option[9:11]                                  # X X 3 4
                            name6 = formattedname_option[0:2] + "-X -" + formattedname_option[6:8] + "-X "                                   # 1 X 3 X
                            name7 = "X -X -" + formattedname_option[6:8] + "-X "                                                             # X X 3 X
                            for name in [name1, name2, name3, name5, name6, name7]:
                                transferrable_dihedral_names.add(name) # add possibility to set (no duplicates)

                    if try_wildcards: # For proper and improper dihedrals, look for wildcards matching (X ... X)                         
                        for line in parent_ff_lines:
                            for transferrable_dihedral_name in transferrable_dihedral_names:
                                if transferrable_dihedral_name in line:
                                    ff_lines_element.append(line)
                                    transferrable_dihedral_name_match = transferrable_dihedral_name # Store name for printing

                        if ff_lines_element: # Check if X ... X helped
                            print(f"Note: no parameters found for {inttypes_names[i]} {formattedname_option} in {parent_ff_file}, but the transferrable dihedral {transferrable_dihedral_name_match} was selected.")
                            noparam_error = False

                    if noparam_error:    
                        print(f"Error: no parameters found for {inttypes_names[i]} \'{formattedname_option}\' in {parent_ff_file}. Check the output file since the entry is written with 999.")
                        ff_lines_element.append("999 999 999 999 999 999 999 999 999 999") # Appending junk numbers to continue writing and not to mess up order. It IS written in the ff.

                ff_lines[i].append(ff_lines_element) # Append the list of lines for one parameter type to the list of relevant lines 

        # Retrieve parameters from the extracted lines
        # Conversion from Amber to Gromacs: https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2015-September/101118.html
        A2nm = 0.1 # Angstrom to nanometer
        kcal2kJ = 4.184
        # BONDS
        for line in ff_lines[0]:
            parameters = []
            for subline in line:
                nowhitespace_subline = re.sub(r"\s*-\s*", "-", subline[:5]) + subline[5:] # Get out the spaces before hyphens, only in the first 5 characters: atom types
                fields = nowhitespace_subline.split()
                force_constant = float(fields[1]) / A2nm**2 * 2 * kcal2kJ
                bond_length = float(fields[2]) * A2nm
                parameters.append([bond_length, force_constant])
            bonds_parameters.append(parameters)
        
        # ANGLES
        for line in ff_lines[1]:
            parameters = []
            for subline in line:
                nowhitespace_subline = re.sub(r"\s*-\s*", "-", subline[:8]) + subline[8:] # Get out the spaces before hyphens, only the first 8 characters: atom types
                fields = nowhitespace_subline.split()
                theta = float(fields[2])
                cth = float(fields[1]) * 2 * kcal2kJ
                parameters.append([theta, cth])
            angles_parameters.append(parameters)
                
        # DIHEDRALS: PROPER AND IMPROPER
        for d, ff_lines_dih in enumerate([ff_lines[2], ff_lines[3]]): # Do both proper and improper dihedrals
            for line in ff_lines_dih:
                parameters = []
                for subline in line:
                    nowhitespace_subline = re.sub(r"\s*-\s*", "-", subline[:11]) + subline[11:] # Get out the spaces before hyphens, only the first 11 characters: atom types
                    fields = nowhitespace_subline.split()

                    if d == 0: # Proper dihedrals have idivf so extra column
                        idivf = float(fields[1])
                        column_shift = 0
                    elif d == 1: # Improper dihedrals do not
                        idivf = 1
                        column_shift = -1
                        
                    phase = float(fields[3 + column_shift])
                    kd = float(fields[2 + column_shift]) * kcal2kJ / idivf 
                    pn = int(float(fields[4 + column_shift]))
                    parameters.append([phase, kd, pn])
                if d == 0: # Proper
                    dihedrals_parameters.append(parameters)
                elif d == 1: # Improper
                    imp_dihedrals_parameters.append(parameters)


    elif ff_type == "gromacs": # GROMACS style of force field

        improper_funcs = [4] # Only func 4 supported now for improper dihedrals
        proper_funcs = [9] # Only func 9 supported now for proper dihedrals

        for i, inttype in enumerate(inttypes): 
            # Gromacs works in columns, so no strings need to be formed for the formatted interaction names
            for atoms_inttype in inttype: # Per interaction, look through the file.
                ff_lines_element = []

                formattedname = "" # Store for printing the interaction name if error
                for atom in atoms_inttype:
                    formattedname = formattedname + "-" + atom
                formattedname = formattedname[1:] # Remove leading - 

                for line in parent_ff_lines:
                    if i == 3:
                        func_idx = i+2-1
                    else:
                        func_idx = i+2

                    line_list = line.split() # Split columns
                    line_types = line_list[:func_idx] # They all have i+2 number of atoms involved, so up to :i+2 are the types (improper i+2-1)
                    types_match = False

                    # First check: see if these are parameters for the right type, if the column next to the types is the func column
                    try:
                        func = line_list[func_idx]
                    except IndexError:
                        func = None # Empty if not a line with func there, not a good parameter line

                    if func and func.isdigit(): # If not, this is an atom type and eg angle parameters instead of bond
                        func = int(func)

                        # Basic requirement for all: all atom types match
                        if Counter(atoms_inttype) == Counter(line_types): 

                            if i == 0: # Bonds: no order necessary, if match it is ok
                                if func == 1: # Only type 1 bonds in this implementation
                                    types_match = True
                                 
                            elif i == 1: # Angles: middle atom needs to match too
                                if func == 1: # Only type 1 angles in this implementation
                                    if line_list[1] == atoms_inttype[1]:
                                        types_match = True
    
                            elif i == 2: # Proper dihedrals: middle atoms need to match and with the correct bordering types 
                                if func in proper_funcs: # Only type 9 proper dihedrals in this implementation 
                                    if Counter(atoms_inttype[1:3]) == Counter(line_list[1:3]): # If two middle ones match
                                        if atoms_inttype[1] == line_list[1]: # If first middle one is on the same place in ff
                                            if atoms_inttype[0] == line_list[0]: # If first border is also the same, then all match
                                                types_match = True
                                        else: # First middle one is second middle in ff
                                            if atoms_inttype[0] == line_list[3]: # If first border is last in ff, all match
                                                types_match = True
    
                            elif i == 3: # Improper dihedrals: matching AMBER-style, putting the center atom as the third X X C X
                                if func in improper_funcs: # Only type 4 improper dihedrals in this implementation 
                                    border_atoms_types = [atoms_inttype[col_nr] for col_nr in [0, 1, 3]]
                                    border_atoms_ff = [line_list[col_nr] for col_nr in [0, 1, 3]]
                                    if Counter(border_atoms_types) == Counter(border_atoms_ff):
                                        types_match = True

                            if types_match:
                                ff_lines_element.append(line)

                if not ff_lines_element: # if empty (i.e. no parameter line found in force field file), report
                    noparam_error = True # Boolean for error at the end, since possibly interfering improper dihedrals make it ok
                    try_wildcards = False # Boolean to see if wildcards were tried

                    if i == 2: # Proper dihedrals: look if involved in an improper dihedral. If not, try wildcards on edge atoms, eg X c3 c3 X

                        # Check if matching atoms with an improper dihedral
                        atom_types_matches = []
                        for imp_dihedrals_type in imp_dihedrals_types:
                            atom_types_matches.append(sum((Counter(imp_dihedrals_type) & Counter(atoms_inttype)).values())) # Count number of matches between dihedral types and the imp dihedral types

                        # From 3 matching atom types on, very probable that an improper dihedral is more adequate. Thus neglect this proper one.
                        if any(atom_types_match >= 3 for atom_types_match in atom_types_matches):
                            print(f"Note: proper dihedral {formattedname} was not found in {parent_ff_file}, but at least three atom types are involved in an improper dihedral. No action required.")
                            ff_lines_element.append("8888 8888 8888 8888 8888 8888 8888 8888 8888 8888") # Appending junk numbers not to mess up order, but not written in ff (8888 detected when writing).
                            noparam_error = False

                        # These atoms are not involved in an improper dihedral. Only checking for proper dihedrals for these general dihedrals X ... X
                        else: 
                            try_wildcards = True
                            # Transferrable dihedral names are not premade now: logic with X and X is made on the spot further

                    elif i == 3: # Improper dihedrals: check for wildcards like eg X -X -c -X
                        try_wildcards = True
                        # Transferrable dihedral names are not premade now: logic with X X C X is made on the spot further

                    if try_wildcards: # For proper and improper dihedrals, look for wildcards matching (X ... X)                         

                        for line in parent_ff_lines: # Look through lines again
                            line_list = line.split() # Split columns
                            line_types = line_list[:4] # All dihedrals have 4 atoms involved, so up to :4 are the types
                            types_match = False
        
                            # First check: see if these are parameters for the right type, if the column next to the types is the func column
                            try:
                                func = line_list[4]
                            except IndexError:
                                func = None

                            if func and func.isdigit(): # If not, these are not dihedral but eg angle parameters
                                func = int(func)

                                # Requirement toned down: not all atom types have to match now. Starting from one.
                                atom_matches = sum((Counter(line_types) & Counter(atoms_inttype)).values())
                                if atom_matches >= 1 and "X" in line_types: # If at least one match and a line with wildcards
                                    
                                    # Proper dihedrals
                                    if i == 2:
                                        same_middle_atoms = Counter(atoms_inttype[1:3]) == Counter(line_types[1:3])
                                        outer_X = line_types[0] == "X" and line_types[3] == "X"
                                        is_proper = func in proper_funcs
                                        if is_proper and same_middle_atoms and outer_X:
                                            types_match = True

                                    # Improper dihedrals
                                    elif i == 3:
                                        same_center = atoms_inttype[2] == line_types[2] # Matching AMBER-style, center needs to be third atom
                                        # no need to check for X A B X, since func 4 makes sure it is not proper dihedral
                                        number_X = Counter(line_types)["X"]
                                        number_match = sum((Counter(atoms_inttype) & Counter(line_types)).values())
                                        no_unmatching_types = number_X + number_match == 4

                                        is_improper = func in improper_funcs
                                        if is_improper and same_center and no_unmatching_types:
                                            types_match = True

                                    if types_match:
                                        ff_lines_element.append(line)
                                        transferrable_dihedral_name_match = f"{line_types[0]} - {line_types[1]} - {line_types[2]} - {line_types[3]}" # Store name for printing

                        if ff_lines_element: # Check if wildcards helped
                            print(f"Note: no parameters found for {inttypes_names[i]} {formattedname} in {parent_ff_file}, but the transferrable dihedral {transferrable_dihedral_name_match} was selected.")
                            noparam_error = False

                    if noparam_error:    
                        print(f"Error: no parameters found for {inttypes_names[i]} \'{formattedname}\' in {parent_ff_file}. Check the output file since the entry is written with 999.")
                        ff_lines_element.append("999 999 999 999 999 999 999 999 999 999") # Appending junk numbers to continue writing and not to mess up order. It IS written in the ff.

                ff_lines[i].append(ff_lines_element) # Append the list of lines for one parameter type to the list of relevant lines  

        # Retrieve parameters from the extracted lines
        # Made for func 1 for bonds, 1 for angles, 9 for proper dihedrals, 4 for improper dihedrals
        # BONDS
        for line in ff_lines[0]:
            parameters = []
            for subline in line:
                line_list = subline.split()
                bond_length = float(line_list[3])
                force_constant = float(line_list[4])
                parameters.append([bond_length, force_constant])
            bonds_parameters.append(parameters)
        
        # ANGLES
        for line in ff_lines[1]:
            parameters = []
            for subline in line:
                line_list = subline.split()
                theta = float(line_list[4])
                cth = float(line_list[5])
                parameters.append([theta, cth])
            angles_parameters.append(parameters)
                
        # DIHEDRALS: PROPER AND IMPROPER
        for d, ff_lines_dih in enumerate([ff_lines[2], ff_lines[3]]): # Do both proper and improper dihedrals
            for line in ff_lines_dih:
                parameters = []
                for subline in line:
                    line_list = subline.split()
                    phase = float(line_list[5])
                    kd = float(line_list[6])
                    pn = int(float(line_list[7]))
                    func = int(line_list[4]) # not used here now
                    parameters.append([phase, kd, pn])
                if d == 0: # Proper
                    dihedrals_parameters.append(parameters)
                elif d == 1: # Improper
                    imp_dihedrals_parameters.append(parameters)

    return bonds_parameters, angles_parameters, dihedrals_parameters, imp_dihedrals_parameters 


def write_atoms_link(atoms_names: list, itp_res1_file: str, itp_res2_file: str, res1_name: str, res2_name: str, link_filename: str):
    """ Write the beginning of the link.ff file and the atoms section. """

    with open(link_filename, "w") as f:
        f.write(f"; Link file created for Polyply using {itp_res1_file} and {itp_res2_file}, by {os.getlogin()} at {datetime.datetime.now()}\n\n") # Opening comment
        f.write("[ link ]\n")
        f.write(write_resnames(res1_name, res2_name)+"\n")
        f.write("[ atoms ]\n")
        for atom in atoms_names:
            f.write(f"{atom:>4s} {{ }}\n") # Curly brackets after

def write_bonded_terms_link(bonds_names: list, angles_names: list, dihedrals_names: list, imp_dihedrals_names: list, bonds_types: list, angles_types: list, dihedrals_types: list, imp_dihedrals_types: list, bonds_parameters: list, angles_parameters: list, dihedrals_parameters: list, imp_dihedrals_parameters: list, res1_name: str, res2_name: str, link_filename: str):
    """ Write the bonds, angles, and dihedrals sections in the link.ff file. """

    with open(link_filename, "a") as f:

        # BONDS
        f.write("\n[ link ]\n")
        f.write(write_resnames(res1_name, res2_name)+"\n")
        f.write("[ bonds ]\n")
        for i, bond_name in enumerate(bonds_names):
            if len(bonds_parameters[i]) > 1:
                print(f"Warning: More than one parameter set for bond name {bond_name[0]}-{bond_name[1]}, i.e. types {bonds_types[i][0]}-{bonds_types[i][1]}. Make sure there is no double counting.")

            for j, bonds_parameter in enumerate(bonds_parameters[i]):
                bond_line = f"{bond_name[0]:>4s}{bond_name[1]:>6s}{1:>4d}{bonds_parameter[0]:>14.4e}{bonds_parameter[1]:>14.4e} ; {bonds_types[i][0]:>4s} - {bonds_types[i][1]:>4s}\n"
                f.write(bond_line)
            
        # ANGLES
        f.write("\n[ link ]\n")
        f.write(write_resnames(res1_name, res2_name)+"\n")
        f.write("[ angles ]\n")
        for i, angle_name in enumerate(angles_names):
            if len(angles_parameters[i]) > 1:
                print(f"Warning: More than one parameter set for angle name {angle_name[0]}-{angle_name[1]}-{angle_name[2]}, i.e. types {angles_types[i][0]}-{angles_types[i][1]}-{angles_types[i][2]}. Make sure there is no double counting.")
            for j, angles_parameter in enumerate(angles_parameters[i]):
                angle_line = f"{angle_name[0]:>4s}{angle_name[1]:>6s}{angle_name[2]:>6s}{1:>4d}{angles_parameter[0]:>8.2f}{angles_parameter[1]:>14.4e} ; {angles_types[i][0]:>4s} - {angles_types[i][1]:>4s} - {angles_types[i][2]:>4s}\n"
                f.write(angle_line)

        # DIHEDRALS: PROPER
        f.write("\n[ link ]\n")
        f.write(write_resnames(res1_name, res2_name)+"\n")
        f.write("[ dihedrals ]\n")
        for i, dihedral_name in enumerate(dihedrals_names):
            for j, dihedrals_parameter in enumerate(dihedrals_parameters[i]): # Multiple functions for one dihedral
                if not dihedrals_parameter[2] == 8888: # If the phase is 8888, this was the discarded proper dihedral due to improper interference
                    dihedral_line = f"{dihedral_name[0]:>4s}{dihedral_name[1]:>6s}{dihedral_name[2]:>6s}{dihedral_name[3]:>6s}{9:>4d}{dihedrals_parameter[0]:>8.2f}{dihedrals_parameter[1]:>10.5f}{dihedrals_parameter[2]:>4d} {{\"version\": \"{j+1}\"}} ; {dihedrals_types[i][0]:>4s} - {dihedrals_types[i][1]:>4s} - {dihedrals_types[i][2]:>4s} - {dihedrals_types[i][3]:>4s}\n"
                    f.write(dihedral_line)

        # DIHEDRALS: IMPROPER
        # (func 4 in Gromacs)
        for i, imp_dihedral_name in enumerate(imp_dihedrals_names):
            if len(imp_dihedrals_parameters[i]) > 1:
                print(f"Warning: More than one parameter set for IMPROPER dihedral name {imp_dihedral_name[0]}-{imp_dihedral_name[1]}-{imp_dihedral_name[2]}-{imp_dihedral_name[3]}, i.e. types {imp_dihedrals_types[i][0]}-{imp_dihedrals_types[i][1]}-{imp_dihedrals_types[i][2]}-{imp_dihedrals_types[i][3]}. This is unusual for IMPROPER dihedrals. Make sure there is no double counting.")
            for j, imp_dihedrals_parameter in enumerate(imp_dihedrals_parameters[i]): 
                imp_dihedral_line = f"{imp_dihedral_name[0]:>4s}{imp_dihedral_name[1]:>6s}{imp_dihedral_name[2]:>6s}{imp_dihedral_name[3]:>6s}{4:>4d}{imp_dihedrals_parameter[0]:>8.2f}{imp_dihedrals_parameter[1]:>10.5f}{imp_dihedrals_parameter[2]:>4d} {{\"version\": \"{j+1}\"}} ; {imp_dihedrals_types[i][0]:>4s} - {imp_dihedrals_types[i][1]:>4s} - {imp_dihedrals_types[i][2]:>4s} - {imp_dihedrals_types[i][3]:>4s}\n"
                f.write(imp_dihedral_line)


def write_pairs_edges_link(pairs_names: list, edges_names: list, res1_name: str, res2_name: str, link_filename: str):
    """ Writes the pairs and edges sections to the link.ff file. No parameters needed."""

    # Pairs
    with open(link_filename, "a") as f:
        f.write("\n[ link ]\n")
        f.write(write_resnames(res1_name, res2_name)+"\n")
        f.write("[ pairs ]\n")
        for pair in pairs_names:
            f.write(f"{pair[0]:>4s}{pair[1]:>6s}{1:>4d}\n")
    
    # Edges
    with open(link_filename, "a") as f:
        f.write("\n[ link ]\n")
        f.write(write_resnames(res1_name, res2_name)+"\n")
        f.write("[ edges ]\n")
        for edge in edges_names:
            f.write(f"{edge[0]:>4s}{edge[1]:>6s}\n")

def write_resnames(res1_name: str, res2_name: str) -> str:
    """ Writes the resname statement under the link section (by passing on the string).
        If res2_name is same or if it equals res1_name, only 'resname "XXXX"' is written.
        If not, 'resname "XXX1|XXX2"' is written."""
    
    if res2_name == "same" or res2_name == res1_name:
        resname_str = f"resname \"{res1_name}\""
    else:
        resname_str = f"resname \"{res1_name}|{res2_name}\""

    return resname_str

def link_filename(res1_name: str, res2_name: str) -> str:
    """ Passes on the string of the file name for the new link file. """

    if res2_name == "same" or res2_name == res1_name:
        resname_str = res1_name + "_" + res1_name
    else:
        resname_str = f"{res1_name}_{res2_name}"
    link_filename = "link_" + resname_str + ".ff"

    return link_filename


if __name__ == "__main__":
    # main: create link.ff file
    itp_res1_file, itp_res2_file, res1_name, res2_name, connection_atomnames_file, parent_ff_file = load_input_command()
    link_filename = link_filename(res1_name, res2_name)

    # Determine type of parent force field file: AMBER (.dat) or GROMACS (.itp)
    ff_type = ff_type(parent_ff_file)

    # Determine connectivity within residues
    atoms_block_int_res1, bonds_block_int_res1, pairs_block_int_res1, namestoidx_res1, atom_types_res1, atom_names_res1 = read_itp_atomsbondspairs(itp_res1_file)
    connectivity_res1 =  determine_connectivity(itp_res1_file, atoms_block_int_res1, bonds_block_int_res1)

    if itp_res1_file == itp_res2_file: # if same, don't go through it again
        connectivity_res2, atoms_block_int_res2, bonds_block_int_res2, pairs_block_int_res2, namestoidx_res2, atom_types_res2, atom_names_res2 = connectivity_res1, atoms_block_int_res1, bonds_block_int_res1, pairs_block_int_res1, namestoidx_res1, atom_types_res1, atom_names_res1
    else:
        atoms_block_int_res2, bonds_block_int_res2, pairs_block_int_res2, namestoidx_res2, atom_types_res2, atom_names_res2 = read_itp_atomsbondspairs(itp_res2_file)
        connectivity_res2 =  determine_connectivity(itp_res2_file, atoms_block_int_res2, bonds_block_int_res2)

    # In this force field, determine for which bond count they are added to pairs. Done for itp1.
    pair_thresholds = determine_pair_thresholds(pairs_block_int_res1, connectivity_res1)

    # Read in connected atom names and determine atom names for atoms, bonds, angles, dihedrals, pairs, and edges between residues
    connections = read_connection_file(connection_atomnames_file, atom_names_res1, atom_names_res2)
    atoms_names, bonds_names, angles_names, dihedrals_names, imp_dihedrals_names, pairs_names, edges_names = link_atom_names(connections, namestoidx_res1, namestoidx_res2, connectivity_res1, connectivity_res2, atom_names_res1, atom_names_res2, bonds_block_int_res1, bonds_block_int_res2, pair_thresholds)

    # For bonds, angles, and dihedrals, retrieve parameters from the parent force field
    bonds_types, angles_types, dihedrals_types, imp_dihedrals_types = get_atom_types_from_names(bonds_names, angles_names, dihedrals_names, imp_dihedrals_names, atom_types_res1, atom_types_res2, namestoidx_res1, namestoidx_res2)
    bonds_parameters, angles_parameters, dihedrals_parameters, imp_dihedrals_parameters = retrieve_bonded_parameters(bonds_types, angles_types, dihedrals_types, imp_dihedrals_types, parent_ff_file, ff_type)
    
    # Write the link file
    write_atoms_link(atoms_names, itp_res1_file, itp_res2_file, res1_name, res2_name, link_filename)
    write_bonded_terms_link(bonds_names, angles_names, dihedrals_names, imp_dihedrals_names, bonds_types, angles_types, dihedrals_types, imp_dihedrals_types, bonds_parameters, angles_parameters, dihedrals_parameters, imp_dihedrals_parameters, res1_name, res2_name, link_filename)
    write_pairs_edges_link(pairs_names, edges_names, res1_name, res2_name, link_filename)
