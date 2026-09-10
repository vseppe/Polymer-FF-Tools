#!/usr/bin/env python3
# title             : remove_caps_ff.py
# description       : Removes specified atom indices from .ff file to remove those caps. 
# use               : python remove_caps_ff.py $itp_res.ff $capping_idx.txt 
# input file        : $itp_res.ff
# output file       : itp_res_uncapped.ff

import numpy as np
import sys
from pathlib import Path
import datetime
import os

def load_input_command():
    """ Process the arguments given in the command line. """
    # Inputs
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <itp_res.ff> <capping_idx.txt>")
        sys.exit(1)
    
    itp_res_file = Path(sys.argv[1])
    capping_idx_file = Path(sys.argv[2])
    
    if not itp_res_file.is_file():
        print(f"Error: itp_res.ff input file '{itp_res_file}' does not exist.")
        sys.exit(1)
    
    if not capping_idx_file.is_file():
        print(f"Error: capping_idx.txt input file '{capping_idx_file}' does not exist.")
        sys.exit(1)
    return itp_res_file, capping_idx_file

# Functions
def section_block(itp_lines: list, sections: dict, section_name: str) -> list:
    """ Returns all the lines of the given itp file for the given section. """

    beginning_line = sections[section_name]
    next_section = min((v for v in sections.values() if v > beginning_line), default=len(itp_lines)) # Next section starts with the first line number that is higher
    block = itp_lines[beginning_line:next_section]
    return block

def find_nonint_lines(block: list) -> list:
    """ Returns line numbers of lines that do not start with an integer. E.g. whitespace, ";" comment. """
    
    nonint_lines = []

    for line_i, line in enumerate(block): # cycle through all lines
        line = line.lstrip()
        if not line: # empty lines
            nonint_lines.append(line_i)
        elif line[0] == ";": # comment line
            nonint_lines.append(line_i)
        elif line[0] == "[": # section line
            nonint_lines.append(line_i)
        else:
            try:
                int(line.split()[0]) # all nonintegers will fail at this
            except ValueError:
                nonint_lines.append(line_i)

    return nonint_lines

def read_itp_res(itp_res_file: str, section_names: list):
    """ Reads the itp_res.ff file, notes down the line numbers of the sections and returns them. """

    # Read itp_res into a list of lines
    with open(itp_res_file, "r") as itp_f:
        itp_lines = itp_f.readlines() 

    # Read and note down line numbers of sections
    section_linenrs = []
    for section_i, section in enumerate(section_names): # Cycle through file for every section
        section = f"[ {section} ]"
        for line_i, line in enumerate(itp_lines):
            if section in line and line.lstrip()[0] != ";": # If this is the section and not commented out by ";"
                section_linenrs.append(line_i) 

        try: # Stop if section not found
            section_linenrs[section_i]
        except IndexError:
            print(f"No section '{ section }' found in '{itp_res_file}'.")
            print("Make sure it is not commented out. If needed specify it, but empty.")
            print("Exiting.")
            sys.exit(1)

    # Dictionary of sections
    sections = dict(zip(section_names, section_linenrs))

    return itp_lines, section_linenrs, sections

def read_capping_idx(capping_idx_file: str) -> list:
    """ Reads the file with atom indices of capping groups, each separated by a newline. Returns list of indices. 1-based indices, error if 0 detected."""

    capping_idx = []

    with open(capping_idx_file, "r") as f:
        capping_lines = f.readlines()

    # Cycle through lines and get index
    for line in capping_lines:
        if not line.lstrip():
            # empty line, do not do anything
            pass
        else:
            line_firstcol = line.split()[0]
            if line_firstcol[0] == "#":
                # comment, do not read
                pass
            else:
                try:
                    int(line_firstcol) # only integers can be indices
                    capping_idx.append(int(line_firstcol)) 
                except ValueError: # nonintegers give error
                    print(f"Error. Noninteger was supplied as index in capping group file '{capping_idx_file}', namely '{line_firstcol}'.")
                    print("Exiting.")
                    sys.exit(1)

    # Errors if not good input
    if 0 in capping_idx:
        print(f"Error. 0 detected in capping group file '{capping_idx_file}'. The atoms are 1-indexed, not 0-indexed.")
        print("Exiting.")
        sys.exit(1)
    elif any(idx < 0 for idx in capping_idx):
        print(f"Error. Negative integer detected in capping group file '{capping_idx_file}'. Reminder: the atoms are 1-indexed, not 0-indexed.")
        print("Exiting.")
        sys.exit(1)
 
    capping_idx = sorted(set(capping_idx)) # sort and remove duplicates

    return capping_idx


def renumber_atoms(itp_lines: list, sections: dict, capping_idx: list):
    """ Renumbers the atoms to remove the capping groups, given by index in a list. Returns new atoms block, new indices in same vector size (e.g. 0 0 0 1 2 3 0 0 0). """

    atoms_block = section_block(itp_lines, sections, "atoms")
    nonint_lines = find_nonint_lines(atoms_block) # lines that need to be added after
    new_nonint_lines = np.array(nonint_lines)

    # Retrieve number of atoms to make vector of new indices
    n_atoms = len(atoms_block) - len(nonint_lines)
    new_atomidx = np.zeros(n_atoms, dtype=int)

    # Cycle through lines in atoms block and check whether capping group. Append new index if not and change atom number in line.
    new_idx = 1
    atoms_block_int = [ x for i, x in enumerate(atoms_block) if i not in nonint_lines ] # indices of list match atom indices-1
    new_atoms_block_int = []
    for line_i, line in enumerate(atoms_block_int): 
        line_list = line.split()
        if int(line_list[0]) in capping_idx: # capping group
            # subtract 1 from nonint lines coming after, so that their order is preserved
            original_line_i = atoms_block.index(line) # get the line number in the whole block to compare
            new_nonint_lines[np.array(nonint_lines) > original_line_i] -= 1 # refer to old order and apply to new indices

        else: # no capping group: append new idx to vector and change atom number
            new_atomidx[line_i] = new_idx
            # Change atom number in line. Formatted string right aligned >
            line_list[0], line_list[5] = new_idx, new_idx
            new_line = write_itp_entry(line_list, "atoms")
            new_atoms_block_int.append(new_line)

            new_idx += 1

    # Print new combined block, also with comments
    new_atoms_block = list(new_atoms_block_int)
    for new_nonint_linenr, nonint_linenr in zip(new_nonint_lines, nonint_lines):
        new_atoms_block[new_nonint_linenr:new_nonint_linenr] = [atoms_block[nonint_linenr]] # insert comments at idx. Because adding elements, do not do it backwards.

    return new_atomidx, new_atoms_block

def renumber_interactions(itp_lines: list, sections: dict, new_atomidx: np.array):
    """ Replace the old atom indices by the new ones, in all sections other than [ atoms ].
        Uses renumber_specific_interaction()."""

    new_bonds = renumber_specific_interaction(itp_lines, sections, "bonds", new_atomidx)
    new_pairs = renumber_specific_interaction(itp_lines, sections, "pairs", new_atomidx)
    new_angles = renumber_specific_interaction(itp_lines, sections, "angles", new_atomidx)
    new_dihedrals = renumber_specific_interaction(itp_lines, sections, "dihedrals", new_atomidx)

    return new_bonds, new_pairs, new_angles, new_dihedrals

def renumber_specific_interaction(itp_lines: list, sections: dict, interaction_type: str, new_atomidx: np.array):
    """ 'Worker' function for renumber_interactions. Renumbers the atoms for each specific interaction."""    

    # Define number of atoms involved
    if interaction_type == "bonds" or interaction_type == "pairs":
        n_atom_nrs = 2
    elif interaction_type == "angles":
        n_atom_nrs = 3
    elif interaction_type == "dihedrals": 
        n_atom_nrs = 4

    # Define the block of the interaction
    old_block = section_block(itp_lines, sections, interaction_type)
    nonint_lines = find_nonint_lines(old_block) # lines that need to be added after
    new_nonint_lines = np.array(nonint_lines)

    # Cycle through each line. Check whether atoms in this line were capping (zero in new_atomidx). If yes, delete line. If not, replace numbers.
    old_block_int = [ x for i, x in enumerate(old_block) if i not in nonint_lines ] # the block with only the actual interaction
    new_block_int = []

    for line_i, line in enumerate(old_block_int):
        line_list = line.split()
        atom_nrs = list(map(int, line_list[:n_atom_nrs])) # get all the atom indices of this interaction. Make integers
        new_atom_nrs = [ new_atomidx[atom_nr - 1] for atom_nr in atom_nrs ] # and the new ones
        if any(new_atom_nr == 0 for new_atom_nr in new_atom_nrs): # discard this entry: at least one atom was capping atom
            # subtract 1 from nonint lines coming after, so that their order is preserved
            original_line_i = old_block.index(line) # get the line number in the whole block to compare
            new_nonint_lines[np.array(nonint_lines) > original_line_i] -= 1 # refer to old order and apply to new indices
        else: # renumber if all still present
            for atom_nr_i, new_atom_nr in enumerate(new_atom_nrs):
                line_list[atom_nr_i] = new_atom_nr
            new_line = write_itp_entry(line_list, interaction_type)
            new_block_int.append(new_line)

    # Print new combined block, also with comments
    new_block = list(new_block_int)
    for new_nonint_linenr, nonint_linenr in zip(new_nonint_lines, nonint_lines):
        new_block[new_nonint_linenr:new_nonint_linenr] = [old_block[nonint_linenr]] # insert comments at idx. Because adding elements, do not do it backwards.

    return new_block

def write_itp_entry(line_list: list, interaction_type: str):
    """ Writes a formatted string for the itp entry line, given a list version of the correct line. Adds correct spacing. If extra text at the end, appended."""

    if interaction_type == "atoms":
        formatted_itp_entry = f"{line_list[0]:>5d}{line_list[1]:>5s}{int(line_list[2]):>6d}{line_list[3]:>6s}{line_list[4]:>6s}{line_list[5]:>5d}{float(line_list[6]):>13.6f}{float(line_list[7]):>13.5f} {" ".join(str(s) for s in line_list[8:])}\n"
    elif interaction_type == "bonds":
        formatted_itp_entry = f"{line_list[0]:>5d}{line_list[1]:>7d}{int(line_list[2]):>4d}{float(line_list[3]):>14.4e}{float(line_list[4]):>14.4e}{line_list[5]:>2s}{line_list[6]:>7s}{line_list[7]:^3s}{line_list[8]:<6} {" ".join(str(s) for s in line_list[9:])}\n"
    elif interaction_type == "pairs":
        formatted_itp_entry = f"{line_list[0]:>5d}{line_list[1]:>7d}{int(line_list[2]):>7d}{line_list[3]:>2s}{line_list[4]:>7s}{line_list[5]:^3s}{line_list[6]:<6} {" ".join(str(s) for s in line_list[7:])}\n"
    elif interaction_type == "angles":
        formatted_itp_entry = f"{line_list[0]:>5d}{line_list[1]:>7d}{line_list[2]:>7d}{int(line_list[3]):>7d}{float(line_list[4]):>14.4e}{float(line_list[5]):>14.4e}{line_list[6]:>2s}{line_list[7]:>7s}{line_list[8]:^3s}{line_list[9]:<7s}{line_list[10]:<2s}{line_list[11]:<6s} {" ".join(str(s) for s in line_list[12:])}\n"
    elif interaction_type == "dihedrals": 
        formatted_itp_entry = f"{line_list[0]:>5d}{line_list[1]:>7d}{line_list[2]:>7d}{line_list[3]:>7d}{int(line_list[4]):>7d}{float(line_list[5]):>9.2f}{float(line_list[6]):>10.5f}{int(line_list[7]):>4d}{line_list[8]:>2s}{line_list[9]:>8s}{line_list[10]:>7s}{line_list[11]:>7s}{line_list[12]:>6s} {" ".join(str(s) for s in line_list[13:])}\n"

    return formatted_itp_entry

def write_itp_res_uncapped(itp_lines: list, section_linenrs: list, new_atoms: list, new_bonds: list, new_pairs: list, new_angles: list, new_dihedrals: list, itp_res_file: str):
    """ Using the new atom numbering, write a new itp file for the uncapped residue. """

    # New file name
    itp_res_uncapped_file = str(itp_res_file)[:-3] + "_uncapped.ff"

    # Combining sections
    notice = [f"; {itp_res_uncapped_file} reworked from {itp_res_file} to delete cap atoms, by {os.getlogin()} at {datetime.datetime.now()}\n"]
    itp_beginning = itp_lines[:min(section_linenrs)]   # the comments etc from original ff file
    moleculetype = section_block(itp_lines, sections, "moleculetype")
    new_itp_lines = notice + itp_beginning + moleculetype + new_atoms + new_bonds + new_pairs + new_angles + new_dihedrals

    # Writing file
    with open(itp_res_uncapped_file, "w") as f:
        for line in new_itp_lines:
            f.write(line) # No newline \n needed, as it is already in each element 

if __name__ == "__main__":
    # Main: rewriting the .ff file for the fragment on itself, deleting caps
    itp_res_file, capping_idx_file = load_input_command()

    section_names = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals"]
    itp_lines, section_linenrs, sections = read_itp_res(itp_res_file, section_names)
    
    capping_idx = read_capping_idx(capping_idx_file) # Atoms are 1-indexed, not 0 !
    new_atomidx, new_atoms = renumber_atoms(itp_lines, sections, capping_idx)
    
    new_bonds, new_pairs, new_angles, new_dihedrals = renumber_interactions(itp_lines, sections, new_atomidx)
    
    write_itp_res_uncapped(itp_lines, section_linenrs, new_atoms, new_bonds, new_pairs, new_angles, new_dihedrals, itp_res_file)
    
    # next: use link_residues.py to make the link.ff file for polyply
