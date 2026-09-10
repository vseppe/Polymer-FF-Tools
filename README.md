# Polymer-FF-Tools
Helper scripts to remove capping groups in GROMACS topologies and create linking files for PolyPly input.

## Installation
Requires Python 3.12 or later.

Requires the NumPy package.

## Quick start
### Removing caps

After obtaining a GROMACS-type of force field file (`.itp`, or renamed to Vermouth `.ff`) from e.g. ACPYPE, the capping group atoms and their corresponding bonded parameters can be removed. The atoms will be renumbered.
The atom indices (with a 1-based index) of the capping atoms need to be specified in a file, each on a new line. Comments are allowed. The name is free to choose, e.g. `capping_idx_POL.txt`.

In the case of a polyol residue with the force field file called `POL.ff`, this will be

```
python remove_caps_ff.py POL.ff capping_idx_POL.txt 
```
This produces the file `POL_uncapped.ff`, which is a functional force field file for the fragment.

### Interactions of linked fragments
PolyPly (Grünewald _et al._) requires linking files between fragments. 
To create the connection between the end of fragment 1 and the beginning of fragment 2, one must specify the atom names (from the `_uncapped.ff` file) that will be linked. If `O3` (fragment 1) will be linked to `C1` (fragment 2), a text file needs to specify both in that order (comments are allowed). Multiple connections can in principle be specified, but are not expected nor tested. The name is free to choose, e.g. `connection_CAR_POL.txt`.
```
# Connection between FRAG1 and FRAG2
O3 C1
```
Next for this stage, the parent force field needs to be supplied from which the bonded parameters can be retrieved. This can either be a `.dat` file, which will be read in as an Amber-type force field file (e.g. from GAFF2), or an `.itp` file as a GROMACS-type force field file (e.g. `ffbonded.itp` from Amber ff14SB). If some parameters cannot be found, the output link file is still written, with a warning printed in the terminal. The parameter that is not found is replaced by 999 and can be handled manually. Alternatively, one can alter the parent force field file after and run again.

For two different residues (e.g. CAR and POL), the command to generate `link_CAR_POL.ff` is
```
python link_interactions_ff.py CAR_uncapped.ff POL_uncapped.ff CAR POL connection_CAR_POL.txt gaff2.dat
```

To link a residue with itself (its end to the beginning of a second one), one can either type the same names or type `"same"` as the following.
```
python link_interactions_ff.py CAR_uncapped.ff same CAR same connection_CAR_CAR.txt gaff2.dat
```

Currently, the fetching of improper dihedrals is only implemented for parent force field files with GROMACS function `4` form (the same form as the proper dihedral), with the central atom as the third in row. Practically, this will work for Amber force fields such as ff14SB or GAFF2.
