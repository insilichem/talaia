# How to use TALAIA on static sutructures

## Inner pore of a Potassium Channel
Here we will represent with TALAIA the residues of the inner pore of a potassiu channel (PDB code: 6eo1).

- Execute command 'open 6eo1' from the command line. The model opens as follows.
- Select the residues from the center of the pore using the initial orientation of the protein. It will select all residues independently of the depth.
- Call 'talaia spec sel' from the comand line. This will represent the TALAIA's figures for every residue selected.
- Changing the orientation of the whole protein you will see the physicochemical characteristics along the pore.

## Binding site of Vitamin D receptor
Another way to use TALAIA might be to represent the character of the binding site of a receptor.
We are now going to visualize the Vitamin D receptor (PDB code: 5xzf).

- Open the 5xzf model in Chimera running the command 'open 5xzf'.
- Selet the ligand
- Run the command 'sel sel zr < 7' to select every residue that is found at less than 7 Arg from the ligand.
- Run TALAIA to visualize the selected residues.

## Histidine representation according to protonation state
TALAIA takes into consideration the protonation state of the Histidine residues. They can be represented with 3 different colors, but for that you need to first make sure the model you are using has Hydrogens.
To see this, we are going to use a Cu/Zn Dismutase (PDB code: ????).

- 


## Non-standard residue treatment in TALAIA
Non-standard residues are diverse and difficult to treat in a unified manner. 
In TALAIA the non-standard residues that are bonded covalently to the backbone of the protein will be represented with a star shape placed on top of the Carbon alpha. 
The rest of the small molecules often present in PDB structures, like ligands, solvents, etc. will be ignored.

