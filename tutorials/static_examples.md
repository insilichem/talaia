# How to use TALAIA on static sutructures

## Inner pore of a Potassium Channel
Here we will represent with TALAIA the residues of the inner pore of a potassiu channel (PDB code: 6eo1).

- Execute command `open 6eo1` from the command line. The model opens as follows.

![static_example1_asopen](https://user-images.githubusercontent.com/63212606/234840722-50d9834d-e1f5-4b11-9b7a-8228276e3bdf.PNG)

- Select the residues from the center of the pore using the initial orientation of the protein. It will select all residues independently of the depth.
To do so, you have to press Shift+CTRL as you drag your mouse clicking with the left botton along all the section you want to select.
Here we select only the pore region, taking advantage of the initial orientation of the model.

![statis_example1_selection_howto_2](https://user-images.githubusercontent.com/63212606/234841986-f94ce56d-175f-4b64-a052-8114ec0a6c2e.png)


- Call `talaia spec sel` from the command line. This will represent the TALAIA's figures for every residue selected.

![statis_example1_selection_talaia](https://user-images.githubusercontent.com/63212606/234841297-eaf85eae-48a2-4567-bb73-420097a41cc7.PNG)

- Changing the orientation of the whole protein you will see the physicochemical characteristics along the pore.

![static_model1_final_result](https://user-images.githubusercontent.com/63212606/234841342-4b1b2a1a-a915-4461-a952-ce4ab764df88.PNG)



## Binding site of Vitamin D receptor
Another way to use TALAIA might be to represent the character of the binding site of a receptor.
We are now going to visualize the Vitamin D receptor (PDB code: 5xzf).

- Open the 5xzf model in Chimera running the command `open 5xzf`.

![static_example2_asopen](https://user-images.githubusercontent.com/63212606/234842474-65d108be-5846-4899-826c-31df4ed90d0a.PNG)

- We want to change the background color to white, for this we run `background solid white`.
- We also would like to change the representation of the ligand. To do so we select the ligand with `sel :8j0` and go to the menu 'Actions -> Atoms/Bonds -> sphere'.
- With the ligand still selected, go to the command line again and run `focus sel`. Reorient the model with your mouse to visualize the ligand as you see fit.
The model looks like this now:

![static_example2_new_situation](https://user-images.githubusercontent.com/63212606/234844568-e8dc7133-676a-410c-b932-35647f8d47ee.PNG)

- With the ligand still selected, run the command 'sel sel zr < 5' to select every residue that is found at less than 5 Arg from the ligand.
- Run TALAIA to visualize the selected residues.

![static_example2_talaia](https://user-images.githubusercontent.com/63212606/234845393-aa014e50-bc02-47db-bad6-3e9e070a0801.PNG)

- If you want you can hide the ribbon representation. This way it is very easy to see that the binding site physicochemical characteristics complement perfectly the ligand structure.

![static_example2_talaia_noribbon](https://user-images.githubusercontent.com/63212606/234845647-d078eeb7-9885-41a0-96b0-9b2331d2eb0e.PNG)


## Histidine representation according to protonation state
TALAIA takes into consideration the protonation state of the Histidine residues. They can be represented with 3 different colors, but for that you need to first make sure the model you are using has Hydrogens.
To see this, we are going to use a Cu/Zn Dismutase (PDB code: ????).

- 


## Non-standard residue treatment in TALAIA
Non-standard residues are diverse and difficult to treat in a unified manner. 
In TALAIA the non-standard residues that are bonded covalently to the backbone of the protein will be represented with a star shape placed on top of the Carbon alpha. 
The rest of the small molecules often present in PDB structures, like ligands, solvents, etc. will be ignored.

