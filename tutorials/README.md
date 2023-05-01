# TALAIA: a visual dictionary for proteins 3D representation

![preview](https://user-images.githubusercontent.com/33349331/195623980-5d61567d-ab7f-41f6-a381-201c72f03744.png)

Originally developed by Mercè Alemany and J.-D. Maréchal.


## Installation


Intalling TALAIA in UCSF Chimera is very easy and quick. You will need to follow these steps:

1. You need UCSF Chimera installed if not download it from https://www.cgl.ucsf.edu/chimera/download.html. 
If it is your first time using UCSF Chimera you can find information here: https://www.cgl.ucsf.edu/Outreach/Tutorials/GettingStarted.html. There is also plenty of tutorials.

Beware, that TALAIA has been developed under UCSF Chimera (which stands on Python 2). Deployments for ChimeraX (Python 3) are not yet finalized.

2. Download TALAIA from the insilichem's Github repository. You can either clone it or download it as ZIP.

![install_download_talaia_fro github_marcat](https://user-images.githubusercontent.com/63212606/234835874-5e6f2524-4ad3-41fb-ae4e-2794bbbf46a4.png)


3. To install TALAIA for the first time, open UCSF Chimera and go to the menu Favorites -> Add to Favorites/Toolbar...
![menu_favourites](https://user-images.githubusercontent.com/63212606/234820955-c83f2c35-323d-4cec-b83c-eda7264dc63e.PNG)

4. Add a new third-party plugin location. Here you have to make sure to select the directory where Talaia's directory is located (not the contents of the directory themselves). Save.

![install_add_location_marcat_2](https://user-images.githubusercontent.com/63212606/234833893-e21aafaa-caad-4f70-8db9-de7b057933d6.PNG)

5. As TALAIA only works through command line, you will need to activate it if it does not appear at the bottom of the Chimera window. In order to do so you have to go to Favorites -> Command line.
![install_step5_command_line](https://user-images.githubusercontent.com/63212606/234821003-8ff69e87-90b1-4834-b593-8617575e319b.PNG)

Now you are ready to use TALAIA!

## Usage


At the moment, Talaia can only be used in its command-line version. Make sure you have UCSF Chimera command line activated, if not, go to menu `Favorites/Command Line`.

Talaia expects a Chimera selection to work with. If none is provided, by default it will select ligands and every residue within 8A.

For custom selections the word spec can be used.
```
# Default behaviour (equivalent to talaia spec ligand zr < 8)

talaia


# Representation of the entire system

talaia spec all


# Representation of given selections

sel :HIS zr < 8
talaia spec sel


# Change transparency (by defaut, opaque; 0.0)

talaia transparency 0.5
```

To disable all Talaia's depictions, use `~talaia`.

Below you will find a few examples on TALAIA's usage.

Example 1 - Inner pore of a Potassium Channel
-----

TALAIA is used to depict the inner residues of the pore of a potassium channel as presented in the corresponding manuscript (PDB code: 6eo1).

- In the command line, execute  `open 6eo1`. The pdb file will be downloaded and model depicted.

![static_example1_asopen](https://user-images.githubusercontent.com/63212606/234840722-50d9834d-e1f5-4b11-9b7a-8228276e3bdf.PNG)

- Select the residues from the center of the pore using the initial orientation of the protein. It will select all residues independently of the depth.
To do so, you have to press Shift+CTRL as you drag your mouse clicking with the left botton along all the section you want to select.
Here we select only the pore region, taking advantage of the initial orientation of the model.

![statis_example1_selection_howto_2](https://user-images.githubusercontent.com/63212606/234841986-f94ce56d-175f-4b64-a052-8114ec0a6c2e.png)


- run `talaia spec sel` from the command line. This will represent the TALAIA's figures for every residue selected.

![statis_example1_selection_talaia](https://user-images.githubusercontent.com/63212606/234841297-eaf85eae-48a2-4567-bb73-420097a41cc7.PNG)

- Changing the orientation of the whole protein you will see the physicochemical characteristics along the pore.

![static_model1_final_result](https://user-images.githubusercontent.com/63212606/234841342-4b1b2a1a-a915-4461-a952-ce4ab764df88.PNG)



Example 2 - Binding site of Vitamin D receptor
-----

We are now going to visualize a Vitamin D receptor (PDB code: 5xzf).

- Open the 5xzf model in Chimera running the command `open 5xzf`.

![static_example2_asopen](https://user-images.githubusercontent.com/63212606/234842474-65d108be-5846-4899-826c-31df4ed90d0a.PNG)

- We want to change the background color to white, for this we run `background solid white`.
- We also would like to change the representation of the ligand. To do so we select the ligand with `sel :8j0` and go to the menu 'Actions -> Atoms/Bonds -> sphere'.
- With the ligand still selected, go to the command line again and run `focus sel`. Reorient the model with your mouse to visualize the ligand as you see fit.
The model looks like this now:

![static_example2_new_situation](https://user-images.githubusercontent.com/63212606/234844568-e8dc7133-676a-410c-b932-35647f8d47ee.PNG)

- With the ligand still selected, run the command `sel sel zr < 5` to select every residue that is found at less than 5 Arg from the ligand.
- Run TALAIA to visualize the selected residues.

![static_example2_talaia](https://user-images.githubusercontent.com/63212606/234845393-aa014e50-bc02-47db-bad6-3e9e070a0801.PNG)

- If you want you can hide the ribbon representation. This way it is very easy to see that the binding site physicochemical characteristics complement perfectly the ligand structure.

![static_example2_talaia_noribbon](https://user-images.githubusercontent.com/63212606/234845647-d078eeb7-9885-41a0-96b0-9b2331d2eb0e.PNG)



Example 3 - Histidine representation according to protonation state
-----

TALAIA takes into consideration the protonation state of the Histidine residues. They can be represented with 3 different colors, but for that you need to first make sure the model you are using has Hydrogens.
To see this, we are going to use a Cu/Zn Dismutase (PDB code: 1to5).

- Open the protein model with the command `open 1to5`.
- Select the residues 721 and 722 with the command `sel :711,712` and focus with `focus sel`.
- Select the Histidine residues using `sel :His` and represent their TALAIA models with `talaia spec sel transparency 0.5`. The resulting figures will be all of the same aquamarine color associated to polar residues without charge.

![static_exemple3_histidine_representations](https://user-images.githubusercontent.com/63212606/234997489-25ebe5b1-f904-427d-b00a-8ffcd5f5bfc6.PNG)

To represent the Histidines according to their protonation state first we need to add the hydrogens. To do so we follow the next steps:
- Run `~talaia` to remove all the current TALAIA representations.
- Run `addh` to add the hydrogen atoms to the whole model.
- Run `talaia spec sel transparency 0.5` again, when you have Histidines selected.
The result of this will look similar to:

![static_exemple3_histidine_representations_withH](https://user-images.githubusercontent.com/63212606/234997528-91327fda-0307-4f14-83be-0bdd3a1fb721.PNG)



Example 4 - Non-standard residue treatment in TALAIA
-----

Non-standard residues are diverse and difficult to treat in a unified manner. 
In TALAIA the non-standard residues that are bonded covalently to the backbone of the protein will be represented with a star shape placed on top of the Carbon alpha. 
The rest of the small molecules often present in PDB structures, like ligands, solvents, etc. will be ignored.
To exemplify the treatment of non-standard residues we are going to use as example the PDB structure of allophycocyanin. (PDB code: 6yx7)

- Run the command `open 6xy7` to open the model in Chimera.
- Select the residue MEN, a modified residue derived from ASN, with `sel :men`.
- To see at atom level only the side chain of the MEN residues run while you have them selected the command `show sel`. This will hide every other atomic structure of the model.
- Run `sel sel zr < 5` to select all the residues in proximity of 5 Arg from the original selection.
- Run `talaia spec sel` to represent every residue selected with TALAIA. It will look similiar to the picture below.


![static_example4_target_men_sphere5_talaia](https://user-images.githubusercontent.com/63212606/234998728-f21d33b2-7ef1-4823-8774-41e41472779e.PNG)

In this protein structure there are several non-standard residues, but only MEN and CYC will be recognized and graphically represented by TALAIA as these are the only ones tethered to the protein chain.



# Talaia's representation for MD trajectories

TALAIA can be used for MD trajectories aswell. To exemplify how to do it we are going to use a trajectory where insulin is interacting with oxaliplatin [1]. Due to this interaction one of the disulphide bridges of the insulin is broken, thus allowing structural changes to the protein structure. With TALAIA we can visualize the forces driving the conformational change along the trajectory.

One way to use TALAIA for MD trajectories is to:

- Load the MD trajectory using UCSF Chimera option from the menu Tools -> MD ensable/Analysis -> MDMovie

![movie_example5_menu_mdmovie2](https://user-images.githubusercontent.com/63212606/235442837-7971d28e-2743-4f21-b7b0-885a7f4d85fd.PNG)

- Load the trajectory selecting the specific format you are using. In this example we have a prmtop and a dcd files.

![movie_example5_mdmoview_file_options](https://user-images.githubusercontent.com/63212606/235442979-085966b9-19dd-4f85-811d-9b5116bf11ec.PNG)

- Once you have the trajectory loaded, open the Per-Frame menu and select Define script...

Here you have a templeate of a script you can use. 
In here you must specify the selection you want to represent. You can also specify transparency value if you wish.
```
# Remove pre-existing TALAIA figures
~talaia
# Make the selection you desire to represent
sel <define your selection in chimera's language>
# Call TALAIA on your current selection
talaia spec sel
```

![movie_example5_perframe_script](https://user-images.githubusercontent.com/63212606/235443498-5058f0e8-45cb-49f1-b610-baf46d8f8ae1.PNG)

As the trajectory advances we can see how the apolar moiety of oxaliplatin tends to establish hydrophobic interactions with non-polar residues by exposing the hydrophobic core of the insulin.

![test_frames3](https://user-images.githubusercontent.com/63212606/235450212-389f2179-d733-4354-921c-67f34bbf4f98.png)






-----

[1] Giuseppe Sciortino, José-Emilio Sánchez-Aparicio, Jaime Rodríguez-Guerra Pedregal, Eugenio Garribba, Jean-Didier Maréchal, Computational insight into the interaction of oxaliplatin with insulin, Metallomics, Volume 11, Issue 4, April 2019, Pages 765–773, https://doi.org/10.1039/c8mt00341f
