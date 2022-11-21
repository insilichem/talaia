# TALAIA: a visual dictionary for proteins 3D representation

![preview](https://user-images.githubusercontent.com/33349331/195623980-5d61567d-ab7f-41f6-a381-201c72f03744.png)

Originally developed by Mercè Alemany and J.-D. Maréchal.


Installation
-----

To use Talaia, you will first need to have UCSF Chimera installed.
Chimera can be downloaded free of charge for academic, government, non-profit and personal use at their website https://www.cgl.ucsf.edu/chimera/download.html.

Once UCSF Chimera is installed, add Talaia as a third-party plugin. To do so, in Chimera, go to the `Favorites/Add to Favorites` menu and add the Talaia's project location in your computer.


Usage
-----

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

# Talaia's representation for MD trajectories

To make Talaia's work on a MD trajectory for all frames, one procedure is:
```
1. load the MD trajectory using the UCSF Chimera option (Tools-->MD ensemble/Analysis --> MD movie)
2. in the MD movie widget, go to "per frame"
3. copy the script:
~talaia
sel (define your selection in chimera's language)
talaia spec sel
```

