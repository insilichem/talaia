# TALAIA visual dictionary for proteins

Simplistic 3D dictionary of amino acids using geometric shapes.
![preview](https://user-images.githubusercontent.com/33349331/195623980-5d61567d-ab7f-41f6-a381-201c72f03744.png)

Originally developed by Mercè Alemany.


Installation
-----

To use Talaia you will first need to have UCSF Chimera installed.
Chimera can be downloaded free of charge for academic, government, non-profit and personal use at their website https://www.cgl.ucsf.edu/chimera/download.html.

Once UCSF Chimera is installed add Talaia as a third-party plugin. To do so go to Chimera'm menu Favorites/Add to Favorites and add the Talaia's project location in your computer.


Usage
-----

Talaia can only be used from the command-line for the moment.

Talaia expects a Chimera selection to work with. If none is provided, by default it will select ligands and every residue within 8A.

For custom selections the word spec can be used.
```
# Default behaviour (ligand zr < 8)

talaia

# Custom selections

talaia spec all

# Use sel for complex selections or quote it
sel :HIS zr < 8
talaia spec sel


# Change transparency (by defaut, opaque; 0.0)

talaia transparency 0.5
```

To disable depictions, use `~talaia`.
