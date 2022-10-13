# TALAIA


Simplistic 3D dictionary of amino acids using geometric shapes.

Originally developed by Mercè Alemany.


Installation
-----

To install Talaia you will need to have UCSF Chimera installed. Chimera can be downloaded free of charge for academic, government, non-profit and personal use at their website https://www.cgl.ucsf.edu/chimera/download.html.

Once UCSF Chimera is installed add Talaia as a third-party plugin. To do so go to Chimera'm menu Favorites/Add to Favorites and add the Talaia's project location in your computer.


Usage
-----

Talaia can only be used from the command-line for the moment.

Talaia expects a Chimera selection to work with. If none is provided it will select ligands and every residue within 8A by default.

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
