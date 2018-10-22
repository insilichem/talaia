# TALAIA


Simplistic 3D dictionary of amino acids using geometric shapes.

Originally developed by Mercè Alemany.

Usage
-----

From the UCSF Chimera command-line (enable only):

Talaia expects a Chimera selection to work with. By default (none
provided), it will select ligands and everything within 8A. You
can change this with `spec`:

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