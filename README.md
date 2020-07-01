# rmsdmat
Superimpose a set of protein structures and report a RSMD matrix, in CSV and Mega-compatible formats.

## Rationale
Given a set of protein structures the pairwise RMSD values between equivalent CA atoms are calculated. For
the RMSD calculation each pair is superimposed, as implemented in Julia's
[BioStructures](https://biojulia.net/BioStructures.jl/stable/) module, which for proteins with
differences in sequences it first imply a sequence alignment and then the corresponding residues are
superimposed following a [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm).

It is up to the user to select sets of structures that make sense to compare. Remember that for
proteins with non-identical sequences the reported RMSD value is for equivalent residues only, any
non-conserved sequence is excluded, as there is not equivalent residue to compare with.

No Pymol-like RMSD optimization cycles are done (where more different atoms are progressively excluded to improve
the fitting), all equivalent CA are compared, so the RMSD values reported by this programs are
generally larger that the ones reported by Pymol and equivalent software (but perhaps make more
sense?)

## Usage
```
rmsdmat [-o OUTPUT] [--version] [-h] structure...

Superimpose a set of protein structures and report a RSMD matrix, in
CSV and Mega-compatible formats.

positional arguments:
  structure            PDB, MMCIF or MMTF structural file

optional arguments:
  -o, --output OUTPUT  output files base name (default: "rmsd_matrix")
  --version            show version information and exit
  -h, --help           show this help message and exit
```

## Installation
This is not a Julia package, just a standalone script. Just download it and put it in your PATH. You need a working Julia environment and to install the dependencies.

## Dependencies
This script depends on the following Julia packages:

* ArgParse
* Printf
* BioStructures

