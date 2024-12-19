# Solvate a protein from PDB

This assumes that you have a TinkerXYZ file for your polymer. If you don't, export the structure from PyMOL as a PDB, and then use Tinker's `pdbxyz` executable to convert it to TinkerXYZ file format. (You need to have a `.key` file for that, which can be the same as that mentioned in [parameterization-for-amoeba.md](./parameterization-for-amoeba.md) - you'll have to change the box size (mentioned below.


### System charge and box size
The charge of the polymer is computed by adding the charge of all monomers. Since you've already parameterized the monomers (or have the parameters from somewhere), you know the charge of each monomer (do not forget to add the charges at the ends of the polymer, if they exist!!). Add them up to get a total charge of the system.

For the box size, open the structure of the polymer you made in PyMOL, and get an approximate idea of the size of the polymer by calculating the distance between the farthest atoms across each of the X,Y,Z axes. For that, there is a nice PyMOL tutorial [here[(https://www.compchems.com/how-to-measure-distances-and-display-interactions-in-pymol/#distances).

We will need to know how big the polymer is to make our solvation box for the next step. A good rule of thumb is to give your protein about 10 angstroms buffer space between each of the sides of the box. It is generally better to slightly overestimate the box size but you want to make sure that the polymer is never interacting with its mirror images across periodic boundaries (unless that is specifically something you're trying to simulate). ***This is especially tricky for flexible molecules!*** The downside to a larger volume box be the longer simulation time that is required. 

### Solvate the polymer on ARC/Local
There are a number of routines available to solvate a protein (Pacmol, GROMACS etc). Here we will use Tinker's `xyzedit` executable.
