# Create Polymer with PyMOL

We will be using PyMOL builder to make the polymer, and using beta-glucose as an example.

You can download PyMOL [here](https://pymol.org/), and get an academic license if you do not have the paid version. 
There are a few extra features in the paid version, but none of them are used here in this tutorial,

The first step is to make a beta-blucose monomer. Open PyMOL, click `Builder` on the top-right side, and build the glucose ring.
You can use the 6-membered ring to start with, that way, you can get the hexagonal ring in one go.
If you think the bonds/angles does not look the way they are supposed to look, click the `Sculpt` button on the builder window to get an approximation of the lowest energy structure.
The beta-glucose monomer looks like this:
![monbuilder](https://github.com/user-attachments/assets/adc46ff0-94bc-482e-aef0-9cbd6ca0741a)

Once you're satisfied with the structure, Go to `File` -> `Export` and export it as a PDB.

Now you have the structure of the monomer. To make a dimer out of this, copy the monomer to a new object, and name it glucose (or anything else you want).
![copyobj](https://github.com/user-attachments/assets/78f4fe45-b538-4608-a3f6-66032f779c2b)

![objname](https://github.com/user-attachments/assets/c6670671-d845-4a00-97ac-bd2308cfe885)

For heteropolymers, it will be worthwhile to create objects for each kind of monomer, so that you can grab them as and when needed.

We only see one ring at this time because the two objects are superimposed onto each other. 
We will drag one of the rings to a position such that the two atoms that are supposed to bond are close to each other.
For example, if the C6 is supposed to bond with the O at C1, we try to superimpose the O bonded to C1 and the O bonded to C6 (that will keep the bond distances close to real values).
The distance we make in this input file should be as close to the real values as possible, that way minimization takes less time. 
Although if you feel you cannot do so, just create the bond as shown in the next step, and minimize the molecule later.

Now we will drag the second monomer to the position we want. To know the drag buttons and other mouse modes, check the bottom right corner of the window (3-Button Editing) and click the arrow next to it.

![buttons](https://github.com/user-attachments/assets/5bdafb35-fade-4671-aaec-b1e399551a6e)

Our system now looks like this.

![two-mer](https://github.com/user-attachments/assets/dc416fe5-a977-42c4-b468-da43780bc012)

Now copy the glucose object to obj01 (which you can rename to pol_progress or something like that if you want).

![copytoold](https://github.com/user-attachments/assets/766b2f96-06c4-4a4b-8035-49a3d9933e01)

NOTE - The order of the atoms in a PDB (or XYZ) file will be the same order you use to add atoms to a structure. 
That is why, we prefer making a monomer, and then making multiple copies of it, adding them to the structure one at a time.
Then the PDB will have atom ordering that shows one block for monomer 1, then monomer 2 and so on. This will help us change atom types, should we need to do so manually.



For heteropolymers, you'll probably have to come up with your own shortcuts to make it quickly, although I will be mentioning some general tips.