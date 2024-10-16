Steps to follow to parameterize a polymer.

For the parameterization, we will be using the Poltype2.
To install and use Poltype2 - follow these links: [Poltype Installation Guide](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Install.md) and [Poltype Usage](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Usage.md)

Here we will parameterize a beta-blucose dimer, to be able to parameterize a polymer with beta-glucose subunits. The molecule looks like this:
![bgdimer](https://github.com/user-attachments/assets/e2cfb451-6813-42a1-8ef6-5f6e8d051ecf)

We will first make the molecule using PyMOL builder. You can download PyMOL [here](https://pymol.org/), and get an academic license if you need it.

The first step is to make a beta-blucose monomer and paramterizing it with Poltype2. Open PyMOL, click `Builder` on the top-right side, and build the glucose ring.
HINT - You can use the 6-membered ring to start with, that way, you can get the hexagonal ring in one go.
NOTE - If you think the bonda/angles does not look the way they is supposed to look, click the `Sculpt` button on the builder window to get an approximation of the lowest energy structure.
The beta-glucose monomer looks like this:
![monbuilder](https://github.com/user-attachments/assets/adc46ff0-94bc-482e-aef0-9cbd6ca0741a)
Once you're satisfied with the structure, Go to `File` -> `Export` and export it as a PDB.
NOTE - Since we will use `.sdf` files as inputs to Poltype2, export another copy as a SDF file for later.

However, you see that the `OH` group on the `C1` is not how it occurs in the dimer. This oxygen atom in the dimer (and similarly a polymer, if need be) is connected to a carbon atom, and not a hydrogen.
Since we would need to model the behaviour of that oxygen as it bonds to the next ring. For that we need the beta-glucose dimer.

Build the beta-glucose dimer, which looks like this, after clicking on `Sculpt`.
![dimbuilder](https://github.com/user-attachments/assets/915c1a51-5e11-49f3-b4c8-e984f4fd7e02)

Save it as a PDB and an SDF file.
