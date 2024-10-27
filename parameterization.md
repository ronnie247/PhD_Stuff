# Steps to follow to parameterize a polymer.

#### Poltype2
For the parameterization, we will be using the Poltype2.

To install and use Poltype2 - follow these links: [Poltype Installation Guide](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Install.md) and [Poltype Usage](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Usage.md)

## Tutorial
Here we will parameterize a beta-blucose dimer, to be able to parameterize a polymer with beta-glucose subunits. It looks like this:
![bgdimer](https://github.com/user-attachments/assets/e2cfb451-6813-42a1-8ef6-5f6e8d051ecf)

### Making the molecule

We will first make the molecule using PyMOL builder. You can download PyMOL [here](https://pymol.org/), and get an academic license if you need it.

The first step is to make a beta-blucose monomer and paramterizing it with Poltype2. Open PyMOL, click `Builder` on the top-right side, and build the glucose ring.

** HINT - You can use the 6-membered ring to start with, that way, you can get the hexagonal ring in one go. **

NOTE - If you think the bonda/angles does not look the way they is supposed to look, click the `Sculpt` button on the builder window to get an approximation of the lowest energy structure.
The beta-glucose monomer looks like this:
![monbuilder](https://github.com/user-attachments/assets/adc46ff0-94bc-482e-aef0-9cbd6ca0741a)
Once you're satisfied with the structure, Go to `File` -> `Export` and export it as a PDB.

NOTE - Since we will use `.sdf` files as inputs to Poltype2, export another copy as a SDF file for later.

However, you see that the `OH` group on the `C1` is not how it occurs in the dimer. This oxygen atom in the dimer (and similarly a polymer, if need be) is connected to a carbon atom, and not a hydrogen.
Since we would need to model the behaviour of that oxygen as it bonds to the next ring. Ideally, we would run poltype for a beta-glucose dimer, and get the parameters for the linkage from there. But, running poltype for the dimer requires too much time and memory to run correctly on ARC, so we improvise. We add a `CH3` group to the oxygen, and since that carbon is supposed to be tertiary, we add one `CH3` group on either side, to account for that environment. The molecule we end up parameterizing is this (check the functional group attached to the OH:
![bgm_oc](https://github.com/user-attachments/assets/8d04e930-14f5-40e4-9d70-686d3e4cca0c)


Save it as a PDB and an SDF file. (Here I've named them `BGM_OC`, and will be referring to this molecule as the "functionalized monomer".)

NOTE - Also make the beta-glucose dimer, and save it as a PDB file. We will use it later, when we check out parameters.

### Structure input files
So far, we have the following SDF files we will use as inputs:
The monomer file `BGM.sdf` looks like:
````sh
BGM
  PyMOL3.0          3D                             0

 24 24  0  0  0  0  0  0  0  0999 V2000
   -0.5700    1.1470   -2.6280 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9670    0.1070   -0.9140 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3920    0.2200   -1.4120 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6560   -0.6770    0.6530 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2830    0.6960   -0.2530 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5580   -1.0810    2.5410 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0860   -0.1860    0.9900 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8900    0.3020    2.0670 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3920   -0.2200    1.4120 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6560    0.6770   -0.6530 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2830   -0.6960    0.2530 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6910    1.7230   -2.9790 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0100   -1.5690    0.6260 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8800   -1.9360    2.2480 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9960    1.2530    1.9840 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0090   -0.2090   -0.5440 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9550    0.5720   -3.4710 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9440    2.3680   -2.3150 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2750    1.9410   -2.3790 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7110   -0.8010   -1.7310 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0320    1.7550   -0.0020 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4340   -1.2240    0.7690 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7110    0.8010    1.7310 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0320   -1.7550    0.0020 H   0  0  0  0  0  0  0  0  0  0  0  0
  1 12  1  0  0  0  0
  1 17  1  0  0  0  0
  1 19  1  0  0  0  0
  1  3  1  0  0  0  0
  2  3  1  0  0  0  0
  3  5  1  0  0  0  0
  3 20  1  0  0  0  0
  4 13  1  0  0  0  0
  5  7  1  0  0  0  0
  5 10  1  0  0  0  0
  5 21  1  0  0  0  0
  6 14  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
  7 22  1  0  0  0  0
  8 15  1  0  0  0  0
  6  9  1  0  0  0  0
  9 11  1  0  0  0  0
  9 23  1  0  0  0  0
 10 16  1  0  0  0  0
  2 11  1  0  0  0  0
  4 11  1  0  0  0  0
 11 24  1  0  0  0  0
 12 18  1  0  0  0  0
M  END
>  <pdb_header>


$$$$

````

The new file for the functionalized monomer `BGM_OC.sdf` looks like:
````sh
bgm_oc
  PyMOL3.0          3D                             0

 33 33  0  0  0  0  0  0  0  0999 V2000
   -0.1991   -0.5323   -0.8443 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0904   -0.0563    0.3141 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8937   -0.9381    1.5577 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5841   -0.9724    1.9791 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4755   -1.4484    0.8208 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3761   -2.7577    1.1795 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3778    0.3947   -2.0610 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7375   -2.7909    1.8987 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4014   -3.7204    1.8829 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7622   -0.1805   -2.9033 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1366    1.6161   -1.7479 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0823    1.1883   -1.8120 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0726   -2.6885    2.8155 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8031    0.5004    2.5510 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8167   -0.9611    0.0237 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5063   -3.0654    0.1420 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2165   -3.3728    2.8993 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1166   -3.8127    1.9133 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.6178   -2.4325    2.9211 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.4428   -2.1498    1.3698 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8355   -4.7198    1.9134 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4608   -3.7501    1.3329 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1600   -0.6452   -0.3466 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8832    0.9706   -2.4118 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7501   -1.8332    3.1089 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6977   -0.4501    2.6349 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4631   -0.0754   -0.0861 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8483   -1.4293    1.2209 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5184   -1.5537   -1.1634 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8392    1.0023    0.5655 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2412   -1.9759    1.3361 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9034    0.0489    2.2982 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2243   -2.5070    0.5693 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  7  1  0  0  0  0
  1 29  1  0  0  0  0
  2  3  1  0  0  0  0
  2 27  1  0  0  0  0
  2 30  1  0  0  0  0
  3  4  1  0  0  0  0
  3 26  1  0  0  0  0
  3 31  1  0  0  0  0
  4  5  1  0  0  0  0
  4 25  1  0  0  0  0
  4 32  1  0  0  0  0
  5 28  1  0  0  0  0
  5 33  1  0  0  0  0
  6  8  1  0  0  0  0
  6  9  1  0  0  0  0
  6 16  1  0  0  0  0
  7 10  1  0  0  0  0
  7 12  1  0  0  0  0
  7 24  1  0  0  0  0
  8 18  1  0  0  0  0
  8 19  1  0  0  0  0
  8 20  1  0  0  0  0
  9 17  1  0  0  0  0
  9 21  1  0  0  0  0
  9 22  1  0  0  0  0
  1 23  1  0  0  0  0
  5 23  1  0  0  0  0
 11 24  1  0  0  0  0
 13 25  1  0  0  0  0
 14 26  1  0  0  0  0
 15 27  1  0  0  0  0
  6 28  1  0  0  0  0
M  END
$$$$

````

### Poltype input files
Now we move on to using Poltype.

REMINDER - To install and use Poltype2 - follow these links: [Poltype Installation Guide](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Install.md) and [Poltype Usage](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Usage.md)

From the Usage guidelines, we know that we need three more files other than the input sdf file. To parameterize the `BGM_OC.sdf` structure, I have the following three input files:

#### paths.sh
The first is the `paths.sh` file, which looks like this for me:
````sh
export PATH=/projects/welbornlab/Poltype2/TinkerEx/:$PATH
export GDMADIR=/projects/welbornlab/Poltype2/bin/
export PATH=/projects/welbornlab/Poltype2/bin/:$PATH
export PATH=/home/mondal/miniconda3/envs/xtbenv/bin/:$PATH
export PSI_SCRATCH=/localscratch/   

````
The only line you need to change is the `export PATH=/home/mondal/miniconda3/envs/xtbenv/bin/:$PATH` line, which will be the path to your `xtbenv/bin`.

#### poltype.ini
The second is the `poltype.ini` file, which looks like this for me:
````sh
structure=BGM_OC.sdf
atmidx=400
new_gdma=True
gdmacommand_Radius_S=0.80
prmmodfile=dma4_hfe2023
````
Check the first line - it is to have the name of the input `sdf` file. The second line is the first value of the atom types. Poltype will assign atom types here starting from 400. I have used 400 for functionalized beta-glucose, and 500 for the beta-glucose monomer but you are free to choose differently, if you have an existing parameter file that overlaps with this numbering.

#### run-poltype.sh
Finally, we have a submit script for Tinkercliffs on ARC, named `run-poltype.sh` which looks like this:
````sh
#!/bin/bash
#SBATCH -J bgm_oc
#SBATCH -A welbornlab
#SBATCH -p welborn_q
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --time=5-23:00:00
#SBATCH --mem=128G
# Run the example
echo "-------- Starting Poltype2: `date` -------"
source activate poltype_psi417
source paths.sh
python /projects/welbornlab/Poltype2/master/PoltypeModules/poltype.py
echo "------- Poltype2 has exited: `date` --------"
````
The details you should know from the Poltype Usage file.

#### Submitting the Job
Upload the four files on ARC (to a folder) and submit the job using the command `sbatch run-poltype.sh`. This will run Poltype for the dimer.
Assuming everything goes as planned, you should see a `final.xyz` and a `final.key` file in that folder. Download the contents of the folder.

TASK - Now repeat this step in a different folder for the monomer, and run Poltype on that one as well.

### Identifying the atom types
Now that we have the parameters for both the monomer and the functionalized monomer. Lets talk about the monomer first. Open the `final.xyz` file of the monomer, and look at the contents.
My `final.xyz` file for the monomer looks like this:
````sh
    24
     1  C      2.591438   -0.511873   -0.335508      500    3    12    17    19
     2  O      0.376997   -1.307845   -0.080451      501    3    11
     3  C      1.131526   -0.162896   -0.476801      502    1     2     5    20
     4  O     -1.657686   -2.326217    0.028653      503   11    13
     5  C      0.726455    1.026047    0.383543      504    3     7    10    21
     6  O     -2.915971    0.253632    0.176162      505    9    14
     7  C     -0.769743    1.256992    0.293519      506    5     8     9    22
     8  O     -1.106868    2.240439    1.273274      507    7    15
     9  C     -1.556303   -0.016723    0.515279      508    6     7    11    23
    10  O      1.421236    2.183361   -0.077005     509    5    16
    11  C     -1.010998   -1.142317   -0.343434     510    2     4     9    24
    12  O      2.900673   -1.467836   -1.344478     511    1    18
    13  H     -1.848650   -2.836579   -0.774144     512    4
    14  H     -3.468190   -0.282456    0.768040     513    6
    15  H     -1.843396    2.760216    0.912316     514    8
    16  H      1.473920    2.791107    0.679336     515   10
    17  H      3.180285    0.406509   -0.447465     516    1
    18  H      3.812661   -1.766067   -1.194504     517   12
    19  H      2.750134   -0.919583    0.672821     516    1
    20  H      0.935553    0.080759   -1.534662     518    3
    21  H      0.953252    0.820114    1.440731     519    5
    22  H     -1.013476    1.615998   -0.717499     520    7
    23  H     -1.454259   -0.354853    1.555096     521    9
    24  H     -1.143683   -0.927498   -1.417403     522   11
````
This is the format of a Tinker xyz file. Column 1 is the atom index, Column 2 is the atom name, Columns 3-5 are the x, y and z coordinates. Column 6 is the atom type, which we need, and the rest of the columns are connectivities. These are the atom indices of the atoms this one is connected to. Which means, atom 1 which is carbon of type 500, is connected to atoms 3, 12, 17 and 19.
Draw the structure of the monomer on paper, that will help you assign the atom types.

Check the bonding of the different atoms in the structure, and the xyz file, and find out which atom is assigned which number as the atom type.

HINT - There is only one carbon that is bonded to two hydrogens (which will be equivalent, i.e. having the same atom type), and there is only one oxygen bonded to two carbons.
This is what my monomer looks like:

![monomer](https://github.com/user-attachments/assets/4c0c7188-8b56-494c-843b-d8edff136147)


TASK - Do the same for the functionalized monomer. For that molecule, my `final.xyz` file looks like this:
````sh
    33 
     1  C     -0.960980    1.206418   -0.263471      400    2     7    23    29
     2  C     -1.970961    0.073744   -0.382255      401    1     3    27    30
     3  C     -1.306140   -1.231243   -0.823281      402    2     4    26    31
     4  C     -0.070436   -1.546024   -0.008426      403    3     5    25    32
     5  C      0.822548   -0.319565    0.062707      404    4    23    28    33
     6  C      2.948836    0.330620    0.907843      405    8     9    16    28
     7  C     -1.483568    2.481025    0.373418      406    1    10    12    24
     8  C      3.917535   -0.017958   -0.208521      407    6    18    19    20
     9  C      3.584998    0.267916    2.279976      407    6    17    21    22
    10  H     -2.047367    2.224056    1.282091     414    7
    11  H     -2.657055    3.948035   -0.126940     415   24
    12  H     -0.608261    3.071664    0.677573     414    7
    13  H      1.141022   -3.057654    0.082234     416   25
    14  H     -2.238731   -2.692437    0.073054     417   26
    15  H     -3.757624   -0.080162   -1.115076     418   27
    16  H      2.523042    1.327392    0.744315     419    6
    17  H      2.834970    0.461742    3.049275     420    9
    18  H      3.430215    0.012603   -1.186066     420    8
    19  H      4.752357    0.689294   -0.225868     420    8
    20  H      4.313389   -1.024579   -0.050850     420    8
    21  H      4.005001   -0.727406    2.449673     420    9
    22  H      4.388166    1.004609    2.367676     420    9
    23  O      0.097664    0.779383    0.604809     408    1     5
    24  O     -2.282116    3.172382   -0.575499     409    7    11
    25  O      0.623721   -2.625116   -0.617532     410    4    13
    26  O     -2.272409   -2.282879   -0.808904     411    3    14
    27  O     -2.958624    0.424489   -1.339580     412    2    15
    28  O      1.863732   -0.625926    0.930784     413    5     6
    29  H     -0.572032    1.448889   -1.265810     421    1
    30  H     -2.417946   -0.117889    0.607715     422    2
    31  H     -1.011139   -1.117505   -1.873023     423    3
    32  H     -0.325822   -1.763749    1.042349     424    4
    33  H      1.188047   -0.057090   -0.946331     425    5

````
and the atom types look like:
![bgm_oc](https://github.com/user-attachments/assets/a46f4343-f058-459e-bacc-8606057dae88)


### Building a parameter file
The file `final.key` contains the parameters for the regular and functionalized monomer. We need to copy them to a `.prm` file, which can be read by Tinker when we run and MD using the AMOEBA Force Field. But the `final.key` file also has a lot of comments, which start with a `#` symbol, which I prefer removing. To do so, I have a python script `copyfile.py`.
````sh
def copy_file_excluding_comments(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                outfile.write(line)
input_file = 'final.key'
output_file = 'beta_glucose_test.prm' #change this line if you need a different filename (like for the BGM_OC)
copy_file_excluding_comments(input_file, output_file)
print(f"Contents of {input_file} have been copied to {output_file} excluding lines that start with #.")
````
You can create this file in the folder that you downloaded for the `BGM` and `BGM_OC`. You can do so by using `cd` to go to that folder, then running the command `vi copyfile.py`, then pressing `I` to insert, and then pasting these lines. Hit `Esc`, then `:wq` to save this file. You can then run it in the same folder using the command `python copyfile.py`.

Now I have a `beta_glucose_test.prm` file in the folder for the monomer that I will be using and adding to, to parameterize the dimer.
This file should have the atom types described in the beginning which look something like this (only showing the first few lines):
````sh
atom          500    500    C     "BGM                 "         6    12.011    4
atom          501    501    O     "BGM                 "         8    15.999    2
atom          502    502    C     "BGM                 "         6    12.011    4
atom          503    503    O     "BGM                 "         8    15.999    2
atom          504    504    C     "BGM                 "         6    12.011    4
````
Copy all the atom definitions, and add them to the list of atom definitions in the `amoebabio18.prm` file (after the atom types for amino acids, NAs, Ions etc, and before the `SOLUTE` parameters). 
For the other parameters in the `beta_glucose_test.prm` file (torsion, SOLUTE, polarize, mutipole, bond, angle, strbend, vdw - the order that I have in my test file), simply add all those lines at the bottom of the `amoebabio18.prm` file. Save this file, and rename it to `beta_glucose_params.prm`.
We will use this file to add new parameters, and run MD in Tinker.

### Adding Parameters for the Linkage
