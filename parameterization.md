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

The dimer file `BGD.sdf` looks like:
````sh
obj01
  PyMOL3.0          3D                             0

 43 44  0  0  0  0  0  0  0  0999 V2000
    0.9081    0.3455   -4.4794 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2457    0.6466   -3.2019 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7280   -0.2821   -2.0808 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8176   -0.2287   -1.9787 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4559   -0.5837   -3.3616 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4780    0.0220   -5.8727 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9759   -0.2734   -3.4958 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0881    0.9854   -0.5808 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2096    0.7625   -1.7447 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3692   -0.9049   -3.9814 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0297    0.6852   -6.6079 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2585   -1.0183   -6.1102 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5523    0.1819   -5.8992 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5517    0.2616   -4.4288 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6677    0.6072   -3.2681 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171    0.0888   -0.8379 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8706   -0.4756   -3.3385 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2160    1.3868   -4.2420 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9627    1.6998   -2.9725 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0384   -1.3291   -2.3033 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1988   -1.6422   -3.5859 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.7572   -1.8705   -0.6388 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3085   -1.2301   -0.3247 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2654   -2.0563    0.9810 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6644   -2.0880    1.6420 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7084   -2.6905    0.6713 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7669   -2.4690   -1.6317 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5153   -0.6090    2.1292 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3687   -3.7548    2.6998 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.9870   -3.2858    2.0672 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.7689   -1.8726   -2.5455 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.4809   -3.4953   -1.8616 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.7685   -2.4587   -1.2068 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3953   -1.7876   -1.1436 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0617   -1.1993   -0.9337 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3050   -1.5122    1.8787 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6237   -2.8409    2.8511 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9901   -2.7234    1.2924 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1002   -0.8377   -0.4054 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5816   -0.1718   -0.0790 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9450   -3.0981    0.7521 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.9679   -1.0514    1.9161 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.4305   -3.7441    0.4363 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1 14  1  0  0  0  0
  1 18  1  0  0  0  0
  2  3  1  0  0  0  0
  2 15  1  0  0  0  0
  2 19  1  0  0  0  0
  3  4  1  0  0  0  0
  3 16  1  0  0  0  0
  3 20  1  0  0  0  0
  4  5  1  0  0  0  0
  4  9  1  0  0  0  0
  5 17  1  0  0  0  0
  5 21  1  0  0  0  0
  6 11  1  0  0  0  0
  6 12  1  0  0  0  0
  6 13  1  0  0  0  0
  2 14  1  0  0  0  0
  7 15  1  0  0  0  0
  8 16  1  0  0  0  0
 10 17  1  0  0  0  0
 22 26  1  0  0  0  0
 22 27  1  0  0  0  0
 22 34  1  0  0  0  0
 22 39  1  0  0  0  0
 23 24  1  0  0  0  0
 23 35  1  0  0  0  0
 23 40  1  0  0  0  0
 24 25  1  0  0  0  0
 24 36  1  0  0  0  0
 24 41  1  0  0  0  0
 25 26  1  0  0  0  0
 25 37  1  0  0  0  0
 25 42  1  0  0  0  0
 26 38  1  0  0  0  0
 26 43  1  0  0  0  0
 27 31  1  0  0  0  0
 27 32  1  0  0  0  0
 27 33  1  0  0  0  0
 23 34  1  0  0  0  0
  4 35  1  0  0  0  0
 28 36  1  0  0  0  0
 29 37  1  0  0  0  0
 30 38  1  0  0  0  0
M  END
$$$$

````

### Poltype input files
Now we move on to using Poltype. 
REMINDER - To install and use Poltype2 - follow these links: [Poltype Installation Guide](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Install.md) and [Poltype Usage](https://github.com/WelbornGroup/Documentation/blob/Workflow_update/Poltype_Usage.md)

From the Usage guidelines, we know that we need three more files other than the input sdf file. To parameterize the `BGD.sdf` structure, I have the following three input files:

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
structure=BGD.sdf
atmidx=400
new_gdma=True
gdmacommand_Radius_S=0.80
prmmodfile=dma4_hfe2023
````
Check the first line - it is to have the name of the input `sdf` file. The second line is the first value of the atom types. Poltype will assign atom types here starting from 400. I have used 400 for beta-glucose, and 500 for the files where I have alpha-glucose, but you are free to choose differently, if you have an existing parameter file that overlaps with this numbering.

#### run-poltype.sh
Finally, we have a submit script for Tinkercliffs on ARC, named `run-poltype.sh` which looks like this:
````sh
#!/bin/bash
#SBATCH -J bgm
#SBATCH -A welbornlab
#SBATCH -p normal_q
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --time=23:00:00

# Run the example
echo "-------- Starting Poltype2: `date` -------"
source activate poltype_psi417
source paths.sh
python /projects/welbornlab/Poltype2/master/PoltypeModules/poltype.py
echo "------- Poltype2 has exited: `date` --------"
````
The details you should know from the Poltype Usage file.

#### Submitting the Job
Upload the four files on ARC (to a folder) and submit the job using the command `sbatch run-poltype.sh`.
Assuming everything goes as planned, you should see a `final.xyz` and a `final.key` file in that folder. Download the contents of the folder.

### Building a parameter file
