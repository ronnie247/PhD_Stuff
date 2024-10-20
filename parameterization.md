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
![BGM_OC](https://github.com/user-attachments/assets/73a42e08-81f5-4194-b796-f9187b05b1b4)

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
BGM_OC
  PyMOL3.0          3D                             0

 33 33  0  0  0  0  0  0  0  0999 V2000
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
    1.0379   -2.3553    2.1043 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4270   -2.6168    2.7156 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0620   -3.4560    2.5600 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1129   -2.3628    1.0169 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1394   -3.3454    3.6256 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8050   -3.5764    2.3629 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3489   -2.6337    3.8026 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.1114   -1.8241    2.4134 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5040   -4.4346    2.3728 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8709   -3.3669    2.0034 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0100   -1.5690    0.6260 H   0  0  0  0  0  0  0  0  0  0  0  0
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
  1 26  1  0  0  0  0
  1 28  1  0  0  0  0
  1  3  1  0  0  0  0
  2  3  1  0  0  0  0
  3  5  1  0  0  0  0
  3 29  1  0  0  0  0
  4 23  1  0  0  0  0
  5  7  1  0  0  0  0
  5 10  1  0  0  0  0
  5 30  1  0  0  0  0
  6 13  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
  7 31  1  0  0  0  0
  8 24  1  0  0  0  0
  6  9  1  0  0  0  0
  9 11  1  0  0  0  0
  9 32  1  0  0  0  0
 10 25  1  0  0  0  0
  2 11  1  0  0  0  0
  4 11  1  0  0  0  0
 11 33  1  0  0  0  0
 12 27  1  0  0  0  0
 13 14  1  0  0  0  0
 13 15  1  0  0  0  0
 13 16  1  0  0  0  0
 14 18  1  0  0  0  0
 14 19  1  0  0  0  0
 14 20  1  0  0  0  0
 15 17  1  0  0  0  0
 15 21  1  0  0  0  0
 15 22  1  0  0  0  0
M  END
>  <pdb_header>


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
Now that we have the parameters for both the monomer and the functionalized monomer. Lets talk about the monomer first.


### Building a parameter file

The file `final.key` contains the parameters for the regular and functionalized monomer. We need to copy them to a `.prm` file, which can be read by Tinker when we run and MD using the AMOEBA Force Field. But the `final.key` file also has a lot of comments, which start with a `#` symbol, which I prefer removing. To do so, I have a python script `copyfile.py`.
````sh
def copy_file_excluding_comments(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                outfile.write(line)
input_file = 'final.key'
output_file = 'beta_glucose.prm' #change this line if you need a different filename (like for the BGM_OC)
copy_file_excluding_comments(input_file, output_file)
print(f"Contents of {input_file} have been copied to {output_file} excluding lines that start with #.")
````
You can create this file in the folder that you downloaded for the `BGM` and `BGM_OC`. You can do so by using `cd` to go to that folder, then running the command `vi copyfile.py`, then pressing `I` to insert, and then pasting these lines. Hit `Esc`, then `:wq` to save this file. You can then run it in the same folder using the command `python copyfile.py`.

Now I have a `beta_glucose.prm` file in the folder for the monomer that I will be using and adding to, to parameterize the dimer.
The file should have first four lines that look something like this:
````sh
parameters /projects/welbornlab/Poltype2/master/ParameterFiles/amoebabio18_header.prm
OPENMP-THREADS 102
digits 8
RESP-WEIGHT 1
````
Remove these lines as they are not to be there in the prm file. You can also remove any literature references if you wish to, it makes copying easy later when we append this file to the `amoebabio18.prm` file. (We will talk about this later).
