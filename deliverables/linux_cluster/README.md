# PSE Molecular Dynamics WS23/24 - Group F

## Members
- Alp Kaan Aksu
- Berke Saylan
- Feryal Ezgi Aşkın

## Code
- Link:     https://github.com/alpkaanaksu/MolSim
- Branch:   master
- Revision: fffd16dfeeacc8750dbddf7859a0e536864e53e9

---

## Performance

It takes **28 seconds and 878 milliseconds** to run the first 1000 iterations of the Rayleigh-Taylor Instability simulation from task 2 with output disabled. (**346273.517 molecule updates  per second**)

You can find the output of the program in the file `molsimf.3922564.i23r05c05s01.out`.

## Running our Program on the Linux Cluster - Step by Step

### Step 1: Login to the cluster

### Step 2: Clone the repository

```bash
git clone https://github.com/alpkaanaksu/MolSim
```

- Go into the directory after cloning the repository

```bash
cd MolSim
```

### Step 3: Create a build directory

```bash
mkdir build
cd build
```

### Step 4: Load the modules

```bash
module load cmake
module load boost
module load xerces-c
module load gcc
```

### Step 5: Run CMake

```bash
ccmake ..
```

- Press `c` to configure
- Press `g` to generate

### Step 6: Compile the program

```bash
make
```

### Step 7: Create a job script

- Run the following command **AFTER** adding your email address to the script

```bash
echo '#!/bin/bash
        
#SBATCH -J molsimf

#SBATCH -o ./%x.%j.%N.out

#SBATCH -D ./

#SBATCH --clusters=serial

#SBATCH --partition=serial_std  

#SBATCH --get-user-env

#SBATCH --mail-type=end

#SBATCH --mem=800mb

#SBATCH --mail-user=alpkaan.aksu@tum.de

#SBATCH --export=NONE

#SBATCH --time=02:00:00

module load slurm_setup

module load gcc
module load boost

./MolSim ../input/performance/contest1.json' > molsimf.cmd
```

### Step 8: Submit the job

```bash
sbatch molsimf.cmd
```

---

- You can check the status of your job with `squeue --clusters serial`
- You will get an email when your job is finished
- The output is saved to a file with the following pattern in its name `molsimf.<jobid>.<node>.out`