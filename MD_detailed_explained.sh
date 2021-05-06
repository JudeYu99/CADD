#! /bin/bash

# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Address: Department of Bioinformatics, Medical College, Soochow University

## Special Attention: All the execution procedures may cause various errors, please check feedback and make adjustments to the following commands.
##					This workflow didn't involve NPT process!!!
## All required .mdp file will be uploaded together with this file.

############################################################################################################

#@ Step One: Prepare the Topology

# Now that the crystal waters are gone and we have verified that all the necessary atoms are present, the PDB file should contain only protein atoms, and is ready to be input into the first GROMACS module, pdb2gmx. 
# The purpose of pdb2gmx is to generate three files:
#	The topology for the molecule.
#	A position restraint file.
#	A post-processed structure file.
# The topology (topol.top by default) contains all the information necessary to define the molecule within a simulation. 
# This information includes nonbonded parameters (atom types and charges) as well as bonded parameters (bonds, angles, and dihedrals). 
# We will take a more detailed look at the topology once it has been generated.

# Execute pdb2gmx by issuing the following command:

nohup gmx pdb2gmx -f CCNB1.pdb -ff amber99sb -water tip3p -o protein.gro -p protein.top -ignh  > pdb2gmx.log 2>&1

# The structure will be processed by pdb2gmx, and you will be prompted to choose a force field.

# -ignh: Ignore H atoms in the PDB file; 
#	   especially useful for NMR structures. 
#	   Otherwise, if H atoms are present, they must be in the named exactly how the force fields in GROMACS expect them to be. 
#	   If you need to preserve the initial H coordinates, but renaming is required, then the Linux sed command is your friend.

# You have now generated three new files:  
#		protein.gro is a GROMACS-formatted structure file that contains all the atoms defined within the force field (i.e., H atoms have been added to the amino acids in the protein). 
#		The protein.top file is the system topology (more on this in a minute). 
#		The posre.itp file contains information used to restrain the positions of heavy atoms (more on this later).

# One final note: 
#		GROMACS can handle many different file formats, with .gro simply being the default for commands that write coordinate files. 
#		It is a very compact format, but it has limited precision. 
#		If you prefer to use, for instance, PDB format, all you need to do is to specify an appropriate file name with .pdb extension as your output. 
#		The purpose of pdb2gmx is to produce a force field-compliant topology; 
#		the output structure is largely a side effect of this purpose and is intended for user convenience. 
#		The format can be just about anything you like (see the GROMACS manual for different formats).


############################################################################################################

#@ Step Two: Defining the Unit Cell

# There are two steps to defining the box and filling it with solvent:
# Define the box dimensions using the editconf module.
# Fill the box with water using the solvate module (formerly called genbox).
# You are now presented with a choice as to how to treat the unit cell. 
#	For the purpose of this tutorial, we will use a simple cubic box as the unit cell. 
#	As you become more comfortable with periodic boundary conditions and box types, 
#	Highly recommend the rhombic dodecahedron, 
#	as its volume is ~71% of the cubic box of the same periodic distance, 
#	thus saving on the number of water molecules that need to be added to solvate the protein.

# Define the box using editconf:

nohup gmx editconf -f protein.gro -o protein-PBC.gro -c -d 0.8 -bt cubic > editconf.log 2>&1

# The above command centers the protein in the box (-c), and places it at least 0.8 nm from the box edge (-d 0.8). 
# The box type is defined as a cube (-bt cubic). 
# The distance to the edge of the box is an important parameter since we will be using periodic boundary conditions, we must satisfy the minimum image convention. 
# That is, a protein should never see its periodic image, otherwise the forces calculated will be spurious. 
# Specifying a solute-box distance of 0.8 nm will mean that there are at least 1.6 nm between any two periodic images of a protein. 
# This distance will be sufficient for just about any cutoff scheme commonly used in simulations.

# gmx editconf -f protein.gro -o protein-triclinic.gro -d 0.8 -bt triclinic

# gmx editconf -f protein.gro -o protein-dodecahedron.gro -d 0.8 -bt dodecahedron

# tail -1 protein-PBC.gro


############################################################################################################

#@ Step Three: Adding Solvent

# Now that we have defined a box, we can fill it with solvent (water). 
# Solvation is accomplished using solvate:

nohup gmx solvate -cp protein-PBC.gro -cs spc216.gro -p protein.top -o protein-water.pdb > solvate.log 2>&1
nohup gmx solvate -cp protein-PBC.gro -cs spc216.gro -p protein.top -o protein-water.gro > solvate.log 2>&1


# The configuration of the protein (-cp) is contained in the output of the previous editconf step, and the configuration of the solvent (-cs) is part of the standard GROMACS installation. 
# We are using spc216.gro, which is a generic equilibrated 3-point solvent model. 
# You can use spc216.gro as the solvent configuration for SPC, SPC/E, or TIP3P water, since they are all three-point water models. 
# The output is called 1AKI_solv.gro, and we tell solvate the name of the topology file (protein.top) so it can be modified. 
# Note the changes to the [ molecules ] directive of protein.top:
#	What solvate has done is keep track of how many water molecules it has added, which it then writes to your topology to reflect the changes that have been made. 
#	Note that if you use any other (non-water) solvent, solvate will not make these changes to your topology! 
#	Its compatibility with updating water molecules is hard-coded.


############################################################################################################

#@ Step Four: Adding Ions

# The tool for adding ions within GROMACS is called genion. 
# What genion does is read through the topology and replace water molecules with the ions that the user specifies. 
# The input is called a run input file, which has an extension of .tpr; 
# this file is produced by the GROMACS grompp module (GROMACS pre-processor), which will also be used later when we run our first simulation. 
# What grompp does is process the coordinate file and topology (which describes the molecules) to generate an atomic-level input (.tpr). 
# The .tpr file contains all the parameters for all of the atoms in the system.
# In my understanding, tpr = top + gro + mdp + ndx

# To produce a .tpr file with grompp, we will need an additional input file, with the extension .mdp (molecular dynamics parameter file); 
# grompp will assemble the parameters specified in the .mdp file with the coordinates and topology information to generate a .tpr file.

# Adding ions is accomplished using grompp and genion:

nohup gmx grompp -v -f ions.mdp -c protein-water.pdb -p protein.top -o protein-water.tpr > ions.log 2>&1
nohup gmx grompp -v -f ions.mdp -c protein-water.gro -p protein.top -o protein-water.tpr > ions.log 2>&1

# For unbalanced system, adding proper ions is a must.
# gmx genion -s protein-water.tpr -o protein-solvated.gro -conc 0.15 -neutral -pname NA -nname CL -p protein.top


############################################################################################################

#@ Step Five: Energy Minimization

# The solvated, electroneutral system is now assembled. 
# Before we can begin dynamics, we must ensure that the system has no steric clashes or inappropriate geometry. 
# The structure is relaxed through a process called energy minimization (EM).

# The process for EM is much like the addition of ions. 
# We are once again going to use grompp to assemble the structure, topology, and simulation parameters into a binary input file (.tpr).
# However, this time, instead of passing the .tpr to genion, we will run the energy minimization through the GROMACS MD engine, mdrun.

# Assemble the binary input using grompp using this input parameter file:

nohup gmx grompp -f minim.mdp -c protein-water.pdb -p protein.top -o protein-EM.tpr > minim.log 2>&1
nohup gmx grompp -f minim.mdp -c protein-water.gro -p protein.top -o protein-EM.tpr > minim.log 2>&1

# Output for .pdb file
nohup gmx mdrun -v -s protein-EM.tpr -deffnm protein-EM-solvated -c protein-EM-solvated.pdb > minim_mdrun.log 2>&1
nohup gmx mdrun -v -s protein-EM.tpr -deffnm protein-EM-solvated > minim_mdrun.log 2>&1

# The protein-EM-solvated.edr file contains all of the energy terms that GROMACS collects during EM. 
# You can analyze any .edr file using the GROMACS energy module:
gmx energy -f protein-EM-solvated.edr -o potential-energy-EM.xvg
# This data in .xvg format can be plotted via numerous plotting tools.


############################################################################################################

#@ Step Six: NVT Equilibration

# EM ensured that we have a reasonable starting structure, in terms of geometry and solvent orientation. 
# To begin real dynamics, we must equilibrate the solvent and ions around the protein. 
# If we were to attempt unrestrained dynamics at this point, the system may collapse. 
# The reason is that the solvent is mostly optimized within itself, and not necessarily with the solute. 
# It needs to be brought to the temperature we wish to simulate and establish the proper orientation about the solute (the protein). 
# After we arrive at the correct temperature (based on kinetic energies), we will apply pressure to the system until it reaches the proper density.

# The purpose of posre.itp is to apply a position restraining force on the heavy atoms of the protein (anything that is not a hydrogen). 
# Movement is permitted, but only after overcoming a substantial energy penalty. 
# The utility of position restraints is that they allow us to equilibrate our solvent around our protein, without the added variable of structural changes in the protein. 
# The origin of the position restraints (the coordinates at which the restraint potential is zero) is provided via a coordinate file passed to the -r option of grompp.

# Equilibration is often conducted in two phases. 
# The first phase is conducted under an NVT ensemble (constant Number of particles, Volume, and Temperature). 
# This ensemble is also referred to as "isothermal-isochoric" or "canonical." 
# The timeframe for such a procedure is dependent upon the contents of the system, but in NVT, the temperature of the system should reach a plateau at the desired value. 
# If the temperature has not yet stabilized, additional time will be required. 
# Depending on your machine, this may take a while (just under an hour if run in parallel on 16 cores or so).

# NOTE: In Gromcs new version, parameters setting for -c -r are necessary. 
# We will call grompp and mdrun just as we did at the EM step:

nohup gmx grompp -f nvt.mdp -c protein-EM-solvated.pdb -r protein-EM-solvated.pdb -p protein.top -o protein-nvt.tpr > nvt.log 2>&1
nohup gmx grompp -f nvt.mdp -c protein-EM-solvated.gro -r protein-EM-solvated.gro -p protein.top -o protein-nvt.tpr > nvt.log 2>&1

nohup gmx mdrun -v -s protein-nvt.tpr -deffnm protein-nvt > nvt_mdrun.log 2>&1

# The protein-nvt.edr file contains all of the energy terms that GROMACS collects during upgrading. 
# You can analyze any .edr file using the GROMACS energy module:
gmx energy -f protein-nvt.edr -o temperature.xvg
# This data in .xvg format can be plotted via numerous plotting tools.

# NOTE: In this script, we didn't take NPT equilibration into consideration.

############################################################################################################

#@ Step Seven: Production MD

# Upon completion of the two equilibration phases, the system is now well-equilibrated at the desired temperature and pressure. 
# We are now ready to release the position restraints and run production MD for data collection. 
# The process is just like we have seen before, as we will make use of the checkpoint file (which in this case now contains preserve pressure coupling information) to grompp.
# Then, execute mdrun for MD sampling:
nohup  gmx grompp -f md.mdp -c protein-nvt.gro -t protein-nvt.cpt -p protein.top -o md_result.tpr > md.log 2>&1
nohup gmx mdrun -deffnm md_result -v  &

# This step took the longest time.
# For a cubic box, the optimal setup will have a PME load of 0.25 (3:1 PP:PME - we're very close to optimal!); 
# for a dodecahedral box, the optimal PME load is 0.33 (2:1 PP:PME). 
# When executing mdrun, the program should automatically determine the best number of processors to assign for the PP and PME calculations. 
# Thus, make sure you indicate an appropriate number of threads/cores for your calculation (the value of -nt X), so that you can get the best performance.


############################################################################################################

#@ Step Eight: Analysis

# Now that we have simulated our protein, we should run some analysis on the system. 
# For this script, a few basic tools will be introduced.

# The first is trjconv, which is used as a post-processing tool to strip out coordinates, correct for periodicity, or manually alter the trajectory (time units, frame frequency, etc). 
# For this exercise, we will use trjconv to account for any periodicity in the system. 
# The protein will diffuse through the unit cell, and may appear "broken" or may "jump" across to the other side of the box. 
# To account for such behaviors, issue the following:

gmx trjconv -s md_result.tpr -f md_result.xtc -o md_result_noPBC.xtc -pbc mol -center

# Then, continue for the following analysis:

# Root mean square deviation (RMSD)
gmx rms -s md_result.tpr -f md_result_noPBC.xtc -o rmsd.xvg -tu ns

# Root mean square fluctuation (RMSF)
gmx rmsf -s md_result.tpr -f md_result_noPBC.xtc -o rmsf.xvg

# Radius of gyration (Rgyr)
gmx gyrate -s md_result.tpr -f md_result_noPBC.xtc -o gyrate.xvg

# Comparison between final conformation and original one.
gmx confrms -f1 protein.gro -f2 md_result.gro -o fit.pdb

# Extract trajectory for a movie, which can be visualized by VMD/PyMol/etc.
gmx trjconv -s md_result.tpr -f md_result_noPBC.xtc -o movie.pdb

