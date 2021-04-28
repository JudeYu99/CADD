#! /bin/bash

# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Address: Department of Bioinformatics, Medical College, Soochow University

## Special Attention: All the execution procedures may cause various errors, please check feedback and make adjustments to the following commands.
##					This workflow didn't involve NPT process!!!
## All required .mdp file will be uploaded together with this file.

## ===================================================================================================== ##

# @ Step One: Prepare the Topology
#
## Here, Amber99SB Force Field and TIP3P water model were used and the purpose of pdb2gmx is to produce a force field-compliant topology.
## Execute pdb2gmx to generate three files prepared for Gromacs by issuing the following command:
#
gmx pdb2gmx -f CCNB1.pdb -ff amber99sb -water tip3p -o protein.gro -p protein.top
#
## Three new files were generated: protein.gro, protein.top, and posre.itp. 



# @ Step Two: Defining the Unit Cell
#
## Here, cubic box type was used via -bt option. The protein was centered in the box (-c), and was placed at least 0.8 nm from the box edge (-d 0.8). 
## Commands to define the box using editconf:
#
gmx editconf -f protein.gro -o protein-PBC.gro -c -d 0.8 -bt cubic
#
## Only one file were generated: protein-PBC.gro



# @ Step Three: Adding Solvent
#
## Here, a empty unit cell requires filling of solvents. As a generic equilibrated 3-point solvent model, spc216.gro was used
## Commands for solvation is accomplished using solvate:
#
gmx solvate -cp protein-PBC.gro -cs spc216.gro -p protein.top -o protein-water.gro
#
## One file generated: protein-water.gro and one file updated: protein.top



# @ Step Four: Adding Ions
#
## Check the whole system whether a equilibrated one.
## ions.mdp file will be uploaded.
## Commands for checking:
#
gmx grompp -v -f ions.mdp -c protein-water.gro -p protein.top -o protein-water.tpr
#
## Two files were generated: mdout.mdp and protein-water.tpr
#
## == Attention ==##
## If the system is positively or negatively charged, extra ions should be added via the next command:
## The number of ions to be added depends on the feedback.
# gmx genion -s protein-water.tpr -o protein-solvated.gro -conc 0.15 -neutral -pname NA -nname CL -p protein.top



# @ Step Five: Energy Minimization
#
## Minimize the energy to improve the entire system. 
## minim.mdp file will be uploaded.
## Commands for minimizing the energy:
#
gmx grompp -f minim.mdp -c protein-water.gro -p protein.top -o protein-EM.tpr
gmx mdrun -v -s protein-EM.tpr -deffnm protein-EM-solvated
#
## One output file for part 1 and a series of protein-EM-solvated suffix files were generated. 
#
## Commands for analysising changes in potential. (Execute and press 10 after feedback shows up)
gmx energy -f protein-EM-solvated.edr -o potential-energy-EM.xvg
#
## One .xgv file will was generated and can be plotted via Xmgrace or Gnuplot.



# @ Step Six: NVT Equilibration
#
## Upgrade the temparature of the whole system to 300K.
## nvt.mdp file will be uploaded.
## Special Attention: Both -c and -r option are a must for Gromacs recent versions.
#
gmx grompp -f nvt.mdp -c protein-EM-solvated.gro -r protein-EM-solvated.gro -p protein.top -o protein-nvt.tpr
gmx mdrun -v -s protein-nvt.tpr -deffnm protein-nvt
#
## One output file for part 1 and a series of protein-nvt suffix files were generated. 
#
gmx energy -f protein-nvt.edr -o temperature.xvg
#
## One .xgv file will was generated and can be plotted via Xmgrace or Gnuplot.



# @ Step Seven: Production MD
#
## Perform the final step of the workflow which is considered as the most siganificant part of MD.
#
gmx grompp -f md.mdp -c protein-nvt.gro -t protein-nvt.cpt -p protein.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 -v
#
## One output file for part 1 and a series of md_0_1 suffix files were generated. 
## This step takes the longest time in the whole workflow and command nohup as background execution is preferred.



# @ Step Eight: Analysis
#
## With the output files from last step, following analysis can be realized.
## Here, RMSD and RMSF analysis are taken for example.
#
## First, .xtc file should be generated for following analysis.
#
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
## Press 1 and 0 in order after feedback shows up.
#
#
## Calculate RMSD for C-alpha:
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
## Press 0 and 3 in order after feedback shows up.
## One .xgv file will was generated and can be plotted via Xmgrace or Gnuplot.
#
#
## Calculate RMSF for C-alpha:
gmx rmsf -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsf.xvg
## Press 3 after feedback shows up.
## One .xgv file will was generated and can be plotted via Xmgrace or Gnuplot.
#
#
## Generate final pdb file for making MD movie.
#
gmx trjconv -s md_0_1.tpr -f md_0_1_noPBC.xtc -o movie.pdb
## VMD programme can be used to generate the movie.


