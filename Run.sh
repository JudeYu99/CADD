#! /bin/bash

# @ Step One: Prepare the Topology
nohup gmx pdb2gmx -f *.pdb -ff amber99sb -water tip3p -o protein.gro -p protein.top -ignh  > pdb2gmx.log 2>&1

# @ Step Two: Defining the Unit Cell
nohup gmx editconf -f protein.gro -o protein-PBC.gro -c -d 0.8 -bt cubic > editconf.log 2>&1

# @ Step Three: Adding Solvent
nohup gmx solvate -cp protein-PBC.gro -cs spc216.gro -p protein.top -o protein-water.gro > solvate.log 2>&1

# @ Step Four: Adding Ions
nohup gmx grompp -v -f ions.mdp -c protein-water.gro -p protein.top -o protein-water.tpr > ions.log 2>&1

# @ Step Five: Energy Minimization
nohup gmx grompp -f minim.mdp -c protein-water.gro -p protein.top -o protein-EM.tpr > minim.log 2>&1
nohup gmx mdrun -v -s protein-EM.tpr -deffnm protein-EM-solvated > minim_mdrun.log 2>&1

python potential-energy-EM-xvg.py
python potential-energy-EM-xvg-plot.py

# @ Step Six: NVT Equilibration
nohup gmx grompp -f nvt.mdp -c protein-EM-solvated.gro -r protein-EM-solvated.gro -p protein.top -o protein-nvt.tpr > nvt.log 2>&1
nohup gmx mdrun -v -s protein-nvt.tpr -deffnm protein-nvt > nvt_mdrun.log 2>&1

python temperature.py
python temperature-plot.py

# Step Seven: Production MD
nohup gmx grompp -f md.mdp -c protein-nvt.gro -t protein-nvt.cpt -p protein.top -o md_result.tpr > md_sample.log 2>&1
nohup gmx mdrun -deffnm md_result -v  > md_sample_mdrun.log 2>&1 &

# Step Eight: Analysis
python md_result_process.py

python rmsd.py
python rmsd-plot.py

python rmsf.py
python rmsf-plot.py

python energy.py
python energy-plot.py

python gyrate.py
python gyrate-plot.py

python confrms.py

python trjconv.py
