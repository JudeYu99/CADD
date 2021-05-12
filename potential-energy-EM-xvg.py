#! /usr/bin

import subprocess

minim_mdrun = subprocess.Popen("gmx energy -f protein-EM-solvated.edr -o potential-energy-EM.xvg", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press Potential
minim_mdrun.stdin.write(b"10\n")

minim_mdrun.stdin.close()
minim_mdrun.stdout.close()

