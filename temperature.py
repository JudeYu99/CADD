#! /usr/bin

import subprocess

# Upgrade the temperature
nvt_mdrun = subprocess.Popen("gmx energy -f protein-nvt.edr -o temperature.xvg", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press temperature
nvt_mdrun.stdin.write(b"16\n")

nvt_mdrun.stdin.close()
nvt_mdrun.stdout.close()

