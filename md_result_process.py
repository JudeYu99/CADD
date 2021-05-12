#! /usr/bin

import subprocess

noPBC = subprocess.Popen("gmx trjconv -s md_result.tpr -f md_result.xtc -o md_result_noPBC.xtc -pbc mol -center", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press Protein
noPBC.stdin.write(b"1\n")
# Press System
noPBC.stdin.write(b"0\n")

noPBC.stdin.close()
noPBC.stdout.close()

