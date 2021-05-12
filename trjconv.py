#! /usr/bin

import subprocess

trjconv = subprocess.Popen("gmx trjconv -s md_result.tpr -f md_result_noPBC.xtc -o movie.pdb", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press Protein
trjconv.stdin.write(b"1\n")

trjconv.stdin.close()
trjconv.stdout.close()
