#! /usr/bin

import subprocess

rmsd = subprocess.Popen("gmx rms -s md_result.tpr -f md_result_noPBC.xtc -o rmsd.xvg -tu ns", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press System
rmsd.stdin.write(b"0\n")
# Press C-alpha
rmsd.stdin.write(b"3\n")

rmsd.stdin.close()
rmsd.stdout.close()

