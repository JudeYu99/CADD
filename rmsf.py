#! /usr/bin

import subprocess

rmsf = subprocess.Popen("gmx rmsf -s md_result.tpr -f md_result_noPBC.xtc -o rmsf.xvg", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press C-alpha
rmsf.stdin.write(b"3\n")

rmsf.stdin.close()
rmsf.stdout.close()
