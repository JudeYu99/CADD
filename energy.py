#! /usr/bin

import subprocess

energy = subprocess.Popen("gmx energy -f md_result.edr -o md_result_potential.xvg", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press 
energy.stdin.write(b"11\n")

energy.stdin.close()
energy.stdout.close()
