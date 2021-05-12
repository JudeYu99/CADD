#! /usr/bin

import subprocess

gyrate = subprocess.Popen("gmx gyrate -s md_result.tpr -f md_result_noPBC.xtc -o gyrate.xvg", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press Protein
gyrate.stdin.write(b"1\n")

gyrate.stdin.close()
gyrate.stdout.close()
