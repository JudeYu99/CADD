#! /usr/bin

import subprocess

hbond_Pro_Pro = subprocess.Popen("gmx hbond -f md_result_noPBC.xtc -s md_result.tpr -num hbnum_Pro-Pro.xvg", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press Protein Twice
hbond_Pro_Pro.stdin.write(b"1\n")
hbond_Pro_Pro.stdin.write(b"1\n")
hbond_Pro_Pro.stdin.close()
hbond_Pro_Pro.stdout.close()
