#! /usr/bin

import subprocess

hbond_Pro_Sol = subprocess.Popen("gmx hbond -f md_result_noPBC.xtc -s md_result.tpr -num hbnum_Pro-Sol.xvg", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press Protein and SOL
hbond_Pro_Sol.stdin.write(b"1\n")
hbond_Pro_Sol.stdin.write(b"13\n")
hbond_Pro_Sol.stdin.close()
hbond_Pro_Sol.stdout.close()
