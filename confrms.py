#! /usr/bin

import subprocess

confrms = subprocess.Popen("gmx confrms -f1 protein.gro -f2 md_result.gro -o fit.pdb", stdout = subprocess.PIPE, stdin = subprocess.PIPE, shell = True)

# Press Protein
confrms.stdin.write(b"1\n")

# Press Protein
confrms.stdin.write(b"1\n")

confrms.stdin.close()
confrms.stdout.close()
