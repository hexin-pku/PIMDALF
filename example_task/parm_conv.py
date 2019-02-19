import numpy as np
import sys

parm_file = str(sys.argv[1])
parm_file_flexible = parm_file + 'flx'

fo1 = open(parm_file, 'r')
fo2 = open(parm_file_flexible, 'w')

lines = fo1.readlines()
if not lines:
    print('Null')
    exit()

n = len(lines)//2
for i in range(n):
    flags = lines[2*i][1:-1].split(',')
    vals  = lines[2*i+1].split()
    for j in range(len(flags)):
        fo2.write('%s  %s\n'%(flags[j],vals[j]))

fo1.close()
fo2.close()

