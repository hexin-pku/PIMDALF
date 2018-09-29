#!/usr/bin/python3
# Filename: rvs.py

import pandas as pd
import sys

fflg=str(sys.argv[1])
tag=str(sys.argv[2])
val=str(sys.argv[3])
    
a=pd.read_csv(fflg,sep='\s+',header=None,dtype=str)

for i in range(len(a[0])):
    if(tag==str(a[0][i])):
        a[1][i]=val
        a.to_csv(fflg,sep='\t',index=False,header=False)

