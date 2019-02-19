#!/usr/bin/python3
# Filename: mespro.py

import numpy as np
import pandas as pd
import sys
import os
import loadfile as lf

narg = len(sys.argv)
if(narg!=4):
    print("mesprocess arguments number error, only recieved for 1 argument!") 
a = pd.read_csv(str(sys.argv[1])+'.ana', header=None, sep='\s+')
a = a.values.T

print("The processing ID is %s"%sys.argv[1])


beta =  4511.066
nfree = 1
nbead = 32
    
pfd = a[1,:]
pff = a[2,:]
v = a[3,:]
kprim = a[4,:]
kvir = a[5,:]
eprim = a[6,:]
evir = a[7,:]
coh = a[8,:]
t1 = a[9,:]
t2 = a[10,:]
t3 = a[11,:]
t4 = a[12,:]

esti_pfd = pfd.mean()
esti_pff = pff.mean()
esti_v = v.mean()
esti_kprim = kprim.mean()
esti_kvir = kvir.mean()
esti_eprim = eprim.mean()
esti_evir = evir.mean()
esti_coh = coh.mean()

cprim = -0.5*nfree*nbead + beta**2/esti_pff * ( 2/beta * esti_kprim +
kprim**2/pff + 2*kprim*v/pff - (esti_eprim**2/esti_pff) ) + (
beta**2/esti_pff * (t1+t2)/pfd )

cvir = beta**2/esti_pff * (kprim*kvir/pff + kprim*v/pff 
- (esti_eprim*esti_evir/esti_pff) + 0.5*nfree/beta * esti_v 
+ esti_kvir/beta ) + (beta**2/esti_pff * ( t1+t2+t3+t4 )/pfd )

esti_cprim = cprim.mean() / 3.157746455E+5
esti_cvir = cvir.mean() / 3.157746455E+5



print("The estimated pfd                   : %.8e"%esti_pfd)
print("The estimated pff                   : %.8e"%esti_pff)
print("The estimated potential energy      : %.8e"%esti_v)
print("The estimated kinetic energy (prim) : %.8e"%esti_kprim)
print("The estimated kinetic energy (vir)  : %.8e"%esti_kvir)
print("The estimated total energy (prim)   : %.8e"%esti_eprim)
print("The estimated total energy (vir)    : %.8e"%esti_evir)
print("The estimated heat capacity (prim)  : %.8e"%esti_cprim)
print("The estimated heat capacity (vir)   : %.8e"%esti_cvir)
print("The estimated cohenrence length     : %.8e"%esti_coh)



