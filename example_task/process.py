#!/usr/bin/python3
# Filename: mespro.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os


narg = len(sys.argv)
if(narg!=4):
    print("mesprocess arguments number error, only recieved for 1 argument!") 
a = pd.read_csv(str(sys.argv[1])+'.ana', header=None, sep='\s+')
a = a.values.T

#print("The processing ID is %s"%sys.argv[1])


beta =  4511.066
nfree = int(sys.argv[3])
nbead = int(sys.argv[2])
    
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

mean_pfd = pfd.mean()
mean_pff = pff.mean()
mean_v = v.mean() 
mean_kprim = kprim.mean() 
mean_kvir = kvir.mean() 
mean_eprim = eprim.mean() 
mean_evir = evir.mean() 
mean_coh = coh.mean() 

esti_v = mean_v / mean_pff
esti_kprim = mean_kprim / mean_pff
esti_kvir = mean_kvir / mean_pff
esti_eprim = mean_eprim / mean_pff
esti_evir = mean_evir / mean_pff
esti_coh = mean_coh / mean_pff 

cprim = -0.5*nfree*nbead + beta**2/mean_pff * ( 
2/beta * kprim +
kprim**2/pff + 2*kprim*v/pff
+ (t1+t2)/pfd  ) - beta**2 * esti_eprim**2

cvir1 = beta**2/mean_pff * ( 
mean_kvir/beta + kprim*kvir/pff + (t1+t2)/pfd
+ ( kprim*v/pff ) * 2
+ ( 0.5*nfree/beta * mean_v + ( t3+t4 )/pfd ) * 0
) - beta**2*(esti_evir*esti_eprim)

cvir2 = beta**2/mean_pff * ( 
mean_kvir/beta + kprim*kvir/pff + (t1+t2)/pfd
+ ( kprim*v/pff ) * 1
+ ( 0.5*nfree/beta * mean_v + ( t3+t4 )/pfd ) * 1
) - beta**2*(esti_evir*esti_eprim)

cvir3 = beta**2/mean_pff * ( 
mean_kvir/beta + kprim*kvir/pff + (t1+t2)/pfd
+ ( kprim*v/pff ) * 0
+ ( 0.5*nfree/beta * mean_v + ( t3+t4 )/pfd ) * 2
) - beta**2*(esti_evir*esti_eprim)
 

# convert to au/K
# K is not a.u. unit

cprim /= 3.157746455E+5
cvir1 /= 3.157746455E+5
cvir2 /= 3.157746455E+5
cvir3 /= 3.157746455E+5

esti_cprim = cprim.mean() 
esti_cvir1 = cvir1.mean()
esti_cvir2 = cvir2.mean() 
esti_cvir3 = cvir3.mean() 

print( mean_pfd, mean_pff, esti_v, esti_kprim, esti_kvir,
esti_eprim, esti_evir,
esti_cprim, esti_cvir1, esti_cvir2, esti_cvir3,
esti_coh )
#print("The estimated pfd                   : %.8e"%mean_pfd)
#print("The estimated pff                   : %.8e"%mean_pff)
#print("The estimated potential energy      : %.8e"%esti_v)
#print("The estimated kinetic energy (prim) : %.8e"%esti_kprim)
#print("The estimated kinetic energy (vir)  : %.8e"%esti_kvir)
#print("The estimated total energy (prim)   : %.8e"%esti_eprim)
#print("The estimated total energy (vir)    : %.8e"%esti_evir)
#print("The estimated heat capacity (prim)  : %.8e"%esti_cprim)
#print("The estimated heat capacity (vir1)  : %.8e"%esti_cvir1)
#print("The estimated heat capacity (vir2)  : %.8e"%esti_cvir2)
#print("The estimated heat capacity (vir3)  : %.8e"%esti_cvir3)
#print("The estimated cohenrence length     : %.8e"%esti_coh)


#plt.hist(cprim, bins=200,normed=True, histtype='step', label='prim')
#plt.hist(cvir1, bins=200,normed=True, histtype='step', label='v1')
#plt.hist(cvir2, bins=200,normed=True, histtype='step', label='v2')
#plt.hist(cvir3, bins=200,normed=True, histtype='step', label ='v3')
#plt.show()



