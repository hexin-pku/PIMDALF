#!/usr/bin/python3
# Filename: cct.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import loadfile as lf

# get arguments from terminal
narg=len(sys.argv)-1
fflg=str(sys.argv[1])  #the filename to plot
clmn1=int(sys.argv[2]) #the column to plot
clmn2=int(sys.argv[3])

myinfo = lf.Load('info.now')
myinfo.add('put.rc')

os.system('\nhead -n 1 %s'%fflg)

a=pd.read_csv(fflg,sep='\s+')
name='null'
if not setinfo.is_number(a.values[0,0]):
    a=pd.read_csv(fflg,sep='\s+')
    name=str(a.columns.values[clmn1])
fflg=(fflg.spilt('.',1))[0]

mtx=a.values.T
x=mtx[0,]
if(clmn1>0):
    y1=mtx[clmn1:clmn1+1,]
    y2=mtx[clmn2:clmn2+1,]
else:
	exit()


rangcct=100 # the max length for stat
u0ui=np.zeros((rangcct,len(y1)))
I1=np.ones(rangcct)
for i in range(1,rangcct+1):
    z1=y1.T[:-i,]
    z2=y2.T[i:,]
    cov=np.multiply(z1,z2)
    u0ui[i-1]=cov.mean(axis=0)
um1=np.outer(I1,y1.mean(axis=1))
um2=np.outer(I1,y2.mean(axis=1))
umum=np.multiply(um1,um2)
yy=np.multiply(y1,y2)
u2m=np.outer(I1,yy.mean(axis=1))
cct=(u0ui-umum)/(u2m-umum)
mcct=cct.T.mean(axis=0)
#

scl=myinfo.args['nstep']//myinfo.args['nsmp'] # adjust the plot range
if( scl > rangcct):
    scl=rangcct
dn=np.arange(0,scl+1)*float(myinfo.args['dtime'])*myinfo.args['nsmp']
#
addend=scl//2

plt.plot(dn[:scl],mcct[:scl],'r--')
rslt=mcct[:addend].sum()*myinfo.args['dtime']


plt.ylabel('Character Corelation Function of\n Potential of Standard') #
plt.xlabel('time [a.u.]') #
plt.savefig('cct_'+fflg+str(name)+'.png')
plt.show()

lst=pd.DataFrame({'cct':mcct})
lst.to_csv('.cctdat.tmp')



