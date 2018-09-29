#!/usr/bin/python3
# Filename: cct.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import loadfile as lf

print('''
### < by Shin He > this procedure is mainly do with IR calculation of water
### each line should give a/or some statistic quantity at each step
### only one argument need
### make sure _put.in _plot.in here\n
''')
fflg=str(sys.argv[1]) #the filename to plot

myinfo = lf.Load('info.now')
myinfo.add('put.rc')

os.system('\nhead -n 1 %s'%fflg)

a=pd.read_csv(fflg,sep='\s+',header=None,low_memory=False)
if not setinfo.is_number(a.values[0,0]):
    a=pd.read_csv(fflg,sep='\s+',low_memory=False)
fflg=(fflg.spilt('.',1))[0]

mtx=a.values.T
x=mtx[0,]
y=mtx[1:,]

rangcct=500000 # the max length for stat
u0ui=np.zeros((rangcct,len(y)))
I1=np.ones(rangcct)
for i in range(1,rangcct+1):
    z1=y.T[:-i,]
    z2=y.T[i:,]
    cov=np.multiply(z1,z2)
    u0ui[i-1]=cov.mean(axis=0)
um=np.outer(I1,y.mean(axis=1))
um2=np.multiply(um,um)
yy=np.multiply(y,y)
u2m=np.outer(I1,yy.mean(axis=1))
cct=u0ui-um2
mcct=cct.T.mean(axis=0)

scl=myinfo.args['nstep']//myinfo.args['nsmp'] # adjust the plot range
if( scl > rangcct):
    scl=rangcct
dn=np.arange(0,scl+1)*float(myinfo.args['dtime'])*myinfo.args['nsmp']

plt.plot(dn[:scl],mcct[:scl],'r--')
plt.ylabel('Character Corelation Function of\n Potential of Standard') #
plt.xlabel('time [a.u.]') #
plt.savefig('IR_'+fflg+'.png')
plt.show()

print(mcct[:scl])
lst=pd.DataFrame({'cct':mcct[:scl]})
lst.to_csv('.irdat.tmp')



