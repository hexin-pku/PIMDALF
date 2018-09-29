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
fflg=str(sys.argv[1]) #the filename to plot
clmn=int(sys.argv[2]) #the column to plot

myinfo = lf.Load('info.now')
myinfo.add('put.rc')
print(myinfo.args)
os.system('\nhead -n 1 %s'%fflg)

a=pd.read_csv(fflg,sep='\s+',header=None)
statname='null'
if not setinfo.is_number(a.values[0,0]):
    a=pd.read_csv(fflg,sep='\s+')
    statname=str(a.columns.values[clmn])

print('statistics on ',statname)

mtx=a.values.T
x=mtx[0,]
if(clmn>0):
    y=mtx[clmn:clmn+1,]
elif(clmn==0):
    y=mtx[1:,]
else:
    pass


rangcct=3000 # the max length for stat
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
cct=(u0ui-um2)/(u2m-um2)
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

plt.title(r'quartic poential, $U=\frac{1}{4}x^4$')
plt.ylabel('character correlation fuction of x') #
plt.xlabel('time [a.u.]') #
txt = plt.text(7.5,0.85,'bead = %d\n'%myinfo.args['bead']
+'dt = %.1e\n'%myinfo.args['dtime']+'temp = %.2e'%myinfo.args['temp'] #+ 'gamma = %.1e\n'%myinfo.args['gammaAD']
#+'\n\n\nmean=%f'%ymean
, bbox=dict(facecolor='yellow', alpha=0.2))
txt.xycoords = 'figure pixels'


plt.savefig('cct_'+fflg+str(statname)+'.png')
plt.show()

lst=pd.DataFrame({'cct':mcct})
lst.to_csv(fflg+str(statname)+'.cctdat.tmp')



