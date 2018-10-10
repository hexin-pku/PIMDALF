#!/usr/bin/python3
# Filename: histo.py

import numpy as np
# np.seterr(divide='ignore', invalid='ignore')
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import loadfile as lf

# get arguments from terminal
narg=len(sys.argv)-1
fflg=str(sys.argv[1]) # the filename to plot
clmn=int(sys.argv[2]) # the column to plot

# load configuration file
icwd=os.getcwd()
print('current working directory is %s\n'%icwd)
if os.path.isfile('./put.rc'):
    myinfo = lf.Load('put.rc')
else:
    print('Warning: put.rc is need in current dir!')
    exit()
if os.path.isfile('./histo.rc'):
    myinfo.Add('histo.rc')
else:
    print('Warning: histo.rc is need in current dir!')
    exit()
if os.path.isfile('./info.now'):
    myinfo.Add('info.now')
else:
    print('Warning: info.now is need in current dir!')
    exit()

# read the data file
if not os.path.isfile(fflg):
    print("%s doesn't exist!")
    exit()
a=pd.read_csv(fflg,sep='\s+',header=None)
name='Nan'
if not lf.is_number(a.values[0,0]):
    a=pd.read_csv(fflg,sep='\s+')
    name=str(a.columns.values[clmn])
fflg=(fflg.split('.',1))[0]

# processing the data, histogram
mtx=a.values.T
if(clmn>0):
    y=mtx[clmn:clmn+1,]
elif(clmn==0):
    y_all=mtx[1:,]
    t=y_all.reshape(1,-1)
else:
    pass
y=y[0]

if 'histo-x1' in myinfo.args and 'histo-x1' in myinfo.args :
    x1=int(myinfo.args['histo-x1'])
    x2=int(myinfo.args['histo-x2'])
else :
    x1 = y.min()
    x2 = y.max()
plt.hist(y,bins=200,normed=True,range=(x1,x2), histtype='step')
y1,y2=plt.ylim()
plt.ylabel(myinfo.args['histo-ylabel'])
plt.xlabel(name)
plt.title(myinfo.args['histo-title'])

ymean=y.mean()
plt.axvline(x=ymean,color='r',ls='--')
print('the mean value of '+name+'\n', ymean)

txt = plt.text(x1+0.75*(x2-x1),0.7*(y2-y1),
'beta=%.2e\n'%myinfo.args['beta']
+'bead=%d\n'%myinfo.args['nbead']
+'dt=%.1e\n'%myinfo.args['dtime']
#+'gamma=%.1e\n'%myinfo.args['gammaAD']
+'\nmean=%f'%ymean,
bbox=dict(facecolor='grey', alpha=0.05))
txt.xycoords = 'figure pixels'

plt.savefig(fflg+'_'+name+'.png')

if(narg==3 and sys.argv[3]=='Y'):
	plt.show()
os.system('echo %.10e > .mth.tmp'%ymean)


