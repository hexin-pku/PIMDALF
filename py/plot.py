#!/usr/bin/env python
# coding=utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import loadfile as lf
import os

# get arguments from terminal
narg=len(sys.argv)-1
fflg=str(sys.argv[1]) # the filename to plot
clmn=int(sys.argv[2]) # the column to plot

# load the configuration file
if os.path.isfile('./put.rc'):
    myinfo = lf.Load('put.rc')
else:
    print('Warning: put.rc is need in current dir!')
    exit()
if os.path.isfile('./plot.rc'):
    myinfo.Add('plot.rc')
else:
    print('Warning: plot.rc is need in current dir!')
    exit()
    
# read the data file
a=pd.read_csv(fflg,sep='\s+',header=None)
name='Nan'
if not lf.is_number(a.values[0,0]):
    a=pd.read_csv(fflg,sep='\s+')
    name=str(a.columns.values[clmn])
fflg=(fflg.split('.',1))[0]

# processing
xs=a.values[:,0]*myinfo.args['dtime']
if(clmn==0):
    ys=a.values[:,1:]
else:
    ys=a.values[:,clmn]

if(clmn==0):
    nclmn = len(ys[0,:])
    for i in range(nclmn):
        plt.plot(xs,ys[:,i],ls='--',c=(i/(nclmn-1),0,1-i/(nclmn-1)), label="%d-th state"%(i+1))
    plt.xlabel(myinfo.args['plot-xlabel'])
    plt.ylabel(myinfo.args['plot-ylabel-all'])
    plt.legend(loc=1)
else:
    plt.plot(xs,ys,'r--')
    plt.xlabel(myinfo.args['plot-xlabel'])
    plt.ylabel(myinfo.args['plot-ylabel'])
plt.title(myinfo.args['plot-title'])
plt.savefig(fflg+'_c'+str(clmn)+'.png')

if(narg==3 and sys.argv[3]=='Y'):
	plt.show()

