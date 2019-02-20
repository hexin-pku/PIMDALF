#!/usr/bin/python3
# Filename: histo.py

import numpy as np
# np.seterr(divide='ignore', invalid='ignore')
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False



# get arguments from terminal
# argv[0]: histo.py
# argv[1]: filename
# argv[2]: column-idx

# atgv[n]: tag values
#--------:
#   xlim : 0:2
#   ylim : 0:1
#   xlabel :

narg=len(sys.argv)-1
fflg=str(sys.argv[1]) # the filename to plot
clmn=int(sys.argv[2]) # the column to plot

# load configuration file
icwd=os.getcwd()
print('current working directory is %s\n'%icwd)



# read the data file
if not os.path.isfile(fflg):
    print("%s doesn't exist!")
    exit()
a=pd.read_csv(fflg,sep='\s+',header=None)
name='Nan'
if not is_number(a.values[0,0]):
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

ymean = y.mean()
print(ymean)

plt.hist(y,bins=200,normed=True, histtype='step')

plt.savefig(fflg+'_'+name+'.png')

if(narg==3 and sys.argv[3]=='Y'):
	plt.show()
os.system('echo %.10e > .mth.tmp'%ymean)


