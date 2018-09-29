import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loadfile as lf
import sys

istate = int(sys.argv[1])

# load argument file
myload = lf.Load('map.rc')

# optimized arguments
timemax = int(1600/myload.args['mapdtime'])
Ndivide = 1000
sigma = 600

a = pd.read_csv('pop.dat',sep='\s+',header=None)
a=a.values.T
if(timemax > len(a[0])):
    timemax = len(a[0])

time = a[0,:timemax]*myload.args['mapdtime']
# and add Gauss smooth
popu = a[istate,:timemax]
popmean = popu.mean()
term = ( popu - popmean ) * np.exp( -time*time/(2*sigma**2) )


W = np.linspace(-0.5,2.5,Ndivide+1)


Wt = np.outer(W,time)
expiWt_r = np.cos(Wt)
expiWt_i = np.sin(Wt)

realpart = np.dot(expiWt_r,term)
imagpart = np.dot(expiWt_i,term)

I = realpart**2 + imagpart**2

plt.plot(W,realpart,'b--',W,imagpart,'g--')
plt.xlabel('frequency [a.u.]')
plt.ylabel('[Arbitrary]')
plt.savefig('freq.png')
plt.show()


