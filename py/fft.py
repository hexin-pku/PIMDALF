#!/usr/bin/python3
# Filename: fft.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loadfile as lf

myinfo = lf.Load('_plot.in')
myinfo.add('_put.in')

val=pd.read_csv('.irdat.tmp')
t=val.values.T[0]
cct=val.values.T[1]
t=np.arange(0,len(cct))
print(t,cct)

N=len(cct)
T=len(t)

cw = np.fft.fft(cct)
freq = np.fft.fftfreq(N,d=myinfo.args['dtime'])
print(len(freq),len(cw))
plt.plot(t,cct,'b--',label='cct plot')
#plt.plot(freq,np.abs(cw),'r-',label='fft plot')
plt.legend()
plt.show()
plt.savefig('IR.png')
#f=np.arange(N)/T
#tmp=(1-np.exp(-8*f[0:cut]))/(8*f[0:cut])*cw[0:cut]
#tmp[0]=cw[0]
#print(cw)
#
#kw=np.arange(N)
#kw[0:cut]=tmp.copy()
#kw[cut1:N]=tmp[::-1]

#exit()

#kcct=np.fft.ifft(kw)

#plt.plot(t,cct,'r-',t[0:cut],kcct[0:cut],'b--')
#plt.show()
