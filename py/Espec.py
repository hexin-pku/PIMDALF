import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loadfile as lf

# load argument file
myload = lf.Load('map.rc')

# optimized arguments
timemax = int(1600/myload.args['mapdtime'])
Ndivide = 1000
sigma = 600

a = pd.read_csv('pj2.dat',sep='\s+',header=None)
a=a.values.T
if(timemax > len(a[0])):
    timemax = len(a[0])

time = a[0,:timemax]*myload.args['mapdtime']
# and add Gauss smooth
real = a[1,:timemax]*np.exp(-time*time/(2*sigma**2)) # real part of <phi|U|Psi>
imag = a[2,:timemax]*np.exp(-time*time/(2*sigma**2)) # imaginary part of <phi|U|Psi>

E = np.linspace(-0.03,-0.01,Ndivide+1)


Et = np.outer(E,time)
expiEt_r = np.cos(Et)
expiEt_i = np.sin(Et)

realpart = ( np.dot(expiEt_r,real) - np.dot(expiEt_i,imag) ) * myload.args['mapdtime']
imagpart = ( np.dot(expiEt_r,imag) + np.dot(expiEt_i,real) ) * myload.args['mapdtime']

plt.xlabel('E [a.u.]')
plt.ylabel(' [Im Arbitrary] ')
plt.plot(E,imagpart,'b-',linewidth=0.5)
plt.savefig('espec.png')
plt.show()


