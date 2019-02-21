import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dist = np.zeros(10000)

for i in range(1,101):
    a=pd.read_csv('distribution%d.out'%i,header=None,sep='\s+')
    a=a.values[:,1]
    dist[0:len(a)] = dist[0:len(a)] + a
    
    print('%d'%i)
    
p = np.linspace(10,25,len(a))

plt.plot(p,dist[0:len(p)]/100,'b--',label=r'$10^6$ trajs MQCIVR')
plt.xlabel(r'$P_{final}$ [a.u.]')
plt.ylabel('Distribution')

exact = pd.read_csv('tul1_mqcivr.dat',sep='\s+').values
plt.plot(exact[:,0],exact[:,1],'k-',label='exact')

plt.legend(loc=1)

plt.savefig('00.png')
plt.show()



b = np.zeros((len(p),2))
b[:,0] = p
b[:,1] = dist[0:len(p)]
b=pd.DataFrame(b)
b.to_csv('0dist.dat')


