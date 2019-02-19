#!/usr/bin/python3
# Filename: mespro.py

import numpy as np
import pandas as pd

a = pd.read_csv('list.dat', header=None, sep='\s+')
a = a.values.T

print(np.mean(a,axis=1))
print(np.std(a,axis=1)/np.mean(a,axis=1))
