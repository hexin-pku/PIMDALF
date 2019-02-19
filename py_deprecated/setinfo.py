#!/usr/bin/python3
# Filename: setinfo.py

import pandas as pd

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

class Infolog:
    dk={}
    def __init__(self,filename):
        a=pd.read_csv(filename,sep='\s+',header=None)
        print(a)
        for two in a.values:
            if(is_number(two[1])):
                if(float(two[1])%1 == 0.0):
                    self.dk[two[0]] = int(two[1])
                else:
                    self.dk[two[0]] = float(two[1])
            else:
                self.dk[two[0]] = two[1].replace(',',' ')
    def add(self,filename):
        a=pd.read_csv(filename,sep='\s+',header=None)
        print(a)
        for two in a.values:
            if two[0] not in self.dk:
                if(is_number(two[1])):
                    if(float(two[1])%1==0.0):
                        self.dk[two[0]] = int(two[1])
                    else:
                        self.dk[two[0]] = float(two[1])
                else:
                    self.dk[two[0]] = two[1].replace(',',' ')
    
        
