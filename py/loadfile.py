#!/usr/bin/python3
# Filename: loadfile.py

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

def is_integer(s):
    try:
        int(s)
        return False
    except ValueError:
        return True
        # for you can't int() a integer representation


class Load:
    args={}
    def __init__(self,filename):
        a=pd.read_csv(filename,sep='\s+',header=None)
        print(a)
        for two in a.values:
            if(is_number(two[1])):
                if( float(two[1])%1==0.0 and not is_integer(two[1])):
                    self.args[two[0]] = int(two[1])
                else:
                    self.args[two[0]] = float(two[1])
            else:
                self.args[two[0]] = two[1].replace(',',' ')
    def Add(self,filename):
        a=pd.read_csv(filename,sep='\s+',header=None)
        print(a)
        for two in a.values:
            if two[0] not in self.args:
                if(is_number(two[1])):
                    if( float(two[1])%1==0.0 and not is_integer(two[1])):
                        self.args[two[0]] = int(two[1])
                    else:
                        self.args[two[0]] = float(two[1])
                else:
                	# original blankspace should written is comma
                    self.args[two[0]] = two[1].replace(',',' ')
    
        
