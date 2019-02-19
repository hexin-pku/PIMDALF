# MES-PIMD
Multi-Eelctronic States Path Integral Molecule Dynamics (2019)  

## Compile  
in the root directory, type `make`  
mespimd.dia.run:  is the diabatic version for mes-pimd  
mespimd.adia.run: is the adiabatic version for mes-pimd  

## Usage  
```
    mespimd.xxx.run -p parameter\_file \
    -s [0|1|e|r|f]  
    -f restart\_configuration\_file
    -o output\_name
```

the parameter file with format like  
```
    # parmN, parmNs, parmdt
    5000   1  100
    # beta, bead
    4511.066    16
    # thermo_gamma
    0.001
    # scheme\_flg, thermo\_flg
    0   1
```

## Scripts  
some bash and python scripts can be found in `example\_task`,
can give statistics result of the analyzation files.  

## Acknowledgement  
Thank Prof. Liu & Xinzijian Liu and other fellows 

# to do  
[ ] Add more thermostat (NHC)
[ ] Add normal mode transform
[ ] Non-adiabatic calculation
[ ] Else  
