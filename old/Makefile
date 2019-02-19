DEB=base.f90 MES_Models.f90 md_pimd.f90
SRC=base.f90 linalgebra.f90 MES_Models.f90 md_pimd.f90 # debug.f90
MAP=base.f90 linalgebra.f90 MES_Models.f90 map.f90 map_pimd.f90 # debug.f90
EXE=mespimd
LIB=-L../lib -L./lib -llapack -lrefblas
default:
	gfortran $(SRC) $(LIB) -o $(EXE)
map:
	gfortran $(MAP) $(LIB) -o mappimd
debug:
	gfortran $(DEB) $(LIB) -o itest
clean:
	rm *.mod
sweep:
	rm *.mod $(EXE) 
remove:
	rm *.trj *.ana *.frm *.rst
