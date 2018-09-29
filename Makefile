SRC=base.f90 MES_Models.f90 md_pimd.f90
EXE=mespimd

default:
	gfortran $(SRC) -o $(EXE)
clean:
	rm *.mod
sweep:
	rm *.mod $(EXE) 
remove:
	rm *.trj *.ana *.frm *.rst
