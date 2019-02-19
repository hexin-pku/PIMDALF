utils=const.f90 random.f90 linalgebra.f90 myobj.f90
model=mes7smors.f90

SRC1=${utils} pisimul.f90 ${model} MES_Models_dia.f90 staging.f90 md_pimd.f90 mespimd.f90
SRC2=${utils} pisimul.f90 ${model} MES_Models_adia.f90 staging.f90 md_pimd.f90 mespimd.f90

EXE=mespimd
LIB=-llapack

dia:
	gfortran $(SRC1) $(LIB) -o $(EXE).dia.run
adia:
	gfortran $(SRC2) $(LIB) -o $(EXE).adia.run
test:
	./mespimd -p s.par -s e -o 0
lin:
	gfortran const.f90 random.f90 linalgebra.f90 latest.f90  $(LIB) -o lintest
clean:
	rm *.mod
