FC = gfortran
FFLAGS = -O2 -llapack

MDL_FILES = const.f90 random.f90 linalgebra.f90
MOD_FILES = ${MDL_FILES:.f90=.o}

all:$(MOD_FILES)

$(MOD_FILES):%.o:%.f90
	echo $@
	$(FC) -c $(FFLAGS) $^ -o $@


.PHONY:clean cleanobj
clean:
	-rm *.mod *.o
cleanobj:
	-rm *.mod

