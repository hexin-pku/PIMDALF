bindir = ../bin/
srcdir = ./src
FC = mpif90 -O2 -g
LIB = -llapack -lblas
.SUFFIXES: .f90 .o

#--------------------------------------------------
MOD = MonteCarlo.f90 \
      traj.f90 \
      supply.f90 \
      params.o 

SRC = timecorr.f90

OBJ = $(MOD) $(SRC)

# rule to make a .o file from a .f90 file
.f90.o:
	$(FC) -c $<
.f.o:
	$(FC) -c $<

dyn.x  :  $(OBJ)
	$(FC) $(OPTION) -o $@ $(OBJ) $(LIB)

clean   :
	-rm *.o 
	-rm *.mod
	-rm fort.*
	-rm mqcivr.*
