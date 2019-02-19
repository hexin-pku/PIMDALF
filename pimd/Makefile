FC = gfortran
FFLAGS = -O2 -llapack

COMMON_DIR = ../common

INC = -I./${COMMON_DIR}

MDL_FILES = pisimul.f90 myobj.f90 staging.f90 def_model.f90 md_pimd.f90

MOD_FILES = ${MDL_FILES:.f90=.o}

COMMON_SRCS = $(wildcard $(COMMON_DIR)/*.f90)
COMMON_OBJS = $(addprefix $(COMMON_DIR)/,$(patsubst %.f90,%.o,$(notdir $(COMMON_SRCS))))

all:$(MOD_FILES)

$(MOD_FILES):%.o:%.f90
	echo $@
	$(FC) $(INC) -c $(FFLAGS) $^ -o $@
	
.PHONY:clean cleanobj
clean:
	-rm *.mod *.o
cleanobj:
	-rm *.mod

