FC = gfortran
FFLAGS = -O2 -llapack

COMMON_DIR = ../common

INC = -I./${COMMON_DIR}

MDL_FILES = simul.f90 myobj.f90 nhc_plus_md.f90 model.f90 md.f90

MOD_FILES = ${MDL_FILES:.f90=.o}

COMMON_SRCS = $(wildcard $(COMMON_DIR)/*.f90)
COMMON_OBJS = $(addprefix $(COMMON_DIR)/,$(patsubst %.f90,%.o,$(notdir $(COMMON_SRCS))))

.PHONY:all
all:$(MOD_FILES) md

$(MOD_FILES):%.o:%.f90
	echo $@
	$(FC) $(INC) -c $(FFLAGS) $^ -o $@

.PHONY:md

md_ho:
	echo $@
	${FC}  $(INC) $(COMMON_OBJS) $(MOD_FIELS) $(FFLAGS) -o mdrun

.PHONY:clean cleanobj
clean:
	-rm -rf *.mod *.o
cleanobj:
	-rm *.mod

