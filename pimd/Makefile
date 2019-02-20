FC = gfortran
FFLAGS = -O2 -llapack

COMMON_DIR = ../common

INC = -I./${COMMON_DIR}

MDL_FILES = pisimul.f90 myobj.f90 trans_plus.f90 nhc_plus.f90 model.f90 md_pimd.f90

MOD_FILES = ${MDL_FILES:.f90=.o}

COMMON_SRCS = $(wildcard $(COMMON_DIR)/*.f90)
COMMON_OBJS = $(addprefix $(COMMON_DIR)/,$(patsubst %.f90,%.o,$(notdir $(COMMON_SRCS))))

.PHONY:all
all:$(MOD_FILES) pimd

$(MOD_FILES):%.o:%.f90
	echo $@
	$(FC) $(INC) -c $(FFLAGS) $^ -o $@

.PHONY:pimd
pimd:pimd_ho pimd_qpw

pimd_ho:
	$(MAKE) -C HO
	
pimd_qpw:
	$(MAKE) -C QPW

pimd_h2o:
	$(MAKE) -C H2O

.PHONY:clean cleanobj
clean:
	$(MAKE) -C HO clean
	$(MAKE) -C QPW clean
	#$(MAKE) -C H2O clean
	-rm -rf *.mod *.o
cleanobj:
	-rm *.mod

