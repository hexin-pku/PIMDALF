FC = gfortran
FFLAGS = -O2 -llapack -fcheck=all -Wall

COMMON_DIR = ../../common
PIMD_DIR = ..

INC = -I./${COMMON_DIR} -I./$(PIMD_DIR)

MDL_FILES = h2opes.f90

MOD_FILES = model.o

COMMON_SRCS = $(wildcard $(COMMON_DIR)/*.f90)
COMMON_OBJS = $(addprefix $(COMMON_DIR)/,$(patsubst %.f90,%.o,$(notdir $(COMMON_SRCS))))

PIMD_SRCS = pisimul.f90 myobj.f90 trans_plus.f90 nhc_plus_pimd.f90
PIMD_OBJS = $(addprefix $(PIMD_DIR)/,$(patsubst %.f90,%.o,$(notdir $(PIMD_SRCS))))

MAIN_SRC = h2o_pimd.f90 ../main.f90

all:model pimd_h2o

model: $(MDL_FILES)
	echo $@
	$(FC) $(INC) -c $(FFLAGS) $(MDL_FILES) -o $(MOD_FILES)
	# mv ho_model.mod model.mod

pimd_h2o:model.o $(MAIN_SRC)
	echo $@
	${FC}  $(INC) $(COMMON_OBJS) $(PIMD_OBJS) \
	model.o $(MAIN_SRC) $(FFLAGS) -o H2O.run

test:
	echo $@
	$(FC) $(FFlAGS) h2opes.f90 testpot.f90 -o _testpot
	
.PHONY:clean cleanobj
clean:
	-rm *.mod *.o
cleanobj:
	-rm *.mod

