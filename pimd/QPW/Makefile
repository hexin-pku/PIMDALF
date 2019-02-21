FC = gfortran
FFLAGS = -O2 -llapack

COMMON_DIR = ../../common
PIMD_DIR = ..

INC = -I./${COMMON_DIR} -I./$(PIMD_DIR)

MDL_FILES = QPW_model.f90

MOD_FILES = model.o

COMMON_SRCS = $(wildcard $(COMMON_DIR)/*.f90)
COMMON_OBJS = $(addprefix $(COMMON_DIR)/,$(patsubst %.f90,%.o,$(notdir $(COMMON_SRCS))))

PIMD_SRCS = pisimul.f90 myobj.f90 trans_plus.f90 nhc_plus_pimd.f90
PIMD_OBJS = $(addprefix $(PIMD_DIR)/,$(patsubst %.f90,%.o,$(notdir $(PIMD_SRCS))))

MAIN_SRC = ../md_pimd.f90 ../main.f90

all:model pimd_qpw

model: $(MDL_FILES)
	echo $@
	$(FC) $(INC) -c $(FFLAGS) $(MDL_FILES) -o $(MOD_FILES)
	# mv ho_model.mod model.mod

pimd_qpw:model.o $(MAIN_SRC)
	echo $@
	${FC} $(INC) $(COMMON_OBJS) $(PIMD_OBJS) \
	model.o $(MAIN_SRC) $(FFLAGS) -o QPW.run

.PHONY:clean cleanobj
clean:
	-rm *.mod *.o
cleanobj:
	-rm *.mod

