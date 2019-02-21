FC = gfortran
FFLAGS = -O2 -llapack

# revise the directory path of lapack
LAPACK_DIR = ../lib

COMMON_DIR = ../common
MODEL_DIR = ../model
PIMD_DIR = ../pimd
BIN_DIR  = ../bin

INC = -I./${COMMON_DIR} -I./$(MODEL_DIR) -I./$(PIMD_DIR)

COMMON_SRCS = $(wildcard $(COMMON_DIR)/*.f90)
COMMON_OBJS = $(addprefix $(COMMON_DIR)/,$(patsubst %.f90,%.o,$(notdir $(COMMON_SRCS))))

PIMD_SRCS = pisimul.f90 myobj.f90 trans_plus.f90 nhc_plus_pimd.f90
PIMD_OBJS = $(addprefix $(PIMD_DIR)/,$(patsubst %.f90,%.o,$(notdir $(PIMD_SRCS))))

SRC_DIA = MES_Models_dia.f90 mes_pimd.f90 mes_main.f90
SRC_ADIA = MES_Models_adia.f90 mes_pimd.f90 mes_main.f90

.PHONY:mes	
mes: mes_dia mes_adia

mes_dia:
	${FC}  $(INC) $(COMMON_OBJS) \
	$(MODEL_DIR)/mes7smors.o $(PIMD_OBJS) $(SRC_DIA) \
	$(FFLAGS) -o $(BIN_DIR)/mespimd.dia.run

mes_adia:
	${FC}  $(INC) $(COMMON_OBJS) \
	$(MODEL_DIR)/mes7smors.o $(PIMD_OBJS) $(SRC_ADIA) \
	-L$(LAPACK_DIR) $(FFLAGS) -o $(BIN_DIR)/mespimd.adia.run
	
.PHONY:clean cleanobj
clean:
	-rm *.mod *.o
cleanobj:
	-rm *.mod
