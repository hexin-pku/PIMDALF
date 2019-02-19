FC = gfortran

COMMON_DIR = common
MODEL_DIR = model
BIN_DIR = bin

.PHONY:all
all: build_common build_model build_pimd build_mespimd
	cp bin/* ./

build_common:
	$(MAKE) -C common
	
build_model:
	$(MAKE) -C model
	
build_pimd: 
	$(MAKE) -C pimd

build_mespimd:
	$(MAKE) -C mespimd
	
.PHONY:clean cleanobj
clean:
	$(MAKE) -C common clean
	$(MAKE) -C model clean
	$(MAKE) -C pimd clean
	$(MAKE) -C mespimd clean
	-rm *.mod *.o
cleanobj:
	-rm *.mod
	-rm mespimd.dia.run mespimd.adia.run bin/*
