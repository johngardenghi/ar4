FC := gfortran
FFLAGS := -O3

ALGENCAN := $(CURDIR)/algencan
BIN := $(CURDIR)/bin
DRV := $(CURDIR)/drv
INC := $(CURDIR)/inc
LIB := $(CURDIR)/lib
OBJ := $(CURDIR)/obj
MGH := $(CURDIR)/mgh
SRC := $(CURDIR)/src

# Please add -lhsl if you included HSL in algencan/sources/hsl folder
ALGENCAN_FLAGS := -lalgencan -lblas -llapack -lhsl

export

all: mghr packr

libr:
	$(MAKE) -C $(CURDIR)/obj all

mghr: libr
	mkdir -p $(BIN)
	$(MAKE) -C $(ALGENCAN)
	$(MAKE) -C $(MGH)
	$(FC) $(FFLAGS) -I$(MGH) -I$(INC) -L$(ALGENCAN)/lib -L$(LIB) -o $(BIN)/mgh $(DRV)/mgh.f08 $(MGH)/mgh.o $(MGH)/set_precision.o -larp $(ALGENCAN_FLAGS)

packr: libr
	mkdir -p $(BIN)
	$(MAKE) -C $(ALGENCAN)
	$(FC) $(FFLAGS) -I$(INC) -L$(LIB) -L$(ALGENCAN)/lib -o $(BIN)/pack $(DRV)/pack.f08 -larp $(ALGENCAN_FLAGS)

clean:
	$(MAKE) -C $(CURDIR)/obj clean
	rm -Rf $(BIN)/mgh $(BIN)/pack

cleanall: clean
	$(MAKE) -C $(ALGENCAN) distclean
	$(MAKE) -C $(MGH) clean
