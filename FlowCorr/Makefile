CXX := $(shell root-config --cxx)
CXX_FLAGS := $(shell root-config --cflags)
CXX_FLAGS += -Wall -Wextra -ggdb
LIBRARIES := $(shell root-config --glibs)

BIN := bin
SRC := src
INCLUDE := include
LIB := lib

CORR2D_D0 := corr2D_trg_d0 
CORR2D_REF := corr2D_trg_ref
C24_D0 := c24_d0
C24_REF := c24_ref

all: $(BIN)/$(CORR2D_D0) $(BIN)/$(CORR2D_REF) $(BIN)/$(C24_D0) $(BIN)/$(C24_REF) $(BIN)/$(MASS_D0)
#all: $(BIN)/$(CORR2D_D0) $(BIN)/$(CORR2D_REF) $(BIN)/$(MASS_D0)

$(BIN)/$(CORR2D_D0): $(SRC)/corr2D_trg_d0.cc $(SRC)/*.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_REF): $(SRC)/corr2D_trg_ref.cc $(SRC)/*.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(C24_D0): $(SRC)/c24_d0.cc $(SRC)/*.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(C24_REF): $(SRC)/c24_ref.cc $(SRC)/*.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*
