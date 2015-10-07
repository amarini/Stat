GCC=g++
CXXFLAGS=`root-config --libs --cflags` -O2 -fPIC  -I./ -std=c++11
## to use RooUnfold
lxplus=$(findstring lxplus, $(shell hostname -f) )
###################### determine if you are on lxplus or not
ifeq ($(strip $(lxplus)),)
$(info You are Not on lxplus)
else
$(info You are on lxplus)
endif

CXXFLAGS += -L$(ROOUNFOLD)/  -ggdb
CXXFLAGS += -I$(ROOUNFOLD)/src/
DICTFLAGS+= -I$(ROOUNFOLD)/src/
SOFLAGS=-shared

SRCDIR=src
BINDIR=bin
HPPDIR=interface

SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJ=$(patsubst $(SRCDIR)/%.cpp, $(BINDIR)/%.o , $(SRC)  )
HPPLINKDEF=$(patsubst $(SRCDIR)/%.cpp, ../interface/%.hpp , $(SRC)  )

.PHONY: all
all: libStat.so
	$(info, "--- Full compilation --- ")	
	$(MAKE) libStat.so


# check if CMSSW is defined
ifndef CMSSW_BASE
$(info No CMSSSW !!!!)
else
$(info CMSSW found: $(CMSSW_BASE) )
endif

## if on mac add the -p to the  DICTFLAGS
UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
$(info You are compiling on mac)
DICTFLAGS += -p
else 
$(info Your are on a linux machine)
endif

# check if Combine is present and compiled 

libStat.so: $(OBJ) Dict | $(BINDIR)
	$(GCC) $(CXXFLAGS) $(SOFLAGS) -o $(BINDIR)/$@ $(OBJ) $(BINDIR)/dict.o

$(OBJ) : $(BINDIR)/%.o : $(SRCDIR)/%.cpp interface/%.hpp | $(BINDIR)
	$(GCC) $(CXXFLAGS) -c -o $(BINDIR)/$*.o $<

.PHONY: Dict
Dict: $(BINDIR)/dict.o

$(BINDIR)/dict.o: $(SRC) | $(BINDIR)
	cd $(BINDIR) && rootcint -v4 -f dict.cc -c -I./ -I../ $(DICTFLAGS)  $(HPPLINKDEF)  ../interface/LinkDef.hpp 
	cd $(BINDIR) && $(GCC) -c -o dict.o $(CXXFLAGS) -I../../ -I../ dict.cc

$(BINDIR):
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	-rm $(OBJ)
	-rm $(BINDIR)/dict*
	-rm $(BINDIR)/libStat.so
	-rmdir $(BINDIR)
