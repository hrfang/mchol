# compiler
CXX     ?= g++
CC      ?= gcc
CFLAGS  ?= -Wall -Wconversion -O3 -fPIC

# mchol library
MCHOL_DIR = ../source
MCHOL_LIB = $(MCHOL_DIR)/libmchol.a

# the wanted
all: gmw se


# for gmw
gmw_driver.o: gmw_driver.cpp
	$(CXX) $(CFLAGS) -I$(MCHOL_DIR) -c gmw_driver.cpp

gmw: gmw_driver.o $(MCHOL_DIR)/libmchol.a
	$(CXX) -o gmw gmw_driver.o $(MCHOL_DIR)/libmchol.a

# for se
se_driver.o: se_driver.cpp
	$(CXX) $(CFLAGS) -I$(MCHOL_DIR) -c se_driver.cpp

se: se_driver.o $(MCHOL_DIR)/libmchol.a
	$(CXX) -o se se_driver.o $(MCHOL_LIB)

# for the mchol library
$(MCHOL_LIB): $(MCHOL_DIR)/*.cpp $(MCHOL_DIR)/*.h
	make -C $(MCHOL_DIR) CFLAGS='$(CFLAGS)' CXX='$(CXX)';


clean:
	rm -f *~ *.o gmw se
