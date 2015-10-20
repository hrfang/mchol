# compiler
CXX     ?= g++
CC      ?= gcc
CFLAGS  ?= -Wall -Wconversion -O3 -fPIC

# archive
AR      ?= ar
ARFLAGS ?= cr

# the wanted
all: libmchol.a


mtx_io.o: mtx_io.cpp mtx_io.h
	$(CXX) $(CFLAGS) -c mtx_io.cpp

choltool.o: choltool.cpp choltool.h
	$(CXX) $(CFLAGS) -c choltool.cpp

gmw.o: gmw.cpp gmw.h choltool.h
	$(CXX) $(CFLAGS) -c gmw.cpp

se.o: se.cpp se.h choltool.h
	$(CXX) $(CFLAGS) -c se.cpp

OBJS = mtx_io.o choltool.o gmw.o se.o

libmchol.a: $(OBJS)
	$(AR) $(ARFLAGS) libmchol.a $(OBJS)


clean:
	rm -f *~ *.o libmchol.a
