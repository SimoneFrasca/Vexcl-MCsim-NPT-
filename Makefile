CXX=g++
CC=gcc
CXXFLAGS=-Wall -std=c++17 -O3
CFLAGS=-Wall -O3
LDFLAGS=-lm

all: vexcl.exe mcsim.exe

vexcl.exe: vexcl.o
	$(CXX) $(LDFLAGS) vexcl.o -o vexcl.exe 

mcsim.exe: mcsim.o
	$(CXX) $(LDFLAGS) mcsim.o -o mcsim.exe 

vexcl.o: vexcl.cpp mcsim.hpp  # <--- dipendenze
	$(CXX) $(CXXFLAGS) -c -o vexcl.o vexcl.cpp 

mcsim.o: mcsim.cpp mcsim.hpp  # <--- dipendenze
	$(CXX) $(CXXFLAGS) -c -o mcsim.o mcsim.cpp 

clean:
	rm -fr *.o vexcl.exe mcsim.exe
