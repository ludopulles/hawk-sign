# This Makefile compiles the implementation in this directory.
.POSIX:

LIBS = -lpthread -lm
CXX = g++
CXXFLAGS = --std=c++17 -W -Wall -Wshadow -O2 -mavx2 -fdiagnostics-color=always -DHAWK_RECOVER_CHECK

OBJ = ../build/common.o ../build/codec.o ../build/fft.o ../build/ffo.o ../build/fpr.o ../build/keygen.o ../build/ntt.o ../build/rng.o ../build/shake.o ../build/sign.o ../build/vrfy.o
PROGS = main

HEAD = fpr.h inner.h

all: $(PROGS)

clean:
	-rm -f main

# Binaries:

# C++
main: main.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o main $(OBJ) main.cpp $(LIBS)
