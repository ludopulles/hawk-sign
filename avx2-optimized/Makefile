# This Makefile compiles the implementation in this directory.

.POSIX:

CC = c99
CFLAGS = -W -Wall -Wshadow -O2 -mavx2
LD = c99
LDFLAGS = 
LIBS = -lpthread
CPP = g++
CPPFLAGS = --std=c++11 $(CFLAGS)

OBJ = build/fft.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o

HEAD = fpr.h inner.h

all: build build/main build/test_keygen build/test_sign

build:
	-mkdir build

clean:
	-rm -f $(OBJ) build/main.o build/main build/test_keygen build/test_keygen.o


build/main: main.c build/main.o $(OBJ)
	$(CC) $(CFLAGS) -o build/main build/main.o $(OBJ) $(LIBS)

build/test_keygen: test_keygen.c build/test_keygen.o $(OBJ)
	$(CC) $(CFLAGS) -o build/test_keygen build/test_keygen.o $(filter-out build/keygen.o,$(OBJ)) $(LIBS)

build/test_sign: test_sign.cpp build/test_sign.o $(OBJ)
	$(CPP) $(CPPFLAGS) -o build/test_sign build/test_sign.o $(OBJ) $(LIBS)


build/fft.o: fft.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/fft.o fft.c
build/fpr.o: fpr.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/fpr.o fpr.c
build/keygen.o: keygen.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/keygen.o keygen.c
build/rng.o: rng.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/rng.o rng.c
build/sampler.o: sampler.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/sampler.o sampler.c
build/shake.o: shake.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/shake.o shake.c
build/sign.o: sign.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/sign.o sign.c
build/vrfy.o: vrfy.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/vrfy.o vrfy.c


build/main.o: main.c
	$(CC) $(CFLAGS) -c -o build/main.o main.c

build/test_keygen.o: test_keygen.c
	$(CC) $(CFLAGS) -c -o build/test_keygen.o test_keygen.c

build/test_sign.o: test_sign.cpp
	$(CPP) $(CPPFLAGS) -c -o build/test_sign.o test_sign.cpp
