# This Makefile compiles the implementation in this directory.

.POSIX:

CC = c99
CFLAGS = -W -Wall -Wshadow -O2
LD = c99
LDFLAGS = 
LIBS = 

OBJ = build/fft.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o

HEAD = fpr.h inner.h

all: build build/main

build:
	-mkdir build

clean:
	-rm -f $(OBJ) build/main.o build/main


build/main: main.c build/main.o $(OBJ)
	$(CC) $(CFLAGS) -o build/main build/main.o $(OBJ) $(LIBS)


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


build/main.o: main.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/main.o main.c
