# This Makefile compiles the implementation in this directory.

.POSIX:

CC = c99
CFLAGS = -W -Wall -Wshadow -O2
LD = c99
LDFLAGS = -W -Wall -Wshadow -O2 -lm
LIBS = 

OBJ = build/compress.o build/fft.o build/ffo.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o

HEAD = fpr.h inner.h

all: build build/main build/generate

build:
	-mkdir build

clean:
	-rm -f $(OBJ) build/main.o build/main build/generate.o build/generate


build/main: build/main.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/main build/main.o $(OBJ) $(LIBS)

build/generate: build/generate.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/generate build/generate.o $(OBJ) $(LIBS)

build/compress.o: compress.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/compress.o compress.c
build/fft.o: fft.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/fft.o fft.c
build/ffo.o: ffo.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/ffo.o ffo.c
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

build/generate.o: generate.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/generate.o generate.c
