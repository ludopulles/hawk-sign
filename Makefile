# This Makefile compiles the implementation in this directory.

.POSIX:

CC = c99
CFLAGS = -W -Wall -Wshadow -O2
LIBS = -lm

OBJ = build/compress.o build/fft.o build/ffo.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o
PROGS = bin/main bin/generate

HEAD = fpr.h inner.h

all: build bin $(PROGS)

build:
	-mkdir -p build
bin:
	-mkdir -p bin

clean:
	-rm -f $(OBJ) $(PROGS)

# Binaries:
bin/main: main.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/main $(OBJ) main.c $(LIBS)

bin/generate: generate.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/generate $(OBJ) generate.c $(LIBS)

# Object files:
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
