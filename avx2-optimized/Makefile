# This Makefile compiles the implementation in this directory along with
# the known answer tests generator located in the
# ../../../KAT/generator/ directory. The output is an executable program
# in the build/ subdirectory, whose name starts with 'kat', followed by
# the implementation name (e.g. 'kat512int' for the 'falcon512int'
# implementation). This program, when executed, generates the .req and
# .rsp files in the expected NIST format.

.POSIX:

CC = c99
CFLAGS = -W -Wall -Wshadow -O2 -mavx2
LD = c99
LDFLAGS = 
LIBS = 

OBJ1 = build/codec.o build/common.o build/fft.o build/fpr.o build/keygen.o build/nist.o build/rng.o build/shake.o build/sign.o build/vrfy.o

HEAD1 = api.h fpr.h inner.h

all: build

build:
	-mkdir build

clean:
	-rm $(OBJ1)


a:
	$(CC) $(CFLAGS) -o a.out lilipu.c

b:
	$(CC) -g -Wall -Wshadow -fsanitize=address,undefined -mavx2 -o b.out lilipu.c

cpp:
	g++ $(CFLAGS) -o a.out main.cpp

build/codec.o: codec.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/codec.o codec.c

build/common.o: common.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/common.o common.c

build/fft.o: fft.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/fft.o fft.c

build/fpr.o: fpr.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/fpr.o fpr.c

build/keygen.o: keygen.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/keygen.o keygen.c

build/nist.o: nist.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/nist.o nist.c

build/rng.o: rng.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/rng.o rng.c

build/shake.o: shake.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/shake.o shake.c

build/sign.o: sign.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/sign.o sign.c

build/vrfy.o: vrfy.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o build/vrfy.o vrfy.c
