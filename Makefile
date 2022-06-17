# This Makefile compiles the implementation in this directory.

.POSIX:

CC = c99
CFLAGS = -W -Wall -Wshadow -O2 -fdiagnostics-color=always -DHAWK_RECOVER_CHECK
#CFLAGS = -W -Wall -Wshadow -g -fsanitize=address,undefined -fdiagnostics-color=always -DHAWK_RECOVER_CHECK
LIBS = -lm

OBJ = build/common.o build/codec.o build/fft.o build/ffo.o build/fpr.o build/keygen.o build/ntt.o build/rng.o build/shake.o build/sign.o build/vrfy.o
PROGS = bin/test_forge bin/speed bin/test_codec bin/test_sampler

HEAD = fpr.h inner.h

all: build bin $(PROGS)

build:
	-mkdir -p build
bin:
	-mkdir -p bin

clean:
	-rm -f $(OBJ) build/hawk.o $(PROGS)

# Binaries:
bin/test_forge: tests/test_forge.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/test_forge $(OBJ) tests/test_forge.c $(LIBS)

bin/speed: tests/speed.c build/hawk.o $(OBJ)
	$(CC) $(CFLAGS) -o bin/speed tests/speed.c build/hawk.o $(OBJ) $(LIBS)

bin/test_codec: tests/test_codec.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/test_codec tests/test_codec.c $(OBJ) $(LIBS)

bin/test_sampler: tests/test_sampler.c build/rng.o build/shake.o
	$(CC) $(CFLAGS) -o bin/test_sampler tests/test_sampler.c $(LIBS)

# Object files:
build/common.o: common.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/common.o common.c
build/codec.o: codec.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/codec.o codec.c
build/fft.o: fft.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/fft.o fft.c
build/ffo.o: ffo.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/ffo.o ffo.c
build/fpr.o: fpr.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/fpr.o fpr.c
build/keygen.o: keygen.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/keygen.o keygen.c
build/ntt.o: ntt.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/ntt.o ntt.c
build/rng.o: rng.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/rng.o rng.c
build/shake.o: shake.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/shake.o shake.c
build/sign.o: sign.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/sign.o sign.c
build/vrfy.o: vrfy.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/vrfy.o vrfy.c

build/hawk.o: hawk.c hawk.h $(HEAD)
	$(CC) $(CFLAGS) -c -o build/hawk.o hawk.c

