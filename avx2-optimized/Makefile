# This Makefile compiles the implementation in this directory.

.POSIX:



CC = c99
CFLAGS = -W -Wall -Wshadow -O3 -march=native
# CFLAGS = -W -Wall -Wshadow -O2 -mavx2
# CFLAGS = -W -Wall -Wshadow -g -fsanitize=address,undefined
LD = $(CC)
LDFLAGS = $(CFLAGS)
LIBS = -lpthread -lm
CPP = g++
CPPFLAGS = --std=c++11 $(CFLAGS)

OBJ = build/fft.o build/ffo.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o

HEAD = fpr.h inner.h

all: build build/main build/test_keygen build/test_sign build/config_fg build/config_FG build/test_babai build/test_enum build/test_forge

build:
	-mkdir build

clean:
	-rm -f $(OBJ) build/main.o build/main build/test_keygen build/test_keygen.o build/test_sign build/test_sign.o build/config_fg build/config_fg.o build/config_FG build/config_FG.o build/test_babai build/test_babai.o build/test_enum build/test_enum.o build/test_forge build/test_forge.o


build/main: main.c build/main.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/main build/main.o $(OBJ) $(LIBS)

build/test_keygen: test_keygen.c build/test_keygen.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/test_keygen build/test_keygen.o $(filter-out build/keygen.o,$(OBJ)) $(LIBS)

build/test_sign: test_sign.cpp build/test_sign.o $(OBJ)
	$(CPP) $(CPPFLAGS) -o build/test_sign build/test_sign.o $(OBJ) $(LIBS)

build/config_fg: config_fg.c build/config_fg.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/config_fg build/config_fg.o $(filter-out build/keygen.o,$(OBJ)) $(LIBS)

build/config_FG: config_FG.c build/config_FG.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/config_FG build/config_FG.o $(filter-out build/keygen.o,$(OBJ)) $(LIBS)

build/test_babai: test_babai.cpp build/test_babai.o $(OBJ)
	$(CPP) $(CPPFLAGS) -o build/test_babai build/test_babai.o $(OBJ) $(LIBS)

build/test_enum: test_enum.cpp build/test_enum.o $(OBJ)
	$(CPP) $(CPPFLAGS) -o build/test_enum build/test_enum.o $(OBJ) $(LIBS)

build/test_forge: test_forge.cpp build/test_forge.o $(OBJ)
	$(CPP) $(CPPFLAGS) -o build/test_forge build/test_forge.o $(OBJ) $(LIBS)

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


build/main.o: main.c
	$(CC) $(CFLAGS) -c -o build/main.o main.c

build/test_keygen.o: test_keygen.c
	$(CC) $(CFLAGS) -c -o build/test_keygen.o test_keygen.c

build/test_sign.o: test_sign.cpp
	$(CPP) $(CPPFLAGS) -c -o build/test_sign.o test_sign.cpp

build/config_fg.o: config_fg.c
	$(CC) $(CFLAGS) -c -o build/config_fg.o config_fg.c

build/config_FG.o: config_FG.c
	$(CC) $(CFLAGS) -c -o build/config_FG.o config_FG.c
	
build/test_babai.o: test_babai.cpp
	$(CPP) $(CFLAGS) -c -o build/test_babai.o test_babai.cpp

build/test_enum.o: test_enum.cpp
	$(CPP) $(CFLAGS) -c -o build/test_enum.o test_enum.cpp

build/test_forge.o: test_forge.cpp
	$(CPP) $(CFLAGS) -c -o build/test_forge.o test_forge.cpp
