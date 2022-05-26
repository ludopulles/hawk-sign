# This Makefile compiles the implementation in this directory.
.POSIX:

CC = c99
#CFLAGS = -W -Wall -Wshadow -O3 -march=native
CFLAGS = -W -Wall -Wshadow -O2 -mavx2 -fdiagnostics-color=always
#CFLAGS = -W -Wall -Wshadow -g -fsanitize=address,undefined
LIBS = -lpthread -lm
CXX = g++
CXXFLAGS = --std=c++17 $(CFLAGS)

OBJ = build/common.o build/codec.o build/fft.o build/ffo.o build/fpr.o build/keygen.o build/ntt.o build/rng.o build/shake.o build/sign.o build/vrfy.o
PROGS = bin/test_keygen bin/test_sign bin/config_fg bin/config_FG bin/test_babai bin/test_enum bin/signature_size bin/speed bin/test_decompress_sig bin/skpk_size

HEAD = fpr.h inner.h

all: build bin $(PROGS)

build:
	-mkdir -p build
bin:
	-mkdir -p bin

clean:
	-rm -f $(OBJ) build/hawk.o $(PROGS)

# Binaries:

# C
bin/test_keygen: tests/test_keygen.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/test_keygen $(filter-out build/keygen.o,$(OBJ)) tests/test_keygen.c $(LIBS)
bin/config_fg: tests/config_fg.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/config_fg $(filter-out build/keygen.o,$(OBJ)) tests/config_fg.c $(LIBS)
bin/config_FG: tests/config_FG.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/config_FG $(filter-out build/keygen.o,$(OBJ)) tests/config_FG.c $(LIBS)
bin/speed: tests/speed.c build/hawk.o $(OBJ)
	$(CC) $(CFLAGS) -o bin/speed build/hawk.o $(OBJ) tests/speed.c $(LIBS)
bin/skpk_size: tests/skpk_size.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/skpk_size $(OBJ) tests/skpk_size.c $(LIBS)

# C++
bin/test_sign: tests/test_sign.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/test_sign $(OBJ) tests/test_sign.cpp $(LIBS)
bin/test_babai: tests/test_babai.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/test_babai $(OBJ) tests/test_babai.cpp $(LIBS)
bin/test_enum: tests/test_enum.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/test_enum $(OBJ) tests/test_enum.cpp $(LIBS)
bin/signature_size: tests/signature_size.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/signature_size $(OBJ) tests/signature_size.cpp $(LIBS)
bin/test_decompress_sig: tests/test_decompress_sig.cpp $(OBJ)
	$(CXX) $(CFLAGS) -o bin/test_decompress_sig tests/test_decompress_sig.cpp $(OBJ) $(LIBS)

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

