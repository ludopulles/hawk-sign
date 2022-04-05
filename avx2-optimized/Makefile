# This Makefile compiles the implementation in this directory.
.POSIX:

CC = c99
#CFLAGS = -W -Wall -Wshadow -O3 -march=native
CFLAGS = -W -Wall -Wshadow -O2 -mavx2 -fdiagnostics-color=always
#CFLAGS = -W -Wall -Wshadow -g -fsanitize=address,undefined
LIBS = -lpthread -lm
CXX = g++
CXXFLAGS = --std=c++17 $(CFLAGS)

OBJ = build/common.o build/compress.o build/fft.o build/ffo.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o
PROGS = bin/main bin/test_keygen bin/test_sign bin/config_fg bin/config_FG bin/test_babai bin/test_enum bin/test_compress_sig bin/speed

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
bin/main: main.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/main $(OBJ) main.c $(LIBS)
bin/test_keygen: test_keygen.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/test_keygen $(filter-out build/keygen.o,$(OBJ)) test_keygen.c $(LIBS)
bin/config_fg: config_fg.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/config_fg $(filter-out build/keygen.o,$(OBJ)) config_fg.c $(LIBS)
bin/config_FG: config_FG.c $(OBJ)
	$(CC) $(CFLAGS) -o bin/config_FG $(filter-out build/keygen.o,$(OBJ)) config_FG.c $(LIBS)
bin/speed: speed.c build/hawk.o $(OBJ)
	$(CC) $(CFLAGS) -o bin/speed speed.c build/hawk.o $(OBJ) $(LIBS)

# C++
bin/test_sign: test_sign.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/test_sign $(OBJ) test_sign.cpp $(LIBS)
bin/test_babai: test_babai.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/test_babai $(OBJ) test_babai.cpp $(LIBS)
bin/test_enum: test_enum.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/test_enum $(OBJ) test_enum.cpp $(LIBS)
bin/test_compress_sig: test_compress_sig.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o bin/test_compress_sig $(OBJ) test_compress_sig.cpp $(LIBS)

# Object files:
build/common.o: common.c $(HEAD)
	$(CC) $(CFLAGS) -c -o build/common.o common.c
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

build/hawk.o: hawk.c hawk.h $(HEAD)
	$(CC) $(CFLAGS) -c -o build/hawk.o hawk.c

