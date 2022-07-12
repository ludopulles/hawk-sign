gcc -O2 -mavx2 -pg ../tests/test_keygen.c ../build/common.o ../build/codec.o ../build/fft.o ../build/ffo.o ../build/fpr.o ../build/ntt.o ../build/rng.o ../build/shake.o ../build/sign.o ../build/vrfy.o -o keygen -lm
time ./keygen
gprof keygen | gprof2dot | dot -Tpng -o performance_keygen.png
