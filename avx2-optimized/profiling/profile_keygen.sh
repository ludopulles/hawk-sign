gcc -O2 -mavx2 -pg ../test_keygen.c ../build/compress.o ../build/fft.o ../build/ffo.o ../build/fpr.o ../build/rng.o ../build/sampler.o ../build/shake.o ../build/sign.o ../build/vrfy.o -o profile -lm
time ./profile
gprof profile | gprof2dot | dot -Tpng -o performance_keygen.png
