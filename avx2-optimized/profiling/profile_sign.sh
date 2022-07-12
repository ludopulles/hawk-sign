gcc -pg -W -Wall -Wshadow -mavx2 ../tests/speed.c ../build/* -o sign
time ./sign $1
gprof sign | gprof2dot | dot -Tpng -o performance_sign.png
