g++ -pg -W -Wall -Wshadow -O2 -mavx2 ../test_sign.cpp ../build/* -o sign
time ./sign
gprof sign | gprof2dot | dot -Tpng -o performance_sign.png
