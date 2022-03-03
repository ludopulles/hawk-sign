g++ -O2 -mavx2 -pg ../test_sign.cpp ../build/* -o profile
time ./profile
gprof profile | gprof2dot | dot -Tpng -o performance_sign.png
