#!/usr/bin/python3
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) <= 1:
    print("Usage: python plot.py datafile")
    sys.exit(1)

#logplot = True
logplot = False

a = []
b = []
with open(sys.argv[1]) as f:
    for line in f.readlines():
        cstQ, fails = line.split()
        a.append(float(cstQ))
        b.append(int(fails))

def estimate(cstQ00):
    # n = 512, sigma_sig = 1.292
    sigma_ei = math.sqrt(cstQ00) * 1.292

    # Perhaps we need a union bound when erf ~ 1
    return 1.0 - math.erf(0.5 / sigma_ei / math.sqrt(2))**512

# the histogram of the data
ys = np.array(b) / 1024
if logplot:
    ys = list(map(math.log, ys))
plt.scatter(a, ys)

sorteda = list(filter(lambda x: estimate(x) > 1e-20, sorted(a)))
# find renormalization factor for discrete gaussian
ys = list(map(estimate, sorteda))
if logplot:
    ys = list(map(math.log, ys))

plt.xlabel('cst(1/Q_{00})')
plt.ylabel('Failure rate (log scale)')
plt.title('Failure probability for a given public key')

# plt.plot(sorteda, ys)
plt.semilogy(sorteda, ys)

plt.show()

