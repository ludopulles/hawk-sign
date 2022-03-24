#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt

xs = []
with open("error_correlation") as f:
    for line in f.readlines():
        xs.append(float(line))

# 40335 / 400000 failed
# Probability: 0.10083750000000000000

plt.hist(xs, 1000, density=True, facecolor='g')
# plt.scatter(a, b, s=0.1)
# plt.plot([-1,-1,1,1,-1],[-1,1,1,-1,-1])

plt.xlabel('error of fraction[0]')
plt.ylabel('error of fraction[n/2]')
plt.title('Correlation between different error-coefficients')
plt.show()

