import math
import numpy as np
import matplotlib.pyplot as plt

data = [] # [ (value, frequency), ... ]

a,b = list(zip(*data))
a = list(a)
b = list(b)

sum_freq = sum(b)
max_freq = max(b)

# normal_samps = 0 + 55 * np.random.randn(1000 * 1000)

# scale = 0.95
scale = 0.008

def normal_prob(x):
    sigma = 58
    return scale * math.exp(-x * x / (2 * sigma * sigma))

# x = []
# for (v,f) in data:
#     for i in range(f):
#         x.append(v)

# the histogram of the data
plt.plot(a, np.array(b) / max_freq)
# plt.hist(x, 50, density=True, facecolor='g')
# plt.plot(a, list(map(normal_prob, a)))

plt.xlabel('s1[i]')
plt.ylabel('Frequency')
plt.title('Histogram of signature values')
plt.show()

