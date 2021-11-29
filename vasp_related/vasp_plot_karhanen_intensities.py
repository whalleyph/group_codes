#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

i, freq, intensity = np.loadtxt('exact.reform.res.txt', unpack=True)
plt.bar(freq, intensity)
plt.show()
