

import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("../recombination_just_saha.txt")


x = data1[:,0]
Xe = data1[:,1]


plt.plot(x, Xe)
plt.show()