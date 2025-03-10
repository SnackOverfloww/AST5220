
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../recombination.txt")

x = data[:,0]
Xe = data[:,1] 
ne = data[:,2]
tau = data[:,3]


# Plot Xe of x
plt.plot(x,tau)
plt.xlabel("x")
plt.ylabel("$X_e$")
plt.yscale("log")
plt.show()