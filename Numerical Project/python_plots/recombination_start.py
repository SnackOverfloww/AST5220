
import numpy as np
import matplotlib.pyplot as plt

#Load data for recombination history
data = np.loadtxt("../recombination.txt")
# data2 = np.loadtxt("../g_tilde_values.txt")

x = data[:,0]
Xe = data[:,1] 
ne = data[:,2]
tau = data[:,3]
tau_deriv = data[:,4]
tau_double_deriv = data[:,5]
g_tilde = data[:,6]


# # Plot Xe of x
# plt.plot(x,Xe)
# plt.xlabel("x")
# plt.ylabel("$X_e$")
# plt.yscale("log")
# plt.show()

# #Plot tau of x and its derivatives
# plt.plot(x,tau, label = "$\\tau$")
# plt.plot(x, -tau_deriv, label = "$-\\tau'(x)$")
# plt.plot(x, tau_double_deriv, label = "$\\tau''(x)$")
# plt.xlabel("x")
# plt.ylabel("$\\tau$")
# plt.xlim(-12, 0)
# plt.ylim(10**-8, 10**8)
# plt.yscale("log")
# plt.legend(fontsize = 16)
# plt.show()

plt.plot(x,g_tilde)
plt.xlabel("x")
plt.ylabel("g")
# plt.yscale("log")
plt.show()
