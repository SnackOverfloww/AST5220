
import numpy as np
import matplotlib.pyplot as plt

#Load data for recombination history
data_no_reionisation = np.loadtxt("../recombination.txt")
data_with_reionisation = np.loadtxt("../recombination_with_reionisation.txt")
# data2 = np.loadtxt("../g_tilde_values.txt")

#Data without reionization
x_no_reionisation = data_no_reionisation[:,0]
Xe_no_reionisation = data_no_reionisation[:,1] 
ne_no_reionisation = data_no_reionisation[:,2]
tau_no_reionisation = data_no_reionisation[:,3]
tau_deriv_no_reionisation = data_no_reionisation[:,4]
tau_double_deriv_no_reionisation = data_no_reionisation[:,5]
g_tilde_no_reionisation = data_no_reionisation[:,6]
g_tilde_deriv_no_reionisation = data_no_reionisation[:,7]
g_tilde_double_deriv_no_reionisation = data_no_reionisation[:,8]

x_with_reionisation = data_with_reionisation[:,0]
Xe_with_reionisation = data_with_reionisation[:,1] 
ne_with_reionisation = data_with_reionisation[:,2]
tau_with_reionisation = data_with_reionisation[:,3]
tau_deriv_with_reionisation = data_with_reionisation[:,4]
tau_double_deriv_with_reionisation = data_with_reionisation[:,5]
g_tilde_with_reionisation = data_with_reionisation[:,6]
g_tilde_deriv_with_reionisation = data_with_reionisation[:,7]
g_tilde_double_deriv_with_reionisation = data_with_reionisation[:,8]

# # Plot Xe of x
# plt.plot(x_with_reionisation,Xe_with_reionisation, linestyle=(0, (2, 2)), linewidth=2, label = "$X_e(x)$ with reionisation")
# plt.plot(x_no_reionisation,Xe_no_reionisation, label = "$X_e(x)$")
# plt.xlabel("x")
# plt.ylabel("$X_e$")
# plt.yscale("log")
# plt.legend()
# plt.show()

# #Plot tau of x and its derivatives (WITH REIONISATION)
# plt.plot(x_with_reionisation,tau_with_reionisation, label = "$\\tau(x)$ with reionisation")
# plt.plot(x_with_reionisation, -tau_deriv_with_reionisation, label = "$-\\tau\\,'(x)$ with reionisation")
# plt.plot(x_with_reionisation, tau_double_deriv_with_reionisation, label = "$\\tau\\,''(x)$ with reionisation")
# plt.xlabel("x")
# # plt.ylabel("$\\tau$")
# plt.xlim(-12, 0)
# plt.ylim(10**-8, 10**8)
# plt.yscale("log")
# plt.legend(fontsize = 16)
# plt.show()

# #Plot tau of x and its derivatives (NO REIONISATION)
# plt.plot(x_no_reionisation,tau_no_reionisation, label = "$\\tau$")
# plt.plot(x_no_reionisation, -tau_deriv_no_reionisation, label = "$-\\tau\\,'(x)$")
# plt.plot(x_no_reionisation, tau_double_deriv_no_reionisation, label = "$\\tau\\,''(x)$")
# plt.xlabel("x")
# # plt.ylabel("$\\tau$")
# plt.xlim(-12, 0)
# plt.ylim(10**-8, 10**8)
# plt.yscale("log")
# plt.legend(fontsize = 16)
# plt.show()

# #Plot g_tilde of x
# plt.plot(x_with_reionisation,g_tilde_with_reionisation, linestyle=(0, (2, 2)), linewidth=2, label = "$\\tilde{g}(x)$ with reionisation")
# plt.plot(x_no_reionisation,g_tilde_no_reionisation, label = "$\\tilde{g}(x)$")
# plt.xlabel("x")
# plt.ylabel("$\\tilde{g}$")
# plt.legend()
# plt.show()

# #Plot g_tilde of x first derivative
# plt.plot(x_with_reionisation,g_tilde_deriv_with_reionisation, linestyle=(0, (2, 2)), linewidth=2, label = "$\\tilde{g}\\,'(x)$ with reionisation")
# plt.plot(x_no_reionisation,g_tilde_deriv_no_reionisation, label = "$\\tilde{g}\\,'(x)$")
# plt.xlabel("x")
# plt.ylabel("$\\tilde{g}\\,'(x)$")
# plt.legend()
# plt.show()

#Plot g_tilde of x second derivative
plt.plot(x_with_reionisation,g_tilde_double_deriv_with_reionisation, linestyle=(0, (2, 2)), linewidth=2, label = "$\\tilde{g}\\,''(x)$ with reionisation")
plt.plot(x_no_reionisation,g_tilde_double_deriv_no_reionisation, label = "$\\tilde{g}\\,''(x)$" )
plt.xlabel("x")
plt.ylabel("$\\tilde{g}\\,''(x)$")
plt.legend()
plt.show()

