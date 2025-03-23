
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u

plt.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize": 15})

data = np.loadtxt("../cosmology.txt")

# xmin = -20
# xmax = 5
# npts = 1000

x = data[:,0]
eta_of_x = data[:,1]
Hp_of_x = data[:,2]
dHpdx_of_x = data[:,3]
ddHpddx_of_x = data[:,4]
OmegaB = data[:,5]
OmegaCDM = data[:,6]
OmegaLambda = data[:,7]
OmegaR = data[:,8]
OmegaNu = data[:,9]
OmegaK = data[:,10]
t_of_x = data[:,11]
eta_of_x_derivative = data[:,12]


total_density_parameters = OmegaB + OmegaCDM + OmegaLambda + OmegaR + OmegaNu + OmegaK
average_density_parameters = np.mean(total_density_parameters)

#print(total_density_parameters)
print(average_density_parameters)



# #-----------------------------------------------------------------------------------------------------------------------
# # for plotting Hp_of_x
# Hp_of_x_adjusted = (Hp_of_x) *(3.086*(10**19)/100)
# print(Hp_of_x_adjusted)

# plt.figure(figsize = (10,7))
# plt.plot(x, Hp_of_x_adjusted, color = "#fe89c8")
# plt.xlabel("x", fontsize = 22)
# plt.ylabel(r"${\mathcal{H}(x)}$", fontsize = 22)
# plt.xlim(-12,0)
# plt.ylim(10**-1, 10**3)
# plt.yscale("log")
# # plt.savefig("plots/Hp_of_x_plot.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------------------

# #-----------------------------------------------------------------------------------------------------------------------
# # # for plotting eta_of_x
# eta_of_x_adjusted = eta_of_x / (3.08*(10**22)) 
# plt.plot(x, eta_of_x_adjusted, color = "#fe89c8")
# plt.xlabel("x", fontsize = 13)
# plt.ylabel("$\eta(x)$", fontsize = 14)
# plt.yscale("log")
# plt.xlim(-12,0)
# plt.ylim(10**0, 3*(10**4))
# # plt.savefig("plots/eta_plot.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------------------


# #-----------------------------------------------------------------------------------------------------------------------
# # # for plotting (eta_of_x * Hp_of_x) / c
# c = const.c
# x_array_for_1 = np.linspace(-20,-10,100) 
# y_array_1 = [1]*len(x_array_for_1)

# plt.figure(figsize = (8,6))
# plt.plot(x, ((eta_of_x * Hp_of_x)/c), color = "#fe89c8", label = "Numerically calculated")
# plt.plot(x_array_for_1, y_array_1, color = "#ff0000", linestyle = "dashed", label = "Radiation dominated")
# plt.xlabel("x", fontsize = 16)
# plt.ylabel(r"$\frac{\eta(x) \mathcal{H}(x)}{c}$", fontsize = 20)
# plt.xlim(-13,0)
# plt.ylim(0.75, 3)
# plt.legend(fontsize = 16)
# plt.savefig("plots/eta_times_Hp_divided_by_c.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------------------


# #-----------------------------------------------------------------------------------------------------------------------
# # #for plotting dHpdx_of_x / Hp_of_x

# x_array_for_minus_1 = np.linspace(-20,-7.5,100) 
# y_array_minus_1 = [-1]*len(x_array_for_minus_1)

# x_array_for_minus_0_5 = np.linspace(-7.5,-0,100) 
# y_array_minus_0_5 = [-0.5] * len(x_array_for_minus_0_5)

# x_array_for_1 = np.linspace(0,5,100) 
# y_array_1 = [1] * len(x_array_for_1)

# plt.figure(figsize=(10, 7))
# plt.vlines(x = -7.5, ymin = -1, ymax = -0.5, color = "#808080", linestyle = "dashed")
# plt.vlines(x = 0, ymin = -0.5, ymax = 1, color = "#808080", linestyle = "dashed")

# plt.plot(x, (dHpdx_of_x / Hp_of_x), color = "#fe89c8", label = "Calculated")
# plt.plot(x_array_for_minus_1, y_array_minus_1, color = "#D55E00", linestyle = "dashed", label = "Radiation dominated")
# plt.plot(x_array_for_minus_0_5, y_array_minus_0_5, color = "#009E73", linestyle = "dashed", label = "Matter dominated")
# plt.plot(x_array_for_1, y_array_1, color = "#882255", linestyle = "dashed", label = "Dark Energy dominated")
# plt.xlabel("x", fontsize = 14)
# plt.ylabel("$\\frac{d\mathcal{H}(x)}{dx}\\frac{1}{\mathcal{H}(x)}$", fontsize = 20)
# plt.legend(fontsize = 18)
# plt.savefig("plots/dHpdx_of_x_divided_by_Hp.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------------------


# #-----------------------------------------------------------------------------------------------------------------------
# # for plotting ddHpddx_of_x / Hp_of_x

# x_array_for_1_4 = np.linspace(-8,-1,100) 
# y_array_1_4 = [1/4]*len(x_array_for_1_4)

# x_array_for_1_first = np.linspace(-20,-8,100)
# y_array_1_first = [1]*len(x_array_for_1_first)

# x_array_for_1_second = np.linspace(-1,5,100)
# y_array_1_second = [1]*len(x_array_for_1_second)

# plt.figure(figsize=(12, 8))
# plt.vlines(x = -8, ymin = 1/4, ymax = 1, color = "#808080", linestyle = "dashed")
# plt.vlines(x = -1, ymin = 1/4, ymax = 1, color = "#808080", linestyle = "dashed")

# plt.plot(x, (ddHpddx_of_x / Hp_of_x), color = "#fe89c8", label = "Calculated") 
# plt.plot(x_array_for_1_first, y_array_1_first, color = "#ff0000", linestyle = "dashed", label = "Radiation dominated")
# plt.plot(x_array_for_1_4, y_array_1_4, color = "#009E73", linestyle = "dashed", label = "Matter dominated")
# plt.plot(x_array_for_1_second, y_array_1_second, color = "#882255", linestyle = "dashed", label = "Dark Energy dominated")
# plt.xlabel("x", fontsize = 20)
# plt.ylabel("$\\frac{d^2 \mathcal{H}(x)}{dx^2}$ $\\frac{1}{\mathcal{H}(x)}$", fontsize = 20)
# plt.legend(fontsize = 16)
# plt.savefig("plots/double_derivative_test.pdf")
# plt.show()

# #-----------------------------------------------------------------------------------------------------------------------

# #-----------------------------------------------------------------------------------------------------------------------
# #for plotting total composition
# plt.figure(figsize=(9,7))
# plt.plot(x, OmegaR+OmegaNu, color = "#b5338a", label = "$\Omega_{Relativistic}$ = $\Omega_{\gamma}$ + $\Omega_{\\nu}$")
# plt.plot(x, OmegaB+OmegaCDM, color = "#dab1da", label = "$\Omega_{Matter}$ = $\Omega_{b}$ + $\Omega_{CDM}$")
# plt.plot(x, OmegaLambda, color = "#4b006e", label = "$\Omega_{\Lambda}$")
# plt.vlines(x = -8.13814, ymin = 0, ymax = 1, color = "k", linestyle = "dashed", label = "Matter - radiation \n equality")
# plt.vlines(x = -0.255255, ymin = 0, ymax = 1, color = "#949494", linestyle = "dashed", label = "Matter - dark energy \n equality")
# plt.xlabel("$x$", fontsize = 15)
# plt.ylabel("Fraction of composition", fontsize = 15)
# plt.ylim(0)
# plt.xlim(-20,5)
# plt.legend(fontsize = 13)
# # plt.savefig("plots/density_domination.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------------------
# # #NOT USED IN REPORT: 
# # For plotting t_of_x and eta_of_x

# t_of_x_scaled = t_of_x / (60*60*365)
# eta_of_x_scaled = (eta_of_x / const.c) / (60*60*365)

# plt.plot(x, t_of_x_scaled)
# plt.plot(x, (eta_of_x_scaled))
# plt.vlines(x = 0, ymin = 0, ymax = 10**13, color = "#808080", linestyle = "dashed")
# plt.xlim(-14,5)
# plt.ylim(10**-0.5, 10**13)
# plt.yscale("log")
# plt.show()
# # -----------------------------------------------------------------------------------------------------------------------
