
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.stats import norm

plt.rcParams.update({"xtick.labelsize": 18, "ytick.labelsize": 18})

bins = 30
sigma1 = 3.53
sigma2 = 6.18 

data = np.loadtxt("../results.txt", skiprows=201)

number_of_points = [np.arange(0,10001)]
# print(number_of_points)

chi_values = data[:,0]
h_values = data[:,1]
OmegaM = data[:,2]
OmegaK = data[:,3] 
OmegaLambda = 1 - OmegaM - OmegaK
print("len chi values is", len(chi_values))


minimum_chi_value_index = np.argmin(chi_values) #this returns the INDEX 
minimum_chi_value = chi_values[minimum_chi_value_index] #this returns the value of the indext from previous line 
print("The minumum chi squared value is", minimum_chi_value)

print(type(OmegaM))
print("Len of OmegaM is",len(OmegaM), "Len of OmegaK is", len(OmegaK), "Len of OmegaLambda is", len(OmegaLambda))

accepted_OmegaM_1_sigma = OmegaM[chi_values < minimum_chi_value + sigma1]
accepted_OmegaK_1_sigma = OmegaK[chi_values < minimum_chi_value + sigma1]
accepted_OmegaLambda_1_sigma = OmegaLambda[chi_values < minimum_chi_value + sigma1]
accepted_h_values_1_sigma = h_values[chi_values < minimum_chi_value + sigma1] 

accepted_OmegaM_2_sigma = OmegaM[chi_values < minimum_chi_value + sigma2]
accepted_OmegaK_2_sigma = OmegaK[chi_values < minimum_chi_value + sigma2]
accepted_OmegaLambda_2_sigma = OmegaLambda[chi_values < minimum_chi_value + sigma2]
accepted_h_values_2_sigma = h_values[chi_values < minimum_chi_value + sigma2] 


print(len(accepted_OmegaK_1_sigma))
print(len(accepted_OmegaM_1_sigma))
print(len(accepted_OmegaLambda_1_sigma))

std_OmegaM_1_sigma = np.std(accepted_OmegaM_1_sigma)
print("The standard deviation for OmegaM",std_OmegaM_1_sigma)
std_OmegaK_1_sigma = np.std(accepted_OmegaK_1_sigma)
print("The standard deviation for OmegaK", std_OmegaK_1_sigma)
std_OmegaLambda_1_sigma = np.std(accepted_OmegaLambda_1_sigma)
print("The standard deviation for OmegaLambda", std_OmegaLambda_1_sigma)
std_h_value_1_sigma = np.std(accepted_h_values_1_sigma)
print("The standard deviation for h", std_h_value_1_sigma)

flat_universe_omega_m = np.linspace(0,1,100)
flat_universe_omega_lambda = 1-flat_universe_omega_m

# plt.figure(figsize = (8,6))
# plt.scatter(accepted_OmegaM_1_sigma, accepted_OmegaLambda_1_sigma, color = "#d36a7e", label = r"$1\sigma$", zorder = 2)
# plt.scatter(accepted_OmegaM_2_sigma, accepted_OmegaLambda_2_sigma, color = "#b13745", label = r"2$\sigma$", zorder = 1)
# plt.plot(flat_universe_omega_m, flat_universe_omega_lambda, color = "#808080", linestyle = "dashed", label = "Flat universe", zorder = 3)
# plt.ylim(0,1.5)
# plt.xlim(0,1)
# plt.xlabel("$\Omega_M$", fontsize = 20)
# plt.ylabel("$\Omega_\Lambda$", fontsize = 20)
# plt.legend(fontsize = 16)
# # plt.savefig("plots/scatter_supernova.pdf")
# plt.show()

print("The minimum value for the chi2 is", minimum_chi_value, "found in the index", minimum_chi_value_index)
print("The best fit for h is", data[minimum_chi_value_index, 1])
print("The best fit for OmegaM is", data[minimum_chi_value_index, 2])
print("The best fit for OmegaK is", data[minimum_chi_value_index, 3]) 
print("the best fit for OmegaLambda is", OmegaLambda[minimum_chi_value_index])



# #-----------------------------------------------------------------------------------------------------------
# #For plotting OmegaM histrogram

# mu_OmegaM = np.mean(accepted_OmegaM_1_sigma)
# sigma_OmegaM = std_OmegaM_1_sigma
# C_OmegaM = 1/(np.sqrt(2*np.pi)*sigma_OmegaM)

# gaussian_x_values = np.linspace(0, 2*mu_OmegaM, 1000)
# gaussian = C_OmegaM * np.exp(-(1/2)*(pow(gaussian_x_values-mu_OmegaM,2)/pow(sigma_OmegaM,2)))

# plt.figure(figsize = (10,7))
# plt.hist(accepted_OmegaM_1_sigma, bins, density = True, color = "#ffa6c9")
# plt.plot(gaussian_x_values, gaussian, color = "#00035b", label = "Gaussian fit")
# plt.xlabel("$\Omega_M$", fontsize = 18)
# plt.axvline(x = data[minimum_chi_value_index, 2] , color = "#ff0000", linestyle = "dashed", label = "Best fit parameter")
# plt.axvline(x = 0.315 , color = "#30d5c8", linestyle = "dashed", label = "Planck best-fit parameter")
# plt.axvline(x = mu_OmegaM, color = "#008000", linestyle = "dashed", label = "Mean")
# plt.legend(fontsize = 12)
# # plt.savefig("plots/omegaM_histogram.pdf")
# plt.show()
# # #-----------------------------------------------------------------------------------------------------------

# #-----------------------------------------------------------------------------------------------------------
# #For plotting OmegaK histrogram

# mu_OmegaK = np.mean(accepted_OmegaK_1_sigma)
# print("mu for OmegaK is", mu_OmegaK)
# sigma_OmegaK = std_OmegaK_1_sigma 
# C_OmegaK = (1/(np.sqrt(2*np.pi)*sigma_OmegaK))

# gaussian_x_values = np.linspace(mu_OmegaK  - (4*sigma_OmegaK), mu_OmegaK  + (4*sigma_OmegaK), 1000)
# gaussian = C_OmegaK * np.exp(-(1/2)*(pow(gaussian_x_values-mu_OmegaK ,2)/pow(sigma_OmegaK,2)))

# plt.figure(figsize = (10,7))
# plt.hist(accepted_OmegaK_1_sigma, bins, density = True, color = "#ffa6c9")
# plt.plot(gaussian_x_values, gaussian, color = "#00035b", label = "Gaussian fit")
# plt.xlabel("$\Omega_\kappa$", fontsize = 18)
# plt.axvline(x = data[minimum_chi_value_index, 3] , color = "#ff0000", linestyle = "dashed", label = "Best fit parameter")
# plt.axvline(x = 0.001 , color = "#30d5c8", linestyle = "dashed", label = "Planck best-fit parameter")
# plt.axvline(x = mu_OmegaK , color = "#008000", linestyle = "dashed", label = "Mean")
# plt.legend(fontsize = 12)
# plt.savefig("plots/omegaK_histogram.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------

# #-----------------------------------------------------------------------------------------------------------
# # #For plotting OmegaLambda histrogram

# mu_OmegaLambda = np.mean(accepted_OmegaLambda_1_sigma)
# sigma_OmegaLambda = std_OmegaLambda_1_sigma
# C_OmegaLambda = 1/(np.sqrt(2*np.pi)*sigma_OmegaLambda)

# gaussian_x_values = np.linspace(mu_OmegaLambda  - (4*sigma_OmegaLambda), mu_OmegaLambda  + (4*sigma_OmegaLambda), 1000)
# gaussian = C_OmegaLambda * np.exp(-(1/2)*(pow(gaussian_x_values-mu_OmegaLambda ,2)/pow(sigma_OmegaLambda,2)))

# plt.figure(figsize = (10,7))
# histo = plt.hist(accepted_OmegaLambda_1_sigma, bins, density = True, color = "#ffa6c9")
# plt.plot(gaussian_x_values, gaussian, color = "#00035b", label = "Gaussian fit")
# plt.xlabel("$\Omega_\Lambda$", fontsize = 18)
# plt.axvline(x = OmegaLambda[minimum_chi_value_index] , color = "#ff0000", linestyle = "dashed", label = "Best fit parameter")
# plt.axvline(x = 0.684 , color = "#30d5c8", linestyle = "dashed", label = "Planck best-fit parameter")
# plt.axvline(x = mu_OmegaLambda , color = "#008000", linestyle = "dashed", label = "Mean")
# plt.legend(fontsize = 12)
# # plt.savefig("plots/omegaLambda_histogram.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------

# #-----------------------------------------------------------------------------------------------------------
# #For plotting H0 histrogram

# H0_values = h_values * 100
# accepted_H0_values_1_sigma = accepted_h_values_1_sigma * 100
# std_H0_value_1_sigma = np.std(accepted_H0_values_1_sigma)

# mu_H0 = np.mean(accepted_H0_values_1_sigma)
# sigma_H0 = std_H0_value_1_sigma
# C_hubble0 = 1/(np.sqrt(2*np.pi)*sigma_H0)

# gaussian_x_values = np.linspace(mu_H0 - (4*sigma_H0), mu_H0 + (4*sigma_H0), 1000)
# gaussian = C_hubble0 * np.exp(-(1/2)*(pow(gaussian_x_values-mu_H0,2)/pow(sigma_H0,2)))

# plt.figure(figsize = (10,7))
# histo = plt.hist(accepted_H0_values_1_sigma, bins, density = True, color = "#ffa6c9")
# plt.plot(gaussian_x_values, gaussian, color = "#00035b", label = "Gaussian fit")
# plt.xlabel("$H_{0}$", fontsize = 18)
# plt.axvline(x = H0_values[minimum_chi_value_index] , color = "#ff0000", linestyle = "dashed", label = "Best fit parameter")
# plt.axvline(x = 67.4, color = "#30d5c8", linestyle = "dashed", label = "Planck best-fit parameter")
# plt.axvline(x = mu_H0 , color = "#008000", linestyle = "dashed", label = "Mean")
# plt.legend(fontsize = 12)
# # plt.savefig("plots/H0_histogram.pdf")
# plt.show()
# #-----------------------------------------------------------------------------------------------------------