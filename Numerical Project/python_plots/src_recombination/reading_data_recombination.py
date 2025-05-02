
import numpy as np

#Load data for recombination history
data_no_reionisation = np.genfromtxt("recombination_no_reionisation.txt")
data_with_reionisation = np.genfromtxt("recombination_with_reionisation.txt")
data_just_saha = np.genfromtxt("recombination_just_saha.txt")
print(data_just_saha)

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
sound_horizon = data_no_reionisation[:,9]

#Data with reionisation
x_with_reionisation = data_with_reionisation[:,0]
Xe_with_reionisation = data_with_reionisation[:,1] 
ne_with_reionisation = data_with_reionisation[:,2]
tau_with_reionisation = data_with_reionisation[:,3]
tau_deriv_with_reionisation = data_with_reionisation[:,4]
tau_double_deriv_with_reionisation = data_with_reionisation[:,5]
g_tilde_with_reionisation = data_with_reionisation[:,6]
g_tilde_deriv_with_reionisation = data_with_reionisation[:,7]
g_tilde_double_deriv_with_reionisation = data_with_reionisation[:,8]

x_just_saha = data_just_saha[:,0]
Xe_just_saha = data_just_saha[:,1]
ne_just_saha = data_just_saha[:,2]
tau_just_saha = data_just_saha[:,3]
tau_deriv_just_saha = data_just_saha[:,4]
tau_double_deriv_just_saha = data_just_saha[:,5]
g_tilde_just_saha = data_just_saha[:,6]
g_tilde_deriv_just_saha = data_just_saha[:,7]
g_tilde_double_deriv_just_saha = data_just_saha[:,8]