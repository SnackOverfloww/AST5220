
import numpy as np
import matplotlib.pyplot as plt

exec(open("python_plots/src_recombination/reading_data_recombination.py").read())

const_m_to_Mpc = 1 / (3.08567758 * pow(10, 16) * pow(10, 6)) 

#==============================================================
#Function to calculate x-value from redshift
def x_from_z(z):
    a = 1/(z + 1)
    x = np.log(a)
    return x
#==============================================================

#==============================================================
#Function to calculate x-value from redshift
def z_from_x(x):
    a = np.exp(x)
    z = (1/a) - 1
    return z
#==============================================================

# #==============================================================
# #Function to normalise data
def normalising(data_set):
    sum_of_data = np.sum(data_set)
    normalised_data = data_set / sum_of_data
    return normalised_data
# #===============================================================


#======================================================================================================================================
# There are two ways to calculate where the decoupling happens via the last scattering: either using where tau = 1, or using where
# the visibility function peaks: 

# Using the visibility function peak to find redshift for last scattering:
numbers = np.arange(0,len(x_no_reionisation))

print("Using the peak of the visibility function:")
# print("The peak value of the visibility function is", np.max(g_tilde_no_reionisation))
x_index_for_peak_visibility_function = np.where(g_tilde_no_reionisation == np.max(g_tilde_no_reionisation))
print("x index is ", x_index_for_peak_visibility_function)
print("The index for the x_value where the peak of the vicibility function is", x_index_for_peak_visibility_function[0][0])
x_value_no_reionisation_visibility_peak = x_no_reionisation[x_index_for_peak_visibility_function][0] 
print("The photon decopling happens at an x-value", x_value_no_reionisation_visibility_peak )
z = (1/np.exp(x_no_reionisation[x_index_for_peak_visibility_function])) - 1 
print("The photon decopling happens at a redshift", str.format('{0:.0f}', z[0]), "\n")

#===========================================================================================================================
# Using where tau = 1 to find redshift for last scattering: 
i = 0
while tau_no_reionisation[i] > 1:
    i+=1
bigger_than_1 = tau_no_reionisation[i]
    
j = len(tau_no_reionisation)
while tau_no_reionisation[i] < 1:
    i-=1
smaller_than_1 = tau_no_reionisation[i]

index_smaller = (np.where(tau_no_reionisation==smaller_than_1)[0][0]) 
index_bigger = (np.where(tau_no_reionisation==bigger_than_1)[0][0]) 
x_value_tau_1 = (x_no_reionisation[index_smaller] + x_no_reionisation[index_bigger])/2
x_value_uncertainty_tau_1 = np.abs(x_no_reionisation[index_smaller] - x_no_reionisation[index_bigger])/2 

print("Using the point where tau = 1: ")
print("The value for x where tau is 1 is ", 
    str.format('{0:.4f}', x_value_tau_1), "+-", 
    str.format('{0:.4f}', x_value_uncertainty_tau_1))
print("This corresponds to a z-value ", str.format('{0:.0f}', z_from_x(x_value_tau_1)), "+-", 
    str.format('{0:.0f}', z_from_x(x_value_tau_1)*x_value_uncertainty_tau_1), "\n")
#========================================================================================================================

#=======================================================================================================================================
# To find where recombination happens, we use the point where X_e = 0.1:
i = 0
while Xe_no_reionisation[i] > 0.1:
    i+=1
bigger_than_1 = Xe_no_reionisation[i]
    
j = len(Xe_no_reionisation)
while Xe_no_reionisation[i] < 0.1:
    i-=1
smaller_than_1 = Xe_no_reionisation[i]


index_smaller = (np.where(Xe_no_reionisation==smaller_than_1)[0][0])

index_bigger = (np.where(Xe_no_reionisation==bigger_than_1)[0][0]) 
x_value_no_reionisation_recombination = (x_no_reionisation[index_smaller] + x_no_reionisation[index_bigger])/2
x_value_uncertainty_no_reionisation_recombination = np.abs(x_no_reionisation[index_smaller] - x_no_reionisation[index_bigger])/2 


print("Using the point where Xe = 0.1: ")
print("The value for x where Xe is 0.1 is ", 
    str.format('{0:.4f}', x_value_no_reionisation_recombination), "+-", 
    str.format('{0:.4f}', x_value_uncertainty_no_reionisation_recombination))
print("This corresponds to a z-value ", str.format('{0:.0f}', z_from_x(x_value_no_reionisation_recombination)), "+-", 
    str.format('{0:.0f}', z_from_x(x_value_no_reionisation_recombination)*x_value_uncertainty_no_reionisation_recombination), "\n")

# #======================================================================================================================================
# print("The freeze-out abundance of free elctrons today is", np.exp(ne_no_reionisation[-1]))
# print(np.min(Xe_no_reionisation))

# #=======================================================================
sound_horizon_in_m = sound_horizon[x_index_for_peak_visibility_function[0][0]] 
sound_horizon_in_Mpc = sound_horizon_in_m * const_m_to_Mpc 
print("The sound horizon at recombination is calculated to be", str.format('{0:.2f}'
    , sound_horizon_in_Mpc), "Mpc. This is when we only use the numerical results (the 5000 points taken from spline), not the spline")
# #========================================================================


#===========================================================================================
print("The freeze-out abundance of free electrons today is", Xe_no_reionisation[-1])
