import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u

import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u

plt.rcParams.update({"xtick.labelsize": 13, "ytick.labelsize": 13})

data = np.loadtxt("../cosmology_density_parameters.txt")
data2 = np.loadtxt("../cosmology.txt")

x_values = data[:,0]
matter = data[:,1] + data[:,2]
relativistic = data[:,4] + data[:,5]
dark_energy = data[:,3]

x = data2[:,0]
eta_of_x = data2[:,1]
t_of_x  = data2[:,11]

t_of_x_scaled = t_of_x / (60*60*24*365)
eta_of_x_scaled = (eta_of_x * u.m / const.c) / (60*60*24*365)

def from_x_to_a(x_index):
    a = np.exp(x_index)
    return a

def from_x_to_z(x_index):
    z = (1/from_x_to_a(x_index)) - 1
    return z

numbers = np.arange(1,1000)


#Density parameters
omega_R = (3.81093*pow(10,-5)) + (5.50896*pow(10,-5)) 
omega_M = (0.05 + 0.267)
omega_Lambda = 0.682907

#Constants
H0 = 2.17132 *pow(10,-18)
t0 = 13.78475393*(10**9)*60*60*24*365
s_to_year = 60*60*24*365
steps = 1000
x_start = np.log(10**-8)


#Function to calculate matter-radiation equality time
def matter_radiation_x_value(omegaR, omegaM):
    a = omegaR/omegaM
    x = np.log(a) 
    return x

#Function to calculate matter-dark energy equality time
def matter_dark_energy_x_value(omegaM, omegaLambda):
    a = pow((omegaM/omegaLambda), (1/3))
    x = np.log(a) 
    return x

#Function to calculate beginning of acceleration time
def beginning_of_acceleration_x_value(omegaM, omegaLambda):
    a = pow((omegaM/(2*omegaLambda)), (1/3))
    x = np.log(a) 
    return x


x_scale_rad = np.linspace(np.log(10**-8), matter_radiation_x_value(omega_R, omega_M), steps)
x_scale_mat = np.linspace(matter_radiation_x_value(omega_R, omega_M), matter_dark_energy_x_value(omega_M, omega_Lambda), steps)
x_scale_lam = np.linspace(matter_dark_energy_x_value(omega_M, omega_Lambda), 5, steps)


#Function for expression t for radiation dominated era
def radiation_dominated_time(H0, omega_R, x):
    t_R = (((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*x))) - ((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*x_start)))) / s_to_year # CORRECT!!!
    return t_R

def matter_dominated_time(H0, omegaM, omegaR, x):
    t_M = radiation_dominated_time(H0, omega_R, matter_radiation_x_value(omegaR, omegaM)) + (((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*x))) 
    - ((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*matter_radiation_x_value(omega_M, omega_R))))) / s_to_year # CORRECT!!! 
    return t_M

def dark_energy_dominated_time(H0, omegaM, omegaR, omegaLambda, x):
    t_Lambda = matter_dominiated_time(H0, omega_M, omega_R, matter_dark_energy_x_value(omegaM, omegaLambda)) 
    + ((1/(H0*np.sqrt(omega_Lambda)) * x) - (matter_dark_energy_x_value(omega_M, omega_Lambda)/(H0*np.sqrt(omega_Lambda))))  / s_to_year #CORRECT FOR NOW!!!
    return t_Lambda


print(radiation_dominated_time(H0, omega_R, matter_radiation_x_value(omega_R, omega_M)))

# eta_R= (((1 / (H0*np.sqrt(omega_R))) * np.exp(x_scale_rad)) - ((1 / (H0*np.sqrt(omega_R))) * np.exp(np.log(10**-8)))) / s_to_year #CORRECT!!!


# eta_M = (51064) + (((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*x_scale_mat)) - ((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*((-8.131921133813))))) / s_to_year #CORRECT I THINK!


# eta_Lambda = (10.38 * 10**9) + (((-1/((H0)*np.sqrt(omega_Lambda)))*np.exp(-1*x_scale_lam)) + ((1/(H0*np.sqrt(omega_Lambda)))*np.exp(0.2558189708133))) / s_to_year #Not correct
# eta_Lambda = ((-1/(H0*np.sqrt(omega_Lambda)))*(np.exp(-1 * x_scale_lam) - np.exp(0.255))) / (60*60*24*365) #alternative but same epression

#Plotting the plot for eta(x), t(x), equality times, and beginning of acceleration
plt.figure(figsize = (16,8))
plt.plot(x, t_of_x_scaled, label = "Time t(x)")
plt.plot(x, (eta_of_x_scaled), label = "Conformal time $\eta(x)/c$")
plt.plot(x_scale_rad, radiation_dominated_time(H0, omega_R, x_scale_rad) , label = "Radiation dominated era", linestyle = "dashed")
plt.plot(x_scale_mat, matter_dominated_time(H0, omega_M, omega_R, x_scale_mat), label = "Matter dominated era", linestyle = "dashed")
# plt.plot(x_scale_lam, t_Lambda, label = "Lambda", linestyle = "dashed")
# plt.plot(x_scale_rad, eta_R, linestyle = "dashed")
# plt.plot(x_scale_mat, eta_M, label = "eta M", linestyle = "dashed")
# plt.plot(x_scale_lam, eta_Lambda, label = "eta lambda", linestyle = "dashed")
plt.xlabel("x", fontsize = 15)
plt.xlim(-14,5)
plt.ylim(10**-0.5, 10**13)
plt.yscale("log")
plt.legend(fontsize = "12")
# plt.savefig("plots/equality_times.pdf")
plt.show()
    
    

