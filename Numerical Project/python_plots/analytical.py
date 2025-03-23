import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u
# import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u

# plt.rcParams.update({"xtick.labelsize": 13, "ytick.labelsize": 13})

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
x_start = np.log(10**-8)

matter_radiation_eq = np.log(omega_R / omega_M)
matter_dark_energy_eq = np.log((omega_M / omega_Lambda)**(1/3))
acceleration = np.log((omega_M / (2*omega_Lambda))**(1/3)) 

print(matter_radiation_eq)
print(matter_dark_energy_eq)
print(acceleration)

x_scale_rad = np.linspace(-20, matter_radiation_eq, 1000)
x_scale_mat = np.linspace(matter_radiation_eq, matter_dark_energy_eq, 1000)
x_scale_lam = np.linspace(matter_dark_energy_eq, 5, 1000)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#calculating t
t_R = (((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*x_scale_rad))) - ((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*np.log(10**-8))))) / s_to_year # CORRECT!!!
matter_radiation_eq_time = (((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*matter_radiation_eq))) - ((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*x_start)))) / s_to_year 
print("The analytical t for matter-radiation equality is", matter_radiation_eq_time)

t_M = matter_radiation_eq_time + (((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*x_scale_mat))) 
    - ((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*matter_radiation_eq)))) / s_to_year # CORRECT!!!
matter_dark_energy_eq_time = matter_radiation_eq_time + (((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*matter_dark_energy_eq))) 
    - ((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*matter_radiation_eq)))) / s_to_year 
print("The analytical t for matter-dark energy equality is", matter_dark_energy_eq_time)

t_Lambda = matter_dark_energy_eq_time + ((1/(H0*np.sqrt(omega_Lambda)) * x_scale_lam) 
        - (matter_dark_energy_eq/(H0*np.sqrt(omega_Lambda))))  / s_to_year #CORRECT FOR NOW!!!
#------------------------------------------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#calculating eta
eta_R= (((1 / (H0*np.sqrt(omega_R))) * np.exp(x_scale_rad)) - ((1 / (H0*np.sqrt(omega_R))) * np.exp(x_start))) / s_to_year 
matter_radiation_eq_eta = (((1 / (H0*np.sqrt(omega_R))) * np.exp(matter_radiation_eq)) - ((1 / (H0*np.sqrt(omega_R))) * np.exp(x_start))) / s_to_year  
print("The analytical eta for matter-radiation equality is", matter_radiation_eq_eta)

eta_M = matter_radiation_eq_eta + (((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*x_scale_mat)) 
    - ((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*((matter_radiation_eq))))) / s_to_year 
matter_dark_energy_eq_eta = matter_radiation_eq_eta + (((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*matter_dark_energy_eq)) 
    - ((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*((matter_radiation_eq))))) / s_to_year #CORRECT I THINK! 
print("The analytical eta for matter-dark energy equality is", matter_dark_energy_eq_eta)

eta_Lambda = matter_dark_energy_eq_eta + (((-1/((H0)*np.sqrt(omega_Lambda)))*np.exp(-1*x_scale_lam)) 
    + ((1/(H0*np.sqrt(omega_Lambda)))*np.exp(-matter_dark_energy_eq))) / s_to_year #Not correct


#------------------------------------------------------------------------------------------------------------------------------------------------------------------



#Plotting the plot for eta(x), t(x), equality times, and beginning of acceleration
plt.figure(figsize = (16,8))
plt.plot(x, t_of_x_scaled, label = "Time t(x)")
plt.plot(x, (eta_of_x_scaled), label = "Conformal time $\eta(x)/c$")
plt.plot(x_scale_rad,t_R, label = "Radiation dominated era $t$", linestyle = "dashed")
plt.plot(x_scale_mat,t_M, label = "Matter dominated era $t$", linestyle = "dashed")
plt.plot(x_scale_lam, t_Lambda, label = "Dark energy dominated era $t$", linestyle = "dashed")
plt.plot(x_scale_rad, eta_R, label = "Radiation dominaed era $\\eta$", linestyle = "dashed")
plt.plot(x_scale_mat, eta_M, label = "Matter dominated era $\\eta$", linestyle = "dashed")
plt.plot(x_scale_lam, eta_Lambda, label = "Dark energy dominated era $\\eta$", linestyle = "dashed")
plt.vlines(x = matter_radiation_eq, ymin = 0, ymax = 10**13, color = "#ff69d4", linestyle = "dotted", label = "Matter - radiation \n equality")
plt.vlines(x = matter_dark_energy_eq, ymin = 0, ymax = 10**13, color = "#50c878", linestyle = "dotted", label = "Matter - radiation \n equality")
plt.vlines(x = acceleration, ymin = 0, ymax = 10**13,color = "#45b6fe", linestyle = "dotted", label = "Acceleration begins")

# plt.vlines(x = 0, ymin = 0, ymax = 10**13, color = "#8a00c4", linestyle = "dotted", label = "Today")
plt.axhline(y = 51064, color = "#949494", linestyle = "dashed")
plt.axhline(y = 10.3782 * 1e9, color = "#949494", linestyle = "dashed")
plt.axhline(y = 7.75249 * 1e9, color = "#949494", linestyle = "dashed")
plt.axhline(y = 3.68439 * 1e8, color = "#949494", linestyle = "dashed")
plt.axhline(y = 42.3666 * 1e9, color = "#949494", linestyle = "dashed", linewidth = 1)
plt.axhline(y = 38.5666 * 1e9, color = "#949494", linestyle = "dashed", linewidth = 1)
plt.yscale("log")
plt.xlabel("x", fontsize = 15)
plt.ylabel("Time in years", fontsize = 15)
plt.xlim(-14,5)


# Set only the specific y-ticks
y_ticks = [51064, 10380000000, 7752490000, 368439000, 42366600000, 38566600000]
y_labels = ["51064 yrs", "10.38 Gyrs", "7.75 Gyrs", "368.4 Myrs", "42.4 Myrs", "38.6 Myrs"] 

x_ticks = [-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
x_labels = [-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5] 


plt.yticks(y_ticks)  # Apply to the axis object
plt.xticks(x_ticks)  # Apply to the axis object
# plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_major_formatter(ScalarFormatter())
plt.gca().yaxis.get_offset_text().set_visible(False)
plt.gca().yaxis.set_major_formatter(plt.FixedFormatter(y_labels))

plt.gca().xaxis.set_major_formatter(ScalarFormatter())
plt.gca().xaxis.get_offset_text().set_visible(False)
plt.gca().xaxis.set_major_formatter(plt.FixedFormatter(x_labels))

plt.xlim(-14, 5)
plt.ylim(10**-0.5, 2 * 1e12)
plt.legend(fontsize = "12")
# plt.savefig("plots/equality_times.pdf")
plt.show()
    

