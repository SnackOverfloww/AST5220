import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u


plt.rcParams.update({"xtick.labelsize": 13, "ytick.labelsize": 13})

data = np.loadtxt("cosmology_density_parameters.txt")
data2 = np.loadtxt("cosmology.txt")

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


#For calculating the matter-radiation equality time
for i in numbers:
    if np.abs(matter[i] - relativistic[i]) < 0.2 and (matter[i] < (0.5 + 0.003) and matter[i] > (0.5 - 0.003) ):
        print("For matter-radiation equality:")
        # print("The difference between matter and relativistic particles is", np.abs(matter[i] - relativistic[i]), "and this difference is found at index", i)
        # print("The x-value for the radiation-matter equality is", x_values[i])
        print("In (normal) time this this equality occurs", t_of_x_scaled[i], "years after the Big Bang")
        # print("This happens at a scale factor", from_x_to_a(x_values[i]), "and redshift", from_x_to_z(x_values[i]))
        print("In conformal time this corresponds to", eta_of_x_scaled[i] / u.s, "years after the Big Bang \n")
        # print("This happens at a scale factor", from_x_to_a(x_values[i]), "and redshift", from_x_to_z(x_values[i]), "\n")

#For calculating the matter-dark energy equality time
for i in numbers:
    if np.abs(matter[i] - dark_energy[i]) < 0.2 and (matter[i] < (0.5 + 0.003) and matter[i] > (0.5 - 0.003) ):
        # print("The difference between matter and relativistic particles is", np.abs(matter[i] - relativistic[i]), "and this difference is found at index", i)
        # print("The x-value for the radiation-matter equality is", x_values[i])
        print("In (normal) time this equality occurs", t_of_x_scaled[i], "years after the Big Bang")
        # print("This happens at a scale factor", from_x_to_a(x_values[i]), "and redshift", from_x_to_z(x_values[i]))
        print("In conformal time this corresponds to", eta_of_x_scaled[i] / u.s, "years after the Big Bang \n")
        # print("This happens at a scale factor", from_x_to_a(x_values[i]), "and redshift", from_x_to_z(x_values[i]), "\n")
        
for i in numbers:
    if (matter[i] < (2*dark_energy[i]) + 0.02) and (matter[i] > (2*dark_energy[i]) - 0.02) and (dark_energy[i] > 10**(-5)):
        print("For start of acceleration the universe:")
        # print("The difference between matter and relativistic particles is", np.abs(matter[i] - relativistic[i]), "and this difference is found at index", i)
        # print("The x-value for the radiation-matter equality is", x_values[i])
        print("In (normal) time this equality occurs", t_of_x_scaled[i], "years after the Big Bang")
        # print("This happens at a scale factor", from_x_to_a(x_values[i]), "and redshift", from_x_to_z(x_values[i]))
        print("In conformal time this corresponds to", eta_of_x_scaled[i] / u.s, "years after the Big Bang \n")
        # print("This happens at a scale factor", from_x_to_a(x_values[i]), "and redshift", from_x_to_z(x_values[i]), "\n")


today_index = np.where(x_values == -0.00500501)
print("The horizon of the Universe since the Big Bang (conformal time) is", (eta_of_x[today_index] / (3.08*(10**25))), "Gpc")
print("The age of the Universe since the Big Bang (time) is", t_of_x_scaled[today_index] / (10**9), "Gyrs")
print("The age of the Universe since the Big Bang (conformal time) is", (eta_of_x_scaled[today_index]) / (u.s * 10**9), "Gyrs")
#These acceptance ranges took a lot of trial and error! It would also maybe be acceptable to end up with a small range of 
#possible values for the x at equality time, but this can also be changed later if needed and time.
#I will also calculate them analytically.

H0 = 2.17132 *pow(10,-18)
t0 = 13.78475393*(10**9)*60*60*24*365
omega_R = (3.81093*pow(10,-5)) + (5.50896*pow(10,-5)) 
omega_M = (0.05 + 0.267)
omega_Lambda = 0.682907
s_to_year = 60*60*24*365

x_scale_rad = np.linspace(-20, -8.131921133813, 1000)
x_scale_mat = np.linspace(-8.131921133813, -0.2558189708133, 1000)
x_scale_lam = np.linspace(-0.2558189708133, 5, 1000)

t_R = (((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*x_scale_rad))) - ((1/(2*H0*np.sqrt(omega_R))) * (np.exp(2*np.log(10**-8))))) / s_to_year # CORRECT!!!
t_M = (51064) + (((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*x_scale_mat))) - ((2/(3*H0*np.sqrt(omega_M))) * (np.exp((3/2)*-8.131921133813)))) / s_to_year # CORRECT!!! 
t_Lambda = (10.38 * 10**9) + ((1/(H0*np.sqrt(omega_Lambda)) * x_scale_lam) + (0.2558189708133/(H0*np.sqrt(omega_Lambda))))  / s_to_year #CORRECT FOR NOW!!!

eta_R= (((1 / (H0*np.sqrt(omega_R))) * np.exp(x_scale_rad)) - ((1 / (H0*np.sqrt(omega_R))) * np.exp(np.log(10**-8)))) / s_to_year #CORRECT!!!
eta_M = (51064) + (((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*x_scale_mat)) - ((2 / (H0*np.sqrt(omega_M))) * np.exp((1/2)*((-8.131921133813))))) / s_to_year #CORRECT I THINK!
eta_Lambda = (10.38 * 10**9) + (((-1/((H0)*np.sqrt(omega_Lambda)))*np.exp(-1*x_scale_lam)) + ((1/(H0*np.sqrt(omega_Lambda)))*np.exp(0.2558189708133))) / s_to_year #Not correct
# eta_Lambda = ((-1/(H0*np.sqrt(omega_Lambda)))*(np.exp(-1 * x_scale_lam) - np.exp(0.255))) / (60*60*24*365) #alternative but same epression

#Plotting the plot for eta(x), t(x), equality times, and beginning of acceleration
plt.figure(figsize = (16,8))
plt.plot(x, t_of_x_scaled, label = "Time t(x)")
plt.plot(x_scale_rad,t_R, label = "Radiation dominated era", linestyle = "dashed")
plt.plot(x_scale_mat,t_M, label = "Matter dominated era", linestyle = "dashed")
plt.plot(x_scale_lam, t_Lambda, label = "Lambda", linestyle = "dashed")
plt.plot(x_scale_rad, eta_R, linestyle = "dashed")
plt.plot(x_scale_mat, eta_M, label = "eta M", linestyle = "dashed")
plt.plot(x_scale_lam, eta_Lambda, label = "eta lambda", linestyle = "dashed")
plt.plot(x, (eta_of_x_scaled), label = "Conformal time $\eta(x)/c$")
plt.vlines(x = 0, ymin = 0, ymax = 10**13, color = "#808080", linestyle = "dotted", label = "Today")
plt.vlines(x = -8.13814, ymin = 0, ymax = 10**13, color = "#008000", linestyle = "dotted", label = "Matter - radiation \n equality num")
plt.vlines(x = -0.255255, ymin = 0, ymax = 10**13, color = "#be2ed6", linestyle = "dotted", label = "Matter - dark energy \n equality num")
plt.vlines(x = -0.48048, ymin = 0, ymax = 10**13, color = "red", linestyle = "dotted", label = "Start of accerelation num")
# plt.vlines(x = -8.131921133813, ymin = 0, ymax = 10**13, linestyle = "dotted", label = "Matter - radiation \n equality ana")
plt.ylabel("Time in years", fontsize = 15)
plt.xlabel("x", fontsize = 15)
plt.xlim(-14,5)
plt.ylim(10**-0.5, 10**13)
plt.yscale("log")
plt.legend(fontsize = "12")
# plt.savefig("python_plots/plots/equality_times.pdf")
plt.show()
    
    