import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u

cosmology = np.loadtxt("../cosmology.txt")
supernova = np.loadtxt("../data/supernovadata.txt")

x = cosmology[:,0]
luminosity_distance_num = (cosmology[:,12] / (3.08*(10**9)*(10**16)))

luminosity_distance_scaled = luminosity_distance_num / ((1/np.exp(x)) -1  )

redshift = supernova[:,0]
luminosity_distance_data = supernova[:,1] / redshift
errors = supernova[:,2] / redshift

a_values = 1/(redshift + 1)
x_values = np.log(a_values)

z_values = (1/np.exp(x))-1  

plt.scatter(redshift, luminosity_distance_data, label = "Given supernova data", color = "#7a3e8d")
plt.errorbar(redshift, luminosity_distance_data, yerr = errors, linestyle="None", color = "#7a3e8d")
plt.plot(z_values, luminosity_distance_scaled, label = "Numerical data", color = "#d30e92")
plt.xlabel("z", fontsize = 12)
plt.ylabel("$d_L [Gpc]$", fontsize = 12)
plt.xscale("log")
plt.xlim(2*1e0, 5*1e-3)
plt.ylim(3.4,8.5)
plt.tick_params(labelsize = 12)
plt.legend(fontsize = 14)
# plt.savefig("plots/supernovafits.pdf")
plt.show()
