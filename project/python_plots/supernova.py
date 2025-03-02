import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u

cosmology = np.loadtxt("../cosmology.txt")
supernova = np.loadtxt("../data/supernovadata.txt")

x = cosmology[:,0]
luminosity_distance_num = (cosmology[:,12] / (3.08*(10**9)*(10**16)))

y = luminosity_distance_num / ((1/np.exp(x)) -1  )

redshift = supernova[:,0]
luminosity_distance_data = supernova[:,1]
errors = supernova[:,2]

a_values = 1/(redshift + 1)
x_values = np.log(a_values)

plt.errorbar(x_values, luminosity_distance_data, errors, label = "Given supernova data")
plt.plot(x, luminosity_distance_num, label = "Numerical data")
plt.xlabel("x")
plt.ylabel("$d_L$")
plt.xlim(-1, 0)
plt.ylim(-3,12)
plt.legend(fontsize = 14)
plt.savefig("plots/supernovafits.pdf")
plt.show()
