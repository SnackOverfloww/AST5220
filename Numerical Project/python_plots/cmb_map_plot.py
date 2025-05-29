
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import colormaps as cmaps
from matplotlib.colors import LinearSegmentedColormap

exec(open("python_plots/src_power_spectrum/reading_data_power_spectrum.py").read())

my_colors = ["#1B1E85", "#15AEE3", "#36B240", "#DDE629", "#FA7119", "#DC3D17", "#AD001E"]
my_cmap = LinearSegmentedColormap.from_list('custom_cmap', my_colors)
N = 2048
m = C_l_data_all_sf_terms_student.tolist()
print(m[0])
# print(len(m))

normalised = []
ells = np.arange(2,2001)

for i in ells:
    # print(i)
    # print(m[i-2])
    norm = m[i-2]*((2*np.pi) / (i*(i+1)))
    # print(norm)
    normalised.append(norm)
# print(normalised)


final_normalised_array = [0,0] + normalised
# print(final_normalised_array) 


hey = hp.sphtfunc.synfast(final_normalised_array, N, lmax = 2000)

sigma = np.std(hey)

hp.mollview(hey, rot=[0, 0, 0], cmap = "jet", min = -3*sigma, max = 3*sigma, title = " ", cbar = False)
# hp.mollview(hey, rot=[0, 0, 0], cmap = "jet")
# hp.graticule()
plt.savefig("python_plots/plots/power_spectrum/cmb_map_student.pdf")
plt.show()

plt.plot(np.arange(0,1999), m)
plt.xscale("log")
# plt.show()