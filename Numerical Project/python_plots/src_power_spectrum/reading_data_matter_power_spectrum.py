import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("matter_power_spectrum.txt")
comparison_data_galaxy_survey = np.genfromtxt("python_plots/src_power_spectrum/power_spectrum_data_compare_galaxy_survey.txt")
comparison_data_cmb = np.genfromtxt("python_plots/src_power_spectrum/power_spectrum_data_compare_cmb.txt")
comparison_data_lyalpha = np.genfromtxt("python_plots/src_power_spectrum/power_spectrum_data_compare_lyalpha.txt")


x = data[:,0]
y = data[:,1]

x_comp_galaxy_survey = comparison_data_galaxy_survey[:,0]
y_comp_galaxy_survey= comparison_data_galaxy_survey[:,1]
error_comp_galaxy_survey = comparison_data_galaxy_survey[:,2]


x_comp_cmb = comparison_data_cmb[:,0]
y_comp_cmb = comparison_data_cmb[:,1]
error_comp_cmb = (comparison_data_cmb[:,2] - comparison_data_cmb[:,1])/2

x_comp_lyalpha = comparison_data_lyalpha[:,0]
y_comp_lyalpha = comparison_data_lyalpha[:,1]
error_comp_lyalpha = (comparison_data_lyalpha[:,2] - comparison_data_lyalpha[:,1])/2
