
import numpy as np

data_C_l_hw = np.loadtxt("C_ells_hw.txt")
data_C_l_all_sf_terms_student = np.loadtxt("C_ells_all_sf_terms_student.txt")
data_C_l_first_sf_term_student = np.loadtxt("C_ells_first_sf_term_student.txt")
data_C_l_second_sf_term_student = np.loadtxt("C_ells_second_sf_term_student.txt")
data_C_l_third_sf_term_student = np.loadtxt("C_ells_third_sf_term_student.txt")
data_C_l_fourth_sf_term_student = np.loadtxt("C_ells_fourth_sf_term_student.txt")

x_axis_values = data_C_l_hw[:,0] #same for all data sets, so we can use this for all plots

C_l_data_hw = data_C_l_hw[:,1]
    
C_l_data_first_sf_term_student = data_C_l_first_sf_term_student[:,1]  
C_l_data_second_sf_term_student = data_C_l_second_sf_term_student[:,1]  
C_l_data_third_sf_term_student = data_C_l_third_sf_term_student[:,1]  
C_l_data_fourth_sf_term_student = data_C_l_fourth_sf_term_student[:,1]  
C_l_data_all_sf_terms_student = data_C_l_all_sf_terms_student[:,1] 
TE_data_all_sf_terms_student = data_C_l_all_sf_terms_student[:,2]
EE_data_all_sf_terms_student = data_C_l_all_sf_terms_student[:,3] 




# Comparison planck data
low_TT_comparison_data = np.loadtxt("python_plots/src_power_spectrum/low_TT_data_planck.txt")
low_TT_ell_values = low_TT_comparison_data[:,0]
low_TT_C_l_data = low_TT_comparison_data[:,1]
low_TT_C_l_err_up = low_TT_comparison_data[:,2]
low_TT_C_l_err_down = low_TT_comparison_data[:,3]

high_TT_comparison_data = np.loadtxt("python_plots/src_power_spectrum/high_TT_data_planck.txt")
high_TT_ell_values = high_TT_comparison_data[:,0]
high_TT_C_l_data = high_TT_comparison_data[:,1]
high_TT_C_l_err_down = high_TT_comparison_data[:,2]
high_TT_C_l_err_up = high_TT_comparison_data[:,3]

high_EE_comparison_data = np.loadtxt("python_plots/src_power_spectrum/high_EE_data_planck.txt")
high_EE_ell_values = high_EE_comparison_data[:,0]
high_EE_C_l_data = high_EE_comparison_data[:,1]
high_EE_C_l_err_down = high_EE_comparison_data[:,2]
high_EE_C_l_err_up = high_EE_comparison_data[:,3]

high_TE_comparison_data = np.loadtxt("python_plots/src_power_spectrum/high_TE_data_planck.txt")
high_TE_ell_values = high_TE_comparison_data[:,0]
high_TE_C_l_data = high_TE_comparison_data[:,1]
high_TE_C_l_err_down = high_TE_comparison_data[:,2]
high_TE_C_l_err_up = high_TE_comparison_data[:,3]

