import numpy as np


#Load all textfiles from Hans Winther cosmology
spline_data_0_001_hw = np.loadtxt("perturbations_k0.001_hw.txt")
spline_data_0_01_hw = np.loadtxt("perturbations_k0.01_hw.txt")
spline_data_0_1_hw = np.loadtxt("perturbations_k0.1_hw.txt")

#Load all textfiles from student cosmology
spline_data_0_001_student = np.loadtxt("perturbations_k0.001_student.txt")
spline_data_0_01_student = np.loadtxt("perturbations_k0.01_student.txt")
spline_data_0_1_student = np.loadtxt("perturbations_k0.1_student.txt")
spline_data_0_3_student = np.loadtxt("perturbations_k0.3_student.txt")
spline_data_3_student = np.loadtxt("perturbations_k3_student.txt")
spline_data_10_student = np.loadtxt("perturbations_k10_student.txt")


spline_data_0_001_student_no_n_or_p = np.loadtxt("perturbations_k0.001_student_no_n_or_p.txt")
spline_data_0_01_student_no_n_or_p = np.loadtxt("perturbations_k0.01_student_no_n_or_p.txt")
spline_data_0_1_student_no_n_or_p = np.loadtxt("perturbations_k0.1_student_no_n_or_p.txt")


spline_data_0_1_testing= np.loadtxt("perturbations_k0.1_testing.txt")


# Array of x-values - these are the same for all textfiles
x_values = spline_data_0_001_student[:,0]

#Categorises the data from the textfile for k=0.001/Mpc Hans Winther cosmology
Theta_0_values_0_001_hw = spline_data_0_001_hw[:,1]
Theta_1_values_0_001_hw = spline_data_0_001_hw[:,2]
Theta_2_values_0_001_hw = spline_data_0_001_hw[:,3]
Phi_values_0_001_hw = spline_data_0_001_hw[:,4]
delta_cdm_values_0_001_hw = spline_data_0_001_hw[:,5]
delta_b_values_0_001_hw = spline_data_0_001_hw[:,6]
v_cdm_values_0_001_hw = spline_data_0_001_hw[:,7]
v_b_values_0_001_hw = spline_data_0_001_hw[:,8]
Psi_values_0_001_hw = spline_data_0_001_hw[:,9]
Pi_values_0_001_hw = spline_data_0_001_hw[:,10]

#Categorises the data from the textfile for k=0.01/Mpc from Hans Winther cosmology
Theta_0_values_0_01_hw = spline_data_0_01_hw[:,1]
Theta_1_values_0_01_hw = spline_data_0_01_hw[:,2]
Theta_2_values_0_01_hw = spline_data_0_01_hw[:,3]
Phi_values_0_01_hw = spline_data_0_01_hw[:,4]
delta_cdm_values_0_01_hw = spline_data_0_01_hw[:,5]
delta_b_values_0_01_hw = spline_data_0_01_hw[:,6]
v_cdm_values_0_01_hw = spline_data_0_01_hw[:,7]
v_b_values_0_01_hw = spline_data_0_01_hw[:,8]
Psi_values_0_01_hw = spline_data_0_01_hw[:,9]
Pi_values_0_01_hw = spline_data_0_01_hw[:,10]

#Categorises the data from the textfile for k=0.1/Mpc from Hans Winther cosmology
Theta_0_values_0_1_hw = spline_data_0_1_hw[:,1]
Theta_1_values_0_1_hw = spline_data_0_1_hw[:,2]
Theta_2_values_0_1_hw = spline_data_0_1_hw[:,3]
Phi_values_0_1_hw = spline_data_0_1_hw[:,4]
delta_cdm_values_0_1_hw = spline_data_0_1_hw[:,5]
delta_b_values_0_1_hw = spline_data_0_1_hw[:,6]
v_cdm_values_0_1_hw = spline_data_0_1_hw[:,7]
v_b_values_0_1_hw = spline_data_0_1_hw[:,8]
Psi_values_0_1_hw = spline_data_0_1_hw[:,9]
Pi_values_0_1_hw = spline_data_0_1_hw[:,10]

#Categorises the data from the textfile for k=0.001/Mpc from student cosmology
Theta_0_values_0_001_student = spline_data_0_001_student[:,1]
Theta_1_values_0_001_student = spline_data_0_001_student[:,2]
Theta_2_values_0_001_student = spline_data_0_001_student[:,3]
Phi_values_0_001_student = spline_data_0_001_student[:,4]
delta_cdm_values_0_001_student = spline_data_0_001_student[:,5]
delta_b_values_0_001_student = spline_data_0_001_student[:,6]
v_cdm_values_0_001_student = spline_data_0_001_student[:,7]
v_b_values_0_001_student = spline_data_0_001_student[:,8]
Psi_values_0_001_student = spline_data_0_001_student[:,9]
Pi_values_0_001_student = spline_data_0_001_student[:,10]
Nu_0_values_0_001_student = spline_data_0_001_student[:,11]
Nu_1_values_0_001_student = spline_data_0_001_student[:,12]
Nu_2_values_0_001_student = spline_data_0_001_student[:,13]
Thetap_0_values_0_001_student = spline_data_0_001_student[:,14]
Thetap_1_values_0_001_student = spline_data_0_001_student[:,15]
Thetap_2_values_0_001_student = spline_data_0_001_student[:,16]
eta_values_0_001_student = spline_data_0_001_student[:,17]

#Categorises the data from the textfile for k=0.01/Mpc from student cosmology
Theta_0_values_0_01_student = spline_data_0_01_student[:,1]
Theta_1_values_0_01_student = spline_data_0_01_student[:,2]
Theta_2_values_0_01_student = spline_data_0_01_student[:,3]
Phi_values_0_01_student = spline_data_0_01_student[:,4]
delta_cdm_values_0_01_student = spline_data_0_01_student[:,5]
delta_b_values_0_01_student = spline_data_0_01_student[:,6]
v_cdm_values_0_01_student = spline_data_0_01_student[:,7]
v_b_values_0_01_student = spline_data_0_01_student[:,8]
Psi_values_0_01_student = spline_data_0_01_student[:,9]
Pi_values_0_01_student = spline_data_0_01_student[:,10]
Nu_0_values_0_01_student = spline_data_0_01_student[:,11]
Nu_1_values_0_01_student = spline_data_0_01_student[:,12]
Nu_2_values_0_01_student = spline_data_0_01_student[:,13]
Thetap_0_values_0_01_student = spline_data_0_01_student[:,14]
Thetap_1_values_0_01_student = spline_data_0_01_student[:,15]
Thetap_2_values_0_01_student = spline_data_0_01_student[:,16]
eta_values_0_01_student = spline_data_0_01_student[:,17]

#Categorises the data from the textfile for k=0.1/Mpc from student cosmology
Theta_0_values_0_1_student = spline_data_0_1_student[:,1]
Theta_1_values_0_1_student = spline_data_0_1_student[:,2]
Theta_2_values_0_1_student = spline_data_0_1_student[:,3]
Phi_values_0_1_student = spline_data_0_1_student[:,4]
delta_cdm_values_0_1_student = spline_data_0_1_student[:,5]
delta_b_values_0_1_student = spline_data_0_1_student[:,6]
v_cdm_values_0_1_student = spline_data_0_1_student[:,7]
v_b_values_0_1_student = spline_data_0_1_student[:,8]
Psi_values_0_1_student = spline_data_0_1_student[:,9]
Pi_values_0_1_student = spline_data_0_1_student[:,10]
Nu_0_values_0_1_student = spline_data_0_1_student[:,11]
Nu_1_values_0_1_student = spline_data_0_1_student[:,12]
Nu_2_values_0_1_student = spline_data_0_1_student[:,13]
Thetap_0_values_0_1_student = spline_data_0_1_student[:,14]
Thetap_1_values_0_1_student = spline_data_0_1_student[:,15]
Thetap_2_values_0_1_student = spline_data_0_1_student[:,16]
eta_values_0_1_student = spline_data_0_1_student[:,17]


# Categorises the data from the textfile for k=0.001/Mpc from student cosmology no neutrinos or polarization
Theta_0_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,1]
Theta_1_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,2]
Theta_2_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,3]
Phi_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,4]
delta_cdm_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,5]
delta_b_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,6]
v_cdm_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,7]
v_b_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,8]
Psi_values_0_001_student_no_n_or_p = spline_data_0_001_student_no_n_or_p[:,9]
Pi_values_0_001_student_no_n_or_p= spline_data_0_001_student_no_n_or_p[:,10]


# Categorises the data from the textfile for k=0.001/Mpc from student cosmology no neutrinos or polarization
Theta_0_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,1]
Theta_1_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,2]
Theta_2_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,3]
Phi_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,4]
delta_cdm_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,5]
delta_b_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,6]
v_cdm_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,7]
v_b_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,8]
Psi_values_0_01_student_no_n_or_p = spline_data_0_01_student_no_n_or_p[:,9]
Pi_values_0_01_student_no_n_or_p= spline_data_0_01_student_no_n_or_p[:,10]


# Categorises the data from the textfile for k=0.001/Mpc from student cosmology no neutrinos or polarization
Theta_0_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,1]
Theta_1_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,2]
Theta_2_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,3]
Phi_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,4]
delta_cdm_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,5]
delta_b_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,6]
v_cdm_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,7]
v_b_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,8]
Psi_values_0_1_student_no_n_or_p = spline_data_0_1_student_no_n_or_p[:,9]
Pi_values_0_1_student_no_n_or_p= spline_data_0_1_student_no_n_or_p[:,10]



#===========================================================================================
# Extra for testing
# Categorises the data from the textfile for k=0.3/Mpc from student cosmology
Theta_0_values_0_3_student = spline_data_0_3_student[:,1]
Theta_1_values_0_3_student = spline_data_0_3_student[:,2]
Theta_2_values_0_3_student = spline_data_0_3_student[:,3]
Phi_values_0_3_student = spline_data_0_3_student[:,4]
delta_cdm_values_0_3_student = spline_data_0_3_student[:,5]
delta_b_values_0_3_student = spline_data_0_3_student[:,6]
v_cdm_values_0_3_student = spline_data_0_3_student[:,7]
v_b_values_0_3_student = spline_data_0_3_student[:,8]
Psi_values_0_3_student = spline_data_0_3_student[:,9]
Pi_values_0_3_student = spline_data_0_3_student[:,10]
Nu_0_values_0_3_student = spline_data_0_3_student[:,11]
Nu_1_values_0_3_student = spline_data_0_3_student[:,12]
Nu_2_values_0_3_student = spline_data_0_3_student[:,13]
eta_values_0_3_student = spline_data_0_3_student[:,14]

# Categorises the data from the textfile for k=0.3/Mpc from student cosmology
Theta_0_values_3_student = spline_data_3_student[:,1]
Theta_1_values_3_student = spline_data_3_student[:,2]
Theta_2_values_3_student = spline_data_3_student[:,3]
Phi_values_3_student = spline_data_3_student[:,4]
delta_cdm_values_3_student = spline_data_3_student[:,5]
delta_b_values_3_student = spline_data_3_student[:,6]
v_cdm_values_3_student = spline_data_3_student[:,7]
v_b_values_3_student = spline_data_3_student[:,8]
Psi_values_3_student = spline_data_3_student[:,9]
Pi_values_3_student = spline_data_3_student[:,10]
eta_values_3_student = spline_data_3_student[:,14]

Theta_0_values_10_student = spline_data_10_student[:,1]
Theta_1_values_10_student = spline_data_10_student[:,2]
Theta_2_values_10_student = spline_data_10_student[:,3]
Phi_values_10_student = spline_data_10_student[:,4]
delta_cdm_values_10_student = spline_data_10_student[:,5]
delta_b_values_10_student = spline_data_10_student[:,6]
v_cdm_values_10_student = spline_data_10_student[:,7]
v_b_values_10_student = spline_data_10_student[:,8]
Psi_values_10_student = spline_data_10_student[:,9]
Pi_values_10_student = spline_data_10_student[:,10]
eta_values_10_student = spline_data_10_student[:,14]




#=======================================================================
Theta_0_values_0_1_testing = spline_data_0_1_testing[:,1]
Theta_1_values_0_1_testing = spline_data_0_1_testing[:,2]
Theta_2_values_0_1_testing = spline_data_0_1_testing[:,3]
Phi_values_0_1_testing = spline_data_0_1_testing[:,4]
delta_cdm_values_0_1_testing = spline_data_0_1_testing[:,5]
delta_b_values_0_1_testing = spline_data_0_1_testing[:,6]
v_cdm_values_0_1_testing = spline_data_0_1_testing[:,7]
v_b_values_0_1_testing = spline_data_0_1_testing[:,8]
Psi_values_0_1_testing = spline_data_0_1_testing[:,9]
Pi_values_0_1_testing = spline_data_0_1_testing[:,10]
Nu_0_values_0_1_testing = spline_data_0_1_testing[:,11]
Nu_1_values_0_1_testing = spline_data_0_1_testing[:,12]
Nu_2_values_0_1_testing = spline_data_0_1_testing[:,13]
Thetap_0_values_0_1_testing = spline_data_0_1_testing[:,14]
Thetap_1_values_0_1_testing = spline_data_0_1_testing[:,15]
Thetap_2_values_0_1_testing = spline_data_0_1_testing[:,16]
eta_values_0_1_testing = spline_data_0_1_testing[:,17]
source_values_0_1_testing = spline_data_0_1_testing[:,18]