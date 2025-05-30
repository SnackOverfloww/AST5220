
import numpy as np
import matplotlib.pyplot as plt

#Imports the categorized data from the text file
exec(open("python_plots/src_perturbations/reading_data_perturbations.py").read())
exec(open("python_plots/src_perturbations/calculating_horizon.py").read())


#=======================================================================================
# # Plots for Hans Winther cosmology 

# # Theta_0
# plt.plot(x_values, Theta_0_values_0_001_hw, label = "k = 0.001/Mpc", color = "#4b006e", zorder = 3)
# plt.plot(x_values, Theta_0_values_0_01_hw, label = "k = 0.01/Mpc", color = "#d9328a")
# plt.plot(x_values, Theta_0_values_0_1_hw, label = "k = 0.1/Mpc", color = "#f1b2e1", zorder = 0)
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "upper left")
# # plt.savefig("python_plots/plots/perturbations/theta_0_plot_hw.pdf")
# plt.show()

# #Theta_1
# plt.plot(x_values, Theta_1_values_0_001_hw, label = "k = 0.001/Mpc", color = "#4b006e", zorder = 3)
# plt.plot(x_values, Theta_1_values_0_01_hw, label = "k = 0.01/Mpc", color = "#d9328a")
# plt.plot(x_values, Theta_1_values_0_1_hw, label = "k = 0.1/Mpc", color = "#f1b2e1", zorder = 0)
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "upper left")
# plt.savefig("python_plots/plots/perturbations/theta_1_plot_hw.pdf")
# plt.show()

# #Phi
# plt.plot(x_values, Phi_values_0_001_hw, label = "k = 0.001/Mpc", color = "#4b006e", zorder = 3)
# plt.plot(x_values, Phi_values_0_01_hw, label = "k = 0.01/Mpc", color = "#d9328a")
# plt.plot(x_values, Phi_values_0_1_hw, label = "k = 0.1/Mpc", color = "#f1b2e1", zorder = 0)
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "lower left")
# plt.savefig("python_plots/plots/perturbations/phi_plot_hw.pdf")
# plt.show()

# #delta_cdm and delta_b
# plt.plot(x_values, np.abs(delta_cdm_values_0_001_hw), label = "k = 0.001/Mpc", color = "#4b006e", zorder = 0)
# plt.plot(x_values, np.abs(delta_cdm_values_0_01_hw), label = "k = 0.01/Mpc", color = "#d9328a", zorder = 1)
# plt.plot(x_values, np.abs(delta_cdm_values_0_1_hw), label = "k = 0.1/Mpc", color = "#f1b2e1", zorder = 2)
# plt.plot(x_values, np.abs(delta_b_values_0_001_hw), color = "#ab61cd", linestyle = "dashed", zorder = 3)
# plt.plot(x_values, np.abs(delta_b_values_0_01_hw), color = "#b7377b", linestyle = "dashed", zorder = 4)
# plt.plot(x_values, np.abs(delta_b_values_0_1_hw), color = "#d0a0c4", linestyle = "dashed", zorder = 5)
# plt.yscale("log")
# plt.ylim(1e-1, 1e5)
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "upper left")
# plt.savefig("python_plots/plots/perturbations/delta_cdm_delta_b_plot_hw.pdf")
# plt.show()

# #v_cdm and v_b
# plt.plot(x_values, np.abs(v_cdm_values_0_001_hw), label = "k = 0.001/Mpc", color = "#4b006e", zorder = 0)
# plt.plot(x_values, np.abs(v_cdm_values_0_01_hw), label = "k = 0.01/Mpc", color = "#d9328a", zorder = 1)
# plt.plot(x_values, np.abs(v_cdm_values_0_1_hw), label = "k = 0.1/Mpc", color = "#f1b2e1", zorder = 2)
# plt.plot(x_values, np.abs(v_b_values_0_001_hw), color = "#ab61cd", linestyle = "dashed", zorder = 3)
# plt.plot(x_values, np.abs(v_b_values_0_01_hw), color = "#b7377b", linestyle = "dashed", zorder = 4)
# plt.plot(x_values, np.abs(v_b_values_0_1_hw), color = "#d0a0c4", linestyle = "dashed", zorder = 5)
# plt.yscale("log")
# plt.ylim(1e-6, 1e2)
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "upper left")
# plt.savefig("python_plots/plots/perturbations/v_cdm_v_b_plot_hw.pdf")
# plt.show()
#=========================================================================================================================

#=========================================================================================================================
#Plots for student cosmology 

# delta_cdm, delta_b and delta_photons k=0.1, k = 0.01, k = 0.001
# fig, axs = plt.subplots(3, 1, figsize = (6,15))
# axs[0].plot(x_values, np.abs(4*Theta_0_values_0_1_student), label = "$\delta_\gamma$", color = "#f1b2e1")
# axs[0].plot(x_values, np.abs(delta_cdm_values_0_1_student), label = "$\delta_{cdm}$", color = "#4b006e")
# axs[0].plot(x_values, np.abs(delta_b_values_0_1_student), label = "$\delta_\mathrm{b}$", linestyle = "dashed", color = "#d9328a")
# axs[0].plot(x_values, np.abs(4*Nu_0_values_0_1_student), label = "$\delta_\\nu$", color = "#addfff", linestyle = "dashed")
# axs[0].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing", zorder = 0)
# axs[0].text(-18, 2e-1, "$k = 0.1$", fontsize = 11, color='black')

# axs[1].plot(x_values, np.abs(4*Theta_0_values_0_01_student), label = "$\delta_\gamma$", color = "#f1b2e1")
# axs[1].plot(x_values, np.abs(delta_cdm_values_0_01_student), label = "$\delta_\mathrm{CDM}$", color = "#4b006e")
# axs[1].plot(x_values, np.abs(delta_b_values_0_01_student), label = "$\delta_\mathrm{b}$", linestyle = "dashed", color = "#d9328a")
# axs[1].plot(x_values, np.abs(4*Nu_0_values_0_01_student), label = "$\delta_\\nu$", color = "#addfff", linestyle = "dashed")
# axs[1].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing", zorder = 0)
# axs[1].text(-18, 1.6e0, "$k = 0.01$", fontsize = 11, color='black')

# axs[2].plot(x_values, np.abs(4*Theta_0_values_0_001_student), label = "$\delta_\gamma$", color = "#f1b2e1")
# axs[2].plot(x_values, np.abs(delta_cdm_values_0_001_student), label = "$\delta_\mathrm{CDM}$", color = "#4b006e")
# axs[2].plot(x_values, np.abs(delta_b_values_0_001_student), label = "$\delta_\mathrm{b}$", linestyle = "dashed", color = "#d9328a")
# axs[2].plot(x_values, np.abs(4*Nu_0_values_0_001_student), label = "$\delta_\\nu$", color = "#addfff", linestyle = "dashed")
# axs[2].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing", zorder = 0)
# axs[2].text(-18, 1.3e0, "$k = 0.001$", fontsize = 11, color='black')



# for ax in axs:
#     ax.axvline(x = -6.9854, color = "black", linestyle = "dotted", label = "Photon decoupling")
#     ax.axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
#     ax.set_yscale("log")
#     ax.set_xlabel("x", fontsize = 13)
#     ax.set_xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
#     ax.tick_params(axis='both', labelsize=11)
#     ax.legend(fontsize = 9)

# plt.savefig("python_plots/plots/perturbations/delta_cdm_delta_b_delta_gamma_delta_nu_plot_student.pdf")
# plt.show()



# # v_cdm, v_b and v_photons k=0.1, k = 0.01, k = 0.001 
fig, axs = plt.subplots(3, 1, figsize = (6,15))
axs[0].plot(x_values, -3*Theta_1_values_0_1_student, label = "$v_\gamma$", color = "#f1b2e1")
axs[0].plot(x_values, np.abs(v_cdm_values_0_1_student), label = "$v_{cdm}$", color = "#4b006e")
axs[0].plot(x_values, np.abs(v_b_values_0_1_student), label = "$v_\mathrm{b}$", linestyle = "dashed", color = "#d9328a")
axs[0].plot(x_values, -3*Nu_1_values_0_1_student, label = "$v_\\nu$", color = "#addfff", linestyle = "dashed")
axs[0].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing", zorder = 0)
axs[0].text(-18, 1, "$k = 0.1$", fontsize = 11, color='black')


axs[1].plot(x_values, -3*Theta_1_values_0_01_student, label = "$v_\gamma$", color = "#f1b2e1")
axs[1].plot(x_values, np.abs(v_cdm_values_0_01_student), label = "$v_\mathrm{CDM}$", color = "#4b006e")
axs[1].plot(x_values, np.abs(v_b_values_0_01_student), label = "$v_\mathrm{b}$", linestyle = "dashed", color = "#d9328a")
axs[1].plot(x_values, -3*Nu_1_values_0_01_student, label = "$v_\\nu$", color = "#addfff", linestyle = "dashed")
axs[1].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing", zorder = 0)
axs[1].text(-18, 0.6, "$k = 0.01$", fontsize = 11, color='black')

axs[2].plot(x_values, -3*Theta_1_values_0_001_student, label = "$v_\gamma$", color = "#f1b2e1")
axs[2].plot(x_values, np.abs(v_cdm_values_0_001_student), label = "$v_\mathrm{CDM}$", color = "#4b006e")
axs[2].plot(x_values, np.abs(v_b_values_0_001_student), label = "$v_\mathrm{b}$", linestyle = "dashed", color = "#d9328a")
axs[2].plot(x_values, -3*Nu_1_values_0_001_student, label = "$v_\\nu$", color = "#addfff", linestyle = "dashed")
axs[2].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing", zorder = 0)
axs[2].text(-18, 0.11, "$k = 0.001$", fontsize = 11, color='black')


for ax in axs:
    ax.axvline(x = -6.9854, color = "black", linestyle = "dotted", label = "Photon decoupling")
    ax.axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
    ax.set_xlabel("x", fontsize = 13)
    ax.set_xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
    ax.tick_params(axis='both', labelsize=11)
    ax.legend(fontsize = 10)

# plt.savefig("python_plots/plots/perturbations/v_cdm_v_b_v_gamma_v_nu_plot_student.pdf")
plt.show()



# Temperature quadrupole and neutrino quadrupole
# fig, axs = plt.subplots(2, 1, figsize = (8,15))
# axs[0].plot(x_values, Theta_2_values_0_001_student, label = "k = 0.001", color = "#4b006e", zorder = 3)
# axs[0].plot(x_values, Theta_2_values_0_01_student, label = "k = 0.01", color = "#d9328a", zorder = 2)
# axs[0].plot(x_values, Theta_2_values_0_1_student, label = "k = 0.1", color = "#f1b2e1", zorder = 1)
# axs[0].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# axs[0].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 0)
# axs[0].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 0)
# axs[0].axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality", zorder = 0)
# axs[0].text(-1, 0.15, "$\Theta_2$", fontsize = 18, color='black')

# axs[1].plot(x_values, Nu_2_values_0_001_student, label = "k = 0.001", color = "#4b006e", zorder = 3)
# axs[1].plot(x_values, Nu_2_values_0_01_student, label = "k = 0.01", color = "#d9328a", zorder = 2)
# axs[1].plot(x_values, Nu_2_values_0_1_student, label = "k = 0.1", color = "#f1b2e1", zorder = 1)
# axs[1].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# axs[1].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 0)
# axs[1].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 0)
# axs[1].axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality", zorder = 0)
# axs[1].text(-1, 0.20, "$\mathcal{N}_2$", fontsize = 18, color='black')


# for ax in axs:
#     ax.set_xlabel("x", fontsize = 17)
#     ax.set_xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
#     ax.tick_params(axis='both', labelsize=15)
#     ax.legend(fontsize = 12)


# plt.savefig("python_plots/plots/perturbations/temperature_quadrupole_and_neutrino_quadrupole_plot_student.pdf")
# plt.show()




# # phi
# plt.plot(x_values, Phi_values_0_001_student, label = "k = 0.001", color = "#4b006e", zorder = 6)
# plt.plot(x_values, Phi_values_0_01_student, label = "k = 0.01", color = "#d9328a")
# plt.plot(x_values, Phi_values_0_1_student, label = "k = 0.1", color = "#f1b2e1", zorder = 4)
# plt.axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# plt.axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 1)
# plt.axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 2)
# plt.axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
# # plt.axvline(x = -0.256, color = "#3f3f3f", linestyle = "dotted", label = "Matter-dark energy equality")
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0], size = 11)
# plt.yticks(size = 11)
# plt.xlabel("x", fontsize = 13)
# plt.xlim(-18, 0)
# plt.legend(loc = "lower left", fontsize = 9)
# plt.savefig("python_plots/plots/perturbations/phi_plot_student.pdf")
# plt.show()



# # Psi + Phi subplots with AND without neutrinos and polarisation
# fig, axs = plt.subplots(2, 1, figsize = (10,14))
# axs[0].plot(x_values, Phi_values_0_001_student + Psi_values_0_001_student, label = "k = 0.001", color = "#4b006e", zorder = 3)
# axs[0].plot(x_values, Phi_values_0_01_student + Psi_values_0_01_student, label = "k = 0.01", color = "#d9328a", zorder = 2)
# axs[0].plot(x_values, Phi_values_0_1_student + Psi_values_0_1_student, label = "k = 0.1", color = "#f1b2e1", zorder = 1)
# axs[0].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# axs[0].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 0)
# axs[0].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 0)
# axs[0].axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
# axs[0].legend(loc = "lower left", fontsize = 14)


# axs[1].plot(x_values, Phi_values_0_001_student_no_n_or_p + Psi_values_0_001_student_no_n_or_p, label = "k = 0.001", color = "#4b006e", zorder = 3)
# axs[1].plot(x_values, Phi_values_0_01_student_no_n_or_p + Psi_values_0_01_student_no_n_or_p, label = "k = 0.01", color = "#d9328a", zorder = 2)
# axs[1].plot(x_values, Phi_values_0_1_student_no_n_or_p + Psi_values_0_1_student_no_n_or_p, label = "k = 0.1", color = "#f1b2e1", zorder = 1)
# axs[1].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# axs[1].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 0)
# axs[1].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 0)
# axs[1].axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
# axs[1].legend(loc = "upper left", fontsize = 14)


# for ax in axs:   
#     ax.set_xlabel("x", fontsize = 19)
#     ax.set_xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
#     ax.tick_params(axis='both', labelsize=17)

# plt.savefig("python_plots/plots/perturbations/phi_plus_psi_plot_student.pdf")
# plt.show()




# For the polarisation multipoles theta_p
# fig, axs = plt.subplots(3, 1, figsize = (6,15))
# axs[0].plot(x_values, Thetap_0_values_0_001_student, label = "k = 0.001", color = "#4b006e", zorder = 3)
# axs[0].plot(x_values, Thetap_0_values_0_01_student, label = "k = 0.01", color = "#d9328a", zorder = 2)
# axs[0].plot(x_values, Thetap_0_values_0_1_student, label = "k = 0.1", color = "#f1b2e1", zorder = 1)
# axs[0].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# axs[0].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 0)
# axs[0].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 0)
# axs[0].axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality", zorder = 0)
# axs[0].text(-18, 0.021, "${\Theta^P}_0$", fontsize=12, color='black')

# axs[1].plot(x_values, Thetap_1_values_0_001_student, label = "k = 0.001", color = "#4b006e", zorder = 3)
# axs[1].plot(x_values, Thetap_1_values_0_01_student, label = "k = 0.01", color = "#d9328a", zorder = 2)
# axs[1].plot(x_values, Thetap_1_values_0_1_student, label = "k = 0.1", color = "#f1b2e1", zorder = 1)
# axs[1].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# axs[1].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 0)
# axs[1].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 0)
# axs[1].axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality", zorder = 0)
# axs[1].text(-18, 0.0075, "${\Theta^P}_1$", fontsize=12, color='black')

# axs[2].plot(x_values, Thetap_2_values_0_001_student, label = "k = 0.001", color = "#4b006e", zorder = 3)
# axs[2].plot(x_values, Thetap_2_values_0_01_student, label = "k = 0.01", color = "#d9328a", zorder = 2)
# axs[2].plot(x_values, Thetap_2_values_0_1_student, label = "k = 0.1", color = "#f1b2e1", zorder = 1)
# axs[2].axvline(x = x_value_horizon_0_001, color = "#4f4f4f", linestyle = "dashed", label = "Horizon crossing k = 0.001", zorder = 0)
# axs[2].axvline(x = x_value_horizon_0_01, color = "#4f4f4f", linestyle = "solid", label = "Horizon crossing k = 0.01", zorder = 0)
# axs[2].axvline(x = x_value_horizon_0_1, color = "#4f4f4f", linestyle = (0, (3, 1, 1, 1, 1, 1)), label = "Horizon crossing k = 0.1", zorder = 0)
# axs[2].axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality", zorder = 0)
# axs[2].text(-18, 0.007, "${\Theta^P}_2$", fontsize=12, color='black')


# for ax in axs:
#     ax.set_xlabel("x", fontsize = 13)
#     ax.set_xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
#     ax.tick_params(axis='both', labelsize=11)
#     ax.legend(loc = "lower left", fontsize = 9)

# plt.savefig("python_plots/plots/perturbations/polarisation_multipoles_plot_student.pdf")
# plt.show()




#======================================================================================================================
# For testing

# # Theta_0 + Psi
# plt.figure(figsize=(10,6))
# plt.plot(x_values, Theta_0_values_0_3_student + Psi_values_0_3_student, label = "$\Theta_0 + \Psi$, $k = 0.3$", color = "#4b006e")
# plt.plot(x_values[:5700], -0.8*np.cos(0.3 * (1./3.08567758e22) * eta_values_0_3_student[:5700] / np.sqrt(3)), label = "$\Theta_0 + \Psi$ analytical, $k = 0.3$",  color = "#f1b2e1", linestyle = "dashed")
# plt.axvline(x = x_value_horizon_0_3, color = "black", linestyle = "dashed", label = "Horizon crossing")
# plt.axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
# plt.text(-17.4,-0.6, "$k=0.3$", fontsize = 18, color='black')
# plt.xticks([-18, -16, -14, -12, -10, -8, -6], size = 17)
# plt.yticks(size = 17)
# plt.xlabel("x", fontsize = 20)
# plt.xlim(-18,-6.5)
# plt.legend(loc = "upper left", fontsize = 15)
# plt.savefig("python_plots/plots/perturbations/testing_theta_0_plot_student.pdf")
# plt.show()


# # Analytical gravitational potential phi
# plt.plot(x_values, Phi_values_0_3_student, label = "$\Phi$", color = "#4b006e")
# y = (0.3 * (1./3.08567758e22) * eta_values_0_3_student) / np.sqrt(3) 
# plt.plot(x_values, Phi_values_0_3_student[0]*((np.sin(y)-(y*np.cos(y)))/((y**3)/3)), label = "$\Phi$ analytical", color = "#f1b2e1", linestyle = "dashed")
# plt.axvline(x = x_value_horizon_0_3, color = "black", linestyle = "dashed", label = "Horizon crossing")
# plt.axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
# plt.text(-3,0.05, "$k=0.3$", fontsize = 14, color='black')
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0], size = 12)
# plt.yticks(size = 12)
# plt.xlabel("x", fontsize = 15)
# plt.legend(loc = "upper right", fontsize = 11)
# plt.savefig("python_plots/plots/perturbations/testing_phi_analytical_plot_student.pdf")
# plt.show()

# # Analytical delta_CDM
# plt.plot(x_values, np.abs(delta_cdm_values_0_3_student), label = "$\delta_\mathrm{CDM}$", color = "#4b006e")
# # plt.plot(x_values[:3334], -0.065*(((x_values[:3334]))), label = "Analytical radiation dominated era", color = "#d9328a", linestyle = "dashed")
# plt.hlines(0.9*1, -18, -12, color = "#d9328a", linestyle = "dashed")
# plt.plot(x_values[6667:], 39000*np.exp(x_values[6667:]), label = "Analytical matter-dominated era", color = "#f1b2e1", linestyle = "dashed")
# plt.axvline(x = x_value_horizon_0_3, color = "black", linestyle = "dashed", label = "Horizon crossing")
# plt.axvline(x = -8.132, color = "#8e8e8e", linestyle = "dotted", label = "Matter-radiation equality")
# plt.text(-16.5,2e0, "$k=0.3$", fontsize = 14, color='black')
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0], size = 12)
# plt.yticks(size = 12)
# plt.xlabel("x", fontsize = 15)
# plt.legend(loc = "upper left", fontsize = 11)
# plt.yscale("log")
# # plt.savefig("python_plots/plots/perturbations/testing_delta_cdm_analytical_plot_student.pdf")
# plt.show()

# plt.plot(x_values, source_values_0_1_testing, color = "red")
# plt.xlim(-8, 0)
# plt.show()






#==============================================================================================================================================
# Might use later
# delta_cdm and delta_b no neutrinos
# plt.plot(x_values, np.abs(delta_cdm_values_0_001_student_no_n_or_p), label = "k = 0.001/Mpc", color = "#4b006e", zorder = 0)
# plt.plot(x_values, np.abs(delta_cdm_values_0_01_student_no_n_or_p), label = "k = 0.01/Mpc", color = "#d9328a", zorder = 1)
# plt.plot(x_values, np.abs(delta_cdm_values_0_1_student_no_n_or_p), label = "k = 0.1/Mpc", color = "#f1b2e1", zorder = 2)
# plt.plot(x_values, np.abs(delta_b_values_0_001_student_no_n_or_p), color = "#ab61cd", linestyle = "dashed", zorder = 3)
# plt.plot(x_values, np.abs(delta_b_values_0_01_student_no_n_or_p), color = "#b7377b", linestyle = "dashed", zorder = 4)
# plt.plot(x_values, np.abs(delta_b_values_0_1_student_no_n_or_p), color = "#d0a0c4", linestyle = "dashed", zorder = 5)
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "upper left")
# plt.yscale("log")
# # plt.savefig("python_plots/plots/perturbations/theta_0_plot.pdf")
# plt.show()


# # neutrino quadrupole
# plt.plot(x_values, Theta_0_values_0_001_student_no_n_or_p, label = "k = 0.001/Mpc", color = "#4b006e", zorder = 3)
# plt.plot(x_values, Theta_0_values_0_01_student_no_n_or_p, label = "k = 0.01/Mpc", color = "#d9328a")
# plt.plot(x_values, Theta_0_values_0_1_student_no_n_or_p, label = "k = 0.1/Mpc", color = "#f1b2e1", zorder = 0)
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "upper left")
# # plt.savefig("python_plots/plots/perturbations/neutrino_quadrupole_plot_student.pdf")
# plt.show()

# plt.plot(x_values, -Psi_values_0_001_student, label = "k = 0.001/Mpc", color = "#4b006e")
# plt.plot(x_values, -Psi_values_0_01_student, label = "k = 0.01/Mpc", color = "#d9328a")
# plt.plot(x_values, -Psi_values_0_1_student, label = "k = 0.1/Mpc", color = "#f1b2e1")
# plt.plot(x_values, Phi_values_0_001_student, label = "k = 0.001/Mpc", color = "#4b006e")
# plt.plot(x_values, Phi_values_0_01_student, label = "k = 0.01/Mpc", color = "#d9328a")
# plt.plot(x_values, Phi_values_0_1_student, label = "k = 0.1/Mpc", color = "#f1b2e1")
# plt.xticks([-18, -16, -14, -12, -10, -8, -6, -4, -2, 0])
# plt.legend(loc = "upper left")
# # plt.savefig("python_plots/plots/perturbations/neutrino_quadrupole_plot_student.pdf")
# plt.show()