
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# plt.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize": 12})

exec(open("src_recombination/calculating_recombination.py").read())

z_reionisation = 7.8 #change this if you want a different redshift for the reionisation

# #=========================================================================================================================================================================
# # Plot Xe of x
# fig, ax = plt.subplots(2,1, figsize = (8,12))
# ax[0].plot(x_no_reionisation,Xe_no_reionisation, label = "$X_e(x)$", color = "#d30e92")
# ax[0].plot(x_with_reionisation,Xe_with_reionisation, linestyle=(0, (2, 2)), linewidth=2, color = "#f6acc5", label = "$X_e(x)$ with reionisation", zorder = 3)
# ax[0].plot(x_just_saha, Xe_just_saha, color = "#cbc3e3", label = "$X_e(x)$ only Saha")
# ax[0].axvline(x = x_value_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[0].axvline(x = x_value_no_reionisation_visibility_peak, linestyle = "dotted", color = "#949494", label = "Decoupling (peak $\\tilde{g}(x)$)")
# ax[0].axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# ax[0].axhline(y = 0.1, color = "#bbbbbb", linestyle = "dashed", zorder = 1)
# ax[0].set_xlabel("x", fontsize = 14)
# ax[0].set_ylabel("$X_e(x)$", fontsize = 14)
# ax[0].set_yscale("log")
# plt.xticks(fontsize = 12)
# plt.yticks(fontsize = 12)
# ax[0].set_xlim(-12, 0)
# ax[0].set_ylim(1e-4, 1e1)
# ax[0].legend(fontsize = 10, loc = "lower left")

# ax[1].plot(x_no_reionisation,Xe_no_reionisation, label = "$X_e(x)$", color = "#d30e92")
# ax[1].plot(x_with_reionisation,Xe_with_reionisation, linestyle=(0, (2, 2)), linewidth=2, color = "#f6acc5", label = "$X_e(x)$ with reionisation", zorder = 3)
# ax[1].plot(x_just_saha, Xe_just_saha, color = "#cbc3e3", label = "$X_e(x)$ only Saha")
# ax[1].axvline(x = x_value_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[1].axvline(x = x_value_no_reionisation_visibility_peak, linestyle = "dotted", color = "#949494", label = "Decoupling (peak $\\tilde{g}(x)$)")
# ax[1].set_xlabel("x", fontsize = 14)
# ax[1].set_ylabel("$X_e(x)$", fontsize = 14)
# ax[1].set_yscale("log")
# plt.xticks(fontsize = 12)
# plt.yticks(fontsize = 12)
# ax[1].set_xlim(-7.01, -6.96)
# ax[1].set_ylim(1e-3, 1)
# ax[1].legend(fontsize = 10, loc = "lower left")
# # plt.savefig("plots/xe_of_x.pdf")
# plt.show()
# # #===================================================================================================================================

# # #======================================================================================================================
# #Plot tau of x and its derivatives (WITH REIONISATION)
# plt.figure(figsize = (12,8))
# plt.plot(x_with_reionisation,tau_with_reionisation, label = "$\\tau(x)$ with reionisation", color = "#d30e92")
# plt.plot(x_with_reionisation, -tau_deriv_with_reionisation, label = "$-\\tau\\,'(x)$ with reionisation", color = "#7a3e8d")
# plt.plot(x_with_reionisation, tau_double_deriv_with_reionisation, label = "$\\tau\\,''(x)$ with reionisation", color = "#f6acc5")
# plt.plot(x_with_reionisation, tau_just_saha, label = "$\\tau(x)$ only Saha", color = "#cbc3e3")
# plt.axvline(x = x_value_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# plt.axvline(x = x_value_no_reionisation_visibility_peak, linestyle = "dotted", color = "#949494", label = "Decoupling (peak $\\tilde{g}(x)$)")
# plt.axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# plt.axhline(y = 1, color = "#bbbbbb", linestyle = "dashed", zorder = 1)
# plt.xlabel("x", fontsize = 20)
# plt.xlim(-12, 0)
# plt.ylim(10**-8, 10**8)
# plt.xticks(fontsize = 16)
# plt.yticks(fontsize = 16)
# plt.yscale("log")
# plt.legend(fontsize = 15)

# # plt.savefig("plots/optical_depth_reionisation.pdf")
# plt.show()
# # #======================================================================================================================


# #======================================================================================================================
# # Subplot g_tilde of x and its derivatives
# fig, ax = plt.subplots(3,1, figsize = (12,20))
# ax[0].plot(x_no_reionisation,g_tilde_no_reionisation, label = "$\\tilde{g}(x)$", color = "#d30e92")
# ax[0].plot(x_with_reionisation,g_tilde_with_reionisation, linestyle=(0, (2, 2)), linewidth=2, color = "#f6acc5", label = "$\\tilde{g}(x)$ with reionisation")
# # ax[0].plot(x_no_reionisation, g_tilde_just_saha, color = "#cbc3e3", zorder = 1)
# ax[0].axvline(x = x_value_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[0].axvline(x = x_value_no_reionisation_visibility_peak, linestyle = "dotted", color = "#949494", label = "Decoupling (peak $\\tilde{g}(x)$)")
# ax[0].axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# ax[0].set_xlabel("x", fontsize = 18)
# ax[0].set_ylabel("$\\tilde{g}(x)$", fontsize = 18)
# ax[0].set_ylim(-0.5, 5)
# ax[0].tick_params(labelsize = 18)
# ax[0].legend(fontsize = 15)


# ax[1].plot(x_no_reionisation,g_tilde_deriv_no_reionisation, label = "$\\tilde{g}\\,'(x)$", color = "#d30e92")
# ax[1].plot(x_with_reionisation,g_tilde_deriv_with_reionisation, linestyle=(0, (2, 2)), label = "$\\tilde{g}\\,'(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# ax[1].plot(x_just_saha, g_tilde_deriv_just_saha, color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}\\,'(x)$ only Saha")
# ax[1].axvline(x = x_value_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[1].axvline(x = x_value_no_reionisation_visibility_peak, linestyle = "dotted", color = "#949494", label = "Decoupling (peak $\\tilde{g}(x)$)")
# ax[1].axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# ax[1].set_xlabel("x", fontsize = 16)
# ax[1].set_ylabel("$\\tilde{g}\\,'(x)$", fontsize = 16)
# ax[1].set_ylim(-32,52)
# ax[1].tick_params(labelsize = 18)
# ax[1].legend(fontsize = 14)

# ax[2].plot(x_no_reionisation,g_tilde_double_deriv_no_reionisation, label = "$\\tilde{g}\\,''(x)$", color = "#d30e92")
# ax[2].plot(x_with_reionisation,g_tilde_double_deriv_with_reionisation, linestyle=(0, (2, 2)), label = "$\\tilde{g}\\,''(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# ax[2].plot(x_just_saha, g_tilde_double_deriv_just_saha, color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}\\,''(x)$ only Saha")
# ax[2].axvline(x = x_value_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[2].axvline(x = x_value_no_reionisation_visibility_peak, linestyle = "dotted", color = "#949494", label = "Decoupling (peak $\\tilde{g}(x)$)")
# ax[2].axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# ax[2].set_xlabel("x", fontsize = 16)
# ax[2].set_ylabel("$\\tilde{g}\\,''(x)$", fontsize = 16)
# ax[2].tick_params(labelsize = 18)
# ax[2].set_ylim(-1200, 1000)
# ax[2].legend(fontsize = 14)

# # plt.savefig("plots/g_tilde_and_derivatives.pdf")
# plt.show()
# #======================================================================================================================

# #======================================================================================================================
# #Plot g_tilde of x smaller peak
# plt.figure()
# plt.plot(x_no_reionisation,g_tilde_no_reionisation, label = "$\\tilde{g}(x)$", color = "#d30e92")
# plt.plot(x_with_reionisation,g_tilde_with_reionisation, linestyle=(0, (2, 2)), label = "$\\tilde{g}(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# plt.plot(x_just_saha, g_tilde_just_saha, color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}(x)$ only Saha")
# plt.axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# plt.xlabel("x", fontsize = 12)
# plt.ylabel("$\\tilde{g}(x)$", fontsize = 12)
# plt.xlim(-3,0)
# plt.ylim(-0.1, 0.3)
# plt.xticks(fontsize = 11)
# plt.yticks(fontsize = 11)
# plt.legend()
# plt.savefig("plots/g_tilde_small_peak.pdf")
# plt.show()
# #======================================================================================================================


# #======================================================================================================================
# #Plot g_tilde of x first derivative smaller peak
# plt.plot(x_no_reionisation,g_tilde_deriv_no_reionisation, label = "$\\tilde{g}\\,'(x)$", color = "#d30e92")
# plt.plot(x_with_reionisation,g_tilde_deriv_with_reionisation, linestyle=(0, (2, 2)), label = "$\\tilde{g}\\,'(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# plt.plot(x_just_saha, g_tilde_deriv_just_saha, color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}\\,'(x)$ only Saha")
# plt.axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# plt.xlabel("x", fontsize = 14)
# plt.ylabel("$\\tilde{g}\\,'(x)$", fontsize = 14)
# plt.xlim(-2.5, -1.8)
# plt.ylim(-4,6)
# plt.xticks(fontsize = 13)
# plt.yticks(fontsize = 13)
# plt.legend(fontsize = 10)
# # plt.savefig("plots/g_tilde_deriv_small_peak.pdf")
# plt.show()
# #======================================================================================================================


# #======================================================================================================================
# #Plot g_tilde of x second derivative smaller peak

# plt.plot(x_no_reionisation,g_tilde_double_deriv_no_reionisation, label = "$\\tilde{g}\\,'(x)$", color = "#d30e92")
# plt.plot(x_with_reionisation,g_tilde_double_deriv_with_reionisation, linestyle=(0, (2, 2)), label = "$\\tilde{g}\\,'(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# plt.plot(x_just_saha, g_tilde_double_deriv_just_saha, color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}\\,'(x)$ only Saha")
# plt.axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# plt.xlabel("x", fontsize = 12)
# plt.ylabel("$\\tilde{g}\\,''(x)$", fontsize = 12)
# plt.xlim(-2.4,-1.95)
# plt.ylim(-100, 100)
# plt.xticks(fontsize = 11)
# plt.yticks(fontsize = 11)
# plt.legend()
# # plt.savefig("plots/g_tilde_double_deriv_small_peak.pdf")
# plt.show()
# #======================================================================================================================

# #=========================================================================================================
# #Plotting sound horizon
# plt.figure(figsize = (9,6))
# plt.plot(x_no_reionisation, sound_horizon * const_m_to_Mpc, color = "#d30e92", label = "Sound horizon")
# plt.axvline(x = x_value_no_reionisation_visibility_peak, linestyle = "dotted", color = "#323232", label = "Decoupling (peak $\\tilde{g}(x)$)") #marks decoupling
# plt.axhline(y = sound_horizon_in_Mpc, color = "#bbbbbb", linestyle = "dashed", zorder = 1) 
# plt.yscale("log")
# plt.xlabel("x", fontsize = 16)
# plt.ylabel("$s(x)$ [Mpc]", fontsize = 16)

# plt.yticks([1, 5, 10, 50, 145.57, 500], labels=[1, 5, 10, 50, (str.format('{0:.2f}', 145.57)), 500])  # Change to desired value
# plt.tick_params(labelsize = 14)
# plt.xlim(-12, 0)
# plt.legend(fontsize = 14)
# # plt.savefig("plots/sound_horizon.pdf")
# plt.show()
# #==========================================================================================================



# #====================================================================================================
# #NOT USED IN REPORT:

# #======================================================================================================================
# #Plot tau of x and its derivatives (NO REIONISATION)
# plt.plot(x_no_reionisation,tau_no_reionisation, label = "$\\tau$")
# plt.plot(x_no_reionisation, -tau_deriv_no_reionisation, label = "$-\\tau\\,'(x)$")
# plt.plot(x_no_reionisation, tau_double_deriv_no_reionisation, label = "$\\tau\\,''(x)$")
# plt.plot(x_just_saha, tau_just_saha)
# plt.xlabel("x")
# # plt.ylabel("$\\tau$")
# plt.xlim(-12, 0)
# plt.ylim(10**-8, 10**8)
# plt.yscale("log")
# plt.legend(fontsize = 16)
# plt.show()
# #======================================================================================================================


#==========================================================================================================================================================================================
# # Subplot g_tilde of x and its derivatives NORMALISED!!!
# fig, ax = plt.subplots(3,1, figsize = (12,20))
# ax[0].plot(x_no_reionisation,normalising(g_tilde_no_reionisation), label = "$\\tilde{g}(x)$", color = "#d30e92")
# ax[0].plot(x_with_reionisation,normalising(g_tilde_with_reionisation), linestyle=(0, (2, 2)), linewidth=2, color = "#f6acc5", label = "$\\tilde{g}(x)$ with reionisation")
# ax[0].plot(x_no_reionisation, normalising(g_tilde_just_saha), color = "#cbc3e3", zorder = 1)
# ax[0].axvline(x = x_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[0].axvline(x = x_no_reionisation_tau_1, linestyle = "dotted", color = "#949494", label = "Decoupling ($\\tau=1$)")
# ax[0].axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# ax[0].set_xlabel("x", fontsize = 18)
# ax[0].set_ylabel("$\\tilde{g}(x)$", fontsize = 18)
# ax[0].set_ylim(-0.001, 0.013)
# ax[0].tick_params(labelsize = 18)
# ax[0].legend(fontsize = 15)


# ax[1].plot(x_no_reionisation,g_tilde_deriv_no_reionisation, label = "$\\tilde{g}\\,'(x)$", color = "#d30e92")
# ax[1].plot(x_with_reionisation,g_tilde_deriv_with_reionisation, linestyle=(0, (2, 2)), label = "$\\tilde{g}\\,'(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# ax[1].plot(x_just_saha, g_tilde_deriv_just_saha, color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}\\,'(x)$ only Saha")
# ax[1].axvline(x = x_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[1].axvline(x = x_no_reionisation_tau_1, linestyle = "dotted", color = "#949494", label = "Decoupling ($\\tau=1$)")
# ax[1].axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# ax[1].set_xlabel("x", fontsize = 16)
# ax[1].set_ylabel("$\\tilde{g}\\,'(x)$", fontsize = 16)
# ax[1].set_ylim(-32,52)
# ax[1].tick_params(labelsize = 18)
# ax[1].legend(fontsize = 14)

# ax[2].plot(x_no_reionisation,g_tilde_double_deriv_no_reionisation, label = "$\\tilde{g}\\,''(x)$", color = "#d30e92")
# ax[2].plot(x_with_reionisation,g_tilde_double_deriv_with_reionisation, linestyle=(0, (2, 2)), label = "$\\tilde{g}\\,''(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# ax[2].plot(x_just_saha, g_tilde_double_deriv_just_saha, color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}\\,''(x)$ only Saha")
# ax[2].axvline(x = x_no_reionisation_recombination, linestyle = "dotted", color = "k", label = "Recombination ($X_e=0.1$)", zorder = 6)
# ax[2].axvline(x = x_no_reionisation_tau_1, linestyle = "dotted", color = "#949494", label = "Decoupling ($\\tau=1$)")
# ax[2].axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# ax[2].set_xlabel("x", fontsize = 16)
# ax[2].set_ylabel("$\\tilde{g}\\,''(x)$", fontsize = 16)
# ax[2].tick_params(labelsize = 18)
# ax[2].set_ylim(-1200, 1000)
# ax[2].legend(fontsize = 14)

# # plt.savefig("plots/g_tilde_and_derivatives_normalised.pdf")
# plt.show()
#========================================================================================================================================================================================

# #======================================================================================================================
# #Plot g_tilde of x smaller peak NORMALISED!!!!
# plt.figure()
# plt.plot(x_no_reionisation,normalising(g_tilde_no_reionisation), label = "$\\tilde{g}(x)$", color = "#d30e92")
# plt.plot(x_with_reionisation,normalising(g_tilde_with_reionisation), linestyle=(0, (2, 2)), label = "$\\tilde{g}(x)$ with reionisation", color = "#f6acc5", zorder = 6)
# plt.plot(x_just_saha, normalising(g_tilde_just_saha), color = "#cbc3e3", zorder = 1, label = "$\\tilde{g}(x)$ only Saha")
# plt.axvline(x = x_from_z(z_reionisation), linestyle=(0, (2, 2)), color = "#545454", label = "Reionisation ($z=7.8$)")
# plt.xlabel("x", fontsize = 12)
# plt.ylabel("$\\tilde{g}(x)$", fontsize = 12)
# plt.xlim(-3,0)
# plt.ylim(-0.0005, 0.001)
# plt.xticks(fontsize = 11)
# plt.yticks(fontsize = 11)
# plt.legend()
# #plt.savefig("plots/g_tilde_small_peak_normalised.pdf")
# plt.show()
# #======================================================================================================================