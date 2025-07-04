#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters for initial background cosmology 
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // // Background parameters from comparison plots (Hans Winter) 
  // double h           = 0.7;
  // double OmegaB      = 0.05;
  // double OmegaCDM    = 0.45;
  // double OmegaK      = 0.0;
  // double Neff        = 0.0;
  // double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0;
  
  // Perturbation parameters
  // double f_nu        = 0;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("cosmology.txt");

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

  // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results.txt");

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  /*Output recombination quantities*/ 
  // rec.output("recombination_with_reionisation.txt");
  rec.output("recombination_no_reionisation.txt"); // change this if you are solving with reionisation, just saha etc.
  // rec.output("recombination_just_saha.txt");


  //=========================================================================
  // Module III
  //=========================================================================
  

  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();


  
  // Output perturbation quantities
  double kvalue = 0.1 / Constants.Mpc; //change this depending on what mode (k) you want to solve for
  pert.output(kvalue, "perturbations_k0.1_student.txt"); //change this depending on what mode (k) you want to solve for
  

  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("C_ells_all_sf_terms_student.txt");
  power.output_transfer("transfer_function.txt", 7, 25, 35, 45, 55); //change this if you want other \ell-values
  // power.output_transfer("transfer_function.txt", 4, 19, 24, 32, 42);
  power.output_pofk("matter_power_spectrum.txt", 0);

  Utils::EndTiming("Everything");

}

