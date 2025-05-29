#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}


//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  // DONE!
   //===================================================================

  double k_start = k_min;
  double k_end = k_max;
  int n_k = 250;
  int n_x = 10000;

  Vector k_array = exp(Utils::linspace(log(k_start), log(k_end), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  Vector x_array_tight;
  Vector x_array_full;
  

  bool neutrinos = Constants.neutrinos;
  bool polarization = Constants.polarization;

// We create vectors for each quantity we want to spline
  Vector deltacdm_array(n_x * n_k);
  Vector deltab_array(n_x * n_k);
  Vector vcdm_array(n_x * n_k);
  Vector vb_array(n_x * n_k);
  Vector Phi_array(n_x * n_k);
  Vector Psi_array(n_x * n_k);
  Vector capital_pi_array(n_x * n_k);

  // Theta vectors
  Vector2D Theta_array(Constants.n_ell_theta);
  for(int i = 0; i < Constants.n_ell_theta; i++){
    Theta_array[i] = Vector(n_x*n_k);
  }

  // Nu vectors
  Vector2D Nu_array(Constants.n_ell_neutrinos);
  for(int i = 0; i < Constants.n_ell_neutrinos; i++){
    Nu_array[i] = Vector(n_x*n_k);
  }

  // Polarization vectors
  Vector2D Thetap_array(Constants.n_ell_thetap);
  for(int i = 0; i < Constants.n_ell_thetap; i++){
    Thetap_array[i] = Vector(n_x*n_k);
  }

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    int ix = 0;
    while(x_array[ix] < x_end_tight){
      ix++;
    }
    x_array_tight = Utils::linspace(x_start, x_array[ix-1], ix);
    x_array_full = Utils::linspace(x_array[ix-1], x_end, n_x-(ix-1));
    // for(int j = 0; j < (ix); j++){
    //   std::cout << "x_arr  " << j << " " << x_array_tight[j] << std::endl;
    // }
    

    //===================================================================
    // Tight coupling integration
    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    
    // The tight coupling ODE system
    ODESolver tight_coup;
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };
    
    // Integrate from x_start -> x_end_tight
    tight_coup.solve(dydx_tight_coupling, x_array_tight, y_tight_coupling_ini);
    
    //==================================================================================================
    // We create arrays for the values from the ODES
    //==================================================================================================
    auto Phi_array_tc = tight_coup.get_data_by_component(Constants.ind_Phi_tc);
    auto deltacdm_array_tc = tight_coup.get_data_by_component(Constants.ind_deltacdm_tc);
    auto deltab_array_tc = tight_coup.get_data_by_component(Constants.ind_deltab_tc);
    auto vcdm_array_tc = tight_coup.get_data_by_component(Constants.ind_vcdm_tc);
    auto vb_array_tc = tight_coup.get_data_by_component(Constants.ind_vb_tc);
  
    
    //For theta we have 2D arrays (array of arrays)
    Vector2D Theta_array_tc(Constants.n_ell_theta_tc);
    for(int i = 0; i < Constants.n_ell_theta_tc; i++){
      Theta_array_tc[i] = tight_coup.get_data_by_component(Constants.ind_start_theta_tc + i);
    }


    //For neutrinos we have 2D arrays (array of arrays)
    Vector2D Nu_array_tc(Constants.n_ell_neutrinos_tc);
    if(neutrinos){
      for(int i = 0; i < Constants.n_ell_neutrinos_tc; i++){
        Nu_array_tc[i] = tight_coup.get_data_by_component(Constants.ind_start_nu_tc + i);
      }
    }   
  
    //====================================================================
    // Full equation integration
    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    Vector y_tight_coupling(Constants.n_ell_tot_tc);
 
    
    for(int i = 0; i < Constants.n_ell_tot_tc; i++){
      auto y_array = tight_coup.get_data_by_component(i);
      y_tight_coupling[i] = y_array[ix-1];
    }
   
    
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_array_full[0], k);

    // ODESolver full_system;
    ODESolver full_system;
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    full_system.solve(dydx_full, x_array_full, y_full_ini);
    

    // We create arrays for the values from the ODES
    auto Phi_array_full = full_system.get_data_by_component(Constants.ind_Phi);
    auto deltacdm_array_full = full_system.get_data_by_component(Constants.ind_deltacdm);
    auto deltab_array_full = full_system.get_data_by_component(Constants.ind_deltab);
    auto vcdm_array_full = full_system.get_data_by_component(Constants.ind_vcdm);
    auto vb_array_full = full_system.get_data_by_component(Constants.ind_vb);
    

     //For theta we have a 2D array (array of arrays)
     Vector2D Theta_array_full(Constants.n_ell_theta);
     for(int i = 0; i < Constants.n_ell_theta; i++){
       Theta_array_full[i] = full_system.get_data_by_component(Constants.ind_start_theta + i);
     }
  

      //For neutrinos we have a 2D array (array of arrays)
     Vector2D Nu_array_full(Constants.n_ell_neutrinos);
     if(neutrinos){
       for(int i = 0; i < Constants.n_ell_neutrinos; i++){
         Nu_array_full[i] = full_system.get_data_by_component(Constants.ind_start_nu + i);
       }
     } 

     //For polarization we have a 2D array (array of arrays)
     Vector2D Thetap_array_full(Constants.n_ell_thetap);
     if(polarization){
       for(int i = 0; i < Constants.n_ell_thetap; i++){
         Thetap_array_full[i] = full_system.get_data_by_component(Constants.ind_start_thetap + i);
       }
     } 
    //=================================================================================================
    // Now we want to fill in the vectors we have created for each quantity. We do this by looping
    // over the arrays we created with values from the ODES  
    //=================================================================================================
    
    double Thetap_0;
    double Thetap_2;
    double Theta_2;
    double H_0                  =cosmo->H_of_x(0);
    double OmegaR_0             =cosmo->get_OmegaR(0);
    double OmegaNu_0            =cosmo->get_OmegaNu(0);

    
    for (int l = 0; l < n_x; l++){
      // std::cout << l << std::endl;
      double Hp_of_x              = cosmo->Hp_of_x(x_array[l]);
      double tau_derivative       = rec->dtaudx_of_x(x_array[l]);
      double a                    = exp(x_array[l]);
      
      if(l < ix-1){
        // Arrays for scalars
        Phi_array[l+(ik*n_x)] = Phi_array_tc[l]; 
        // std::cout << k * Constants.Mpc << " " << l << " " <<  Phi_array_tc[l] << std::endl; 
        deltacdm_array[l+(ik*n_x)] = deltacdm_array_tc[l];   
        deltab_array[l+(ik*n_x)] = deltab_array_tc[l];
        vcdm_array[l+(ik*n_x)] = vcdm_array_tc[l];
        vb_array[l+(ik*n_x)] = vb_array_tc[l];
        

        // Theta arrays
        for(int j = 0; j < Constants.n_ell_theta; j++){
          if(j<Constants.n_ell_theta_tc){
            Theta_array[j][l+(ik*n_x)] = Theta_array_tc[j][l];
          } else{
            Theta_array[j][l+(ik*n_x)] = (-j / (2.0*j + 1.0)) * ((Constants.c*k)/(Hp_of_x*tau_derivative))*Theta_array[j-1][l];
          }
        }
        
        // Nu arrays
        if(neutrinos){
          for(int j = 0; j < Constants.n_ell_neutrinos; j++){
            if(j>Constants.n_ell_neutrinos_tc-1){
              Nu_array[j][l+(ik*n_x)] = ((Constants.c * k) / ((2*j + 1)*Hp_of_x)) * Nu_array[j-1][l];
            } else{
              Nu_array[j][l+(ik*n_x)] = Nu_array_tc[j][l]; 
            }
          }
        }

        if(polarization){
          for(int j = 0; j < Constants.n_ell_thetap; j++){
            Thetap_array[j][l+(ik*n_x)] = 0; 
          }
        } 
        

        if(polarization){
          Theta_2 = (Theta_array_tc[1][l]) * (-((8. * Constants.c * k) / (15. * Hp_of_x * tau_derivative)));
          Thetap_0 = (5./4.) * Theta_2; 
          Thetap_2 = (1./4.) * Theta_2; 
        } else {
          Theta_2 = (Theta_array_tc[1][l]) * (-((20. * Constants.c * k) / (45. * Hp_of_x * tau_derivative)));
          Thetap_0 = 0;
          Thetap_2 = 0;
        }
        

        if(neutrinos){ 
          Psi_array[l+(ik*n_x)] = (- Phi_array_tc[l]) - ((12 * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a, 2)))) * ((OmegaR_0*Theta_2)+(Nu_array_tc[2][l]*OmegaNu_0));
        }else{
          Psi_array[l+(ik*n_x)] = (- Phi_array_tc[l]) - ((12 * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a, 2)))) * (OmegaR_0*Theta_2);
        }
        
        if(polarization){ //problem is here
          capital_pi_array[l+(ik*n_x)] = Theta_2 + Thetap_0 + Thetap_2;
        } else{
          capital_pi_array[l+(ik*n_x)] = Theta_2;
        }
        // std::cout << "fuck you" << std::endl;

      } else{
  
        // // Arrays for scalars
        Phi_array[l+(ik*n_x)] = Phi_array_full[l-(ix-1)]; // problemet er her  
        deltacdm_array[l+(ik*n_x)] = deltacdm_array_full[l-(ix-1)];    
        deltab_array[l+(ik*n_x)] = deltab_array_full[l-(ix-1)];
        vcdm_array[l+(ik*n_x)] = vcdm_array_full[l-(ix-1)];
        vb_array[l+(ik*n_x)] = vb_array_full[l-(ix-1)];
  
        // Theta arrays
        for(int j = 0; j < Constants.n_ell_theta; j++){
          Theta_array[j][l+(ik*n_x)] = Theta_array_full[j][l-(ix-1)]; 
        }
        
        // Nu arrays
        if(neutrinos){
          for(int j = 0; j < Constants.n_ell_neutrinos; j++){
            Nu_array[j][l+(ik*n_x)] = Nu_array_full[j][l-(ix-1)];
          }
        }

        // Thetap arrays
        if(polarization){
          for(int j = 0; j < Constants.n_ell_thetap; j++){
            Thetap_array[j][l+(ik*n_x)] = Thetap_array_full[j][l-(ix-1)];
          }
        }
  
        if(neutrinos){ 
          Psi_array[l+(ik*n_x)] = (- Phi_array_full[l-(ix-1)]) - ((12 * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a, 2)))) * ((OmegaR_0*Theta_array_full[2][l-(ix-1)])+(Nu_array_full[2][l-(ix-1)]*OmegaNu_0));
        }else{
          Psi_array[l+(ik*n_x)] = (- Phi_array_full[l-(ix-1)]) - ((12 * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a, 2)))) * (OmegaR_0*Theta_array_full[2][l-(ix-1)]);
        } 

        if(polarization){
          capital_pi_array[l+(ik*n_x)] = Theta_array_full[2][l-(ix-1)] + Thetap_array_full[0][l-(ix-1)] + Thetap_array_full[2][l-(ix-1)];
        } else{
          capital_pi_array[l+(ik*n_x)] = Theta_array_full[2][l-(ix-1)];
        }

      }
    }
  } 
 
  //=============================================================================
  // We now make all the splines we need
  Phi_spline.create(x_array, k_array, Phi_array);
  Psi_spline.create(x_array, k_array, Psi_array);
  delta_b_spline.create(x_array, k_array, deltab_array);
  delta_cdm_spline.create(x_array, k_array, deltacdm_array);
  v_b_spline.create(x_array, k_array, vb_array);
  v_cdm_spline.create(x_array, k_array, vcdm_array);
  Pi_spline.create(x_array, k_array, capital_pi_array);

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for(int j = 0; j < Constants.n_ell_theta; j++){
    Theta_spline[j].create(x_array, k_array, Theta_array[j]);
  }

  Nu_spline = std::vector<Spline2D>(Constants.n_ell_neutrinos);
  for(int j = 0; j < Constants.n_ell_neutrinos; j++){
    Nu_spline[j].create(x_array, k_array, Nu_array[j]);
  }

  Theta_p_spline = std::vector<Spline2D>(Constants.n_ell_thetap);
  for(int j = 0; j < Constants.n_ell_thetap; j++){
    Theta_p_spline[j].create(x_array, k_array, Thetap_array[j]);
  }
 
  Utils::EndTiming("integrateperturbation");
    
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;
  
  // References to the tight coupling quantities
  double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // Set the initial conditions in the tight coupling regime
  const double a = exp(x);
  const double Hp_of_x            = cosmo->Hp_of_x(x);
  const double H_0                = cosmo->H_of_x(0);
  const double dtaudx_of_x        = rec->dtaudx_of_x(x);
  const double OmegaNu_0          = cosmo->get_OmegaNu(0);
  const double OmegaR_0           = cosmo->get_OmegaR(0);
  const double f_nu               = OmegaNu_0/(OmegaNu_0 + OmegaR_0);


  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  double Psi = -1. / ((3./2.) + ((2.* f_nu)/5.));
  Phi_tc = -(1.+ ((2. *f_nu)/5.)) * Psi;
  delta_cdm_tc = -(3./2.) * Psi;
  delta_b_tc = delta_cdm_tc;
  v_cdm_tc = -((Constants.c * k)/(2. * Hp_of_x))*Psi; //fixed   
  v_b_tc = v_cdm_tc;


  // SET: Photon temperature perturbations (Theta_ell)
  *Theta_tc = -((1./2.)*Psi);
  *(Theta_tc+1) = ((Constants.c * k) / (6. * Hp_of_x)) * Psi;
  //Theta ell?

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    *Nu_tc = -((1./2.)*Psi);
    *(Nu_tc+1) = ((Constants.c * k) / (6. * Hp_of_x)) * Psi;
    *(Nu_tc+2) = -(pow(Constants.c, 2) * pow(k, 2) * pow(a, 2) * (Psi + Phi_tc))
                  /(12 * pow(H_0, 2) * OmegaNu_0); //fixed
    
    for(int ell = 3; ell < n_ell_neutrinos_tc; ell++){
      *(Nu_tc+ell) = *(Nu_tc+ell-1)*((Constants.c * k)/(((2*ell)+1)*Hp_of_x)); 
    }

  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  double *Theta_p         = &y[Constants.ind_start_thetap];
  double *Nu              = &y[Constants.ind_start_nu];

  //Quantities we must fetch from previous milestones
  const double Hp_of_x            = cosmo->Hp_of_x(x);
  const double H_0                = cosmo->H_of_x(0);
  const double OmegaNu_0          = cosmo->get_OmegaNu(0);
  const double OmegaR_0           = cosmo->get_OmegaR(0);
  const double dtaudx_of_x        = rec->dtaudx_of_x(x); 


  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm= v_cdm_tc;   
  v_b = v_b_tc;

 // SET: Photon temperature perturbations (Theta_ell)
 *Theta = *Theta_tc;
 *(Theta+1) = *(Theta_tc+1);
 if(polarization){
   *(Theta+2) = -((8. * Constants.c * k) / (15. * Hp_of_x * dtaudx_of_x)) * *(Theta+1); 
 } else {
   *(Theta+2) = -((20. * Constants.c * k) / (45. * Hp_of_x * dtaudx_of_x)) * *(Theta+1);
 }

 for(int ell = 3; ell < n_ell_theta; ell++){
  *(Theta+ell) = *(Theta+ell-1)*(-ell)/((2*ell) + 1) * (Constants.c * k)/(Hp_of_x * dtaudx_of_x);
 }

 // SET: Neutrino perturbations (N_ell)
 if(neutrinos){
  *Nu = *Nu_tc;
  *(Nu+1) = *(Nu_tc+1);
  *(Nu+2) = *(Nu_tc+2);

  for(int ell = 3; ell < n_ell_neutrinos; ell++){
    *(Nu+ell) = *(Nu_tc+ell);
   }
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    *Theta_p = *(Theta+2) * (5./4.);
    *(Theta_p+1) = *(Theta+2) * -((Constants.c * k) / (4 * Hp_of_x * dtaudx_of_x));
    *(Theta_p+2) = *(Theta+2) * (1./4.);
    
    for(int ell = 3; ell < n_ell_thetap; ell++){
      *(Theta_p+ell) = -(*(Theta_p+ell-1))*((ell)/((2.0*ell) + 1.0)) * ((Constants.c * k)/(Hp_of_x * dtaudx_of_x));
     }
  }

  return y;
}

//====================================================
// Calculate the time when tight coupling end
double Perturbations::get_tight_coupling_time(const double k) const{
  double max_x_value_for_tight_coupling = -8.3; // if the condidion never kicks inn, our maximum allowed x-value will be -8.3

  Vector x_arr = Utils::linspace(x_start, -8.3, 1000); //this is a serach array for where tight coupling ends

  for (auto x: x_arr){
    double tau_derivative = std::abs(rec->dtaudx_of_x(x));
    if (tau_derivative < 10 or tau_derivative < 10 * (Constants.c * k)/(cosmo->Hp_of_x(x))){
      return x;
    }
  }
  return max_x_value_for_tight_coupling;
}

//====================================================
// After integrating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  // for (int l = 0; l < k_array.size(); l++){
  //   std::cout << k_array[l] << std::endl;
  // }
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ik = 0; ik < k_array.size(); ik++){
    const double k = k_array[ik];
    for(auto ix = 0; ix < x_array.size(); ix++){
      const double x = x_array[ix];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + (n_x * ik);

      // Fetch all the things we need
      bool neutrinos               = Constants.neutrinos;
      bool polarization            = Constants.polarization;
   
      double a                                  = exp(x);
      const double Hp                           = cosmo->Hp_of_x(x);
      const double Hp_deriv                     = cosmo->dHpdx_of_x(x);
      const double Hp_double_deriv              = cosmo->ddHpddx_of_x(x);
      
      const double tau                          = rec->tau_of_x(x);
      const double visibility                   = rec->g_tilde_of_x(x);
      const double visibility_deriv             = rec->dgdx_tilde_of_x(x);
      const double visibility_double_deriv      = rec->ddgddx_tilde_of_x(x);
      const double theta_0                      = get_Theta(x, k, 0);
      const double psi                          = get_Psi(x,k);
      const double PI                           = get_Pi(x, k);
      const double PI_deriv                     = get_Pi_derivative_x(x, k);
      const double PI_double_deriv              = get_Pi_double_derivative_x(x, k);
      const double phi                          = get_Phi(x, k);
      const double psi_deriv                    = get_Psi_derivative_x(x, k);
      const double phi_deriv                    = get_Phi_derivative_x(x, k);
      const double v_b                          = get_v_b(x, k);
      const double v_b_deriv                    = get_v_b_derivative_x(x, k);
      const double eta_0                        = cosmo->eta_of_x(0); 
      const double eta                          = cosmo->eta_of_x(x);
     
      double first_term = visibility*(theta_0 + psi + (PI/4.));
      double second_term = exp(-tau) * (psi_deriv - phi_deriv); 
      double third_term = -1./(Constants.c * k) * ((Hp_deriv*visibility*v_b) + (Hp * visibility_deriv * v_b) + (Hp * visibility * v_b_deriv));
      double fourth_term = (3./(4*pow(Constants.c, 2)*pow(k,2)))*((Hp_deriv * Hp_deriv * visibility * PI) + (Hp * Hp_double_deriv * visibility * PI)
                            + (Hp * Hp_deriv * visibility_deriv * PI) + (Hp * Hp_deriv * visibility * PI_deriv) + (Hp_deriv * Hp *visibility_deriv * PI)
                            + (Hp * Hp_deriv * visibility_deriv * PI) + (Hp * Hp * visibility_double_deriv * PI) + (Hp * Hp * visibility_deriv * PI_deriv)
                            + (Hp_deriv * Hp * visibility * PI_deriv) + (Hp * Hp_deriv * visibility * PI_deriv) + (Hp * Hp * visibility_deriv * PI_deriv)
                            + (Hp * Hp * visibility * PI_double_deriv));

      // Temperature source
      ST_array[index] = first_term + second_term + third_term + fourth_term; //looks like this works correctly!
      // ST_array[index] = visibility;

      // Polarization source
      if(Constants.polarization){
        if(x > -3.){
          SE_array[index] = SE_array[index-1];
        }else{
          SE_array[index] = (3. * visibility * PI) / (4 * pow(k, 2) * pow(eta_0 - eta, 2));
        }
        
      }
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;
  const bool polarization       = Constants.polarization;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  double a = exp(x);
  double Psi;
  double Theta_2;

  // Fetch quantities from previous milestones
  const double Hp_of_x                  = cosmo->Hp_of_x(x);
  const double H_prime_derivative       = cosmo->dHpdx_of_x(x);
  const double H_0                      = cosmo->get_H0();
  const double OmegaCDM_0               = cosmo->get_OmegaCDM(0);
  const double OmegaB_0                 = cosmo->get_OmegaB(0);
  const double OmegaR_0                 = cosmo->get_OmegaR(0);
  const double OmegaNu_0                = cosmo->get_OmegaNu(0);
  const double conformal_time           = cosmo->eta_of_x(x);
  const double tau_derivative           = rec->dtaudx_of_x(x);
  const double tau_double_derivative    = rec->ddtauddx_of_x(x);

  
  if(polarization){
    Theta_2 = *(Theta+1) * (-((8. * Constants.c * k) / (15. * Hp_of_x * tau_derivative))); 
  } else {
    Theta_2 = *(Theta+1) * (-((20. * Constants.c * k) / (45. * Hp_of_x * tau_derivative)));
  }

  
  if(neutrinos){ 
    Psi = (- Phi) - (((12. * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a, 2)))) * ((OmegaR_0*Theta_2)+(*(Nu+2)*OmegaNu_0)));
  }else{
    Psi = (- Phi) - (((12. * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a, 2)))) * (OmegaR_0*Theta_2));
  }

  // SET: Scalar quantities
  if(neutrinos){
    dPhidx = Psi - (((pow(Constants.c, 2)*pow(k, 2)) / (3*pow(Hp_of_x, 2)))*Phi)
    + ((pow(H_0, 2) / (2.*pow(Hp_of_x, 2))) * ((OmegaCDM_0*pow(a, -1)*delta_cdm)
    +(OmegaB_0 * pow(a, -1) * delta_b)+((*Theta)*4.*OmegaR_0*pow(a, -2))
    +((*Nu)*4.*OmegaNu_0*pow(a, -2))));
  }else{
    dPhidx = Psi - (((pow(Constants.c, 2)*pow(k, 2)) / (3*pow(Hp_of_x, 2)))*Phi)
    + ((pow(H_0, 2) / (2.*pow(Hp_of_x, 2))) * ((OmegaCDM_0*pow(a, -1)*delta_cdm)
    +(OmegaB_0 * pow(a, -1) * delta_b)+((*Theta)*4.*OmegaR_0*pow(a, -2))));
  }

  *dThetadx = (*(Theta+1) * (-(Constants.c * k)/(Hp_of_x)))-dPhidx;

  //Needed for tight coupling regime
  double R = (4.*OmegaR_0)/(3.*OmegaB_0*a);
  double q = ((-(((1.-R)*tau_derivative) + ((1+R)*tau_double_derivative))*((*(Theta+1)*3)+(v_b)))
              -(((Constants.c * k)/(Hp_of_x))*Psi)
              + ((1-((H_prime_derivative)/(Hp_of_x)))*((Constants.c * k)/(Hp_of_x))*(-*Theta 
              + ((Theta_2)*2)))-(*(dThetadx)*((Constants.c * k)/(Hp_of_x)))) 
              / (((1.+R)*tau_derivative) + ((H_prime_derivative)/(Hp_of_x)) - 1.);
  
  ddelta_cdmdx = (((Constants.c * k)/(Hp_of_x))*(v_cdm)) - (3.*dPhidx);
  
  ddelta_bdx = (((Constants.c * k)/(Hp_of_x))*(v_b)) - (3.*dPhidx);
  
  dv_cdmdx = (-v_cdm) - (((Constants.c * k)/(Hp_of_x))*Psi);
  
  dv_bdx = (1. / (1.+R)) * (-v_b - (((Constants.c * k)/(Hp_of_x))*Psi) 
          + (R*(q+(((Constants.c * k)/Hp_of_x)*(-*Theta+(2.*Theta_2))) 
          - (((Constants.c * k)/Hp_of_x)*Psi))));


  // SET: Photon multipoles 
  *(dThetadx+1) = (1./3.)*(q-dv_bdx);
  
  // SET: Neutrino mutlipoles
  if(neutrinos){
    *dNudx = (*(Nu+1)*(-((Constants.c * k)/(Hp_of_x)))) - dPhidx;
    *(dNudx+1) = (*Nu * ((Constants.c * k)/(3*Hp_of_x))) - (*(Nu+2)*((2 * Constants.c * k)/(3*Hp_of_x)))
    + (((Constants.c * k)/(3* Hp_of_x))*Psi);
    
    for(int ell = 2; ell < (n_ell_neutrinos_tc-1); ell++){
      *(dNudx+ell) = (*(Nu+ell-1)*((ell * Constants.c * k)/(((2*ell)+1)*Hp_of_x)))
      - (*(Nu+ell+1)*(((ell+1)*Constants.c * k)/(((2*ell)+1)*Hp_of_x)));
    }

    *(dNudx+(n_ell_neutrinos_tc-1)) = *(Nu+(n_ell_neutrinos_tc-1)-1)*((Constants.c * k)/(Hp_of_x))
        - *(Nu+(n_ell_neutrinos_tc-1))*(Constants.c)*(((n_ell_neutrinos_tc-1) + 1)/(Hp_of_x * conformal_time));
  }
  
  return GSL_SUCCESS;
}


// ====================================================
// The right hand side of the full ODE
int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  double capital_pi;
  double Psi;
  double a = exp(x);

  // Quantities we must fetch from previous milestones
  double Hp_of_x                      = cosmo->Hp_of_x(x);
  double H_0                          = cosmo->H_of_x(0);
  double OmegaB_0                     = cosmo->get_OmegaB(0);
  double OmegaR_0                     = cosmo->get_OmegaR(0);
  double OmegaCDM_0                   = cosmo->get_OmegaCDM(0);
  double OmegaNu_0                    = cosmo->get_OmegaNu(0);
  const double conformal_time         = cosmo->eta_of_x(x);

  // Recombination variables
  double tau_derivative               = rec->dtaudx_of_x(x);

  //=============================================================================
  // Fill in the expressions for all the derivatives
  if(neutrinos){ 
    Psi = (- Phi) - ((12 * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a,2)))) * ((*(Theta+2)*OmegaR_0)+(*(Nu+2)*OmegaNu_0));
  } else{
    Psi = (- Phi) - ((12 * pow(H_0,2))/(pow(Constants.c, 2) * pow(k, 2) * (pow(a,2)))) * (*(Theta+2)*OmegaR_0);
  }

  // SET: Scalar quantities
  if(neutrinos){
    dPhidx = Psi - (((pow(Constants.c, 2)*pow(k, 2)) / (3*pow(Hp_of_x, 2)))*Phi)
    + ((pow(H_0, 2) / (2.*pow(Hp_of_x, 2))) * ((OmegaCDM_0*pow(a, -1)*delta_cdm)
    +(OmegaB_0 * pow(a, -1) * delta_b)+((*Theta)*4*OmegaR_0*pow(a, -2))
    +((*Nu)*4.*OmegaNu_0*pow(a, -2))));
  }else{
    dPhidx = Psi - (((pow(Constants.c, 2)*pow(k, 2)) / (3*pow(Hp_of_x, 2)))*Phi)
    + ((pow(H_0, 2) / (2*pow(Hp_of_x, 2))) * ((OmegaCDM_0*pow(a, -1)*delta_cdm)
    +(OmegaB_0 * pow(a, -1) * delta_b)+((*Theta)*4*OmegaR_0*pow(a, -2))));
  }

  ddelta_cdmdx = (((Constants.c * k)/(Hp_of_x))*(v_cdm))-(3.*dPhidx);
  
  ddelta_bdx = (((Constants.c * k)/(Hp_of_x))*(v_b)) - (3.*dPhidx);
  
  dv_cdmdx = (-v_cdm) - (((Constants.c) * k/(Hp_of_x))*Psi);

  double R = (4.*OmegaR_0)/(3.*OmegaB_0*a);
  dv_bdx = (-v_b) - (((Constants.c * k)/(Hp_of_x))*Psi) + (tau_derivative*(R*((*(Theta+1)*3)+v_b)));

  // SET: Photon polarization multipoles
  if(polarization){
    capital_pi = *(Theta+2) + *(Theta_p) + *(Theta_p+2);
    *(dTheta_pdx) = -(*(Theta_p+1)*((Constants.c * k)/(Hp_of_x))) + (tau_derivative*(*Theta_p - (capital_pi/2)));
    for(int ell = 1; ell < (n_ell_thetap-1); ell++){
      if(ell==2){
        *(dTheta_pdx+ell) = (*(Theta_p+ell-1)*((ell * Constants.c * k)/(((2*ell)+1.0)*Hp_of_x)))
        -(*(Theta_p+ell+1)*(((ell+1)*Constants.c * k)/(((2*ell)+1.0)*Hp_of_x)))
        + (tau_derivative*(*(Theta_p+ell)-((1./10.)*capital_pi))); 
      }else{
        *(dTheta_pdx+ell) = (*(Theta_p+ell-1)*((ell * Constants.c * k)/(((2*ell)+1.0)*Hp_of_x)))
        -(*(Theta_p+ell+1)*(((ell+1)*Constants.c * k)/(((2*ell)+1.0)*Hp_of_x)))
        + (tau_derivative*(*(Theta_p+ell)));
      }
    }
   
    *(dTheta_pdx + (n_ell_thetap-1)) = (*(Theta_p+(n_ell_thetap-1)-1)*((Constants.c * k)/(Hp_of_x)))
                                    - (*(Theta_p+(n_ell_thetap-1))*(Constants.c*(((n_ell_thetap-1)+1)/(Hp_of_x*conformal_time))))
                                    + (*(Theta_p+(n_ell_thetap-1))*tau_derivative);
 
  } else{
    capital_pi = *(Theta+2); 
  }

  // SET: Photon multipoles 
  *dThetadx = (*(Theta+1)*(-((Constants.c * k) / (Hp_of_x)))) - dPhidx;
  *(dThetadx+1) = *Theta*((Constants.c * k)/(3*Hp_of_x)) - (*(Theta+2)*((2 * Constants.c * k)/(3 * Hp_of_x)))
                + (((Constants.c * k)/(3 * Hp_of_x))*Psi) + (tau_derivative*(*(Theta+1) + (v_b/3.))); 
  
  for(int ell = 2; ell < (n_ell_theta-1); ell++){
    if(ell==2){
      *(dThetadx+ell) = (*(Theta+ell-1)*((ell * Constants.c * k)/(((2*ell)+1)*Hp_of_x)))
      -(*(Theta+ell+1)*(((ell+1)*Constants.c * k)/(((2*ell)+1)*Hp_of_x)))
      + (tau_derivative*(*(Theta+ell)-((1./10.)*capital_pi)));
    } else{
      *(dThetadx+ell) = (*(Theta+ell-1)*((ell * Constants.c * k)/(((2*ell)+1)*Hp_of_x)))
      -(*(Theta+ell+1)*(((ell+1)*Constants.c * k)/(((2*ell)+1)*Hp_of_x)))
      + (tau_derivative*(*(Theta+ell)));
    }
  }
  
  *(dThetadx+(n_ell_theta-1)) = (*(Theta+(n_ell_theta-1)-1)*((Constants.c * k)/(Hp_of_x)))
                          - (*(Theta+(n_ell_theta-1))*(Constants.c*(((n_ell_theta-1)+1)/(Hp_of_x*conformal_time))))
                          + (*(Theta+(n_ell_theta-1))*tau_derivative);

  
  // SET: Neutrino mutlipoles
  if(neutrinos){
    *dNudx = (*(Nu+1)*(-((Constants.c * k)/(Hp_of_x)))) - dPhidx;
    *(dNudx+1) = (*Nu * ((Constants.c * k)/(3*Hp_of_x))) - (*(Nu+2)*((2 * Constants.c * k)/(3*Hp_of_x)))
    + (((Constants.c * k)/(3* Hp_of_x))*Psi);
    
    for(int ell = 2; ell < (n_ell_neutrinos-1); ell++){
      *(dNudx+ell) = (*(Nu+ell-1)*((ell * Constants.c * k)/(((2*ell)+1)*Hp_of_x)))
      - (*(Nu+ell+1)*(((ell+1)*Constants.c * k)/(((2*ell)+1)*Hp_of_x)));
    }

    *(dNudx+(n_ell_neutrinos-1)) = *(Nu+(n_ell_neutrinos-1)-1)*((Constants.c * k)/(Hp_of_x))
        - *(Nu+(n_ell_neutrinos-1))*(Constants.c)*(((n_ell_neutrinos-1) + 1)/(Hp_of_x * conformal_time));
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_v_b_derivative_x(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Phi_derivative_x(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Psi_derivative_x(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Pi_derivative_x(const double x, const double k) const{
  return Pi_spline.deriv_x(x,k);
}
double Perturbations::get_Pi_double_derivative_x(const double x, const double k) const{
  return Pi_spline.deriv_xx(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc                  << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc                    << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc                      << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc                        << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc                       << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc               << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc                   << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc                << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc             << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc                     << "\n";
  std::cout << "eta(0) is:          " << cosmo->eta_of_x(0) / Constants.Mpc                         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  bool neutrinos = Constants.neutrinos;
  bool polarization = Constants.polarization;
  const int npts = 10000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_delta_cdm(x,k)       << " ";
    fp << get_delta_b(x,k)       << " ";
    fp << get_v_cdm(x,k)       << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    if(neutrinos){
      fp << get_Nu(x,k, 0)        << " ";
      fp << get_Nu(x,k, 1)        << " ";
      fp << get_Nu(x,k, 2)        << " ";
    }
    if(polarization){
      fp << get_Theta_p(x,k, 0)        << " ";
      fp << get_Theta_p(x,k, 1)        << " ";
      fp << get_Theta_p(x,k, 2)        << " ";
    }
    fp << cosmo->eta_of_x(x)                                   << " ";
    fp << get_Source_T(x,k)                                    << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";

    if(polarization){
      fp << get_Source_E(x,k)                                    << " "; 
      fp << get_Source_E(x,k) * Utils::j_ell(5,   arg)           << " ";
      fp << get_Source_E(x,k) * Utils::j_ell(50,  arg)           << " ";
      fp << get_Source_E(x,k) * Utils::j_ell(500, arg)           << " ";
      fp << "\n";
    }
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
