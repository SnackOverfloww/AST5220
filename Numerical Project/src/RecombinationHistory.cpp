#include"RecombinationHistory.h"
#include <algorithm>

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
  
  // Compute and spline Xe, ne
  solve_number_density_electrons();

  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();

  calculate_sound_horizon();
}

//=====================================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//=====================================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // Sets up x-array and make arrays to store X_e(x) and n_e(x) on DONE
  
  const double xmin = x_start;
  const double xmax = x_end;
  const int    npts = 1000;
  double z_reion = 8.;
  double delta_z_reion = 0.5;
  
  Vector x_array = Utils::linspace(xmin, xmax, npts);
  Vector Xe_arr(npts);
  Vector Xe_reion_arr(npts);
  Vector ne_arr(npts);
  Vector log_ne_arr(npts);
  

  // Calculate recombination history
  bool just_saha = false;
  bool saha_regime = true;
  // std::cout << Xe_saha_limit << std::endl;
  bool reionisation = true; //change this based on whether you want reionisation or not
  for(int i = 0; i < npts ; i++){

    //==================================================================
    // Get X_e from solving the Saha equation:
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;
    
    if(just_saha){
      //Store the result we got from the Saha equation
      Xe_arr[i] = Xe_current;
      log_ne_arr[i] = log(ne_current);

    } else{
      if(saha_regime){
      
        //=============================================================================
        //Store the result we got from the Saha equation
        Xe_arr[i] = Xe_current;
        log_ne_arr[i] = log(ne_current);
        
  
        //=============================================================================
  
      //==================================================================
  
      } else {
  
        if(i==0) std::cout << "Error: i = 0" << std::endl; // Tells us if we skip the whole Saha regime
  
        // The Peebles ODE equation
        ODESolver peebles_Xe_ode;
        ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
          return rhs_peebles_ode(x, Xe, dXedx);
        };
        
        Vector x_array_peebles{x_array[i-1], x_array[i]};
        Vector initial_condition{Xe_arr[i-1]};
        peebles_Xe_ode.solve(dXedx, x_array_peebles, initial_condition);
        auto Xe_array = peebles_Xe_ode.get_data_by_component(0);
  
        Xe_arr[i] = Xe_array[1];
        
      }

    }

  }

  for (int i = 0; i < 1000; i++) {
    double y = exp((-3.*x_array[i])/2.);
    double fHe = Yp/(4.*(1-Yp));
    double y_reion = pow((1.+z_reion), (3./2.));
    double delta_y_reion = (3./2.) * sqrt(1.+z_reion) * delta_z_reion; 
        
    //=======================================================================================================
    if (reionisation){
      Xe_arr[i] = Xe_arr[i] + (((1.+fHe)/2.) * (1. + tanh((y_reion - y)/(delta_y_reion)))); 
    } else{ 
      Xe_arr[i] = Xe_arr[i]; 
    }
     
    double nH = RecombinationHistory::nH_of_x(x_array[i]);
    double ne = nH * Xe_arr[i];
    // std::cout << "Xe = " << Xe_arr[i] << ", n_e = " << ne << std::endl;
    log_ne_arr[i] = log(ne);
    // std::cout << "Xe = " << Xe_arr[i] << ", n_e = " << ne_arr[i] << std::endl;
  }

  // Make the splines
  Xe_of_x_spline.create(x_array, Xe_arr, "Xe");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");


  Utils::EndTiming("Xe");
}

//====================================================
// Solves the Saha equation to get ne and Xe

std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB(0);
  const double H0          = cosmo->get_H0();
  const double TCMB_0      = 2.7255;

  // Compute critical density
  double critical_density = (3.*pow(H0, 2.))/(8.*M_PI*Constants.G); 

  
  // Compute baryon density n_b
  double baryon_density = (OmegaB * critical_density)/(Constants.m_H * pow(a, 3.));

  //Compute baryon temprature
  double baryon_temperature = TCMB_0/a;
    
  // Compute constant for solving quadratc for X_e
  double b_for_quadratic = (1./(baryon_density * pow(Constants.hbar,3.)))*pow(((Constants.m_e*baryon_temperature*Constants.k_b)/(2.*M_PI)), (3./2.))*exp(-Constants.epsilon_0 / (baryon_temperature * Constants.k_b));
  
  double solving_quadratic = (-b_for_quadratic + sqrt(pow(b_for_quadratic, 2.) + (4*b_for_quadratic)))/2.;
  

  //=============================================================================
  // Compute Xe and ne from the Saha equation
  double Xe;
  double ne;
  if (4/b_for_quadratic > 0.001){
    Xe = solving_quadratic;
    if (Xe < pow(10,-30)){
      Xe = 1e-30;
    }
  } 
  else{
    Xe = 1;
  }
  ne = Xe * baryon_density;
  
  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double fine_structure_constant = 1/137.0359992;

  // Cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB(0);
  const double H0          = cosmo->get_H0();
  const double H_of_x      = cosmo->H_of_x(x);
  const double TCMB_0      = 2.7255; // TODO


  //=============================================================================
  //Expression for dXedx

  // double H_of_x = H0*sqrt(OmegaLambda + (OmegaK*pow(a, -2)) + ((OmegaB + OmegaCDM)*pow(a, -3)) + ((OmegaR+OmegaNu)*pow(a, -4)));
  double baryon_temperature = (TCMB_0 / a) * Constants.K;
  double phi2 = 0.448 * log(Constants.epsilon_0/(baryon_temperature * Constants.k_b));
  double alpha2 = ((64.*M_PI)/sqrt(27.*M_PI)) * (pow(fine_structure_constant, 2)/pow(Constants.m_e, 2)) * sqrt(Constants.epsilon_0 / (baryon_temperature*Constants.k_b))
         * (pow(Constants.hbar, 2) / Constants.c) * phi2;
  double beta = alpha2 * pow((Constants.m_e * baryon_temperature)/(2*M_PI), (3./2.)) * exp(-Constants.epsilon_0 / (baryon_temperature * Constants.k_b))
          * pow(Constants.hbar, -3) * pow(Constants.k_b, (3./2.));
  
  double beta2 = alpha2 * pow((Constants.m_e * baryon_temperature)/(2*M_PI), (3./2.)) * pow(Constants.hbar, -3) * pow(Constants.k_b, (3./2.)) * exp((-1.*Constants.epsilon_0)/(4*baryon_temperature*Constants.k_b)); 
  
  // double beta2 = alpha2 * pow((Constants.m_e * baryon_temperature)/(2*M_PI), (3./2.)) * exp((-1*Constants.epsilon_0)/(4*baryon_temperature*Constants.k_b));
  double n_H = (1.- Yp) * ((3*pow(H0, 2)*OmegaB)/(8.*M_PI*Constants.G*Constants.m_H*pow(a,3)));
  double n_1s = (1.- X_e)*n_H;
  double lambda_alpha = H_of_x * ((pow((3*Constants.epsilon_0), 3)) / (pow((8.*M_PI), 2) * n_1s)) * pow(Constants.c, -3) 
          * pow(Constants.hbar, -3); 
  double lambda_2s_1s = 8.227;
  double C_r = (lambda_2s_1s + lambda_alpha)/(lambda_2s_1s + lambda_alpha + beta2);
    
  
  dXedx[0] = ((C_r)/(H_of_x)) * ((beta*(1-X_e)) - (n_H * alpha2 * pow(X_e, 2)));

  return GSL_SUCCESS;
}
  //=============================================================================

//====================================================
// Solves for the optical depth tau, compute the visibility function and spline the result

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  int xmin = x_start;
  int xmax = 0;
  Vector x_array = Utils::linspace(xmin, xmax, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

  //=============================================================================
  // Expression for dtaudx
  
  const double a           = exp(x);

  // Physical constants in SI units
  const double c           = Constants.c;
  const double sigma_T     = Constants.sigma_T;
  
  // Cosmological parameters
  const double H_of_x      = cosmo->H_of_x(x);

  // Set the derivative for photon optical depth
    
    double ne = ne_of_x(x);
    dtaudx[0] = ( -Constants.c *ne* Constants.sigma_T) / H_of_x;

    return GSL_SUCCESS;
  };


  ODESolver ode;
  Vector initial_condition{0};
  Vector x_array_reverse = Utils::linspace(xmax, xmin, npts); 

  ode.solve(dtaudx, x_array_reverse, initial_condition);
  auto tau_array = ode.get_data_by_component(0);
  std::reverse(tau_array.begin(), tau_array.end()); 

  // Make the spline
  tau_of_x_spline.create(x_array, tau_array, "opticaldepth");

  // Utils::EndTiming("opticaldepth");
  
  
  //============================================================================
  // Compute visibility functions and spline

  Vector g_tilde(npts);
  for (int i = 0; i<npts; i++){
    g_tilde[i] = -1*tau_of_x_spline.deriv_x(x_array[i]) * exp(-1 * tau_of_x_spline(x_array[i]));
  }
  
  g_tilde_of_x_spline.create(x_array, g_tilde, "g"); // the problem is here!

  Utils::EndTiming("opticaldepth");
}


void RecombinationHistory::calculate_sound_horizon(){
  Utils::StartTiming("sound");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  int xmin = x_start;
  int xmax = 0;
  Vector x_array = Utils::linspace(xmin, xmax, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){

  //=============================================================================
  // Expression for dsdx
  
  const double a           = exp(x);

  // Physical constants in SI units
  const double c           = Constants.c;

  
  // Cosmological parameters
  const double Hp_of_x      = cosmo->Hp_of_x(x);
  const double OmegaB      = cosmo->get_OmegaB(0);
  const double H0          = cosmo->get_H0();
  const double OmegaR      = cosmo->get_OmegaR(0);

  double R = (4*OmegaR) / (3*OmegaB*a);
  double c_s = Constants.c * sqrt((R)/(3*(1+R)));

  // Set the derivative for photon optical depth
    
    dsdx[0] = c_s / Hp_of_x;

    return GSL_SUCCESS;
  };

  const double OmegaR      = cosmo->get_OmegaR(0);
  const double OmegaB      = cosmo->get_OmegaB(0);
  
  
  double initial_c_s = Constants.c * sqrt(((4*OmegaR) / (3*OmegaB*(exp(x_start)))) / (3*(1+((4*OmegaR) / (3*OmegaB*exp(x_start))))));
  double initial_Hp_of_x = cosmo->Hp_of_x(x_start);
  double initial = initial_c_s / initial_Hp_of_x;

  ODESolver ode;
  Vector initial_condition{initial};

  ode.solve(dsdx, x_array, initial_condition);
  auto sound_horizon_array = ode.get_data_by_component(0);

  // Make the spline
  sound_horizon_of_x_spline.create(x_array, sound_horizon_array, "sound");
}

//====================================================
// Get methods

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}


double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{  
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::log_ne_of_x(double x) const{
  return log_ne_of_x_spline(x);
}

double RecombinationHistory::nH_of_x(double x) const{

  double a = exp(x);
  const double OmegaB      = cosmo->get_OmegaB(0);
  const double H0          = cosmo->get_H0();
  
  double n_H = (1.- Yp) * ((3*pow(H0, 2)*OmegaB)/(8.*M_PI*Constants.G*Constants.m_H*pow(a,3)));

  return n_H;
}

double RecombinationHistory::get_sound_horizon(double x) const{
  return sound_horizon_of_x_spline(x);
}


double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << "Today we have a value of $X_e$ of " << Xe_of_x(0) << std::endl;
  std::cout << "The freeze-out abundance of free electrons today is " << exp(log_ne_of_x(0)) << std::endl;
  std::cout << "The sound horizon at decoupling is " << sound_horizon_of_x_spline(-6.9854) / (3.08*pow(10,16)*pow(10,6)) << std::endl;
  std::cout << std::endl;
  
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{ 
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = -12;
  const double x_max   = 0;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << log(ne_of_x(x))      << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << sound_horizon_of_x_spline(x) << " ";

    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}