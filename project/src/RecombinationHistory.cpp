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
  
  std::cout << "testing" << std::endl;
  // Compute and spline Xe, ne
  solve_number_density_electrons();

  // std::cout << "testing2" << std::endl;
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//=====================================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//=====================================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on DONE
  //=============================================================================
  
  const double xmin = -12.;
  const double xmax = 0.;
  const int    npts = 1000;
  
  Vector x_array = Utils::linspace(xmin, xmax, npts);
  Vector Xe_arr(npts);
  Vector ne_arr(npts);
  

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < 1000 ; i++){

    //==================================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==================================================================

    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current; 
      ne_arr[i] = log(ne_current);

      std::cout << "Xe = " << Xe_arr[i] << ", n_e = " << ne_arr[i] << std::endl;

    } else {

      if(i==0) std::cout << "Error: i = 0" << std::endl; // Tells us if we skip the whole Saha regime

      //========================================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //======================================================================== Xe / n_h = n_e

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
      double nH = RecombinationHistory::nH_of_x(x_array[i]);
      double ne = nH * Xe_arr[i];
      // std::cout << "Xe = " << Xe_arr[i] << ", n_e = " << ne << std::endl;
      ne_arr[i] = log(ne); 
      std::cout << "Xe = " << Xe_arr[i] << ", n_e = " << ne_arr[i] << std::endl;
      
    }

    //===============================================================================
    // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
    // functions are working
    //===============================================================================
    
    // Make the splines
    Xe_of_x_spline.create(x_array, Xe_arr, "Xe");
    log_ne_of_x_spline.create(x_array, ne_arr, "ne");

  }

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
  // std::cout << "a is" << a << std::endl;
  // std::cout << "x is" << x << std::endl;
 
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
  const double OmegaR      = cosmo->get_OmegaR(0);
  const double TCMB_0      = 2.7255;

  // Compute critical density
  double critical_density = (3.*pow(H0, 2.))/(8.*M_PI*Constants.G); 

  
  // Compute baryon density n_b
  double baryon_density = (OmegaB * critical_density)/(Constants.m_H * pow(a, 3.));
  std::cout << "Baryon density is " << baryon_density << std::endl;

  //Compute baryon temprature
  double baryon_temperature = TCMB_0/a;
    
  // Compute constant for solving quadratc for X_e
  double constant_for_quadratic = (1./(baryon_density * pow(Constants.hbar,3.)))*pow(((Constants.m_e*baryon_temperature*Constants.k_b)/(2.*M_PI)), (3./2.))*exp(-Constants.epsilon_0 / (baryon_temperature * Constants.k_b));
  
  double solving_quadratic = (-constant_for_quadratic + sqrt(pow(constant_for_quadratic, 2.) + (4*constant_for_quadratic)))/2.;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  double Xe;
  double ne;
  if (4/constant_for_quadratic > 0.001){
    
    Xe = solving_quadratic;
    ne = Xe * baryon_density;
  } 
  
  else{
    Xe = 1;
    ne = baryon_density;

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
  const double OmegaR      = cosmo->get_OmegaR(0);
  const double OmegaCDM    = cosmo->get_OmegaCDM(0);
  const double OmegaNu     = cosmo->get_OmegaNu(0);
  const double OmegaLambda = cosmo->get_OmegaLambda(0);
  const double OmegaK      = cosmo->get_OmegaK(0);
  const double H_of_x      = cosmo->H_of_x(x);
  const double TCMB_0      = 2.7255; // TODO


  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================

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

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(-12, 0, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================

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
  const double OmegaR      = cosmo->get_OmegaR(0);
  const double OmegaCDM    = cosmo->get_OmegaCDM(0);
  const double OmegaNu     = cosmo->get_OmegaNu(0);
  const double OmegaLambda = cosmo->get_OmegaLambda(0);
  const double OmegaK      = cosmo->get_OmegaK(0);
  const double TCMB_0      = 2.7255;


    // Set the derivative for photon optical depth
    
    double H_of_x = H0*sqrt(OmegaLambda + (OmegaK*pow(exp(x), -2)) + ((OmegaB + OmegaCDM)*pow(exp(x), -3)) + ((OmegaR+OmegaNu)*pow(exp(x), -4)));
    double ne = exp(ne_of_x(x)); //need to inverse of the ne values here! 
    
    dtaudx[0] = ( -Constants.c *ne* Constants.sigma_T) / H_of_x; // need to inverse the H_of_x_values here!

    return GSL_SUCCESS;
  };


  ODESolver ode;
  Vector initial_condition{0};
  Vector x_array_reverse = Utils::linspace(0, -12, npts); 

  ode.solve(dtaudx, x_array, initial_condition);
  auto tau_array = ode.get_data_by_component(0);

  // Make the spline
  tau_of_x_spline.create(x_array, tau_array, "opticaldepth");

  // Utils::EndTiming("opticaldepth");
  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}


double RecombinationHistory::dtaudx_of_x(double x) const{

  //==================================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //==================================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return log_ne_of_x_spline(x);
}

double RecombinationHistory::nH_of_x(double x) const{

  double a = exp(x);
  const double OmegaB      = cosmo->get_OmegaB(0);
  const double H0          = cosmo->get_H0();
  
  double n_H = (1.- Yp) * ((3*pow(H0, 2)*OmegaB)/(8.*M_PI*Constants.G*Constants.m_H*pow(a,3)));

  return n_H;
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
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    // fp << dtaudx_of_x(x)       << " ";
    // fp << ddtauddx_of_x(x)     << " ";
    // fp << g_tilde_of_x(x)      << " ";
    // fp << dgdx_tilde_of_x(x)   << " ";
    // fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

