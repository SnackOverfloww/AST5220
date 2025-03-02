#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================
  
  // We compute OmegaR, OmegaNu, OmegaLambda, H0 TODAY, i.e. technically OmegaR0, OmegaNu0, OmegaLambda0, and of course H0: 
  H0 = h * Constants.H0_over_h;
  OmegaR = 2*((pow(M_PI,2))/30)*((pow(Constants.k_b*TCMB, 4))/(pow(Constants.hbar, 3)*pow(Constants.c,5)))*((8*M_PI*Constants.G)/(3*pow(H0, 2))); 
  OmegaNu = Neff *(7./8.)*pow((4./11.),(4./3.))*OmegaR;
  OmegaLambda = 1 - (OmegaR+OmegaB+OmegaCDM+OmegaK+OmegaNu);
}


// Do all the solving. Compute eta(x) DONE
// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
 
  // Set the range of x and the number of points for the splines. For this Utils::linspace(x_start, x_end, npts) is useful
  const double xmin = -20.;
  const double xmax = 5.;
  const int    npts = 1000;
  
  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

  // Set the rhs of the detadx ODE
     detadx[0] = Constants.c / Hp_of_x(x) ;

    return GSL_SUCCESS;
  };

  
  // TODO: Set the initial condition, set up the ODE system, solve and make  the spline eta_of_x_spline 
  Vector x_array = Utils::linspace(xmin, xmax, npts);
 
  ODESolver ode;

  Vector initial_condition_eta{Hp_of_x(x_start)};
  
  ode.solve(detadx, x_array, initial_condition_eta);
  auto eta_array = ode.get_data_by_component(0);

  // Make the spline
  eta_of_x_spline.create(x_array, eta_array, "eta");

  Utils::EndTiming("Eta");


  // Now, in the same manner, solve the ODE for time, i.e. dt/dx 
  Utils::StartTiming("Time");
  
  // The ODE for deta/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    dtdx[0] = 1/(H_of_x(x));

    return GSL_SUCCESS;
  };

  // TODO: Set the initial condition, set up the ODE system, solve and make the spline t_of_x_spline 
  Vector initial_condition_time{1/(2*H_of_x(x_start))};
  
  ode.solve(dtdx, x_array, initial_condition_time);
  auto time_array = ode.get_data_by_component(0);

 // Make the spline
 t_of_x_spline.create(x_array, time_array, "time");

  Utils::EndTiming("Time");

}

// Get methods

//Function to calculate scale factor a as a function of x 
double BackgroundCosmology::a(double x) const{

  double scale_factor = exp(x);

  return scale_factor;
}

// Function to calculate for H(x) 
double BackgroundCosmology::H_of_x(double x) const{

  double H_of_x_calculated = H0*sqrt(OmegaLambda + (OmegaK*pow(a(x), -2)) + ((OmegaB + OmegaCDM)*pow(a(x), -3)) + ((OmegaR+OmegaNu)*pow(a(x), -4))); 

  return H_of_x_calculated;
}

// Function to calculate for Hprime(x)
double BackgroundCosmology::Hp_of_x(double x) const{
 
  double Hp_of_x_calculated = a(x)* H_of_x(x);

  return Hp_of_x_calculated;
}

// Function to calculate the first derivative of Hprime(x) 
double BackgroundCosmology::dHpdx_of_x(double x) const{
  
 double dHpdx_of_x_calculated = ((1./2.)*H0*pow((((OmegaB + OmegaCDM)*exp(-x)) + ((OmegaR + OmegaNu)*exp(-2.*x)) + OmegaK + (OmegaLambda*exp(2.*x))), -1./2.))
                                *((-1.*(OmegaB + OmegaCDM)*exp(-x)) - (2.*(OmegaR + OmegaNu)*exp(-2.*x)) + (2.*OmegaLambda*exp(2.*x)));

  return dHpdx_of_x_calculated;
}


// // Function to calculate the second derivative of Hprime(x) 
// double BackgroundCosmology::ddHpddx_of_x(double x) const{
  
// double ddHpddx_of_x_calculated = ( ((1./2.)*(-1./2.)*H0*pow((((OmegaB+OmegaCDM)*exp(-x)) + ((OmegaR + OmegaNu)*exp(-2.*x)) + OmegaK + (OmegaLambda*exp(2.*x))), (-3./2.))) 
//                                  *((-1.*(OmegaB + OmegaCDM)*exp(-x)) - (2.*(OmegaR+OmegaNu)*exp(-2.*x)) + (2.*OmegaLambda*exp(2.*x))) )
                                  
//                                  + ( ((1./2.)*H0*pow(((OmegaB + OmegaCDM)*exp(-x)) + ((OmegaR + OmegaNu)*exp(-2.*x)) + OmegaK + (OmegaLambda*exp(2.*x)), (-1./2.)))
//                                  * (((OmegaB + OmegaCDM)*exp(-x)) + (4.*(OmegaR + OmegaNu)*exp(-2.*x)) + (4.*OmegaLambda*exp(2.*x))) );
//   return ddHpddx_of_x_calculated;
// }

// Function to calculate the second derivative of Hprime(x) 
double BackgroundCosmology::ddHpddx_of_x(double x) const{
  
  double ddHpddx_of_x_calculated = (H0 * (1./2.) * (-1./2.) * pow(((OmegaB + OmegaCDM)*exp(-1.*x)) + ((OmegaR + OmegaNu)*exp(-2.*x)) + OmegaK + (OmegaLambda)*exp(2.*x), (-3./2.))
                                    * ((-1.*(OmegaB + OmegaCDM)*exp(-1.*x)) + (-2.*(OmegaR + OmegaNu)*exp(-2.*x)) + (2.*OmegaLambda*exp(2.*x)))
                                    * ((-1.*(OmegaB + OmegaCDM)*exp(-1.*x)) + (-2.*(OmegaR + OmegaNu)*exp(-2.*x)) + (2.*OmegaLambda*exp(2.*x))))
                                    
                                    + (H0 * (1./2.) * pow(((OmegaB + OmegaCDM)*exp(-1.*x) + ((OmegaR + OmegaNu)*exp(-2.*x)) + OmegaK + (OmegaLambda*exp(2.*x))), (-1./2.))
                                    * (((OmegaB + OmegaCDM)*exp(-1.*x)) + (4.*(OmegaR + OmegaNu)* exp(-2.*x)) + (4*(OmegaLambda)*exp(2.*x))));

  return ddHpddx_of_x_calculated;
  }
  



//Function to calculate OmegaB
double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  double Omega_B_of_x = (OmegaB*pow(H0,2))/(pow(a(x), 3) * pow(H_of_x(x), 2)); 
  
  return Omega_B_of_x;
}

//Function to calculate OmegaR
double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  
  double Omega_R_of_x = (OmegaR*pow(H0,2))/(pow(a(x), 4) * pow(H_of_x(x), 2)); 

  return Omega_R_of_x;
}

//Function to calculate OmegaNu 
double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  double Omega_Nu_of_x = (OmegaNu*pow(H0,2))/(pow(a(x), 4) * pow(H_of_x(x), 2)); 

  return Omega_Nu_of_x;
}

//Function to calculate OmegaCMD
double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  
  double Omega_CDM_of_x = (OmegaCDM*pow(H0,2))/(pow(a(x), 3) * pow(H_of_x(x), 2)); 

  return Omega_CDM_of_x;
}

//Function to calculate OmegaLambda
double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  
  double Omega_Lambda_of_x = (OmegaLambda*pow(H0,2))/(pow(H_of_x(x), 2)); 

  return Omega_Lambda_of_x;
}

//Function to calculate OmegaK
double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;
  
  double Omega_K_of_x = (OmegaK*pow(H0,2))/(pow(a(x), 2) * pow(H_of_x(x), 2)); 

  return Omega_K_of_x;
}
  
//Function to calculate comoving distance (Chi)
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{

  double get_comoving_distance_of_x_calculated = eta_of_x(0.) - eta_of_x(x);

  return get_comoving_distance_of_x_calculated;
}

//Function to calculate r as a function of Chi (ultimately a function of x)
double BackgroundCosmology::r_of_chi(double x) const{
  
  double r;

  if (OmegaK < 0) {
    r = get_comoving_distance_of_x(x) * ((sin(sqrt(abs(OmegaK))*H0*(get_comoving_distance_of_x(x)/Constants.c)))
    /(sqrt(abs(OmegaK))*H0*(get_comoving_distance_of_x(x)/Constants.c))); 
  }

  if (OmegaK > 0) {
    r = get_comoving_distance_of_x(x) * ((sinh(sqrt(abs(OmegaK))*H0*(get_comoving_distance_of_x(x)/Constants.c)))
        /(sqrt(abs(OmegaK))*H0*(get_comoving_distance_of_x(x)/Constants.c))); 
  }

  if (OmegaK == 0) {
    r = get_comoving_distance_of_x(x); 
  }

  return r;
}

// Function to calculate the luminosity distance d_L
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{

  double luminosity_distance_of_x_calculated = r_of_chi(x)/a(x);

  return luminosity_distance_of_x_calculated;
}

//Function to calculate the angular diameter distance d_A
double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{

  double anugular_diameter_distance_of_x_calculated = a(x) * r_of_chi(x);

  return anugular_diameter_distance_of_x_calculated;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << "H0:          " << H0          << "\n";
  std::cout << "The age of the Universe is " << t_of_x(0) / (60*60*24*365*pow(10,9))<< " Gyrs \n";
  std::cout << "Matter-radiation equality occur at a time " << t_of_x(-8.131921133813) / (60*60*24*365) << " yrs \n";
  std::cout << "Matter-dark energy equality occur at a time " << t_of_x(-0.2558189708133) / (60*60*24*365*pow(10,9)) << " Gyrs \n";
  std::cout << "The universe starts to accelerate at a time " << t_of_x(-0.4868680309999) / (60*60*24*365*pow(10,9)) << " Gyrs \n";
  std::cout << "Matter-radiation equality occur at a conformal time " << eta_of_x(-8.131921133813) / (Constants.c * 60*60*24*365) << " yrs \n";
  std::cout << "Matter-dark energy equality occur at a conformal time " << eta_of_x(-0.2558189708133) / (Constants.c * 60*60*24*365*pow(10,9)) << " Gyrs \n";
  std::cout << "The universe starts to accelerate at a conformal time " << eta_of_x(-0.4868680309999) / (Constants.c * 60*60*24*365*pow(10,9)) << " Gyrs \n";
  std::cout << "The conformal time of the Universe is " << eta_of_x(0) / (Constants.c * 60*60*24*365*pow(10,9)) << " Gyrs \n";
  //std::cout << "H0:        " << H0        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -20.0;
  const double x_max =  5.0;
  const int    n_pts =  1000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    // deriv_x
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << t_of_x(x)          << " ";
    fp << get_luminosity_distance_of_x(x) << " ";
    // fp << H0                 << " ";
    // fp << get_comoving_distance_of_x(x) << " ";
    // fp << calculating_r_of_chi(x) << " ";
    fp <<"\n";
  };
  
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
