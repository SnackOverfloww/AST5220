#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){
  double eta_0      = cosmo->eta_of_x(0);
  double k_start    = k_min;
  double k_end      = 0.3 / Constants.Mpc;
  double x_start    = Constants.x_start;
  double x_end      = Constants.x_end;
  int n_k_los        = std::ceil(std::abs(k_end - k_start) / ((2.*M_PI)/(6.*eta_0)));
  // std::cout << "n_k_los is: " << n_k_los << std::endl;
  int n_k_C_l        = std::ceil(std::abs(log(k_end) - log(k_start)) / ((2*M_PI)/(6*eta_0 * k_end)));
  // std::cout << "n_k_C_l is: " << n_k_C_l << std::endl;
 
  


  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector k_array_los = Utils::linspace(k_start, k_end, n_k_los);
  
  Vector log_k_array_C_l = (Utils::linspace(log(k_start), log(k_end), n_k_C_l));
  
  generate_bessel_function_splines();

  line_of_sight_integration(k_array_los);

  auto cell_TT = solve_for_cell(log_k_array_C_l, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  if(Constants.polarization){
    auto cell_TE = solve_for_cell(log_k_array_C_l, thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");

    auto cell_EE = solve_for_cell(log_k_array_C_l, thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_EE_spline.create(ells, cell_EE, "Cell_EE_of_ell");
  }

}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  
  double eta_0      = cosmo->eta_of_x(0);
  double k_end = 0.3 / Constants.Mpc;

  Vector z = Utils::linspace(0, eta_0 * k_end, int((eta_0 * k_end) / ((2*M_PI)/16)));
  Vector j_ell_vec(z.size());

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    for (size_t iz = 0; iz < z.size(); iz++){
        j_ell_vec[iz] = Utils::j_ell(ell, z[iz]); 
    }
    j_ell_splines[i].create(z, j_ell_vec, "besselspline"); 
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array_los, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  double x_start        = Constants.x_start;
  double x_end          = Constants.x_end;

  double delta_x = 0.03;

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array_los.size()));
  Vector x_array = Utils::linspace(x_start, x_end, int((x_end - x_start) / delta_x));

  double eta;
  double eta_0 = cosmo->eta_of_x(0);
  double x;
  double k;
  double integral;
  double function_to_integrate;

  for(size_t il = 0; il < ells.size(); il++){
    for(size_t ik = 0; ik < k_array_los.size(); ik++){
      k = k_array_los[ik];
      
      // Goal: integrate Theta_ell(k) = Int F(x,k) dx
      integral = 0.0;
      for(int ix = 0; ix < x_array.size(); ix++){
        x = x_array[ix];
        eta = cosmo->eta_of_x(x);

        function_to_integrate = j_ell_splines[il](k*(eta_0 - eta)) * source_function(x, k); 
      
        integral += function_to_integrate * delta_x;

      }
      // Store the result for Source_ell(k) in results[ell][ik]
      result[il][ik] = integral; //Temperature multipole (theta) for all ell-values
  }
}
  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
void PowerSpectrum::line_of_sight_integration(Vector & k_array_los){
  const int n_k         = k_array_los.size();
  const int nells       = ells.size();

  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array_los, source_function_T);
  
  // Spline the result and store it in thetaT_ell_of_k_spline
  for(int j = 0; j < nells; j++){
    thetaT_ell_of_k_spline[j].create(k_array_los, thetaT_ell_of_k[j], "ThetaT_spline");
  }

  
  if(Constants.polarization){

    std::function<double(double,double)> source_function_E = [&](double x, double k){
    return pert->get_Source_E(x,k);
    };

    thetaE_ell_of_k_spline = std::vector<Spline>(nells);

    Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array_los, source_function_E);

    for(int j = 0; j < nells; j++){
      thetaE_ell_of_k_spline[j].create(k_array_los, thetaE_ell_of_k[j], "ThetaE_spline");
    }

  }
}



//====================================================
// Compute Cell (could be TT or TE or EE) 
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array_C_l,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size(); 

  // double eta_0 = cosmo->eta_of_x(0);
  double integral;
  double delta_log_k = log_k_array_C_l[1] - log_k_array_C_l[0];
  double k;
  double function_to_integrate;
  Vector result(nells);


  for(size_t il = 0; il < ells.size(); il++){
    integral = 0.0;
    for(size_t ik = 0; ik < log_k_array_C_l.size(); ik++){
      k = exp(log_k_array_C_l[ik]);

      function_to_integrate = 4 * M_PI * primordial_power_spectrum(k) * f_ell_spline[il](k) * g_ell_spline[il](k); 

      integral += function_to_integrate * delta_log_k;
      
    }
  
    result[il] = integral;
  }



  return result;
}

//====================================================
// The primordial power-spectrum
double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  double a = exp(x);
  const double H_0          = cosmo->get_H0();
  const double Phi          = pert->get_Phi(x, k_mpc / Constants.Mpc);
  const double Omega_B      = cosmo->get_OmegaB(0); 
  const double Omega_CDM    = cosmo->get_OmegaCDM(0);
  
  double Omega_M            = Omega_B + Omega_CDM;

  double delta_M = (pow(Constants.c, 2) * pow(k_mpc / Constants.Mpc, 2) * Phi) / ((3./2. * Omega_M * pow(a, -1) * pow(H_0, 2)));
  
  pofk = (pow(delta_M, 2) * primordial_power_spectrum(k_mpc / Constants.Mpc)) * ((2*pow(M_PI, 2)) / pow((k_mpc / Constants.Mpc),3))* pow(cosmo->get_h(), 3) / pow(Constants.Mpc, 3);  
  
  return pofk;
}

//====================================================
// Get methods
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}
double PowerSpectrum::get_transfer_function(const double k, int ell_ix) const{
  return thetaT_ell_of_k_spline[ell_ix](k); 
}


//====================================================
// Output the cells to file
void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorEE = pow(1e6 * cosmo->get_TCMB(), 2) * (((ell+2)*(ell+1)*(ell)*(ell-1)) / ((pow(10, -5))));
    double normfactorTE = pow(1e6 * cosmo->get_TCMB(), 2) * sqrt(((ell+2)*(ell+1)*(ell)*(ell-1))) * (ell * (ell+1)) / (2.0 * M_PI)  ; 
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << get_cell_TT(ell) * normfactor  << " ";
    if(Constants.polarization){
      fp << get_cell_TE(ell) * normfactorTE  << " ";
      fp << get_cell_EE(ell) * normfactorEE  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}



// I know hard coding is bad but i keep getting errors if I dont...
void PowerSpectrum::output_transfer(std::string filename, const int ell_ix_1, 
  const int ell_ix_2, const int ell_ix_3, const int ell_ix_4, const int ell_ix_5) const{

  double k_start = k_min;
  double k_end = k_max;

  Vector k_array = exp(Utils::linspace(log(0.00005 / Constants.Mpc), log(0.3 / Constants.Mpc), 10000));
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double k) {
    fp << k * (cosmo->eta_of_x(0)) << " ";
    fp << get_transfer_function(k, ell_ix_1) << " ";
    fp << get_transfer_function(k, ell_ix_2) << " ";
    fp << get_transfer_function(k, ell_ix_3) << " ";
    fp << get_transfer_function(k, ell_ix_4) << " ";
    fp << get_transfer_function(k, ell_ix_5) << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}



void PowerSpectrum::output_pofk(std::string filename, const double x) const{

  double k_start = k_min * Constants.Mpc;
  double k_end = k_max * Constants.Mpc;

  Vector k_array = exp(Utils::linspace(log(k_start), log(k_end), 1000));
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double k) {
    fp << k / cosmo->get_h() << " ";
    fp << get_matter_power_spectrum(x, k) << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
  
}