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

// Function to di trapezoidal integration 
double PowerSpectrum::integrate_trapezoid(double dx, Vector y_array){
  // Declare and define needed variables 
  double integral_value = 0;
  for(int i=0; i < y_array.size()-1; i++){
    integral_value += 0.5*(y_array[i+1] + y_array[i]);
  }
  return integral_value * dx;
}

// Function to get sufficient amount of points in linspace from given delta
Vector PowerSpectrum::get_linspace(double x_start, double x_end, double dx){
  int npts = (x_end - x_start)/dx;
  return Utils::linspace(x_start, x_end, npts);
}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================

  double dk_theta = 2.*M_PI / (eta0 * n_k);
  Vector k_theta_array = get_linspace(k_min, k_max, dk_theta);

  //=========================================================================
  // Make splines for j_ell. 
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================

  line_of_sight_integration(k_theta_array);

  //=========================================================================
  // Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================

  // Cannot use get_linspace(log(k_min), log(k_max), dk_ps) as this over-estimates npts
  double dk_ps = 2.*M_PI / (eta0 * n_k_ps);
  int npts_ps = (k_max - k_min)/dk_ps;
  Vector log_k_ps_array = Utils::linspace(log(k_min), log(k_max), npts_ps);

  auto cell_TT = solve_for_cell(log_k_ps_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  int const size_ell = ells.size();
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(size_ell);
    
  // Determine argument interval
  double z_min = 0.0;
  double z_max = k_max*eta0;
  double dz = 2.*M_PI/n_bessel;
  Vector z_array = get_linspace(z_min, z_max, dz);

  #pragma omp parallel for schedule(dynamic, 1)
  for(int il=0; il < size_ell; il++){
    const int ell  = ells[il];

    // Make vector to store j_ell
    Vector j_ell_array(z_array.size());

    //  Loop over z-array
    for(int iz=0; iz < z_array.size(); iz++){
      j_ell_array[iz] = Utils::j_ell(ell, z_array[iz]);
    }
    j_ell_splines[il].create(z_array, j_ell_array);
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================


Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  // Create arrays
  double dx = 2.*M_PI/n_x_los;
  Vector x_array = get_linspace(x_start_los, x_end_los, dx);

  #pragma omp parallel for schedule(dynamic, 1)
  for(int ik=0; ik < k_array.size(); ik++){
     // Progress bar...
    if(ik*10/k_array.size() != ((ik + 1)*10) / k_array.size()){
      std::cout << (ik + 1)*100 / k_array.size() << "% " << std::flush;
      if(ik == k_array.size() - 1) std::cout << std::endl;
    }
    double k_value = k_array[ik]; // k-value for each iteration
    for(int il=0; il < ells.size(); il++){
      double ell = ells[il]; // ell-value for each iteration

      Vector integrand(x_array.size());
      for(int i=0; i < x_array.size(); i++){
        integrand[i] = source_function(x_array[i], k_value) * j_ell_splines[il](k_value*(eta0 - cosmo->eta_of_x(x_array[i])));
      }
      result[il][ik] = integrate_trapezoid(dx, integrand);
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int nells = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert -> get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

 for(int il=0; il < ells.size(); il++){
    thetaT_ell_of_k_spline[il].create(k_array, thetaT_ell_of_k[il]);
  }
  
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  // I choose to do the dCell/dlogk
  //============================================================================

  Vector result(nells);
  int N = log_k_array.size();
  double dlogk = (log_k_array[N-1] - log_k_array[0])/N;

  // Loop over and integrate for all ells
  for(int il=0; il < nells; il++){
    double ell = ells[il];
    Vector integrand(log_k_array.size());
    for(int i=0; i < log_k_array.size(); i++){
      double k_value = exp(log_k_array[i]);
      integrand[i] = primordial_power_spectrum(k_value) * abs(f_ell_spline[il](k_value)*g_ell_spline[il](k_value));
    }

    result[il] = 4.*M_PI * integrate_trapezoid(dlogk, integrand);
  }

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow(Constants.Mpc*k/kpivot_mpc, n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  // Variables and constants
  double OmegaM = cosmo -> get_OmegaM(x);
  double Hp     = cosmo -> Hp_of_x(x);
  double Phi    = pert -> get_Phi(x, k_mpc);
  double c      = Constants.c;
  double P_prim = 2.*pow(M_PI, 2.)*primordial_power_spectrum(k_mpc) / pow(k_mpc, 3.);

  // Calculate Delta_M
  double Delta_M = 2.0*pow(c*k_mpc, 2.0)*Phi / (3.0*OmegaM*pow(Hp, 2.));

  // Calculate P(k,x)
  double pofk = pow(Delta_M, 2.) * P_prim ;
  
  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}

double PowerSpectrum::get_bessel_func(const int il_index, const double z) const{
  return j_ell_splines[il_index](z);
}

double PowerSpectrum::get_thetaT_ell_of_k_spline(const int il_index, const double k) const{
  return thetaT_ell_of_k_spline[il_index](k);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell + 1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(0.0), 2.);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline(ell) * normfactor  << " ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

void PowerSpectrum::output_matter_PS(std::string filename) const{
  std::ofstream fp(filename.c_str());
  double Mpc = Constants.Mpc;
  double h = cosmo -> get_h();
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, 10000);
  Vector k_array = exp(log_k_array);

  auto print_data = [&] (const double k) {
    fp << k * Mpc/h                           << " ";
    fp << get_matter_power_spectrum(0.0, k) * pow(h/Mpc, 3.)  << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output_theta(std::string filename) const{
  // Output of theta l for l={6, 100, 200, 500, 1000}
  std::ofstream fp(filename.c_str());
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, 10000); // Linearly spaced logarithmic values
  Vector k_array = exp(log_k_array);
  double c = Constants.c;
  double H0 = cosmo -> get_H0();

  auto print_data = [&] (const double k){
    fp << k*eta0 << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[0], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[1], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[2], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[3], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[4], k) << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output_bessel_function(std::string filename) const{
  // Output bessel functions
  std::ofstream fp(filename.c_str());
  Vector z_array = Utils::linspace(0.0,1000,5000);
  auto print_data = [&] (const double z) {
    fp << z << " ";
    fp << get_bessel_func(test_ell_index[0], z) << " ";
    fp << get_bessel_func(test_ell_index[1], z) << " ";
    fp << get_bessel_func(test_ell_index[2], z) << " ";
    fp << get_bessel_func(test_ell_index[3], z) << " ";
    fp << get_bessel_func(test_ell_index[4], z) << " ";
    fp << "\n";
  };
  std::for_each(z_array.begin(), z_array.end(), print_data);
}

