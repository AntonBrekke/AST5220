#include"RecombinationHistory.h"

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

// Declare relevant constants in class 
const double m_e         = Constants.m_e;
const double epsilon_0   = Constants.epsilon_0;
const double hbar        = Constants.hbar;
const double k_b         = Constants.k_b;
const double m_H         = Constants.m_H;
const double G           = Constants.G;
const double c           = Constants.c;
const double H0_over_h   = Constants.H0_over_h;
const double sigma_T     = Constants.sigma_T;
const double Lambda_2s1s = Constants.lambda_2s1s;

// For Saha equation to not blow up 
const double global_tol = 2e-4;

// Instead of setting the same variables everywhere, I just declare them in class
void RecombinationHistory::set_cosmological_constants(){
  OmegaB = cosmo -> get_OmegaB();
  TCMB = cosmo -> get_TCMB();
  H0 = cosmo -> get_H0();
  rho_c0 = 3.0*H0 * H0 / (8.0 * M_PI * G);   // Value for critical density today
}

void RecombinationHistory::solve(){

  // Set cosmological constants
  set_cosmological_constants();
  std::cout << "Set cosmological constants\n";
  
  // Compute and spline Xe, ne

  solve_number_density_electrons();
  std::cout << "Solved number density electrons\n";
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
  std::cout << "Solved optical depth\n";

  // Compute and spline sound horizon
  solve_for_sound_horizon();
  std::cout << "Solved sound horizon\n";
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");

  // Set up x-array and make arrays to store X_e(x) and n_e(x)
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_arr_saha(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  int i_last_saha = 0;

  // Calculate recombination history
  bool saha_regime = true;
  for (int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;
    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      // Store the result we got from the Saha equation
      Xe_arr[i] = Xe_current;
      Xe_arr_saha[i] = Xe_current;
      ne_arr[i] = ne_current;
    } 

    else {
      // Store last index that was solved with Saha-equation
      i_last_saha = i-1;
      std::cout << "Break loop!\n";
      break;
      }
    }

  // Fill rest of Xe_arr_saha with Saha solution 
  // Since lim K->0 K/2. * (-1 + sqrt(1 + 4./K)) = 0, we can safely return global tol 
  // when Saha-equation no longer work (negative, nan)
  for(int i=i_last_saha + 1; i < npts_rec_arrays; i++){
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    double const Xe_current = Xe_ne_data.first;
    // Check for non-zero and negative values
    double Xe_current_non_neg_zero_nan = Xe_current < global_tol || std::isnan(Xe_current) ? global_tol : Xe_current;
    Xe_arr_saha[i] = Xe_current_non_neg_zero_nan;
  }

  // Make an array with x-values that has not been solved for,
  // to be solved using Peebles equation
  const int npts_ode_array = npts_rec_arrays - i_last_saha;
  Vector x_array_ode(npts_ode_array);
  for (int j=0; j < npts_ode_array; j++){
    x_array_ode[j] = x_array[j + i_last_saha];
  }

  // The Peebles ODE equation
  ODESolver peebles_Xe_ode;
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    return rhs_peebles_ode(x, Xe, dXedx);
  };
    
  // Setting initial conditions 
  double Xe_init_val = Xe_arr[i_last_saha];
  Vector Xe_init_vec{Xe_init_val};
  
  peebles_Xe_ode.solve(dXedx, x_array_ode, Xe_init_vec);
  auto Xe_ode = peebles_Xe_ode.get_data_by_component(0);
  for (int j=i_last_saha; j < npts_rec_arrays; j++){
    Xe_arr[j] = Xe_ode[j - i_last_saha];
    ne_arr[j] = Xe_ode[j - i_last_saha] * nH_of_x(x_array[j]);
    }
 
  // Make splines of the results
  Xe_of_x_spline.create(x_array, Xe_arr, "Xe");
  XeSaha_of_x_spline.create(x_array, Xe_arr_saha, "XeSaha");
  ne_of_x_spline.create(x_array, ne_arr, "ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a = exp(x);
 
  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // Compute Xe and ne from the Saha equation
  //=============================================================================
  double Tb = cosmo -> get_TCMB(x);
  double nb = nb_of_x(x);

  double K = 1. / nb * pow(k_b*m_e*Tb / (2.0*M_PI*pow(hbar, 2)), 1.5) * exp(-epsilon_0/(k_b*Tb));
  if(4.0/K < global_tol)
    Xe = 1.0;
  else{
    Xe = K/2. * (-1 + sqrt(1 + 4./K));
  }

  const double nH = nH_of_x(x);
  ne = Xe*nH;

  return std::pair<double,double>(Xe, ne);
}


int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){
    // Current value of a and X_e
    const double X_e = Xe[0];
    const double a = exp(x);

    // Constants for finding RHS of peebles eq.
    const double Tb = cosmo -> get_TCMB(x);
    const double H = cosmo -> H_of_x(x);
    const double eps_Tb = epsilon_0 / (k_b*Tb);
    const double alpha = c / hbar*sqrt(3.0*sigma_T/(8.0*M_PI))*m_e; // Dimensionless fine structure constant. 
    
    // Write down all terms/constants involved in RHS 
    const double phi2 = 0.448 * log(eps_Tb); // dimensionless
    const double alpha2 = pow(hbar, 2) / c * 64.0*M_PI*pow(alpha, 2)*phi2 / (pow(m_e, 2)) * sqrt(eps_Tb/(27.0*M_PI)); // dimension m^3/s
    const double beta = alpha2* pow(k_b*m_e*Tb/(2.0*M_PI*pow(hbar, 2)), 1.5) * exp(-eps_Tb);
    double beta2; // dimension 1/s

    // Checking for large exponent to avoid overflow
    if(eps_Tb > 200.0){
      beta2 = 0.0;
    }
    else{
      beta2 = beta*exp(eps_Tb*3.0/4.0); // dimension 1/s
    }

    const double nH = rho_c0*OmegaB / (m_H * pow(a, 3)); // dimension 1/m^3. No Helium -> Y_p = 0 
    const double n1s = (1.0 - X_e)*nH; // dimension 1/m^3
    const double Lambda_alpha = H * pow(3*epsilon_0, 3) / (pow(hbar*c, 3)) / (pow(8*M_PI, 2) * n1s); // dimension 1/s
    const double Cr = (Lambda_2s1s + Lambda_alpha) / (Lambda_2s1s + Lambda_alpha + beta2); // dimensionless

    const double rhs = Cr/H * (beta*(1.0 - X_e) - nH*alpha2*pow(X_e, 2));

    dXedx[0] = rhs;

    return GSL_SUCCESS;
};

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  // Want to integrate backwards in time. ODESolver need increasing values, hence minus signs. 
  Vector x_array_tau_rev = Utils::linspace(-x_end, -x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODESolver tau_ode;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Set the derivative for photon optical depth
    dtaudx[0] = c*sigma_T*ne_of_x(-x) / (cosmo -> H_of_x(-x));
    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set up and solve the ODE and make tau splines
  //=============================================================================
  // Initial Xe will be the last added value from the Saha regime. 
  Vector tau_init_vec{0.0};

  tau_ode.solve(dtaudx, x_array_tau_rev, tau_init_vec);

  auto tau_vec_inv = tau_ode.get_data_by_component(0);

  // Reverse the arrays to get correct order in elements 
  Vector x_array_tau(npts);
  Vector tau_vec(npts);
  for(int i=0; i < npts; i++){
    x_array_tau[i] = -x_array_tau_rev[npts-1 - i];
    tau_vec[i] = tau_vec_inv[npts-1 - i];
  }

  // Generate dtaudx
  Vector dtau_vec(npts);
  for(int i=0; i < npts; i++){
    double const x_i = x_array_tau[i];
    dtau_vec[i] = -c*ne_of_x(x_i)*sigma_T / (cosmo->H_of_x(x_i));
  }

  //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================
  tau_of_x_spline.create(x_array_tau, tau_vec, "tau");
  dtaudx_of_x_spline.create(x_array_tau, dtau_vec, "dtaudx");

  Utils::EndTiming("opticaldepth");
  Utils::StartTiming("Visibility");
   // Generate gtilde
  Vector g_tilde(npts);
  for(int i=0; i < npts; i++){
    double const x_i = x_array_tau[i];
    g_tilde[i] = -dtaudx_of_x(x_i) * exp(-tau_of_x(x_i));
  }


  // Generate dgtildedx
  Vector dg_tildedx(npts);
  for(int i=0; i < npts; i++){
    double x_i = x_array_tau[i];
    dg_tildedx[i] = exp(-tau_of_x(x_i)) * (dtaudx_of_x(x_i)*dtaudx_of_x(x_i) - ddtauddx_of_x(x_i));
  }
  
  g_tilde_of_x_spline.create(x_array_tau, g_tilde, "g");
  dg_tildedx_of_x_spline.create(x_array_tau, dg_tildedx, "dgdx");

  Utils::EndTiming("Visibility");
}


//====================================================
// Solve for the sound horizon s and spline the result
//====================================================

void RecombinationHistory::solve_for_sound_horizon(){
  Utils::StartTiming("Sound horizon");

  const int npts = 1000;
  Vector x_array_s = Utils::linspace(x_start, x_end, npts);

  ODESolver sound_horizon_ode;
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    dsdx[0] = cs_of_x(x) / cosmo -> Hp_of_x(x);
    return GSL_SUCCESS;
  };
  double init_val = cs_of_x(x_array_s[0])/cosmo -> Hp_of_x(x_array_s[0]);
  Vector s_init_vec{init_val};
  sound_horizon_ode.solve(dsdx, x_array_s, s_init_vec);
  auto s_vec = sound_horizon_ode.get_data_by_component(0);
  s_of_x_spline.create(x_array_s, s_vec, "s");

  Utils::EndTiming("Sound horizon");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return Xe_of_x_spline(x);
}
double RecombinationHistory::XeSaha_of_x(double x) const{
  return XeSaha_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return ne_of_x_spline(x);
}

double RecombinationHistory::nb_of_x(double x) const{
  return OmegaB * rho_c0 / m_H * exp(-3.0*x);
}

double RecombinationHistory::nH_of_x(double x) const{
  return nb_of_x(x); // We have no helium on project
}

double RecombinationHistory::R_of_x(double x) const{
  return (3.0 * cosmo -> get_OmegaB(x))/(4.0*cosmo -> get_OmegaR(x));
}

double RecombinationHistory::cs_of_x(double x) const{
  double R = R_of_x(x);
  return c*sqrt(1. / (3.*(1. + R)));
}

double RecombinationHistory::s_of_x(double x) const{
  return s_of_x_spline(x);
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
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                      << " ";
    fp << cosmo -> t_of_x(x)     << " ";
    fp << cosmo -> get_z_of_x(x) << " ";
    fp << Xe_of_x(x)             << " ";
    fp << XeSaha_of_x(x)         << " ";
    fp << ne_of_x(x)             << " ";
    fp << tau_of_x(x)            << " ";
    fp << dtaudx_of_x(x)         << " ";
    fp << ddtauddx_of_x(x)       << " ";
    fp << g_tilde_of_x(x)        << " ";
    fp << dgdx_tilde_of_x(x)     << " ";
    fp << ddgddx_tilde_of_x(x)   << " ";
    fp << s_of_x(x)              << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

