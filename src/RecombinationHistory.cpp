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

// Declare relevant constants 
const double m_e         = Constants.m_e;
const double epsilon_0   = Constants.epsilon_0;
const double hbar        = Constants.hbar;
const double k_b         = Constants.k_b;
const double m_H         = Constants.m_H;
const double G           = Constants.G;
const double c           = Constants.c;
const double H0_over_h   = Constants.H0_over_h;
const double sigma_T     = Constants.sigma_T;
const double lambda_2s1s = Constants.lambda_2s1s;

// For Saha equation to not blow up 
const double global_tol = 1e-7;

void RecombinationHistory::set_cosmological_constants(){
  OmegaB = cosmo -> get_OmegaB();
  TCMB = cosmo -> get_TCMB();
  H0 = cosmo -> get_H0();
  rho_c0 = 3.0*H0 * H0 / (8.0 * M_PI * G);   // Value for critical density today
}

void RecombinationHistory::solve(){

  // Set constants
  set_cosmological_constants();
  std::cout << "Set cosmological constants\n";
  
  // Compute and spline Xe, ne
  solve_number_density_electrons();
  std::cout << "Solved number density electrons\n";
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
  std::cout << "Solved optical depth\n";
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_arr_saha(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  int i_last_saha = 0;

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
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
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      Xe_arr_saha[i] = Xe_current;
      ne_arr[i] = ne_current;
    } 

    else {
      i_last_saha = i-1;
      std::cout << "Break loop!\n";
      break;
      }
    }
  //==============================================================
  // TODO: Compute X_e from current time til today by solving 
  // the Peebles equation (NB: if you solve all in one go remember to
  // exit the for-loop!)
  // Implement rhs_peebles_ode
  //==============================================================
  const int npts_ode_array = npts_rec_arrays - i_last_saha;
  Vector x_ode(npts_ode_array);
  for (int j=0; j < npts_ode_array; j++){
    x_ode[j] = x_array[j + i_last_saha];
  }

  // The Peebles ODE equation
  ODESolver peebles_Xe_ode;
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    // Current value of a and X_e
    const double X_e = Xe[0];
    const double a = exp(x);

    // Constants for finding RHS of peebles eq.
    const double Tb = cosmo -> get_TCMB(x);
    const double eps_Tb = epsilon_0 / (k_b*Tb);
    const double alpha = c / hbar*sqrt(3.0*sigma_T/(8.0*M_PI))*m_e; // Dimensionless fine structure constant. 
    const double H = cosmo -> H_of_x(x);
    
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
    const double Cr = (lambda_2s1s + Lambda_alpha) / (lambda_2s1s + Lambda_alpha + beta2); // dimensionless

    const double rhs = Cr/H * (beta*(1.0 - X_e) - nH*alpha2*pow(X_e, 2));

    dXedx[0] = rhs;

    return GSL_SUCCESS;
  };
    
  //=============================================================================
  // TODO: Set up IC, solve the ODE and fetch the result 
  //=============================================================================

  double Xe_init_val = Xe_arr[i_last_saha];
  Vector Xe_init_vec{Xe_init_val};
  
  peebles_Xe_ode.solve(dXedx, x_ode, Xe_init_vec);
  auto Xe_ode = peebles_Xe_ode.get_data_by_component(0);
  for (int j=i_last_saha; j < npts_rec_arrays; j++){
    Xe_arr[j] = Xe_ode[j - i_last_saha];
    ne_arr[j] = Xe_ode[j - i_last_saha] * nH_of_x(x_array[j]);
    }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================

  // Make splines of the result
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
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  double Tb = cosmo -> get_TCMB(x);
  double nb = nb_of_x(x);

  double K = 1./nb * pow(k_b*m_e*Tb/(2.0*M_PI*pow(hbar, 2)), 1.5) * exp(-epsilon_0/(k_b*Tb));
  if(4.0/K < global_tol)
    Xe = 1.0;
  else{
    Xe = K/2. * (-1+sqrt(1 + 4./K));
  }

  const double nH = nH_of_x(x);
  ne = Xe*nH;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  // Want to integrate backward in time. ODESolver need increasing values, hence minus signs. 
  Vector x_array_tau_rev = Utils::linspace(-x_end, -x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODESolver tauOde;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    // Set the derivative for photon optical depth
    dtaudx[0] = c*sigma_T*ne_of_x(-x) / (cosmo -> H_of_x(-x));
    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  // Initial Xe will be the last added value from the Saha regime. 
  Vector tau_init_vec{0.0};

  tauOde.solve(dtaudx, x_array_tau_rev, tau_init_vec);

  auto tau_vec_inv = tauOde.get_data_by_component(0);

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
  // TODO: Compute visibility functions and spline everything
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
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return dg_tildedx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return ne_of_x_spline(x);
}

double RecombinationHistory::nb_of_x(double x) const{
  return OmegaB * rho_c0 / m_H * exp(-3.0*x);
}

double RecombinationHistory::nH_of_x(double x) const{
  return nb_of_x(x); // We have no helium on project
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
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

