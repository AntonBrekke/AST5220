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
  H0 = h * Constants.H0_over_h;
  OmegaR = 2. * pow(Constants.pi, 2) / 30 * pow(Constants.k_b * TCMB, 4) / (pow(Constants.hbar, 3) 
             * pow(Constants.c, 5)) * 8 * Constants.pi * Constants.G / (3 * pow(H0, 2));
  OmegaNu = Neff * 7. / 8 * pow(4 / 11, 4 / 3) * OmegaR;
  OmegaLambda = 1. - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu); 

}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  // Utils::StartTiming("Eta");
  // Utils::StartTiming("Time");  
  //=============================================================================
  // Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  double x_start = Constants.x_start;
  double x_end = Constants.x_end;
  int npts = 1e4;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  //=============================================================================
  // Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline and t_of_x_spline
  //=============================================================================

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    //=============================================================================
    // Set the rhs of the detadx ODE
    //=============================================================================
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };
  
  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    //=============================================================================
    // Set the rhs of the dtdx ODE
    //=============================================================================
    dtdx[0] = 1. / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // Initial condition for eta and t
  double eta_0 = Constants.c / Hp_of_x(x_start);
  double t_0 = 1. / (2*H_of_x(x_start));
  Vector eta_i{eta_0};
  Vector t_i{t_0};

  // Solve using ODESolver 
  ODESolver ODE_eta_of_x; 
  ODE_eta_of_x.solve(detadx, x_array, eta_i);
  Vector eta_array = ODE_eta_of_x.get_data_by_component(0); 
  eta_of_x_spline.create(x_array, eta_array, "eta_of_x");

  ODESolver ODE_t_of_x; 
  ODE_t_of_x.solve(dtdx, x_array, t_i);
  Vector t_array = ODE_t_of_x.get_data_by_component(0);
  t_of_x_spline.create(x_array, t_array, "t_of_x");

  // Utils::EndTiming("Eta");
  // Utils::EndTiming("Time");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  //=============================================================================
  // Using x = log(a) with natural logarithm. 
  //=============================================================================
  double H = H0 * sqrt((OmegaNu + OmegaR)*exp(-4*x) + (OmegaCDM + OmegaB)*exp(-3*x) + OmegaK*exp(-2*x) + OmegaLambda);

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  //=============================================================================
  // Using x = log(a) with natural logarithm.
  //=============================================================================
  double Hp = exp(x) * H_of_x(x);

  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  //=============================================================================
  // Take derivative of the above Hp_of_x.
  //=============================================================================
  double dHpdx = 0.5*pow(H0 / H_of_x(x), 2) * (OmegaLambda*exp(x) - 2*(OmegaR + OmegaNu)*exp(-3*x) - (OmegaB + OmegaCDM)*exp(-2*x));

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  //=============================================================================
  // Ugly derivative, double check answer. 
  //=============================================================================
  double A = OmegaR + OmegaNu;
  double B = OmegaB + OmegaCDM;
  double dHdx = 0.5*pow(H0, 2) / H_of_x(x) * (-4*A*exp(-4*x) - 3*B*exp(-3*x) - 2*OmegaK*exp(-2*x));
  double first_term = 0.5*pow(H0 / H_of_x(x), 2) * (OmegaLambda*exp(x) + 6*A*exp(-3*x) + 2*B*exp(-2*x));
  double second_term = pow(H0, 2) / pow(H_of_x(x), 3) * dHdx*(OmegaLambda*exp(x) - 2*A*exp(-3*x) - B*exp(-2*x));
  double ddHpddx = first_term - second_term; 

  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if (x == 0.0) return OmegaB;
  else return OmegaB * pow(H0, 2) / (exp(3*x) * pow(H_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if (x == 0.0) return OmegaR;
  else return OmegaR * pow(H0, 2) / (exp(4*x) * pow(H_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if (x == 0.0) return OmegaNu;
  else return OmegaNu * pow(H0, 2) / (exp(4*x) * pow(H_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if (x == 0.0) return OmegaCDM;
  else return OmegaCDM * pow(H0, 2) / (exp(3*x) * pow(H_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if (x == 0.0) return OmegaLambda;
  else return OmegaLambda * pow(H0, 2) / (pow(H_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if (x == 0.0) return OmegaK;
  else return OmegaK * pow(H0, 2) / (exp(2*x) * pow(H_of_x(x), 2));
}
    
double BackgroundCosmology::get_r_of_x(double x) const{
  double Chi = get_comoving_distance_of_x(x);
  double K = H0*Chi / Constants.c;
  if (abs(OmegaK) < 1e-5) return Chi;
  
  if (OmegaK < 0) return Chi * sin(sqrt(abs(OmegaK))*K) / (sqrt(OmegaK)*K);

  else return Chi * sinh(sqrt(abs(OmegaK))*K) / (sqrt(OmegaK)*K);
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // Symbol is d_L.
  //=============================================================================
  return get_r_of_x(x) * exp(-x);
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // Symbol is Chi.
  //=============================================================================
  return eta_of_x(0) - eta_of_x(x);
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  //=============================================================================
  // Symbol is d_A.
  //=============================================================================
  return get_r_of_x(x) * exp(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::detadx_of_x(double x) const{
  return eta_of_x_spline.deriv_x(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_z_of_x(double x) const{ 
  return exp(-x) - 1; 
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
  // std::cout << "t0:          " << t_0         << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -20.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << detadx_of_x(x)     << " ";
    fp << t_of_x(x)          << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << get_luminosity_distance_of_x(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

