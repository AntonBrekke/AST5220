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
  // Can make array with any spaving f(x) by f^-1(linspace(f(x0), f(x1), npts)). 
  // Idea is to make a linspace along f-axis, and project it down on x-axis.
  // I choose log since I like it better than square.
  //===================================================================
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);

  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Declare vectors storing values we are going to calculate
  Vector Phi_array(n_x*n_k);
  Vector Psi_array(n_x*n_k);
  Vector delta_cdm_array(n_x*n_k);
  Vector delta_b_array(n_x*n_k);
  Vector v_cdm_array(n_x*n_k);
  Vector v_b_array(n_x*n_k);
  std::vector<Vector> Theta_array(Constants.n_ell_theta, Vector(n_x*n_k));

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik)/n_k != (10*ik + 10)/n_k ) {
      std::cout << (100*ik + 100)/n_k << "% " << std::flush;
      if(ik == n_k - 1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to. I switch from time to index since I find this more intuitive
    int end_tc_index = get_tight_coupling_index(k, x_array);
    double x_end_tc = x_array[end_tc_index];

    // Find number of points in each regime 
    int npts_tc_regime = end_tc_index + 1;
    int npts_full_regime = n_x - npts_tc_regime;
    int first_full_regime_index = npts_tc_regime;

    // Make x_array for tight coupling solution
    Vector x_tc_array(npts_tc_regime);
    for(int i=0; i < npts_tc_regime; i++){
      x_tc_array[i] = x_array[i];
    }

    // Make x_array for full solution
    Vector x_full(npts_full_regime + 1);
    for (int i=0; i <= npts_full_regime; i++){
      x_full[i] = x_array[(npts_tc_regime - 1) + i];
    }

    // The tight coupling ODE system
    ODEFunction dydx_tc = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    ODESolver ode_tight_coupling;
      auto y_tc_init = set_ic(x_start, k);
      ode_tight_coupling.solve(dydx_tc, x_tc_array, y_tc_init);

      auto solution_tight_coupling = ode_tight_coupling.get_data();
      auto derivatives_tight_coupling = ode_tight_coupling.get_derivative_data();

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    ODESolver ode_full;
      Vector y_full_init = solution_tight_coupling[end_tc_index];
      auto y_full_init_vec = set_ic_after_tight_coupling(y_full_init, x_end_tc, k);
      ode_full.solve(dydx_full, x_full, y_full_init_vec);

      auto solution_full = ode_full.get_data();
      auto derivative_full = ode_full.get_derivative_data();

    // Start filling arrays form the tight coupling regime 
    for (int ix=0; ix < npts_tc_regime; ix++){
      int index = ix + n_x*ik;
      auto y_tc = solution_tight_coupling[ix];

      // Get quantities needed for computation
      double x = x_tc_array[ix];
      double c = Constants.c;
      double Hp = cosmo -> Hp_of_x(x);
      double OmegaR = cosmo -> get_OmegaR(x);
      double ckHp = c*k/Hp;
      double dtau = rec -> dtaudx_of_x(x);

      // Copy from below (line "References to the tight coupling quantities")
      double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
      double &delta_b      =  y_tc[Constants.ind_deltab_tc];
      double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
      double &v_b          =  y_tc[Constants.ind_vb_tc];
      double &Phi          =  y_tc[Constants.ind_Phi_tc];
      double *Theta        = &y_tc[Constants.ind_start_theta_tc];

      // Calculate quantities
      Phi_array[index] = Phi;
      Psi_array[index] = -Phi - 12.*pow(Hp/(c*k), 2) * OmegaR * Theta[2];
      delta_cdm_array[index] = delta_cdm;
      delta_b_array[index] = delta_b;
      v_cdm_array[index] = v_cdm;
      v_b_array[index] = v_b;
      Theta_array[0][index] = Theta[0];
      Theta_array[1][index] = Theta[1];
      Theta_array[2][index] = -20.*ckHp/(45.*dtau) * Theta[1];   // No polarization
      // Solve for rest of the l-values
      for(int l=3; l < Constants.n_ell_theta; l++){
        Theta_array[l][index] = -l/(2.*l+1.)*ckHp/dtau * Theta[l-1];
      }
    }

    // Now fill rest of arrays from the full regime 
    for (int ix=npts_tc_regime; ix < n_x; ix++){
      int index = ix + n_x*ik;
      auto y_full = solution_full[ix - npts_tc_regime];

      // Get quantities needed for computation
      double x = x_full[ix - npts_tc_regime];
      double c = Constants.c;
      double Hp = cosmo -> Hp_of_x(x);
      double ckHp = c*k/Hp;
      double OmegaR = cosmo -> get_OmegaR(x);
      double dtau = rec -> dtaudx_of_x(x);

      // Copy from below (line "References to the tight coupling quantities")
      double &delta_cdm    =  y_full[Constants.ind_deltacdm_tc];
      double &delta_b      =  y_full[Constants.ind_deltab_tc];
      double &v_cdm        =  y_full[Constants.ind_vcdm_tc];
      double &v_b          =  y_full[Constants.ind_vb_tc];
      double &Phi          =  y_full[Constants.ind_Phi_tc];
      double *Theta        = &y_full[Constants.ind_start_theta_tc];

      // Calculate quantities
      Phi_array[index] = Phi;
      Psi_array[index] = -Phi - 12.*pow(Hp/(c*k), 2) * OmegaR * Theta[2];
      delta_cdm_array[index] = delta_cdm;
      delta_b_array[index] = delta_b;
      v_cdm_array[index] = v_cdm;
      v_b_array[index] = v_b;
      // Solve for all of the l-values
      for(int l=0; l < Constants.n_ell_theta; l++){
        Theta_array[l][index] = Theta[l];
      }
    }
  }
  Utils::EndTiming("integrateperturbation");

  // Spline results
  Phi_spline.create(x_array, k_array, Phi_array);
  Psi_spline.create(x_array, k_array, Psi_array);
  delta_cdm_spline.create(x_array, k_array, delta_cdm_array);
  delta_b_spline.create(x_array, k_array, delta_b_array);
  v_cdm_spline.create(x_array, k_array, v_cdm_array);
  v_b_spline.create(x_array, k_array, v_b_array);
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for(int l=0; l < Constants.n_ell_theta; l++){
    Theta_spline[l].create(x_array, k_array, Theta_array[l]);
  }
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  // Initial conditions
  double Psi  = -2./3.;
  double Hp   = cosmo -> Hp_of_x(x);
  double ckHp = Constants.c*k/Hp;

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  Phi       = -Psi;
  delta_cdm = -3./2.*Psi;
  delta_b   = -3./2.*Psi;
  v_cdm     = -0.5*ckHp*Psi;
  v_b       = -0.5*ckHp*Psi;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = -0.5*Psi;
  Theta[1] = ckHp/6.*Psi;

  // No neutrinos

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

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
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  double Hp = cosmo -> Hp_of_x(x);
  double dtau = rec -> dtaudx_of_x(x);
  double ckHp = Constants.c*k/Hp;

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = -20.*ckHp/(45.*dtau) * Theta[1];
  for(int l=3; l < n_ell_theta; l++){
    Theta[l] = -l/(2*l+1)*ckHp/dtau * Theta[l-1];
  }

  // No polarizations or neutrinos

  return y;
}

//====================================================
// The index when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_index(const double k, Vector x_array) const{
    // Find index 
    int tight_coupling_index = 0;
    bool tight_coupling = true;
    double x_current;
    double ckHp;
    double dtaudx;
    while(tight_coupling){
      x_current = x_array[tight_coupling_index];
      ckHp = Constants.c*k / cosmo -> Hp_of_x(x_current);
      dtaudx = abs(rec -> dtaudx_of_x(x_current));
      // Check conditions from Callin
      if(dtaudx < 10){tight_coupling = false;}
      else if(x_current > -8.3){tight_coupling = false;}
      else if(dtaudx < 10*ckHp){tight_coupling = false;}
      else{tight_coupling_index = tight_coupling_index + 1;}
    }

    return tight_coupling_index;
  }

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  // Setting up k-array and x-array
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x*ik;

      // Fetch functions from BackgroundCosmology
      const double Hp   = cosmo -> Hp_of_x(x);
      const double dHp  = cosmo -> dHpdx_of_x(x);
      const double ddHp = cosmo -> ddHpddx_of_x(x);

      // Fetch functions from Recombination
      const double tau   = rec -> tau_of_x(x);
      const double dtau  = rec -> dtaudx_of_x(x);
      const double ddtau = rec -> ddtauddx_of_x(x);
      const double g     = rec -> g_tilde_of_x(x);
      const double dg    = rec -> dgdx_tilde_of_x(x);
      const double ddg   = rec -> ddgddx_tilde_of_x(x);

      // Get functions from Perturbation
      double T0  = get_Theta(x,k,0);
      double T1  = get_Theta(x,k,1);
      double T2  = get_Theta(x,k,2);
      double T3  = get_Theta(x,k,3);
      double v_b = get_v_b(x,k);
      double Psi = get_Psi(x,k);

      double dPsi = get_dPsi(x,k);
      double dPhi = get_dPhi(x,k);
      double dv_b = get_dv_b(x,k);
      double dT1  = get_dTheta(x,k,1);
      double dT2  = get_dTheta(x,k,2);
      double dT3  = get_dTheta(x,k,3);

      double c = Constants.c;
      double ckHp = c*k/Hp;

      double ddT2 = 1./5.*ckHp*(2.*dT1 - 3.*dT3) - ckHp*dHp/(5.*Hp)*(T1 - T3) + 9./10.*(ddtau*T2 + dtau*dT2);
      
      // Temperatur source. Remember: k -> ik.
      double term_1 = Hp/c * (g*(T0 + 1./4.*T2) - dPhi*exp(-tau));
      double term_2 = Hp/c * (dPsi*exp(-tau) + Psi*g + 1./(c*k)*(dHp*v_b*g + Hp*dv_b*g + Hp*v_b*dg));
      double term_3 = -3./(4.*pow(k, 2)) * dHp*Hp/pow(c, 3) * (dHp*g*T2 + Hp*dg*T2 + Hp*g*dT2);
      double term_4 = -3./(4.*pow(k, 2)) * pow(Hp, 2)/pow(c, 3) * (ddHp*g*T2 + Hp*ddg*T2 + Hp*g*ddT2 + 2.*(dHp*dg*T2 + dHp*g*dT2 + Hp*dg*dT2));
      ST_array[index] = term_1 + term_2 + term_3 + term_4;
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

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

  //=============================================================================
  // Fill in the expressions for all the derivatives
  //=============================================================================

  // Fetch functions from BackgroundCosmology
  double H0       = cosmo -> get_H0();
  double Hp       = cosmo -> Hp_of_x(x);
  double dHp      = cosmo -> dHpdx_of_x(x);
  double OmegaCDM = cosmo -> get_OmegaCDM(x);
  double OmegaB   = cosmo -> get_OmegaB(x);
  double OmegaR   = cosmo -> get_OmegaR(x);

  // Fetch functions from Recombination
  double R     = rec -> R_of_x(x);
  double dtau  = rec -> dtaudx_of_x(x);
  double ddtau = rec -> ddtauddx_of_x(x);

  double ckHp   = Constants.c*k/Hp;
  double Theta2 = -20.*ckHp/(45.*dtau) * Theta[1];
  double Psi    = -Phi - 12.*pow(Hp/(Constants.c*k), 2) * OmegaR * Theta2;

  // SET: Scalar quantities (Phi, delta, v, ...). Notice that OmegaCDM contains OmegaCDM0*a^{-1}*(H0/H)**2 etc. No neutrinos.
  dPhidx       = Psi - 1./3. * pow(ckHp, 2)*Phi + 0.5*(OmegaCDM*delta_cdm + OmegaB*delta_b + 4.*OmegaR * Theta[0]);
  ddelta_cdmdx = ckHp*v_cdm - 3.*dPhidx;
  ddelta_bdx   = ckHp*v_b - 3.*dPhidx;
  dv_cdmdx     = -v_cdm - ckHp*Psi;

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0]  = -ckHp * Theta[1] - dPhidx;
  double num = -((1. - R)*dtau + (1. + R)*ddtau)*(3.*Theta[1] + v_b) - ckHp*Psi + (1 - dHp/Hp)*ckHp*(2*Theta2 - Theta[0]) - ckHp*dThetadx[0];
  double denom = (1. + R)*dtau + dHp/Hp - 1;
  double q = num / denom;
  dThetadx[1] = 1./3.*(q - dv_bdx);

  // No neutrinos 

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

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

  // Fetch functions from BackgroundCosmology
  double H0       = cosmo -> get_H0();
  double Hp       = cosmo -> Hp_of_x(x);
  double dHp      = cosmo -> dHpdx_of_x(x);
  double OmegaCDM = cosmo -> get_OmegaCDM(x);
  double OmegaB   = cosmo -> get_OmegaB(x);
  double OmegaR   = cosmo -> get_OmegaR(x);

  // Fetch functions from Recombination
  double R     = rec -> R_of_x(x);
  double dtau  = rec -> dtaudx_of_x(x);
  double ddtau = rec -> ddtauddx_of_x(x);

  double ckHp   = Constants.c*k/Hp;
  double Psi    = -Phi - 12.*pow(Hp/(Constants.c*k), 2) * OmegaR * Theta[2];

  //=============================================================================
  // Fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx       = Psi - 1./3. * pow(ckHp, 2)*Phi + 0.5*(OmegaCDM*delta_cdm + OmegaB*delta_b + 4.*OmegaR * Theta[0]);
  ddelta_cdmdx = ckHp*v_cdm - 3.*dPhidx;
  ddelta_bdx   = ckHp*v_b - 3.*dPhidx;
  dv_cdmdx     = -v_cdm - ckHp*Psi;
  dv_bdx       = -v_b - ckHp*Psi + dtau*R*(3*Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -ckHp*Theta[1] - dPhidx;
  dThetadx[1] = ckHp/3.*Theta[0] - 2./3.*ckHp*Theta[2] + ckHp/3.*Psi + dtau*(Theta[1] + 1./3.*v_b);
  dThetadx[2] = 2./5.*ckHp*Theta[1] - 3./5.*ckHp + 9./10.*dtau*Theta[2];

  // Solve for 2 < l < lmax
  for(int l=3; l < n_ell_theta-1; l++){
      dThetadx[l] = l/(2.*l+1.)*ckHp*Theta[l-1] - (l+1.)/(2.*l+1.)*ckHp*Theta[l+1] + dtau*Theta[l];
    }
  // Solve for l=lmax
  int lmax = n_ell_theta-1;
  dThetadx[lmax] = ckHp*Theta[lmax-1] - Constants.c*(lmax+1.)*Theta[lmax]/(Hp*cosmo->eta_of_x(x)) + dtau*Theta[lmax];

  // No neutrinos or polarization

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

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

double Perturbations::get_dv_b(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}

double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}

double Perturbations::get_dPhi(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}

double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}

double Perturbations::get_dPsi(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}

double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
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

double Perturbations::get_dTheta(const double x, const double k, const int ell) const{
  return Theta_spline[ell].deriv_x(x,k);
}

double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}

double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}


//====================================================
// Print some useful info about the class
//====================================================

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
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

