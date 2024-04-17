#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  // Control what should run 
  bool output = true;
  bool supernovafit = false; 

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);

  Utils::StartTiming("Solve");
  cosmo.solve();
  cosmo.info();
  Utils::EndTiming("Solve");
  
  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

  // Output background evolution quantities
  if (output){
    Utils::StartTiming("Output");
    cosmo.output("data/cosmology.txt");

    // Output best fit values
    // Minimum chi^2 found chi^2 = 29.2811 h = 0.70189 OmegaM = 0.25932 (OmegaB + OmegaCDM) OmegaK = 0.0673887
    BackgroundCosmology bestFit(0.702, OmegaB, 0.259 - OmegaB, 0.067, Neff, TCMB);
    // Utils::StartTiming("Solve best params");
    // bestFit.solve();
    // bestFit.info();
    // bestFit.output("data/bestFitBackground.txt");
    // Utils::EndTiming("Solve best params");

    Utils::EndTiming("Output");
  }

  if (supernovafit){
    Utils::StartTiming("SupernovaFit");
    mcmc_fit_to_supernova_data("data/supernovadata.txt", "data/results_supernovafitting.txt");
    Utils::EndTiming("SupernovaFit");
  }

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  Utils::StartTiming("Solve");
  rec.solve();
  rec.info();
  Utils::EndTiming("Solve");

  // Output recombination quantities
  if (output){
    Utils::StartTiming("Output");
    rec.output("data/recombination.txt");
    Utils::EndTiming("Output");
  }

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations.
  // polarization, neutrinos = false in Utils.h
  // Set Neff = 0 since we have no neutrinos, Yp = 0 since no Helium.
  Neff = 0;
  Yp = 0;
  BackgroundCosmology cosmo_no_neutrino(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo_no_neutrino.solve();

  RecombinationHistory rec_no_neutrino(&cosmo_no_neutrino, Yp);
  rec_no_neutrino.solve();

  Perturbations pert(&cosmo_no_neutrino, &rec_no_neutrino);
  Utils::StartTiming("Solve");
  pert.solve();
  pert.info();
  Utils::EndTiming("Solve");
  
  // Output perturbation quantities
  if (output){
    Utils::StartTiming("Output");
    double kvalue1 = 0.001 / Constants.Mpc;
    pert.output(kvalue1, "data/perturbations_k0.001.txt");

    double kvalue2 = 0.01 / Constants.Mpc;
    pert.output(kvalue2, "data/perturbations_k0.01.txt");
    
    double kvalue3 = 0.1 / Constants.Mpc;
    pert.output(kvalue3, "data/perturbations_k0.1.txt");
    Utils::EndTiming("Output");
  }
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
