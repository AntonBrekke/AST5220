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
  bool supernovafit = true; 

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

  if (output){
    // Output background evolution quantities
    Utils::StartTiming("Output");
    cosmo.output("data/cosmology.txt");

    // Output best fit values
    // BackgroundCosmology bestFit(0.702, 0.05, 0.209, 0.067, Neff, TCMB);
    // Utils::StartTiming("Solve best params");
    // bestFit.solve();
    // bestFit.info();
    // bestFit.output("data/bestFitBackground.txt");
    // Utils::EndTiming("Solve best params");
    // Utils::EndTiming("Output background");

    Utils::EndTiming("Output");
  }

  if (supernovafit){
    Utils::StartTiming("SupernovaFit");
    mcmc_fit_to_supernova_data("data/supernovadata.txt", "data/results_supernovafitting.txt");
    Utils::EndTiming("SupernovaFit");
  }

  // return 0;

  

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");
  
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