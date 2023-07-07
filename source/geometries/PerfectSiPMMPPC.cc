#include "PerfectSiPMMPPC.h"
#include "OpticalMaterialProperties.h"

#include <G4SystemOfUnits.hh>

#include <utility>

namespace nexus {

  PerfectSiPMMPPC::PerfectSiPMMPPC():
  SiPMMPPC()
  {
  }

  PerfectSiPMMPPC::~PerfectSiPMMPPC()
  {
  }

  std::pair<G4int, G4double*> PerfectSiPMMPPC::GetSensareaEnergyArray(){
    G4double* energy = new double[this->GetNPoints()]{opticalprops::optPhotMinE_, 
                                                      1.2*eV, 2.2*eV, 9.5*eV, 10.5*eV, 
                                                      opticalprops::optPhotMaxE_};
    return std::make_pair(this->GetNPoints(), energy);
  }
 
  std::pair<G4int, G4double*> PerfectSiPMMPPC::GetSensareaEfficiencyArray(){
    G4double* efficiency = new double[this->GetNPoints()]{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    return std::make_pair(this->GetNPoints(), efficiency);
  }

  std::pair<G4int, G4double*> PerfectSiPMMPPC::GetSensareaReflectivityArray(){
    G4double* reflectivity = new double[this->GetNPoints()]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    return std::make_pair(this->GetNPoints(), reflectivity);
  }    
}