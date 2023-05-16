#include "FbkNuvHdCryoTT.h"
#include "OpticalMaterialProperties.h"

#include <G4SystemOfUnits.hh>

#include <utility>

namespace nexus {

  FbkNuvHdCryoTT::FbkNuvHdCryoTT():
  SiPMMPPC()
  {
  }

  FbkNuvHdCryoTT::~FbkNuvHdCryoTT()
  {
  }
 
  std::pair<G4int, G4double*> FbkNuvHdCryoTT::GetSensareaEnergyArray(){
    G4double* energy = new double[this->GetNPoints()]{opticalprops::optPhotMinE_, h_Planck*c_light/(435*nm),
                                    opticalprops::optPhotMaxE_};
    return std::make_pair(this->GetNPoints(), energy);
  }

  std::pair<G4int, G4double*> FbkNuvHdCryoTT::GetSensareaEfficiencyArray(){
    G4double* efficiency = new double[this->GetNPoints()]{0.000, 0.425, 0.000}; // This is he only 
                                                              // data we've got so far
    return std::make_pair(this->GetNPoints(), efficiency);
  }

  std::pair<G4int, G4double*> FbkNuvHdCryoTT::GetSensareaReflectivityArray(){
    G4double* reflectivity = new double[this->GetNPoints()]{0.0, 0.0, 0.0};
    return std::make_pair(this->GetNPoints(), reflectivity);
  }
}