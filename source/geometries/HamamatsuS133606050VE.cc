#include "HamamatsuS133606050VE.h"
#include "OpticalMaterialProperties.h"

#include <G4SystemOfUnits.hh>

#include <utility>

namespace nexus {

  HamamatsuS133606050VE::HamamatsuS133606050VE():
  SiPMMPPC()
  {
  }

  HamamatsuS133606050VE::~HamamatsuS133606050VE()
  {
  }

  std::pair<G4int, G4double*> HamamatsuS133606050VE::GetSensareaEnergyArray(){
    G4double* energy = new double[this->GetNPoints()]{opticalprops::optPhotMinE_, 1.384*eV, 1.4003*eV, 1.4255*eV, 1.4516*eV, 
                                    1.4787*eV, 1.5068*eV, 1.536*eV, 1.5663*eV, 1.5979*eV, 1.6279*eV, 
                                    1.6563*eV, 1.6856*eV, 1.716*eV, 1.7475*eV, 1.7802*eV, 1.8107*eV, 
                                    1.8387*eV, 1.8676*eV, 1.8974*eV, 1.9281*eV, 1.9559*eV, 1.9803*eV, 
                                    2.0053*eV, 2.031*eV, 2.0574*eV, 2.0844*eV, 2.1122*eV, 2.1407*eV, 
                                    2.17*eV, 2.2001*eV, 2.2311*eV, 2.2629*eV, 2.2956*eV, 2.3293*eV, 
                                    2.3641*eV, 2.4059*eV, 2.462*eV, 2.5341*eV, 2.6177*eV, 2.7071*eV, 
                                    2.8028*eV, 2.8966*eV, 2.969*eV, 3.0256*eV, 3.0845*eV, 3.1353*eV, 
                                    3.1771*eV, 3.2201*eV, 3.2644*eV, 3.3098*eV, 3.3565*eV, 3.3924*eV, 
                                    3.4291*eV, 3.4792*eV, 3.5049*eV, 3.5308*eV, 3.5572*eV, 3.5841*eV, 
                                    3.5976*eV, 3.6112*eV, 3.625*eV, 3.6389*eV, 3.6528*eV, 3.6669*eV, 
                                    3.6907*eV, 3.7195*eV, 3.7488*eV, 3.7836*eV, 3.8138*eV, 3.8447*eV, 
                                    opticalprops::optPhotMaxE_};
    return std::make_pair(this->GetNPoints(), energy);
  }

  std::pair<G4int, G4double*> HamamatsuS133606050VE::GetSensareaEfficiencyArray(){
    G4double* efficiency   = new double[this->GetNPoints()]{0.0, 0.0367, 0.0408, 0.0474, 0.0544, 0.0623, 0.07, 0.0783, 0.0871, 
                                          0.0961, 0.1056, 0.1151, 0.1248, 0.1348, 0.145, 0.1554, 0.166, 0.1758, 
                                          0.1855, 0.1961, 0.2073, 0.2197, 0.2321, 0.2424, 0.2519, 0.2606, 0.2709, 
                                          0.2809, 0.2912, 0.3018, 0.3124, 0.3224, 0.3333, 0.3433, 0.3536, 0.3628, 
                                          0.3738, 0.3849, 0.393, 0.399, 0.4017, 0.3951, 0.3856, 0.3758, 0.3661, 
                                          0.3563, 0.3465, 0.3352, 0.3217, 0.3063, 0.2904, 0.2734, 0.2612, 0.2498, 
                                          0.2319, 0.2201, 0.2103, 0.1993, 0.1859, 0.1741, 0.1656, 0.1546, 0.1461, 
                                          0.1351, 0.127, 0.1107, 0.0928, 0.0752, 0.0586, 0.044, 0.0318, 0.0};
    return std::make_pair(this->GetNPoints(), efficiency);
  }

  std::pair<G4int, G4double*> HamamatsuS133606050VE::GetSensareaReflectivityArray(){
    G4double* reflectivity = new double[this->GetNPoints()]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    return std::make_pair(this->GetNPoints(), reflectivity);
  }
}