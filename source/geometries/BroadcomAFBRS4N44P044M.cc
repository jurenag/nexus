#include "BroadcomAFBRS4N44P044M.h"
#include "OpticalMaterialProperties.h"

#include <G4SystemOfUnits.hh>

#include <utility>

namespace nexus {

  BroadcomAFBRS4N44P044M::BroadcomAFBRS4N44P044M():
  SiPMMPPC()
  {
  }

  BroadcomAFBRS4N44P044M::~BroadcomAFBRS4N44P044M()
  {
  }

  std::pair<G4int, G4double*> BroadcomAFBRS4N44P044M::GetSensareaEnergyArray(){
    G4double* energy = new double[this->GetNPoints()]{opticalprops::optPhotMinE_, 
                                                      h_Planck*c_light/(897.82*nm), h_Planck*c_light/(888.81*nm), h_Planck*c_light/(874.64*nm), 
                                                      h_Planck*c_light/(860.48*nm), h_Planck*c_light/(846.31*nm), h_Planck*c_light/(832.14*nm), 
                                                      h_Planck*c_light/(817.97*nm), h_Planck*c_light/(803.8*nm), h_Planck*c_light/(789.63*nm), 
                                                      h_Planck*c_light/(775.46*nm), h_Planck*c_light/(761.29*nm), h_Planck*c_light/(747.11*nm), 
                                                      h_Planck*c_light/(732.94*nm), h_Planck*c_light/(718.76*nm), h_Planck*c_light/(704.59*nm), 
                                                      h_Planck*c_light/(690.41*nm), h_Planck*c_light/(676.24*nm), h_Planck*c_light/(662.07*nm), 
                                                      h_Planck*c_light/(648.89*nm), h_Planck*c_light/(636.29*nm), h_Planck*c_light/(625.33*nm), 
                                                      h_Planck*c_light/(615.01*nm), h_Planck*c_light/(604.69*nm), h_Planck*c_light/(593.09*nm), 
                                                      h_Planck*c_light/(579.56*nm), h_Planck*c_light/(565.39*nm), h_Planck*c_light/(551.22*nm), 
                                                      h_Planck*c_light/(538.33*nm), h_Planck*c_light/(526.72*nm), h_Planck*c_light/(515.76*nm), 
                                                      h_Planck*c_light/(505.45*nm), h_Planck*c_light/(495.13*nm), h_Planck*c_light/(484.81*nm), 
                                                      h_Planck*c_light/(473.85*nm), h_Planck*c_light/(461.61*nm), h_Planck*c_light/(448.07*nm), 
                                                      h_Planck*c_light/(433.9*nm), h_Planck*c_light/(419.73*nm), h_Planck*c_light/(405.57*nm), 
                                                      h_Planck*c_light/(394.66*nm), h_Planck*c_light/(388.25*nm), h_Planck*c_light/(383.13*nm), 
                                                      h_Planck*c_light/(378.01*nm), h_Planck*c_light/(372.24*nm), h_Planck*c_light/(364.55*nm), 
                                                      h_Planck*c_light/(352.98*nm), h_Planck*c_light/(342.06*nm), h_Planck*c_light/(335.66*nm), 
                                                      h_Planck*c_light/(331.18*nm), h_Planck*c_light/(325.42*nm), h_Planck*c_light/(314.5*nm), 
                                                      h_Planck*c_light/(305.51*nm), h_Planck*c_light/(301.68*nm), h_Planck*c_light/(298.5*nm), 
                                                      h_Planck*c_light/(295.95*nm), h_Planck*c_light/(293.41*nm), h_Planck*c_light/(290.87*nm), 
                                                      h_Planck*c_light/(288.33*nm), h_Planck*c_light/(285.79*nm), h_Planck*c_light/(283.25*nm), 
                                                      h_Planck*c_light/(280.71*nm), h_Planck*c_light/(278.17*nm), h_Planck*c_light/(275.62*nm), 
                                                      h_Planck*c_light/(273.07*nm), h_Planck*c_light/(269.88*nm), h_Planck*c_light/(265.4*nm),
                                                      opticalprops::optPhotMaxE_};
    return std::make_pair(this->GetNPoints(), energy);
  }
 
  std::pair<G4int, G4double*> BroadcomAFBRS4N44P044M::GetSensareaEfficiencyArray(){
    G4double* efficiency = new double[this->GetNPoints()]{0.000, 
                                                          0.05268, 0.0585, 0.06811, 0.07816, 0.08875, 0.10002, 0.11197, 0.12475, 0.13844, 0.15273, 
                                                          0.16785, 0.18404, 0.20113, 0.21814, 0.23546, 0.25218, 0.26798, 0.28432, 0.30067, 0.31906, 
                                                          0.33759, 0.35631, 0.37471, 0.392, 0.40781, 0.4218, 0.43692, 0.4541, 0.47259, 0.49159, 
                                                          0.50979, 0.52809, 0.54587, 0.56335, 0.58199, 0.6002, 0.61608, 0.6281, 0.63022, 0.61433, 
                                                          0.59464, 0.57592, 0.55721, 0.53716, 0.51817, 0.50756, 0.49384, 0.47173, 0.45225, 0.43243, 
                                                          0.42188, 0.40178, 0.38043, 0.35465, 0.33302, 0.30931, 0.28352, 0.25981, 0.23236, 0.20782, 
                                                          0.18287, 0.15957, 0.13795, 0.12006, 0.10065, 0.08221,
                                                          0.000};
    return std::make_pair(this->GetNPoints(), efficiency);
  }

  std::pair<G4int, G4double*> BroadcomAFBRS4N44P044M::GetSensareaReflectivityArray(){
    G4double* reflectivity = new double[this->GetNPoints()]{0.0, 
                                                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                            0.0};
    return std::make_pair(this->GetNPoints(), reflectivity);
  }    
}