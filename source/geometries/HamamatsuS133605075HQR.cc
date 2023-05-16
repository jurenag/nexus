#include "HamamatsuS133605075HQR.h"
#include "OpticalMaterialProperties.h"

#include <G4SystemOfUnits.hh>

#include <utility>

namespace nexus {

  HamamatsuS133605075HQR::HamamatsuS133605075HQR():
  SiPMMPPC()
  {
  }

  HamamatsuS133605075HQR::~HamamatsuS133605075HQR()
  {
  }

  std::pair<G4int, G4double*> HamamatsuS133605075HQR::GetSensareaEnergyArray(){
    G4double* energy = new double[this->GetNPoints()]{opticalprops::optPhotMinE_, 
                                    h_Planck*c_light/(994.74*nm), h_Planck*c_light/(978.11*nm), h_Planck*c_light/(961.47*nm), 
                                    h_Planck*c_light/(944.84*nm), h_Planck*c_light/(928.2*nm), h_Planck*c_light/(911.57*nm), 
                                    h_Planck*c_light/(894.94*nm), h_Planck*c_light/(878.3*nm), h_Planck*c_light/(861.67*nm), 
                                    h_Planck*c_light/(845.03*nm), h_Planck*c_light/(828.4*nm), h_Planck*c_light/(811.76*nm), 
                                    h_Planck*c_light/(795.13*nm), h_Planck*c_light/(779.25*nm), h_Planck*c_light/(763.37*nm), 
                                    h_Planck*c_light/(747.49*nm), h_Planck*c_light/(732.37*nm), h_Planck*c_light/(718.0*nm), 
                                    h_Planck*c_light/(705.15*nm), h_Planck*c_light/(693.05*nm), h_Planck*c_light/(680.95*nm), 
                                    h_Planck*c_light/(670.37*nm), h_Planck*c_light/(660.54*nm), h_Planck*c_light/(649.95*nm), 
                                    h_Planck*c_light/(640.12*nm), h_Planck*c_light/(631.8*nm), h_Planck*c_light/(624.24*nm), 
                                    h_Planck*c_light/(615.93*nm), h_Planck*c_light/(606.85*nm), h_Planck*c_light/(597.78*nm), 
                                    h_Planck*c_light/(588.71*nm), h_Planck*c_light/(579.63*nm), h_Planck*c_light/(571.32*nm), 
                                    h_Planck*c_light/(563.0*nm), h_Planck*c_light/(554.68*nm), h_Planck*c_light/(546.36*nm), 
                                    h_Planck*c_light/(537.29*nm), h_Planck*c_light/(528.22*nm), h_Planck*c_light/(519.14*nm), 
                                    h_Planck*c_light/(508.56*nm), h_Planck*c_light/(495.7*nm), h_Planck*c_light/(480.58*nm), 
                                    h_Planck*c_light/(463.95*nm), h_Planck*c_light/(447.31*nm), h_Planck*c_light/(430.68*nm), 
                                    h_Planck*c_light/(416.31*nm), h_Planck*c_light/(405.72*nm), h_Planck*c_light/(397.41*nm), 
                                    h_Planck*c_light/(390.6*nm), h_Planck*c_light/(386.22*nm), h_Planck*c_light/(381.68*nm), 
                                    h_Planck*c_light/(377.14*nm), h_Planck*c_light/(370.94*nm), h_Planck*c_light/(359.6*nm), 
                                    h_Planck*c_light/(342.97*nm), h_Planck*c_light/(330.87*nm), h_Planck*c_light/(323.31*nm), 
                                    h_Planck*c_light/(315.75*nm), h_Planck*c_light/(308.94*nm), h_Planck*c_light/(301.38*nm), 
                                    h_Planck*c_light/(298.35*nm), h_Planck*c_light/(291.85*nm), h_Planck*c_light/(288.22*nm), 
                                    h_Planck*c_light/(286.58*nm), h_Planck*c_light/(285.25*nm), opticalprops::optPhotMaxE_};
    return std::make_pair(this->GetNPoints(), energy);
  }
 
  std::pair<G4int, G4double*> HamamatsuS133605075HQR::GetSensareaEfficiencyArray(){
    G4double* efficiency = new double[this->GetNPoints()]{0.000, 0.0108, 0.0143, 0.0186, 0.0232, 0.0286, 0.0349, 0.0422, 0.0507, 0.0599, 0.0696, 
                                        0.0795, 0.0905, 0.1026, 0.1158, 0.129, 0.1423, 0.1566, 0.1706, 0.184, 0.1973, 0.2112, 
                                        0.2256, 0.2393, 0.2526, 0.266, 0.2807, 0.2967, 0.3112, 0.3238, 0.3365, 0.3501, 0.3639, 
                                        0.3778, 0.3915, 0.4059, 0.4204, 0.435, 0.4492, 0.463, 0.476, 0.4891, 0.4998, 0.5028, 
                                        0.4997, 0.4911, 0.4784, 0.4634, 0.4494, 0.4352, 0.4201, 0.4048, 0.3909, 0.3785, 0.3618, 
                                        0.3547, 0.3388, 0.325, 0.3108, 0.297, 0.2687, 0.2548, 0.2153, 0.1776, 0.1548, 0.1274, 
                                        0.000};
    return std::make_pair(this->GetNPoints(), efficiency);
  }

  std::pair<G4int, G4double*> HamamatsuS133605075HQR::GetSensareaReflectivityArray(){
    G4double* reflectivity = new double[this->GetNPoints()]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    return std::make_pair(this->GetNPoints(), reflectivity);
  }    
}