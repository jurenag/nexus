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
    G4double energy[]       = {opticalprops::optPhotMinE_, h_Planck*c_light/(285.25*nm), h_Planck*c_light/(286.58*nm), 
                              h_Planck*c_light/(288.22*nm), h_Planck*c_light/(291.85*nm), h_Planck*c_light/(298.35*nm), 
                              h_Planck*c_light/(301.38*nm), h_Planck*c_light/(308.94*nm), h_Planck*c_light/(315.75*nm), 
                              h_Planck*c_light/(323.31*nm), h_Planck*c_light/(330.87*nm), h_Planck*c_light/(342.97*nm), 
                              h_Planck*c_light/(359.6*nm), h_Planck*c_light/(370.94*nm), h_Planck*c_light/(377.14*nm), 
                              h_Planck*c_light/(381.68*nm), h_Planck*c_light/(386.22*nm), h_Planck*c_light/(390.6*nm), 
                              h_Planck*c_light/(397.41*nm), h_Planck*c_light/(405.72*nm), h_Planck*c_light/(416.31*nm), 
                              h_Planck*c_light/(430.68*nm), h_Planck*c_light/(447.31*nm), h_Planck*c_light/(463.95*nm), 
                              h_Planck*c_light/(480.58*nm), h_Planck*c_light/(495.7*nm), h_Planck*c_light/(508.56*nm), 
                              h_Planck*c_light/(519.14*nm), h_Planck*c_light/(528.22*nm), h_Planck*c_light/(537.29*nm), 
                              h_Planck*c_light/(546.36*nm), h_Planck*c_light/(554.68*nm), h_Planck*c_light/(563.0*nm), 
                              h_Planck*c_light/(571.32*nm), h_Planck*c_light/(579.63*nm), h_Planck*c_light/(588.71*nm), 
                              h_Planck*c_light/(597.78*nm), h_Planck*c_light/(606.85*nm), h_Planck*c_light/(615.93*nm), 
                              h_Planck*c_light/(624.24*nm), h_Planck*c_light/(631.8*nm), h_Planck*c_light/(640.12*nm), 
                              h_Planck*c_light/(649.95*nm), h_Planck*c_light/(660.54*nm), h_Planck*c_light/(670.37*nm), 
                              h_Planck*c_light/(680.95*nm), h_Planck*c_light/(693.05*nm), h_Planck*c_light/(705.15*nm), 
                              h_Planck*c_light/(718.0*nm), h_Planck*c_light/(732.37*nm), h_Planck*c_light/(747.49*nm), 
                              h_Planck*c_light/(763.37*nm), h_Planck*c_light/(779.25*nm), h_Planck*c_light/(795.13*nm), 
                              h_Planck*c_light/(811.76*nm), h_Planck*c_light/(828.4*nm), h_Planck*c_light/(845.03*nm), 
                              h_Planck*c_light/(861.67*nm), h_Planck*c_light/(878.3*nm), h_Planck*c_light/(894.94*nm), 
                              h_Planck*c_light/(911.57*nm), h_Planck*c_light/(928.2*nm), h_Planck*c_light/(944.84*nm), 
                              h_Planck*c_light/(961.47*nm), h_Planck*c_light/(978.11*nm), h_Planck*c_light/(994.74*nm),
                              opticalprops::optPhotMaxE_};
    return std::make_pair(67, energy);
  }
 
  std::pair<G4int, G4double*> HamamatsuS133605075HQR::GetSensareaEfficiencyArray(){
    G4double efficiency[]   = {0.000, 0.127, 0.155, 0.178, 0.215, 0.255, 0.269, 0.297, 0.311, 0.325, 0.339, 0.355, 
                              0.362, 0.379, 0.391, 0.405, 0.420, 0.435, 0.449, 0.463, 0.478, 0.491, 0.500, 0.503, 
                              0.500, 0.489, 0.476, 0.463, 0.449, 0.435, 0.420, 0.406, 0.392, 0.378, 0.364, 0.350, 
                              0.336, 0.324, 0.311, 0.297, 0.281, 0.266, 0.253, 0.239, 0.226, 0.211, 0.197, 0.184, 
                              0.171, 0.157, 0.142, 0.129, 0.116, 0.103, 0.090, 0.079, 0.070, 0.060, 0.051, 0.042, 
                              0.035, 0.029, 0.023, 0.019, 0.014, 0.011, 0.000};
    return std::make_pair(67, efficiency);
  }

  std::pair<G4int, G4double*> HamamatsuS133605075HQR::GetSensareaReflectivityArray(){
    G4double reflectivity[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    return std::make_pair(67, reflectivity);
  }    
}