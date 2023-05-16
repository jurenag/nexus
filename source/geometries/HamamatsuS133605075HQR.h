#ifndef HAMAMATSU_S13360_5075_HQR_H
#define HAMAMATSU_S13360_5075_HQR_H

#include "SiPMMPPC.h"

namespace nexus {

  class HamamatsuS133605075HQR: public SiPMMPPC
  {
  public:
    // Constructor for a Hamamatsu S13360 5075 HQR MPPC
    HamamatsuS133605075HQR();

    // Destructor
    ~HamamatsuS133605075HQR();

    G4double const GetTransverseDim()         override  { return 6.4*mm;          }   
    G4double const GetSensareaTransverseDim() override  { return 6.0*mm;          }
    G4double const GetThickness()             override  { return 1.5*mm;          } // This one differs from HamamatsuS133606050VE
    G4double const GetWindowThickness()       override  { return 0.15*mm;         } // This one differs from HamamatsuS133606050VE  
    G4double const GetSensareaThickness()     override  { return 0.1*mm;          }
    G4String const GetModel()                 override  { return "HS133605075HQR";}

  private:
    G4int GetNPoints() override { return 67;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif