#ifndef PERFECT_SIPM_MPPC_H
#define PERFECT_SIPM_MPPC_H

#include "SiPMMPPC.h"

namespace nexus {

  class PerfectSiPMMPPC: public SiPMMPPC
  {
  public:
    // Constructor for a SiPMMPPC with the dimensions of Hamamatsu S13360 5075 HQR 
    // MPPC but with 100% detection efficiency from optPhotMinE_ to optPhotMaxE_
    PerfectSiPMMPPC();

    // Destructor
    ~PerfectSiPMMPPC();

    G4double const GetTransverseDim()         override  { return 6.4*mm;          }   
    G4double const GetSensareaTransverseDim() override  { return 6.0*mm;          }
    G4double const GetThickness()             override  { return 1.5*mm;          }
    G4double const GetWindowThickness()       override  { return 0.15*mm;         }
    G4double const GetSensareaThickness()     override  { return 0.1*mm;          }
    G4String const GetModel()                 override  { return "PERFECTSIPMMPPC";}

  private:
    G4int GetNPoints() override { return 6;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif