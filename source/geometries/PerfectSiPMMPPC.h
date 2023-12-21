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

    G4double GetTransverseDim()           override { return 6.4   *mm         ; }
    G4int GetSiPMsNoAlongX()              override { return 1                 ; }
    G4int GetSiPMsNoAlongZ()              override { return 1                 ; }
    G4double GetOuterFrameWidthAlongX()   override { return 0.2   *mm         ; }
    G4double GetOuterFrameWidthAlongZ()   override { return 0.2   *mm         ; }
    G4double GetInnerFramesWidthAlongX()  override { return 0.0               ; }
    G4double GetInnerFramesWidthAlongZ()  override { return 0.0               ; }
    G4double GetThickness()               override { return 1.5   *mm         ; }
    G4double GetWindowThickness()         override { return 0.15  *mm         ; }
    G4double GetSensareaThickness()       override { return 0.1   *mm         ; }
    G4String GetModel()                   override { return "PERFECTSIPMMPPC" ; }

  private:
    G4int GetNPoints() override { return 6;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif