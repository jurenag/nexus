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

    G4double GetTransverseDim()           override { return 6.4   *mm       ; }
    G4int GetSiPMsNoAlongX()              override { return 1               ; }
    G4int GetSiPMsNoAlongZ()              override { return 1               ; }
    G4double GetOuterFrameWidthAlongX()   override { return 0.2   *mm       ; }
    G4double GetOuterFrameWidthAlongZ()   override { return 0.2   *mm       ; }
    G4double GetInnerFramesWidthAlongX()  override { return 0.0             ; }
    G4double GetInnerFramesWidthAlongZ()  override { return 0.0             ; }
    G4double GetThickness()               override { return 1.5   *mm       ; } ///< This one differs from HamamatsuS133606050VE
    G4double GetWindowThickness()         override { return 0.15  *mm       ; } ///< This one differs from HamamatsuS133606050VE
    G4double GetSensareaThickness()       override { return 0.1   *mm       ; }
    G4String GetModel()                   override { return "HS133605075HQR"; } ///< This one differs from HamamatsuS133606050VE


  private:
    G4int GetNPoints() override { return 67;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif