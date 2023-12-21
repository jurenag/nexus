#ifndef BROADCOM_AFBR_S4N44P044M_H
#define BROADCOM_AFBR_S4N44P044M_H

#include "SiPMMPPC.h"

namespace nexus {

  class BroadcomAFBRS4N44P044M: public SiPMMPPC
  {
  public:
    // Constructor for a Broadcom AFBR S4N44P044M
    // I got the specifications for this SiPMMPPC from the data sheet which can be found in
    // https://www.broadcom.com/products/optical-sensors/silicon-photomultiplier-sipm/afbr-s4n44p044m
    BroadcomAFBRS4N44P044M();

    // Destructor
    ~BroadcomAFBRS4N44P044M();

    G4double GetTransverseDim()           override { return 8.26  *mm       ; }
    G4int GetSiPMsNoAlongX()              override { return 2               ; }
    G4int GetSiPMsNoAlongZ()              override { return 2               ; }
    G4double GetOuterFrameWidthAlongX()   override { return 0.19  *mm       ; }
    G4double GetOuterFrameWidthAlongZ()   override { return 0.19  *mm       ; }
    G4double GetInnerFramesWidthAlongX()  override { return 0.2   *mm       ; }
    G4double GetInnerFramesWidthAlongZ()  override { return 0.4   *mm       ; }
    G4double GetThickness()               override { return 1.23  *mm       ; }
    G4double GetWindowThickness()         override { return 0.20  *mm       ; }   ///< This one was not specified in the datasheet
    G4double GetSensareaThickness()       override { return 0.10  *mm       ; }   ///< The same goes for this one
    G4String GetModel()                   override { return "BroadcomAFBRS4N44P044M"; }


  private:
    G4int GetNPoints() override { return 68;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif