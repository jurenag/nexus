#ifndef FBK_NUV_HD_CRYO_TT_H
#define FBK_NUV_HD_CRYO_TT_H

#include "SiPMMPPC.h"

namespace nexus {

  class FbkNuvHdCryoTT: public SiPMMPPC
  {
  public:
    // Constructor for a FBK NUV HD CRYO TRIPLE TRENCH MPPC.
    FbkNuvHdCryoTT();

    // Destructor
    ~FbkNuvHdCryoTT();

    G4double GetTransverseDim()           override { return 7.35  *mm       ; }
    G4int GetSiPMsNoAlongX()              override { return 1               ; }
    G4int GetSiPMsNoAlongZ()              override { return 1               ; }
    G4double GetOuterFrameWidthAlongX()   override { return 0.625 *mm     ; }
    G4double GetOuterFrameWidthAlongZ()   override { return 0.625 *mm     ; }
    G4double GetInnerFramesWidthAlongX()  override { return 0.0             ; }
    G4double GetInnerFramesWidthAlongZ()  override { return 0.0             ; }
    G4double GetThickness()               override { return 1.8   *mm       ; }
    G4double GetWindowThickness()         override { return 0.6   *mm       ; }
    G4double GetSensareaThickness()       override { return 0.6   *mm       ; }
    G4String GetModel()                   override { return "FbkNuvHdCryoTT"; }

  private:
    G4int GetNPoints() override { return 3;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif