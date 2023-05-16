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

    G4double const GetTransverseDim()         override  { return 7.35*mm;         }
    G4double const GetSensareaTransverseDim() override  { return 6.1*mm;          }
    G4double const GetThickness()             override  { return 1.8*mm;          }
    G4double const GetWindowThickness()       override  { return 0.6*mm;          }
    G4double const GetSensareaThickness()     override  { return 0.6*mm;          }
    G4String const GetModel()                 override  { return "FbkNuvHdCryoTT";}

  private:
    G4int GetNPoints() override { return 3;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif