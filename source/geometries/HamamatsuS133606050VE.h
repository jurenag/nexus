#ifndef HAMAMATSU_S13360_6050VE_H
#define HAMAMATSU_S13360_6050VE_H

#include "SiPMMPPC.h"

namespace nexus {

  class HamamatsuS133606050VE: public SiPMMPPC
  {
  public:
    // Constructor for a Hamamatsu S13360-6050VE MPPC.
    // See hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360-2050ve_etc_kapd1053e.pdf
    HamamatsuS133606050VE();

    // Destructor
    ~HamamatsuS133606050VE();

    G4double const GetTransverseDim()         override  { return 6.4*mm;          }   
    G4double const GetSensareaTransverseDim() override  { return 6.0*mm;          }
    G4double const GetThickness()             override  { return 1.35*mm;         }
    G4double const GetWindowThickness()       override  { return 0.1*mm;          }   
    G4double const GetSensareaThickness()     override  { return 0.1*mm;          }   
    G4String const GetModel()                 override  { return "HS133606050VE"; }

  private:
    G4int GetNPoints() override { return 72;}
    std::pair<G4int, G4double*> GetSensareaEnergyArray() override;
    std::pair<G4int, G4double*> GetSensareaEfficiencyArray() override;
    std::pair<G4int, G4double*> GetSensareaReflectivityArray() override;
  };

} // namespace nexus

#endif