#ifndef SIPM_MPPC_H
#define SIPM_MPPC_H

#include "GeometryBase.h"

template <class T1, class T2>
class pair;

namespace nexus {

  class SiPMMPPC: public GeometryBase
  {
  public:
    // Constructor for a the base (abstract) clase from which to derive the classes
    // for every SiPM MPPC, such as HPK S13360-5075-HD-HQR, HPK S13360-6050-VE or
    // FBK NUV HD CRYO TT.
    SiPMMPPC();

    // Destructor
    ~SiPMMPPC();

    void Construct();
    
    //Set methods
    inline void SetVisibility(G4bool input)                 { visibility_ = input; } 
    inline void SetReflectiveSupports(G4bool input)         { reflective_support_ = input; }
  
    //Get methods

    virtual G4double const GetTransverseDim() =0;         ///< One of the overall transverse dimensions 
                                                          ///  (the MPPC section is square)
    virtual G4double const GetSensareaTransverseDim() =0; ///< One of the transverse dimensions for the 
                                                          ///  sensitive volume (its section is also an square) 
    virtual G4double const GetThickness() =0;             ///< Overall thickness of the MPPC
    virtual G4double const GetWindowThickness() =0;       ///< Thickness of the epoxy window    
    virtual G4double const GetSensareaThickness() =0;     ///< Thickness of the sensitive volume. This is an 
                                                          ///  approximation. This dimension does not 
                                                          ///  matter since the sensarea will be coated with 
                                                          ///  a sensitive optical surface and the only exposed 
                                                          ///  face is the front one.
    virtual G4String const GetModel() =0;                 ///< Name of the SiPM model   

  private:
    virtual G4int GetNPoints()=0;
    virtual std::pair<G4int, G4double*> GetSensareaEnergyArray()=0;
    virtual std::pair<G4int, G4double*> GetSensareaEfficiencyArray()=0;
    virtual std::pair<G4int, G4double*> GetSensareaReflectivityArray()=0;

    G4bool visibility_;                   ///<Whether the MPPC is visible or not
    G4bool reflective_support_;           ///<Whether the FR4 support is vikuiti-coated


  };

} // namespace nexus

#endif