#ifndef SIPM_MPPC_H
#define SIPM_MPPC_H

#include "GeometryBase.h"

template <class T1, class T2>
class pair;

namespace nexus {

  class SiPMMPPC: public GeometryBase
  {
  public:
    // Constructor for a base (abstract) clase from which to derive the classes
    // for every SiPM MPPC, such as HPK S13360-5075-HD-HQR, HPK S13360-6050-VE or
    // FBK NUV HD CRYO TT. This class can also model SiPMMPPC arrays, such as 
    // Broadcom AFBR-S4N44P044M.
    SiPMMPPC();

    // Destructor
    ~SiPMMPPC();

    void Construct();
    
    //Set methods
    inline void SetVisibility(G4bool input)                 { visibility_ = input; } 
    inline void SetReflectiveSupports(G4bool input)         { reflective_support_ = input; }
  
    //Get methods

    virtual G4double GetTransverseDim() =0;           ///< One of the overall transverse dimensions 
                                                      ///< (the MPPC section is square).
    virtual G4int GetSiPMsNoAlongX() =0;              ///< Number of SiPMs along the x axis.
    virtual G4int GetSiPMsNoAlongZ() =0;              ///< Number of SiPMs along the z axis.
    virtual G4double GetOuterFrameWidthAlongX() =0;   ///< Let the 'frame' be the dead space in between
                                                      ///< the different SiPMs of this SiPMMPPC array.
                                                      ///< Then this function returns the width of the
                                                      ///< outer frame along the x axis.
    virtual G4double GetOuterFrameWidthAlongZ() =0;   ///< The width of the outer frame along the z axis.
    virtual G4double GetInnerFramesWidthAlongX() =0;  ///< The width of the inner frames along the x axis.
    virtual G4double GetInnerFramesWidthAlongZ() =0;  ///< The width of the inner frames along the z axis.
    virtual G4double GetThickness() =0;               ///< Overall thickness of the MPPC.
    virtual G4double GetWindowThickness() =0;         ///< Thickness of the epoxy window.   
    virtual G4double GetSensareaThickness() =0;       ///< Thickness of the sensitive volume. This is an 
                                                      ///  approximation. This dimension does not 
                                                      ///  matter since the sensarea will be coated with 
                                                      ///  a sensitive optical surface and the only exposed 
                                                      ///  face is the front one.
    virtual G4String GetModel() =0;                   ///< Name of the SiPM model.

    G4double GetSensareaTransverseDimAlongX();        ///< Transverse dimension of the 
                                                      ///< sensitive areas along the x axis the 
    G4double GetSensareaTransverseDimAlongZ();        ///< Transverse dimension of the 
                                                      ///< sensitive areas along the z axis the 
  private:
    virtual G4int GetNPoints()=0;
    virtual std::pair<G4int, G4double*> GetSensareaEnergyArray()=0;
    virtual std::pair<G4int, G4double*> GetSensareaEfficiencyArray()=0;
    virtual std::pair<G4int, G4double*> GetSensareaReflectivityArray()=0;

    G4bool visibility_;                   ///<Whether the MPPC is visible or not
    G4bool reflective_support_;           ///<Whether the FR4 support is vikuiti-coated
    G4double reflectivity_scale_factor_;  ///<Only used if reflective_support_ is true. In that case, this
                                          ///<is a double which must belong to the [0.0, 1.0] range which
                                          ///<is used to scale factor for the vikuiti reflectivity curve.


  };

} // namespace nexus

#endif