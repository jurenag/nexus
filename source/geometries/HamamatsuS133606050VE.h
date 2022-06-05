#ifndef HAMAMATSU_S13360_6050VE_H
#define HAMAMATSU_S13360_6050VE_H

#include "GeometryBase.h"

namespace nexus {

  class HamamatsuS133606050VE: public GeometryBase
  {
  public:
    // Constructor for a Hamamatsu S13360-6050VE MPPC.
    // See hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360-2050ve_etc_kapd1053e.pdf
    HamamatsuS133606050VE();

    // Destructor
    ~HamamatsuS133606050VE();

    static const G4double transverse_dim_;      ///<One of the overall transverse dimensions (the MPPC section is square)
    static const G4double thickness_;           ///<Overall thickness of the MPPC
    static const G4double window_thickness_;    ///<Thickness of the epoxy window
    static const G4double sensarea_td_;         ///<One of the transverse dimensions for the sensitive volume (its section is also an square)        
    static const G4double sensarea_thickness_;  ///<Thickness of the sensitive volume. This is an approximation. This dimension does not matter since the sensarea will be coated with a sensitive optical surface and the only exposed face is the front one.

    //Get methods
    static G4double GetWidth();
    static G4double GetHeight();
    static G4double GetThickness();

    void Construct();

    //Set methods
    inline void SetVisibility(G4bool input)                 { visibility_ = input; } 
    inline void SetReflectiveSupports(G4bool input)         { reflective_support_ = input; }

  private:

    G4bool visibility_;                         ///<Whether the MPPC is visible or not
    G4bool reflective_support_;                 ///<Whether the FR4 support is vikuiti-coated
  };

} // namespace nexus

#endif