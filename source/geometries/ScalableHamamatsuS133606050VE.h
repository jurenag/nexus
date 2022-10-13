#ifndef SCALABLE_HAMAMATSU_S13360_6050VE_H
#define SCALABLE_HAMAMATSU_S13360_6050VE_H

#include "GeometryBase.h"

namespace nexus {

  // This class aims to model a SiPM with the photon-detection specifications of HamamatsuS133606050VE
  // but with non-fixed transverse dimensions. Such dimensions can be chosen for this class' objects.
  class ScalableHamamatsuS133606050VE: public GeometryBase
  {
  public:
    // Constructor for a Hamamatsu S13360-6050VE MPPC with a given transverse dimensions
    // See hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360-2050ve_etc_kapd1053e.pdf
    ScalableHamamatsuS133606050VE(G4double height = 6.4*mm, G4double width = 6.4*mm);

    // Destructor
    ~ScalableHamamatsuS133606050VE();

    static const G4double thickness_;           ///<Overall thickness of the MPPC
    static const G4double window_thickness_;    ///<Thickness of the epoxy window
    static const G4double sensarea_thickness_;  ///<Thickness of the sensitive volume. This is an approximation. This dimension does not matter since the sensarea will be coated with a sensitive optical surface and the only exposed face is the front one.

    //Get methods
    static G4double GetThickness();
    G4double GetHeight()                                    { return height_; }
    G4double GetWidth()                                     { return width_; }

    void Construct();

    //Set methods
    void SetHeight(G4double input)                          { height_               = input; }
    void SetWidth(G4double input)                           { width_                = input; }
    inline void SetVisibility(G4bool input)                 { visibility_           = input; } 
    inline void SetReflectiveSupports(G4bool input)         { reflective_support_   = input; }

  private:

    G4double height_;                           ///<Height of the MPPC
    G4double width_;                            ///<Width of the MPPC
    G4bool visibility_;                         ///<Whether the MPPC is visible or not
    G4bool reflective_support_;                 ///<Whether the FR4 support is vikuiti-coated
  };

} // namespace nexus

#endif