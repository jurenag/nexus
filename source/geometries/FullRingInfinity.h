#ifndef FULLRINGINF_
#define FULLRINGINF_

#include "BaseGeometry.h"

class G4GenericMessenger;
class G4LogicalVolume;
namespace nexus {
  class SiPMpetFBK;
  class QuadFBK;
  class CylinderPointSampler;
}

namespace nexus {
  class FullRingInfinity : public BaseGeometry {

  public:
    // Constructor
    FullRingInfinity();
    //Destructor
    ~FullRingInfinity();

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    private:
    void Construct();
    void BuildCryostat();
    void BuildQuadSensors();
    void BuildSensors(); 
    void BuildPhantom(); 

    QuadFBK* quad_;
    SiPMpetFBK* sipm_;
    
    G4LogicalVolume* lab_logic_;
    G4LogicalVolume* LXe_logic_;
    G4LogicalVolume* active_logic_;

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_; 

    G4double lat_dimension_cell_;
    G4int n_cells_; ///< number of virtual cells of ~ 5 cm of side I want to fit in the ring
    G4int lin_n_sipm_per_cell_; ///< linear number of sipms sin a cell (the side, not the area)
    G4int lin_n_quad_per_cell_; ///< linear number of quads in a cell (the side, not the area)
    G4double quad_pitch_; ///< space between the centres of two quads in both dimensions
    G4double kapton_thickn_;
    G4double depth_;

    G4double internal_radius_, external_radius_;
    G4double cryo_width_, cryo_thickn_;

    G4double phantom_diam_; 
    G4double phantom_length_;

    G4double max_step_size_;

    G4bool quad_arrangement_;

    CylinderPointSampler* cylindric_gen_;

  };
}
#endif
