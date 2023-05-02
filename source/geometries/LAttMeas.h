#ifndef LATT_MEAS
#define LATT_MEAS

#include "GeometryBase.h"

class G4VPhysicalVolume;
class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4UserLimits;

namespace nexus {

  /// Setup used to measure the attenuation length of a G2p WLS plate

  class LAttMeas: public GeometryBase
  {

  public:
    ///Constructor
    LAttMeas();
    ///Destructor
    ~LAttMeas();

    void Construct();                   ///< Constructs the geometry

  private:

    G4VPhysicalVolume* ConstructBlackBox(G4VPhysicalVolume*) const;         ///< Called by Construct(). Returns the black box physical volume.

    void ConstructWLSPlate(G4VPhysicalVolume*) const;         ///< Called by Construct(). Adds the WLS plate.
    void ConstructPlateAdapter(G4VPhysicalVolume*) const;     ///< Called by Construct(). Adds the plate-PMT adapter.
    void ConstructPMT(G4VPhysicalVolume*) const;              ///< Called by Construct(). Adds the PMT.
    void ConstructPlateHolders(G4VPhysicalVolume*) const;     ///< Called by Construct(). Adds the plate holders.
    
    G4ThreeVector GenerateVertex(const G4String&) const;
    G4bool geometry_is_ill_formed();                          ///<Check that the geometry is feasible.

  private:
    G4int config_code_;                                             ///< Value in (1,2) which labels 
                                                                    ///< 1 - For LED measure 
                                                                    ///< 2 - For Laser measure along long side
                                                                    ///< 3 - For Laser measure along short side

    G4bool with_box_;                           ///< Whether to construct a black box surrounding the plate or not
    G4double black_box_internal_dx_;            ///< Internal dimensions of the black box   
    G4double black_box_internal_dy_;
    G4double black_box_internal_dz_;

    G4double plate_x_pos_;                      ///< Cartesian coordinates of the geometric
    G4double plate_y_pos_;                      ///< center of the WLS plate with respect to
    G4double plate_z_pos_;                      ///< the geometric center of the black box

    G4double plate_dx_, plate_dy_, plate_dz_;   ///< Dimensions of the WLS plate
    G4double wls_attlength_;                    ///< Attenuation length of the WLS plate
    G4double gap_;                              ///< Gap between the WLS plate and the adapter/PMT

    G4bool with_holder_;                                    ///< Whether to construct the holder or not
    G4double holder_x_pos_, holder_y_pos_, holder_z_pos_;   ///< Coordinates of the holder_encasing
    G4double holder_dx_, holder_dy_, holder_dz_;            ///< Dimensions of the holder_encasing

    G4double pmt_case_depth_;                   ///< Depth of the PMT case                    
    G4double pmt_case_diameter_;                ///< Diameter of the PMT case
    G4double pmt_sensarea_diameter_;            ///< Diameter of the PMT sensitive area (window)

    G4bool with_adapter_;                       ///< Whether to place an adapter in between the WLS plate and the PMT

    G4double adp_diameter_;                     ///< Adapter diameter
    G4double adp_thickn_;                       ///< Adapter thickness
    G4double adp_hole_height_;                  ///< Height of the adapter hole
    G4double adp_slot_width_;                   ///< Width of the adapter slot
    G4double adp_slot_thickn_;                  ///< Thickness of the adapter right at the slot  
  
    G4double world_length_;                     ///< Extra thickness for the surrounding box world to wrap the XArapuca

    ///---- Internal attributes ----///         ///< These attributes are internal. They must not be set by the user.
    G4double pmt_x_pos_;                        ///< Cartesian coordinates of the geometric center of the PMT
    G4double pmt_y_pos_;  
    G4double pmt_z_pos_; 
    G4double adapter_x_pos_;                    ///< Cartesian coordinates of the geometric center of the adapter
    G4double adapter_y_pos_; 
    G4double adapter_z_pos_;
    G4double gen_vertex_y_;                     ///< Cartesian coordinates of the generation vertex
    G4double gen_vertex_z_;                     
    //------------------------------------------------------------------------------------------

    G4double gen_vertex_x_;                     ///< Cartesian coordinates of the generation vertex

    G4GenericMessenger* msg_;                   ///< Messenger for the definition of control commands
    G4UserLimits* ul_;                          ///< Useful to set a maximum track length in the plate
  };
}

#endif