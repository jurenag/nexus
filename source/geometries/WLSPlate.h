#ifndef WLS_PLATE
#define WLS_PLATE

#include "GeometryBase.h"

#include <G4ThreeVector.hh>
#include <G4SystemOfUnits.hh>

class G4MaterialPropertiesTable;
class G4GenericMessenger;

namespace nexus {

  /// WLS Plate

  class WLSPlate: public GeometryBase
  {
  public:
    ///Default constructor
    WLSPlate(   G4bool with_LAr = true, 
                G4bool dimples_at_x_plus = false, 
                G4bool dimples_at_x_minus = false,
                G4bool dimples_at_z_plus = false, 
                G4bool dimples_at_z_minus = false, 
                G4String dimple_type = "cylindrical", 
                G4int how_many_dimples = 24, 
                G4double flat_dimple_width = 6.1*mm, 
                G4double flat_dimple_depth = 2.*mm, 
                G4double curvy_dimple_radius = 3.1*mm,
                G4bool cut_plate = true,
                G4double cut_angle = 45.*deg,
                G4double cut_thickness = 1.*mm,
                G4double tunneling_probability = 0.0);
    ///Construct a WLSPlate with given dimensions
    WLSPlate(   G4double, 
                G4double, 
                G4double, 
                G4bool with_LAr = false, 
                G4bool dimples_at_x_plus = false, 
                G4bool dimples_at_x_minus = false, 
                G4bool dimples_at_z_plus = false, 
                G4bool dimples_at_z_minus = false, 
                G4String dimple_type = "cylindrical", 
                G4int how_many_dimples = 24, 
                G4double flat_dimple_width = 6.1*mm, 
                G4double flat_dimple_depth = 2.*mm, 
                G4double curvy_dimple_radius = 3.1*mm,
                G4bool cut_plate = false,
                G4double cut_angle = 45.*deg,
                G4double cut_thickness = 1.*mm,
                G4double tunneling_probability = 0.0);
    ///Construct a WLSPlate with given dimensions and optical properties
    WLSPlate(   G4double, 
                G4double, 
                G4double, 
                G4MaterialPropertiesTable*, 
                G4bool with_LAr = false, 
                G4bool dimples_at_x_plus = false, 
                G4bool dimples_at_x_minus = false, 
                G4bool dimples_at_z_plus = false, 
                G4bool dimples_at_z_minus = false, 
                G4String dimple_type = "cylindrical", 
                G4int how_many_dimples = 24, 
                G4double flat_dimple_width = 6.1*mm, 
                G4double flat_dimple_depth = 2.*mm, 
                G4double curvy_dimple_radius = 3.1*mm,
                G4bool cut_plate = false,
                G4double cut_angle = 45.*deg,
                G4double cut_thickness = 1.*mm,
                G4double tunneling_probability = 0.0);
    ///Destructor
    ~WLSPlate();

    G4ThreeVector GetDimensions();
    void          SetOpticalProperties(G4MaterialPropertiesTable*);
    void Construct();
    void ConstructWLSPlate(G4LogicalVolume*);
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4double dx_, dy_, dz_;                             ///< WLSPlate dimensions 
    G4bool dimples_at_x_plus_;                          ///< Whether to carve dimples at the WLSPlate side which is contained within x>0.0
    G4bool dimples_at_x_minus_;                         ///< Whether to carve dimples at the WLSPlate side which is contained within x<0.0
    G4bool dimples_at_z_plus_;                          ///< Whether to carve dimples at the WLSPlate side which is contained within z>0.0
    G4bool dimples_at_z_minus_;                         ///< Whether to carve dimples at the WLSPlate side which is contained within z<0.0
    G4String dimple_type_;                              ///< Dimple type. Might be 'flat', 'cylindrical' or 'spherical'.
    G4int how_many_dimples_;                            ///< How many dimples to carve on EACH selected side of the plate
    G4double flat_dimple_width_, flat_dimple_depth_;    ///< Used for flat dimples. The width of the dimple (along the board direction) and its depth, perpendicular to the plate surface.
    G4double curvy_dimple_radius_;                      ///< Used for cylindrical or spherical dimples. Radius of the dimple.
    G4bool cut_plate_;                                  ///< Whether to split the WLS plate up into two pieces
                                                        ///< WARNING: Although the application works fine, the visualization may
                                                        ///<          break down when dimples are used and cut_plate_ is true.
    G4double cut_angle_;                                ///< This parameter only makes a difference if cut_plate_ is true. Angle of the cut, with respect to the Z axis (the axis along which the WLS plate is dz_ long)
    G4double cut_thickness_;                            ///< This parameter only makes a difference if cut_plate_ is true. Thickness of the cut/crack that is carved from the plate.
    G4double tunneling_probability_;                    ///< Probability that a photon tunnels through the surface of the WLSPlate (inwards or outwards), meaning that it is 
                                                        ///< straightfoward-ly transmitted (no deflection due to Frensel refraction). The goal of this parameter is to let 
                                                        ///< us model surface imperfections of the WLS plate.
    G4double generation_y_pos_;                         ///< The y-position of the generation vertex of photons.
    G4String generation_mode_;                          ///< Controls the generation vertex of photons. It can take the following values
                                                        ///< "random"       - The x- and z-positions of the generation vertex are randomly sampled 
                                                        ///<                  in the intervals (-dx_/2, dx_/2) and (-dz_/2, dz_/2), respectively.
                                                        ///< any other case - The x- and z-positions of the generation vertex are set to 0.0.
    G4MaterialPropertiesTable* mpt_;                    ///< WLS optical properties
    G4bool with_LAr_env_;                               ///< Whether to build a LAr environment or not 

    G4GenericMessenger* msg_;
  };

  inline G4ThreeVector  WLSPlate::GetDimensions()                                         { return G4ThreeVector(dx_, dy_, dz_); }
  inline void           WLSPlate::SetOpticalProperties(G4MaterialPropertiesTable* input)  { mpt_ = input; }

}

#endif