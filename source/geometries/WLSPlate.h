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
    WLSPlate(   G4bool with_LAr = true, G4bool with_dimples = false, G4String dimple_type = "cylindrical", 
                G4int dimples_no = 24, G4bool along_both_directions = false, G4double flat_dimple_width = 6.1*mm, 
                G4double flat_dimple_depth = 2.*mm, G4double curvy_dimple_radius = 3.1*mm);
    ///Construct a WLSPlate with given dimensions
    WLSPlate(   G4double, G4double, G4double, G4bool with_LAr = false, G4bool with_dimples = false, 
                G4String dimple_type = "cylindrical", G4int dimples_no = 24, G4bool along_both_directions = false, 
                G4double flat_dimple_width = 6.1*mm, G4double flat_dimple_depth = 2.*mm, 
                G4double curvy_dimple_radius = 3.1*mm);
    ///Construct a WLSPlate with given dimensions and optical properties
    WLSPlate(   G4double, G4double, G4double, G4MaterialPropertiesTable*, G4bool with_LAr = false, 
                G4bool with_dimples = false, G4String dimple_type = "cylindrical", G4int dimples_no = 24, 
                G4bool along_both_directions = false, G4double flat_dimple_width = 6.1*mm, 
                G4double flat_dimple_depth = 2.*mm, G4double curvy_dimple_radius = 3.1*mm);
    ///Destructor
    ~WLSPlate();

    G4ThreeVector GetDimensions();
    void          SetOpticalProperties(G4MaterialPropertiesTable*);
    void Construct();
    void ConstructWLSPlate(G4LogicalVolume*);
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4double dx_, dy_, dz_;                             ///< WLSPlate dimensions 
    G4bool with_dimples_;                               ///< Whether the plate has carved dimples
    G4String dimple_type_;                              ///< Dimple type. Might be 'flat', 'cylindrical' or 'spherical'.
    G4int how_many_dimples_;                            ///< How many dimples to carve on EACH selected side (up to along_both_directions_) of the plate
    G4bool along_both_directions_;                      ///< If so, dimples are carved along dx_ and dz_ sides of the plate. If not, they are only carved along the two dx_ sides.
    G4double flat_dimple_width_, flat_dimple_depth_;    ///< Used for flat dimples. The width of the dimple (along the board direction) and its depth, perpendicular to the plate surface.
    G4double curvy_dimple_radius_;                      ///< Used for cylindrical or spherical dimples. Radius of the dimple.
    G4MaterialPropertiesTable* mpt_;                    ///< WLS optical properties
    G4bool with_LAr_env_;                               ///< Whether to build a LAr environment or not 

    G4GenericMessenger* msg_;
  };

  inline G4ThreeVector  WLSPlate::GetDimensions()                                         { return G4ThreeVector(dx_, dy_, dz_); }
  inline void           WLSPlate::SetOpticalProperties(G4MaterialPropertiesTable* input)  { mpt_ = input; }

}

#endif