#ifndef DIMPLE_TEST
#define DIMPLE_TEST

#include "GeometryBase.h"

#include <G4ThreeVector.hh>

class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4VPhysicalVolume;

namespace nexus {

  /// Dimple test

  class DimpleTest: public GeometryBase
  {
  public:
    ///Default constructor
    DimpleTest();
    ///Destructor
    ~DimpleTest();

    void Construct();
    void ConstructWLSPlate(G4VPhysicalVolume*);
    void ConstructSiPM(G4VPhysicalVolume*);
    void ConstructArtificialDetector(G4VPhysicalVolume*);
    void ConstructReflectiveEnclosure(G4VPhysicalVolume*);
    G4bool geometry_is_ill_formed() const;
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4double world_radius_;                             ///< LAr sphere radius.
    G4int how_many_dimples_;                            ///< How many dimples to carve in front of the photosensor. how_many_dimples_==0 matches the case of no dimples.
    G4double gap_in_between_dimples_;                   ///< Gap in between two consecutive dimples
    G4String dimple_type_;                              ///< Dimple type. Must be one among 'cylindrical', 'spherical' and 'flat'
    G4double flat_dimple_width_;                        ///< Width of the flat dimple
    G4double flat_dimple_depth_;                        ///< Depth of the flat dimple
    G4double curvy_dimple_radius_;                      ///< Radius of the curvy dimples (i.e. cylindrical or spherical ones)
    G4double plate_dx_, plate_dy_, plate_dz_;           ///< WLS plate sample dimensions
    G4bool ref_phsensor_support_;                       ///< Whether the SiPM support is reflective (VIKUITI)
    G4double gap_;                                      ///< Gap between the plane surface of the WLS plate and the SiPM sensitive surface
    G4bool with_board_;                                 ///< Whether to add a reflective board behind the SiPM
    G4double board_thickn_;                             ///< Thickness of the SiPM-mounting board or the artificial detector
    G4double board_height_;                             ///< Height of the SiPM-mounting board or the artificial detector
    G4bool with_real_sipm_;                             ///< Whether to place a real SiPM or an artificial detector
    G4bool artificial_board_is_efficient_everywhere_ ;  ///< Whether the artificial detector is 100% efficient everywhere or it has a reflective chunk out of the virtual sipm boundaries
    G4double enclosure_gap_;                            ///< Gap between the reflective enclosure and the span of the geometry that is within the enclosure

    G4MaterialPropertiesTable* mpt_;                    ///< WLS optical properties
    G4GenericMessenger* msg_;                           ///< Messenger for the definition of control commands
  };

}

#endif