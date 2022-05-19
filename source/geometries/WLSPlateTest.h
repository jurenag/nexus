#ifndef WLS_PLATE_TEST
#define WLS_PLATE_TEST

#include "GeometryBase.h"

#include <G4ThreeVector.hh>

class G4MaterialPropertiesTable;

namespace nexus {

  /// WLS Plate test

  class WLSPlateTest: public GeometryBase
  {
  public:
    ///Default constructor
    WLSPlateTest();
    ///Destructor
    ~WLSPlateTest();

    void Construct();
    void ConstructCollectors(G4LogicalVolume*, G4bool, G4bool);
    void ConstructCovers(G4LogicalVolume*);
    void ConstructWLSPlate(G4LogicalVolume*);
    G4bool GeometryIsIllFormed();
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4double world_dx_, world_dy_, world_dz_, wlsp_dx_, wlsp_dy_, wlsp_dz_; ///<World and WLSPlate dimensions
    G4double gap_;                                              ///<Thickness of LAr atmosphere between wls plate and collectors. If covers are set, then this is the thickness of such atmosphere between the covers and the collectors.
    G4bool upper_collector_, lower_collector_;                              ///<Whether to construct the upper and lower collectors
    G4double collector_thickn_;                                             ///<Thickness of the photon collector  
    G4bool upper_cover_, lower_cover_;                     ///<Whether to add a perfect absorption layer over and/or under WLSPlate.
    G4double cover_thickn_;                                                 ///<Thickness of the absorption covers.
    G4MaterialPropertiesTable* mpt_;                                        ///<WLS optical properties
  };
}

#endif