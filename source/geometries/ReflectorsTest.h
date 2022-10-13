#ifndef REFLECTORS_TEST
#define REFLECTORS_TEST

#include "GeometryBase.h"

#include <G4ThreeVector.hh>

class G4MaterialPropertiesTable;

namespace nexus {

  /// Reflectors test

  class ReflectorsTest: public GeometryBase
  {
  public:
    ///Default constructor
    ReflectorsTest();
    ///Destructor
    ~ReflectorsTest();

    void Construct();
    void ConstructReflector(G4LogicalVolume*);
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4double world_dx_, world_dy_, world_dz_, reflector_dx_, reflector_dy_, reflector_dz_; ///<World and WLSPlate dimension
    G4MaterialPropertiesTable* mpt_;                                        ///<WLS optical properties
  };
}

#endif