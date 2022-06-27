#ifndef DICHROIC_TEST
#define DICHROIC_TEST

#include "GeometryBase.h"

#include <G4ThreeVector.hh>

class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4VPhysicalVolume;

namespace nexus {

  /// Dichroic filter test

  class DichroicTest: public GeometryBase
  {
  public:
    ///Default constructor
    DichroicTest();
    ///Destructor
    ~DichroicTest();

    G4ThreeVector GetDimensions();
    void          SetOpticalProperties(G4MaterialPropertiesTable*);
    void Construct();
    void ConstructDichroicFilter(G4VPhysicalVolume*);
    void ConstructDetectors(G4VPhysicalVolume*);
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4double radius_, thickn_;                      ///<Cylindrical filter dimensions
    G4double det_thickn_, det_depth_;  ///<Cylindrical detectors dimensions
    G4MaterialPropertiesTable* mpt_;                ///<Dichroic filter substrate optical properties

    G4GenericMessenger* msg_;
  };

  inline G4ThreeVector  DichroicTest::GetDimensions()                                         { return G4ThreeVector(radius_, radius_, thickn_); }
  inline void           DichroicTest::SetOpticalProperties(G4MaterialPropertiesTable* input)  { mpt_ = input; }

}

#endif