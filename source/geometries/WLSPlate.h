#ifndef WLS_PLATE
#define WLS_PLATE

#include "GeometryBase.h"

#include <G4ThreeVector.hh>

class G4MaterialPropertiesTable;
class G4GenericMessenger;

namespace nexus {

  /// WLS Plate

  class WLSPlate: public GeometryBase
  {
  public:
    ///Default constructor
    WLSPlate(G4bool with_LAr = true);
    ///Construct a WLSPlate with given dimensions
    WLSPlate(G4double, G4double, G4double, G4bool with_LAr = false);
    ///Construct a WLSPlate with given dimensions and optical properties
    WLSPlate(G4double, G4double, G4double, G4MaterialPropertiesTable*, G4bool with_LAr = false);
    ///Destructor
    ~WLSPlate();

    G4ThreeVector GetDimensions();
    void          SetOpticalProperties(G4MaterialPropertiesTable*);
    void Construct();
    void ConstructWLSPlate(G4LogicalVolume*);
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4double dx_, dy_, dz_;             ///<WLSPlate dimensions 
    G4MaterialPropertiesTable* mpt_;    ///<WLS optical properties
    G4bool with_LAr_env_;               ///<Whether to build a LAr environment or not 

    G4GenericMessenger* msg_;
  };

  inline G4ThreeVector  WLSPlate::GetDimensions()                                         { return G4ThreeVector(dx_, dy_, dz_); }
  inline void           WLSPlate::SetOpticalProperties(G4MaterialPropertiesTable* input)  { mpt_ = input; }

}

#endif