#include "LArSphere.h"

#include "MaterialsList.h"
#include "FactoryBase.h"
#include "OpticalMaterialProperties.h"
#include "SpherePointSampler.h"

#include <G4GenericMessenger.hh>
#include <G4Orb.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>

#include <CLHEP/Units/SystemOfUnits.h>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(LArSphere, GeometryBase)

namespace nexus {

  LArSphere::LArSphere():
    GeometryBase(), radius_(1.*m)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/LArSphere/",
      "Control commands of geometry LArSphere.");

    G4GenericMessenger::Command& radius_cmd =
      msg_->DeclareProperty("radius", radius_, "Radius of the liquid argon sphere.");
    radius_cmd.SetUnitCategory("Length");
    radius_cmd.SetParameterName("radius", false);
    radius_cmd.SetRange("radius>0.");

    // Create a vertex generator for a sphere
    sphere_vertex_gen_ = new SpherePointSampler(radius_, 0.);
  }

  LArSphere::LArSphere(G4double radius):
  GeometryBase(), radius_(radius)
  {
    // Create a vertex generator for a sphere
    sphere_vertex_gen_ = new SpherePointSampler(radius_, 0.);
  }

  LArSphere::~LArSphere()
  {
    if(msg_) delete msg_;
    delete sphere_vertex_gen_;
   
  }

  void LArSphere::Construct()
  {
    G4String name = "LAR_SPHERE";

    // Define solid volume as a sphere
    G4Orb* sphere_solid = new G4Orb(name, radius_);

    // Define the material
    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");

    // Set the optical material properties
    lAr->SetMaterialPropertiesTable(opticalprops::LAr());

    // Define the logical volume of the sphere using the material
    // and the solid volume defined above
    G4LogicalVolume* sphere_logic =
    new G4LogicalVolume(sphere_solid, lAr, name);
    GeometryBase::SetLogicalVolume(sphere_logic);
  }



  G4ThreeVector LArSphere::GenerateVertex(const G4String& region) const
  {
    return sphere_vertex_gen_->GenerateVertex(region);
  }


} // end namespace nexus
