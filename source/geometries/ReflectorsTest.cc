#include "ReflectorsTest.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4SubtractionSolid.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4NistManager.hh>
#include <G4GenericMessenger.hh>
#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

using namespace nexus;

REGISTER_CLASS(ReflectorsTest, GeometryBase)

namespace nexus{

  ReflectorsTest::ReflectorsTest():
  GeometryBase(),
  world_dx_(10.*m),
  world_dy_(10.*m),
  world_dz_(10.*m),
  reflector_dx_(10.*cm),
  reflector_dy_(3.5*mm),
  reflector_dz_(10.*cm),
  mpt_(opticalprops::Vikuiti())
  {
  }

  ReflectorsTest::~ReflectorsTest()
  {
  }

  void ReflectorsTest::Construct()
  {

    const G4String world_name = "WORLD";

    // Create the solid volume
    G4Box* world_solid = 
        new G4Box(world_name, world_dx_/2., world_dy_/2., world_dz_/2.);

    // Define the material
    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");

    // Set the optical material properties
    lAr->SetMaterialPropertiesTable(opticalprops::LAr());

    // Create the logical volume
    G4LogicalVolume* world_logic = 
        new G4LogicalVolume(world_solid, lAr, world_name);
    //world_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    this->SetLogicalVolume(world_logic);
    ConstructReflector(world_logic);
    return;
  }

  void ReflectorsTest::ConstructReflector(G4LogicalVolume* world_logic_vol)
  {

    const G4String reflector_name = "REFLECTOR_PLATE";

    G4Box* reflector_solid =
      new G4Box(reflector_name, reflector_dx_/2., reflector_dy_/2., reflector_dz_/2.);

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    G4LogicalVolume* reflector_logic =
      new G4LogicalVolume(reflector_solid, pvt, reflector_name);

    G4OpticalSurface* reflector_opsurf =
      new G4OpticalSurface("REFLECTOR_OPSURF", unified, ground, dielectric_metal);
    reflector_opsurf->SetMaterialPropertiesTable(mpt_);
    new G4LogicalSkinSurface("REFLECTOR_OPSURF", reflector_logic, reflector_opsurf);

    if(world_logic_vol){
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), 
                    reflector_logic, reflector_name, world_logic_vol, false, 0, true);
    }
    else{
        G4Exception("[ReflectorsTest]", "ConstructReflector()", FatalException,
        "World logical volume not found.");
    }
    return;
  }

  G4ThreeVector ReflectorsTest::GenerateVertex(const G4String&) const
  {
    G4double vertex_height_over_reflector_surf = 5.*cm;
    return G4ThreeVector(0., (reflector_dy_/2.)+vertex_height_over_reflector_surf, 0.);
  }
}