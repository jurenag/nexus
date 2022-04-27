#include "WLSPlate.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "LArSphere.h"

#include <G4Box.hh>
#include <G4Orb.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4GenericMessenger.hh>
#include <G4SystemOfUnits.hh>

using namespace nexus;

REGISTER_CLASS(WLSPlate, GeometryBase)

namespace nexus{

  WLSPlate::WLSPlate(G4bool with_LAr):
  GeometryBase(),
  dx_(487.*mm),
  dy_(3.5*mm),
  dz_(93.*mm), //<Default WLS plate dimensions taken from FD1 TDR vol. IX
  mpt_(opticalprops::EJ286()),
  with_LAr_env_(with_LAr)
  {

    msg_ = new G4GenericMessenger(this, "/Geometry/WLSPlate/",
				"Control commands of geometry WLSPlate.");
    
    G4GenericMessenger::Command& dx_cmd =
      msg_->DeclareProperty("depth", dx_,
			    "Depth of the WLS plate.");
    dx_cmd.SetUnitCategory("Length");
    dx_cmd.SetParameterName("depth", false);
    dx_cmd.SetRange("depth>0.");

    G4GenericMessenger::Command& dy_cmd =
      msg_->DeclareProperty("height", dy_,
			    "Height of the WLS plate.");
    dy_cmd.SetUnitCategory("Length");
    dy_cmd.SetParameterName("height", false);
    dy_cmd.SetRange("height>0.");

    G4GenericMessenger::Command& dz_cmd =
      msg_->DeclareProperty("width", dz_,
			    "Width of the WLS plate.");
    dz_cmd.SetUnitCategory("Length");
    dz_cmd.SetParameterName("width", false);
    dz_cmd.SetRange("width>0.");

    msg_->DeclareProperty("LAr", with_LAr_env_,
    "Build a LAr sphere surrounding the WLSPlate.");

  }

  WLSPlate::WLSPlate(G4double dx, G4double dy, G4double dz, G4bool with_LAr):
  GeometryBase(),
  dx_(dx),
  dy_(dy),
  dz_(dz),
  mpt_(nullptr),
  with_LAr_env_(with_LAr)
  {
  }

  WLSPlate::WLSPlate(G4double dx, G4double dy, G4double dz, G4MaterialPropertiesTable* mpt, G4bool with_LAr):
  GeometryBase(),
  dx_(dx),
  dy_(dy),
  dz_(dz),
  mpt_(mpt),
  with_LAr_env_(with_LAr)
  {
  }

  WLSPlate::~WLSPlate()
  {
    if(msg_) delete msg_;
  }

  void WLSPlate::Construct()
  {

    G4LogicalVolume* sphere_logic = nullptr;
    if(with_LAr_env_){
        G4String name = "LAR_SPHERE";

        // Define solid volume as a sphere
        G4Orb* sphere_solid = new G4Orb(name, dx_+dy_+dz_+(5.*m));

        // Define the material
        G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");

        // Set the optical material properties
        lAr->SetMaterialPropertiesTable(opticalprops::paulucci_LAr());

        // Define the logical volume of the sphere using the material
        // and the solid volume defined above
        sphere_logic = new G4LogicalVolume(sphere_solid, lAr, name);
        sphere_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
        this->SetLogicalVolume(sphere_logic);
    }
    ConstructWLSPlate(sphere_logic);
    return;

  }

  void WLSPlate::ConstructWLSPlate(G4LogicalVolume* world_logic_vol)
  {

    const G4String plate_name = "WLS_PLATE";

    G4Box* plate_solid =
      new G4Box(plate_name, dx_/2., dy_/2., dz_/2.);

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    if(mpt_ && pvt){
      pvt->SetMaterialPropertiesTable(mpt_);
    }
    else{
      G4Exception("[WLSPlate]", "Construct()", JustWarning,
      "The optical properties of the WLS Plate were not set");
    }

    G4LogicalVolume* plate_logic =
      new G4LogicalVolume(plate_solid, pvt, plate_name);

    if(world_logic_vol){
        new G4PVPlacement(nullptr, G4ThreeVector{}, plate_logic, plate_name, 
                        world_logic_vol, false, 0, true);
    }
    else{
        this->SetLogicalVolume(plate_logic);
    }
    return;
  }

  G4ThreeVector WLSPlate::GenerateVertex(const G4String&) const
  {
    return G4ThreeVector(0., 100.*mm, 0.);
  }
}