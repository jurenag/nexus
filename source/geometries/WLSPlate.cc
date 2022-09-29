#include "WLSPlate.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "LArSphere.h"

#include <G4VSolid.hh>
#include <G4MultiUnion.hh>
#include <G4SubtractionSolid.hh>
#include <G4Transform3D.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Orb.hh>
#include <G4RotationMatrix.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4GenericMessenger.hh>

using namespace nexus;

REGISTER_CLASS(WLSPlate, GeometryBase)

namespace nexus{

  WLSPlate::WLSPlate(G4bool with_LAr, G4bool with_dimples, G4String dimple_type, 
                    G4int dimples_no, G4bool along_both_directions, G4double flat_dimple_width, 
                    G4double flat_dimple_depth, G4double curvy_dimple_radius):
  GeometryBase(),
  dx_(487.*mm),
  dy_(3.5*mm),
  dz_(93.*mm), //<Default WLS plate dimensions taken from FD1 TDR vol. IX
  with_LAr_env_(with_LAr),
  with_dimples_(with_dimples),
  dimple_type_(dimple_type),
  how_many_dimples_(dimples_no),
  along_both_directions_(along_both_directions),
  flat_dimple_width_(flat_dimple_width),
  flat_dimple_depth_(flat_dimple_depth),
  curvy_dimple_radius_(curvy_dimple_radius),
  mpt_(opticalprops::EJ286())
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

  WLSPlate::WLSPlate(G4double dx, G4double dy, G4double dz, G4bool with_LAr, 
                    G4bool with_dimples, G4String dimple_type, G4int dimples_no,
                    G4bool along_both_directions, G4double flat_dimple_width, 
                    G4double flat_dimple_depth, G4double curvy_dimple_radius):
  GeometryBase(),
  dx_(dx),
  dy_(dy),
  dz_(dz),
  with_LAr_env_(with_LAr),
  with_dimples_(with_dimples),
  dimple_type_(dimple_type),
  how_many_dimples_(dimples_no),
  along_both_directions_(along_both_directions),
  flat_dimple_width_(flat_dimple_width),
  flat_dimple_depth_(flat_dimple_depth),
  curvy_dimple_radius_(curvy_dimple_radius),
  mpt_(opticalprops::EJ286())
  {
  }

  WLSPlate::WLSPlate(G4double dx, G4double dy, G4double dz, G4MaterialPropertiesTable* mpt, 
                    G4bool with_LAr, G4bool with_dimples, G4String dimple_type, G4int dimples_no, 
                    G4bool along_both_directions, G4double flat_dimple_width, 
                    G4double flat_dimple_depth, G4double curvy_dimple_radius):
  GeometryBase(),
  dx_(dx),
  dy_(dy),
  dz_(dz),
  mpt_(mpt),
  with_dimples_(with_dimples),
  dimple_type_(dimple_type),
  how_many_dimples_(dimples_no),
  along_both_directions_(along_both_directions),
  flat_dimple_width_(flat_dimple_width),
  flat_dimple_depth_(flat_dimple_depth),
  curvy_dimple_radius_(curvy_dimple_radius),
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
        lAr->SetMaterialPropertiesTable(opticalprops::LAr());

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

    G4VSolid* geometry_solid = nullptr;
    G4Box* whole_plate_solid = new G4Box(plate_name, dx_/2., dy_/2., dz_/2.);

    if(with_dimples_ && how_many_dimples_>=1)
    {

        G4double tolerance = 0.5*mm; // To avoid boolean subtraction of matching surfaces
        G4VSolid* carving_solid = nullptr;
        if(dimple_type_=="flat"){
            // Yes, using 2*flat_dimple_depth_ as the whole dimension of the carving along the z-axis is ok 
            // (the carvings are subtracted from the very edge of the plate)
            carving_solid = dynamic_cast<G4VSolid*>(new G4Box("AUX", flat_dimple_width_/2., dy_/2. +tolerance, flat_dimple_depth_));
        }
        else if(dimple_type_=="spherical"){
            carving_solid = dynamic_cast<G4VSolid*>(new G4Orb("AUX", curvy_dimple_radius_));
        }
        else{ //Default is cylindrical dimples
            carving_solid = dynamic_cast<G4VSolid*>(
                                        new G4Tubs("AUX", 0., curvy_dimple_radius_, dy_/2. +tolerance, 0., twopi));
        }

        G4RotationMatrix* rot = new G4RotationMatrix();
        if(dimple_type_=="cylindrical"){
            rot->rotateX(+90.*deg);
        }

        G4Transform3D* transform_ptr = nullptr;
        G4MultiUnion* carvings_multiunion_solid = new G4MultiUnion("CARVINGS");
        for(G4int i=0; i<how_many_dimples_; i++){
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector((-1.*dx_/2.)+(1.*(0.5+i)*dx_/(1.*how_many_dimples_)), 0., +dz_/2.));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector((-1.*dx_/2.)+(1.*(0.5+i)*dx_/(1.*how_many_dimples_)), 0., -1.*dz_/2.));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
        }

        if(dimple_type_=="flat"){
            rot->rotateY(+90.*deg);
        }

        if(along_both_directions_){
            for(G4int i=0; i<how_many_dimples_; i++){
                transform_ptr = new G4Transform3D(*rot, G4ThreeVector(+dx_/2., 0., (-1.*dz_/2.)+(1.*(0.5+i)*dz_/(1.*how_many_dimples_))));
                carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
                transform_ptr = new G4Transform3D(*rot, G4ThreeVector(-1.*dx_/2., 0., (-1.*dz_/2.)+(1.*(0.5+i)*dz_/(1.*how_many_dimples_))));
                carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
            }
        }
        carvings_multiunion_solid->Voxelize();

        G4SubtractionSolid* dimpled_plate_solid = new G4SubtractionSolid(plate_name, 
                                                                        whole_plate_solid, carvings_multiunion_solid);
        geometry_solid = dynamic_cast<G4VSolid*>(dimpled_plate_solid);
    }
    else{
        geometry_solid = dynamic_cast<G4VSolid*>(whole_plate_solid);
    }

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    if(mpt_ && pvt){
      pvt->SetMaterialPropertiesTable(mpt_);
    }
    else{
      G4Exception("[WLSPlate]", "Construct()", JustWarning,
      "The optical properties of the WLS Plate were not set");
    }

    G4LogicalVolume* geometry_logic =
      new G4LogicalVolume(geometry_solid, pvt, plate_name);

    if(world_logic_vol){
        new G4PVPlacement(nullptr, G4ThreeVector{}, geometry_logic, plate_name, 
                        world_logic_vol, false, 0, true);
    }
    else{
        this->SetLogicalVolume(geometry_logic);
    }
    return;
  }

  G4ThreeVector WLSPlate::GenerateVertex(const G4String&) const
  {
    return G4ThreeVector(0., 100.*mm, 0.);
  }
}