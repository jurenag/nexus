#include "WLSPlateTest.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "WLSPlate.h"

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

REGISTER_CLASS(WLSPlateTest, GeometryBase)

namespace nexus{

  WLSPlateTest::WLSPlateTest():
  GeometryBase(),
  world_dx_(10.*m),
  world_dy_(12.*mm),
  world_dz_(10.*m),
  //wlsp_dx_(487.*mm),
  wlsp_dx_(10.*cm),
  wlsp_dy_(3.5*mm),
  //wlsp_dz_(93.*mm),
  //wlsp_dz_(50.*mm),   //<Default WLS plate dimensions taken from FD1 TDR vol. IX
  wlsp_dz_(10.*cm),
  gap_(1.*mm),
  //gap_(0.),
  upper_collector_(true),
  lower_collector_(true),
  collector_thickn_(1.*mm),
  upper_cover_(false),
  lower_cover_(false),
  cover_thickn_(1.5*mm),
  mpt_(opticalprops::EJ286())
  {
  }

  WLSPlateTest::~WLSPlateTest()
  {
  }

  void WLSPlateTest::Construct()
  {

    if(GeometryIsIllFormed()){
        G4Exception("[WLSPlateTest]", "Construct()", FatalException,
        "Geometry is ill-formed.");
    }

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
    if(upper_collector_ || lower_collector_) ConstructCollectors(world_logic, upper_collector_, lower_collector_);
    if(upper_cover_ || lower_cover_) ConstructCovers(world_logic);
    ConstructWLSPlate(world_logic);
    return;
  }

  void WLSPlateTest::ConstructCollectors(G4LogicalVolume* world_logic_vol, G4bool upper_collector, G4bool lower_collector)
  {

    const G4String collector_name = "COLLECTOR";

    G4Box* collector_solid =
      new G4Box(collector_name, world_dx_/2., collector_thickn_/2., world_dz_/2.);

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    

    G4LogicalVolume* collector_logic =
      new G4LogicalVolume(collector_solid, pvt, collector_name);

    //  Detection coating
    G4OpticalSurface* collector_opsurf =
      new G4OpticalSurface("COLLECTOR_OPSURF", unified, polished, dielectric_metal);
    collector_opsurf->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());
    new G4LogicalSkinSurface("COLLECTOR_OPSURF", collector_logic, collector_opsurf);

    if(world_logic_vol){
        if(lower_collector){

            G4double lower_y_pos = (-1.*world_dy_/2.)+(collector_thickn_/2.);
            new G4PVPlacement(nullptr, G4ThreeVector(0., lower_y_pos, 0.), 
                            collector_logic, "LOWER_COLLECTOR", world_logic_vol, false, 1, true);        
        }
        if(upper_collector){

            G4double upper_y_pos = -1.*world_dy_/2.;
            if(lower_collector_) upper_y_pos += collector_thickn_;
            upper_y_pos += gap_;
            if(lower_cover_) upper_y_pos += cover_thickn_;
            upper_y_pos +=  wlsp_dy_;
            if(upper_cover_) upper_y_pos += cover_thickn_;
            upper_y_pos += gap_ +(collector_thickn_/2.);

            new G4PVPlacement(nullptr, G4ThreeVector(0., upper_y_pos, 0.), 
                            collector_logic, collector_name, world_logic_vol, false, 0, true);
        }
    }
    else{
      G4Exception("[WLSPlateTest]", "ConstructCollectors()", FatalException,
      "World logical volume not found.");
    }
    return;   
  }

  void WLSPlateTest::ConstructCovers(G4LogicalVolume* world_logic_vol)
  {
      // The aim of this cover is to enforce photons to scape the WLS plate via its lateral
      // faces. If we use the usual generator (throwing photons from the LAr gap that is
      // found above the WLS plate, with a direction that is perpendicular to the WLS plate
      // bigger face), this UpperCover would prevent photons to reach the WLS plate in the
      // first place. To prevent this, there are different solutions. One of them (1) comprises
      // changing the generator vertex. Maybe throwing photons from underneath the WLS plate, from
      // one side, or even generate them within the WLS plate. To parse the WLShifted photons, this
      // this should not be a drawback. However, since I want to compare this results to another 
      // results that I got with the current generator, I want change the conditions as little as
      // possible, so I'm going to stick to the current generator, including its generator vertex.
      // I suggest the following solution. The upper cover is a plate with an inner notch in one of
      // its surfaces. The surface with the notch faces the WLS plate. The cover is placed so that
      // the generator vertex falls within the notch.

    const G4String cover_name = "COVER";
    G4double notch_thickn = 1.25*mm;
    G4double hole_radius = 0.1*mm;

    G4Box* box_solid =
      new G4Box("AUX_1", wlsp_dx_/2., cover_thickn_/2., wlsp_dz_/2.);

    G4Tubs* hole_solid = 
        new G4Tubs("AUX_2", 0., hole_radius, cover_thickn_/2., 0., twopi);

    G4RotationMatrix* aux_rot = new G4RotationMatrix{};
    aux_rot->rotateX(90.*degree);

    G4ThreeVector aux_vec = G4ThreeVector(0., -1.*(cover_thickn_-notch_thickn), 0.);

    G4SubtractionSolid* cover_solid = 
        new G4SubtractionSolid(cover_name, box_solid, hole_solid, 
                                aux_rot, aux_vec);

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_lXe");
    pvt->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());

    G4LogicalVolume* cover_logic =
      new G4LogicalVolume(cover_solid, pvt, cover_name);

    G4OpticalSurface* cover_opsurf =
      new G4OpticalSurface("COVER_OPSURF", unified, polished, dielectric_metal);
    cover_opsurf->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());
    new G4LogicalSkinSurface("COVER_OPSURF", cover_logic, cover_opsurf);


    if(world_logic_vol){
        G4double minigap = 0.05*mm;
        if(lower_cover_){
            G4double lower_y_pos = -1.*world_dy_/2.;
            if(lower_collector_) lower_y_pos += collector_thickn_;
            lower_y_pos += gap_ +(cover_thickn_/2.) -minigap;

            new G4PVPlacement(nullptr, G4ThreeVector(0., lower_y_pos, 0.), 
                        cover_logic, cover_name, world_logic_vol, false, 0, true);
        }
        if(upper_cover_){
            G4double upper_y_pos = -(world_dy_/2.);
            if(lower_collector_) upper_y_pos += collector_thickn_;
            upper_y_pos += gap_;
            if(lower_cover_) upper_y_pos += cover_thickn_;
            upper_y_pos +=  wlsp_dy_ +(cover_thickn_/2.) +minigap;

            new G4PVPlacement(nullptr, G4ThreeVector(0., upper_y_pos, 0.), 
                            cover_logic, cover_name, world_logic_vol, false, 0, true);
        }
    }
    else{
        G4Exception("[WLSPlateTest]", "ConstructCovers()", FatalException,
        "World logical volume not found.");
    }
    return;
  }

  void WLSPlateTest::ConstructWLSPlate(G4LogicalVolume* world_logic_vol)
  {

    const G4String plate_name = "WLS_PLATE";

    WLSPlate* plate = new WLSPlate(wlsp_dx_, wlsp_dy_, wlsp_dz_, false);
    plate->SetOpticalProperties(opticalprops::EJ286());
    plate->Construct();
    G4LogicalVolume* plate_logic = plate->GetLogicalVolume();

    G4double y_pos = -1.*world_dy_/2.;
    if(lower_collector_) y_pos += collector_thickn_;
    y_pos += gap_;
    if(lower_cover_) y_pos += cover_thickn_;
    y_pos += wlsp_dy_/2.;

    if(world_logic_vol && plate_logic){
        new G4PVPlacement(nullptr, G4ThreeVector(0., y_pos, 0.), 
                        plate_logic, plate_name, world_logic_vol, false, 0, true);
    }
    else{
      G4Exception("[WLSPlateTest]", "Construct()", FatalException,
      "World logical volume not found.");
    }
    return;
  }

  G4bool WLSPlateTest::GeometryIsIllFormed()
  {
    G4double minimum_required_world_thickn = wlsp_dy_ + (2.*gap_);
    if(lower_collector_) minimum_required_world_thickn += collector_thickn_;
    if(lower_cover_) minimum_required_world_thickn += cover_thickn_;
    if(upper_cover_) minimum_required_world_thickn += cover_thickn_;
    if(upper_collector_) minimum_required_world_thickn += collector_thickn_;

    return world_dy_<minimum_required_world_thickn ? true : false;
  }

  G4ThreeVector WLSPlateTest::GenerateVertex(const G4String&) const
  {
    G4double arapuca_upper_surf_pos = (-1.*world_dy_/2.);+gap_+wlsp_dy_; 
    if (lower_collector_) arapuca_upper_surf_pos += collector_thickn_;
    arapuca_upper_surf_pos += gap_;
    if (lower_cover_) arapuca_upper_surf_pos += cover_thickn_;
    arapuca_upper_surf_pos += wlsp_dy_;

    G4double vertex_height_over_arapuca_surf = 1.*mm;

    return G4ThreeVector(0., arapuca_upper_surf_pos+vertex_height_over_arapuca_surf, 0.);
  }
}