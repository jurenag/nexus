#include "DimpleTest.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "LArSphere.h"
#include "SensorSD.h"
#include "Visibilities.h"
#include "HamamatsuS133606050VE.h"

#include <G4VSolid.hh>
#include <G4Tubs.hh>
#include <G4Orb.hh>
#include <G4Box.hh>
#include <G4MultiUnion.hh>
#include <G4Transform3D.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4SubtractionSolid.hh>
#include <G4PVPlacement.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SDManager.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4GenericMessenger.hh>
#include <G4SystemOfUnits.hh>

using namespace nexus;

REGISTER_CLASS(DimpleTest, GeometryBase)

namespace nexus{

  DimpleTest::DimpleTest():
  GeometryBase(),
  world_radius_             (1.*m),
  how_many_dimples_         (1),
  gap_in_between_dimples_   (1.*mm),
  dimple_type_              ("cylindrical"),
  flat_dimple_width_        (8.*mm),
  flat_dimple_depth_        (1.2*mm),
  curvy_dimple_radius_      (1.5*mm),
  plate_dx_                 (5.*mm),
  plate_dy_                 (3.5*mm),
  plate_dz_                 (5.*mm),
  ref_phsensor_support_     (true),
  gap_                      (0.5*mm),
  with_board_               (true),
  board_thickn_             (2.*mm),
  with_real_sipm_           (true),
  enclosure_gap_            (0.5*mm),
  mpt_                      (opticalprops::EJ286())
  {

    msg_ = new G4GenericMessenger(this, "/Geometry/DimpleTest/",
				"Control commands of geometry DimpleTest.");

    G4GenericMessenger::Command& wr_cmd =
      msg_->DeclareProperty("world_radius", world_radius_,
			    "LAr sphere radius.");
    wr_cmd.SetUnitCategory("Length");
    wr_cmd.SetParameterName("world_radius", false);
    wr_cmd.SetRange("world_radius>0.");

    G4GenericMessenger::Command& hmd_cmd =
      msg_->DeclareProperty("how_many_dimples", how_many_dimples_,
			    "How many dimples to carve in front of the photosensor.");
    hmd_cmd.SetParameterName("how_many_dimples", false);
    hmd_cmd.SetRange("how_many_dimples>=0");

    G4GenericMessenger::Command& gibd_cmd =
      msg_->DeclareProperty("gap_in_between_dimples", gap_in_between_dimples_,
			    "Gap in between two consecutive dimples.");
    gibd_cmd.SetUnitCategory("Length");
    gibd_cmd.SetParameterName("gap_in_between_dimples", false);
    gibd_cmd.SetRange("gap_in_between_dimples>=0.");
    
    G4GenericMessenger::Command& dty_cmd =
      msg_->DeclareProperty("dimple_type", dimple_type_,
			    "Dimple type. Might be 'flat', 'cylindrical' or 'spherical'.");

    G4GenericMessenger::Command& fdw_cmd =
      msg_->DeclareProperty("flat_dimple_width", flat_dimple_width_,
			    "Width of the flat dimples.");
    fdw_cmd.SetUnitCategory("Length");
    fdw_cmd.SetParameterName("flat_dimple_width", false);
    fdw_cmd.SetRange("flat_dimple_width>0.");

    G4GenericMessenger::Command& fdd_cmd =
      msg_->DeclareProperty("flat_dimple_depth", flat_dimple_depth_,
			    "Depth of the flat dimples.");
    fdd_cmd.SetUnitCategory("Length");
    fdd_cmd.SetParameterName("flat_dimple_depth", false);
    fdd_cmd.SetRange("flat_dimple_depth>0.");

    G4GenericMessenger::Command& cdr_cmd =
      msg_->DeclareProperty("curvy_dimple_radius", curvy_dimple_radius_,
			    "Radius of the cylindrical or spherical dimples.");
    cdr_cmd.SetUnitCategory("Length");
    cdr_cmd.SetParameterName("curvy_dimple_radius", false);
    cdr_cmd.SetRange("curvy_dimple_radius>0.");

    G4GenericMessenger::Command& dx_cmd =
      msg_->DeclareProperty("length", plate_dx_,
			    "Length of the WLS plate.");
    dx_cmd.SetUnitCategory("Length");
    dx_cmd.SetParameterName("length", false);
    dx_cmd.SetRange("length>0.");

    G4GenericMessenger::Command& dy_cmd =
      msg_->DeclareProperty("thickn", plate_dy_,
			    "Thickness of the WLS plate.");
    dy_cmd.SetUnitCategory("Length");
    dy_cmd.SetParameterName("thickn", false);
    dy_cmd.SetRange("thickn>0.");

    G4GenericMessenger::Command& dz_cmd =
      msg_->DeclareProperty("width", plate_dz_,
			    "Width of the WLS plate.");
    dz_cmd.SetUnitCategory("Length");
    dz_cmd.SetParameterName("width", false);
    dz_cmd.SetRange("width>0.");

    G4GenericMessenger::Command& g_cmd =
      msg_->DeclareProperty("gap", gap_,
			    "Gap between the plane surface of the WLS plate and the SiPM.");
    g_cmd.SetUnitCategory("Length");
    //g_cmd.SetParameterName("gap", false);
    //g_cmd.SetRange("gap>0.");    // These are commented so that gap_ can help modelate the immersion of the SiPMs into the flat dimple

    G4GenericMessenger::Command& wd_cmd =
      msg_->DeclareProperty("ref_phsensor_support", ref_phsensor_support_,
			    "Whether support of the SiPM is reflective.");

    G4GenericMessenger::Command& wb_cmd =
      msg_->DeclareProperty("with_board", with_board_,
			    "Whether there is a reflective board behind the SiPM.");

    G4GenericMessenger::Command& bt_cmd =
      msg_->DeclareProperty("board_thickn", board_thickn_,
			    "Thickness of the SiPM board or the artificial detector.");
    bt_cmd.SetUnitCategory("Length");
    bt_cmd.SetParameterName("board_thickn", false);
    bt_cmd.SetRange("board_thickn>0.");

    G4GenericMessenger::Command& bh_cmd =
      msg_->DeclareProperty("board_height", board_height_,
			    "Height of the SiPM board or the artificial detector.");
    bh_cmd.SetUnitCategory("Length");
    bh_cmd.SetParameterName("board_height", false);
    bh_cmd.SetRange("board_height>0.");

    G4GenericMessenger::Command& wrs_cmd =
      msg_->DeclareProperty("with_real_sipm", with_real_sipm_,
			    "Whether to place a real SiPM or an artificial piece-wise detector.");

    G4GenericMessenger::Command& abiee_cmd =
      msg_->DeclareProperty("artificial_board_is_efficient_everywhere", artificial_board_is_efficient_everywhere_,
			    "Whether the artificial detector is perfectly efficient everywhere or it has a reflective chunk out of the virtual sipm boundaries.");

    G4GenericMessenger::Command& eg_cmd =
      msg_->DeclareProperty("enclosure_gap", enclosure_gap_,
			    "Gap between the reflective enclosure and the span of the geometry that is within the enclosure.");
    eg_cmd.SetUnitCategory("Length");
    eg_cmd.SetParameterName("enclosure_gap", false);
    eg_cmd.SetRange("enclosure_gap>0.");

  }

  DimpleTest::~DimpleTest()
  {
    if(msg_) delete msg_;
  }

    // KEEP REVIEWING THE CODE AS OF HERE!!!!!!!!!!!!!!!!!!!!!!!

  void DimpleTest::Construct()
  {

    if(this->geometry_is_ill_formed()){
      G4Exception("[DimpleTest]", "Construct()",
                  FatalException, "Geometry test failed.");
    }

    // VACUUM CAPSULE
    const G4String world_name = "VACUUM_CAPSULE";
    G4Orb* world_solid =
        new G4Orb(world_name, world_radius_+1.*m);
    G4Material* vacuum =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4LogicalVolume* world_logic = 
        new G4LogicalVolume(world_solid, vacuum, world_name, 0, 0, 0, true);
    world_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(world_logic);

    // LAR SPHERE
    G4String sphere_name = "LAR_SPHERE";
    G4Orb* sphere_solid = new G4Orb(sphere_name, world_radius_);
    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    lAr->SetMaterialPropertiesTable(opticalprops::LAr());
    G4LogicalVolume* sphere_logic = new G4LogicalVolume(sphere_solid, lAr, sphere_name);
    sphere_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    // Place the LAr sphere within the vacuum world
    G4VPhysicalVolume* mother_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement (new G4RotationMatrix(), G4ThreeVector(0., 0., 0.), sphere_logic, 
                        sphere_name, world_logic, false, 0, false));

    ConstructWLSPlate(mother_physical);
    if(with_real_sipm_) ConstructSiPM(mother_physical);
    else ConstructArtificialDetector(mother_physical);
    ConstructReflectiveEnclosure(mother_physical);
    
    return;

  }

  void DimpleTest::ConstructWLSPlate(G4VPhysicalVolume* mother_physical)
  {
    const G4String plate_name = "WLS_PLATE";

    G4VSolid* geometry_solid = nullptr;
    G4Box* whole_plate_solid =
      new G4Box(plate_name, plate_dx_/2., plate_dy_/2., plate_dz_/2.);

    G4double tolerance = 0.5*mm; // To avoid boolean subtraction of matching surfaces
    
    if(how_many_dimples_>0){

        G4double dimple_z_span;
        G4RotationMatrix* rot = new G4RotationMatrix();
        G4VSolid* carving_solid = nullptr;
        if(dimple_type_=="flat"){
            // Yes, using 2*flat_dimple_depth_ as the whole dimension of the carving along the z-axis is ok 
            // (the carvings are subtracted from the very edge of the plate)
            carving_solid = dynamic_cast<G4VSolid*>(new G4Box("AUX", flat_dimple_depth_, plate_dy_/2. +tolerance, flat_dimple_width_/2.));
            dimple_z_span = flat_dimple_width_;
        }
        else if(dimple_type_=="cylindrical"){
            carving_solid = dynamic_cast<G4VSolid*>(
                                        new G4Tubs("AUX", 0., curvy_dimple_radius_, plate_dy_/2. +tolerance, 0., twopi));
            dimple_z_span = 2*curvy_dimple_radius_;
            rot->rotateX(+90.*deg);
        }
        else{// Default is spherical dimple
            carving_solid = dynamic_cast<G4VSolid*>(new G4Orb("AUX", curvy_dimple_radius_));
            dimple_z_span = 2*curvy_dimple_radius_;
        }

        G4Transform3D* transform_ptr = nullptr;
        G4double z_pos = (-0.5*((how_many_dimples_*dimple_z_span)+((how_many_dimples_-1)*gap_in_between_dimples_))) +(0.5*dimple_z_span);
        G4MultiUnion* carvings_multiunion_solid = new G4MultiUnion("CARVINGS");
        for(G4int i=0; i<how_many_dimples_; i++){
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector(0., 0., z_pos));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
            z_pos += dimple_z_span + gap_in_between_dimples_;
        }
        carvings_multiunion_solid->Voxelize();

        G4SubtractionSolid* dimpled_plate_solid = new G4SubtractionSolid(   plate_name, whole_plate_solid, carvings_multiunion_solid, 
                                                                            nullptr, G4ThreeVector(plate_dx_/2., 0., 0.));
        geometry_solid = dynamic_cast<G4VSolid*>(dimpled_plate_solid);
    }
    else{
        geometry_solid = dynamic_cast<G4VSolid*>(whole_plate_solid);
    }

    G4Material* plate_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    if(mpt_){
    plate_mat->SetMaterialPropertiesTable(mpt_);
    }
    else{
    G4Exception("[DimpleTest]", "ConstructWLSPlate()", JustWarning,
    "The optical properties of the WLS plate substrate were not set");
    }

    G4LogicalVolume* plate_logic =
      new G4LogicalVolume(geometry_solid, plate_mat, plate_name);

    /*
    G4VisAttributes plate_col = nexus::RedAlpha();
    plate_col.SetForceSolid(true);
    plate_logic->SetVisAttributes(plate_col);
    */

    G4VPhysicalVolume* plate_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0), plate_name, plate_logic,
        mother_physical, false, 0, false));

    return;
  }

  void DimpleTest::ConstructSiPM(G4VPhysicalVolume* mother_physical)
  {

    HamamatsuS133606050VE sipm;
    sipm.SetReflectiveSupports(ref_phsensor_support_);
    sipm.Construct();
    G4double sipm_thickn = sipm.GetThickness();
    G4LogicalVolume* sipm_logic_vol = sipm.GetLogicalVolume();

    G4VisAttributes sipm_col = nexus::Red();
    sipm_col.SetForceSolid(true);
    sipm_logic_vol->SetVisAttributes(sipm_col);

    if (!sipm_logic_vol) {
      G4Exception("[DimpleTest]", "ConstructSiPM()",
                  FatalException, "Null pointer to logical volume.");
    }

    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateZ(-90.*degree);
    G4ThreeVector pos((plate_dx_+sipm_thickn)/2.+gap_, 0., 0.);
    new G4PVPlacement(  rot, pos, "S133606050VE_MPPC", sipm_logic_vol,
                        mother_physical, false, 0, true);

    if(with_board_){
        const G4String board_name = "BOARD";
        G4Box* board_solid = new G4Box(board_name, board_thickn_/2., board_height_/2., plate_dz_/2.);
        G4LogicalVolume* board_logic = new G4LogicalVolume(board_solid, materials::FR4(), board_name);

        // Set its color for visualization purposes
        G4VisAttributes board_col = nexus::WhiteAlpha();
        board_col.SetForceSolid(true);
        board_logic->SetVisAttributes(board_col);

        //Now create the reflective optical surface
        const G4String refsurf_name = "REF_SURFACE";
        G4OpticalSurface* refsurf_opsurf = 
        new G4OpticalSurface(refsurf_name, unified, ground, dielectric_metal, 1);
    
        refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
        new G4LogicalSkinSurface("BOARD_REF_SURFACE", board_logic, refsurf_opsurf);   

        new G4PVPlacement(  nullptr, G4ThreeVector(plate_dx_/2. +gap_ +sipm_thickn +board_thickn_/2., 0., 0.),
                            board_name, board_logic,
                            mother_physical, false, 0, true);
    }
    return;
  }

  void DimpleTest::ConstructArtificialDetector(G4VPhysicalVolume* mother_physical)
  {

    // Say that the upper view of the setup with a real sipm is the following
    //
    //  |                |
    //  |   WLS plate    |
    //  |                |
    //  ------------------
    //        ------
    //  ------+----+------   
    //  |     Board      |
    //
    // Where you can see the WLS plate, the SiPM and the SiPM-mounting board
    // Then, in the configuration with the artificial detector what you have is
    //
    //  |                |
    //  |   WLS plate    |
    //  |                |
    //  ------------------
    //        
    //  ------------------   
    //  |  R  |  E |  R  |
    //
    // Where the "R" chunks are vikuiti coated, and the "E" chunk is coated
    // with a perfectly efficient detector. In this way, the analysis integrates
    // the effect of photons that are reflected in the reflective SiPM-mounting
    // boards for the eventual X-ARAPUCA geometry, but also allowing a flat 
    // screen over which we will study the photon impacts. The transverse dimensions
    // of the "E" chunk are those of the actual SiPM

    HamamatsuS133606050VE sipm;
    G4double sipm_height = sipm.GetTransverseDim();
    G4double sipm_width = sipm.GetTransverseDim();

    const G4String central_chunk_name = "CENTRAL_CHUNK";
    const G4String ext_chunk_name = "EXT_CHUNK";

    G4Box* central_chunk_solid        =   new G4Box(central_chunk_name, board_thickn_/2., sipm_height/2., sipm_width/2.);
    G4Box* central_chunk_subtrahend   =   new G4Box("INNER_AUX", board_thickn_, sipm_height/2., sipm_width/2.);
    G4Box* ext_chunk_minuend      =   new G4Box("OUTTER_AUX", board_thickn_/2., board_height_/2., plate_dz_/2.);
    G4SubtractionSolid* ext_chunk_solid = new G4SubtractionSolid( ext_chunk_name, ext_chunk_minuend, central_chunk_subtrahend, 
                                                                nullptr, G4ThreeVector(0., 0., 0.));

    G4LogicalVolume* central_chunk_logic = new G4LogicalVolume(central_chunk_solid, materials::FR4(), central_chunk_name);
    G4LogicalVolume* ext_chunk_logic = new G4LogicalVolume(ext_chunk_solid, materials::FR4(), ext_chunk_name);

    // Set its color for visualization purposes
    G4VisAttributes central_chunk_col = nexus::LightBlueAlpha();
    central_chunk_col.SetForceSolid(true);
    central_chunk_logic->SetVisAttributes(central_chunk_col);

    G4VisAttributes ext_chunk_col = nexus::LillaAlpha();
    ext_chunk_col.SetForceSolid(true);
    ext_chunk_logic->SetVisAttributes(ext_chunk_col);

    //Create the efficient optical surface
    const G4String e_surf_name = "EFFICIENT_SURFACE";
    G4OpticalSurface* e_surf_opsurf = 
    new G4OpticalSurface(e_surf_name, unified, polishedfrontpainted, dielectric_dielectric, 1);
    e_surf_opsurf->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());
    new G4LogicalSkinSurface(e_surf_name, central_chunk_logic, e_surf_opsurf);   

    if(!artificial_board_is_efficient_everywhere_){
        //Create the reflective optical surface
        const G4String r_surf_name = "REFLECTIVE_SURFACE";
        G4OpticalSurface* r_surf_opsurf = 
        new G4OpticalSurface(r_surf_name, unified, ground, dielectric_metal, 1);
        r_surf_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
        new G4LogicalSkinSurface(r_surf_name, ext_chunk_logic, r_surf_opsurf);  
    }
    else{
        new G4LogicalSkinSurface(e_surf_name, ext_chunk_logic, e_surf_opsurf);  
    }

    // Place both chunks
    new G4PVPlacement(  nullptr, G4ThreeVector(plate_dx_/2. +gap_ +board_thickn_/2., 0., 0.),
                        central_chunk_name, central_chunk_logic,
                        mother_physical, false, 0, true);

    new G4PVPlacement(  nullptr, G4ThreeVector(plate_dx_/2. +gap_ +board_thickn_/2., 0., 0.),
                        ext_chunk_name, ext_chunk_logic,
                        mother_physical, false, 0, true);
    
    return;
  }

  void DimpleTest::ConstructReflectiveEnclosure(G4VPhysicalVolume* mother_physical)
  {
    // The reflective enclosure is designed so as to tightly fit the lateral faces of the WLS 
    // plate. One can think of this unrealistic reflection on the lateral faces of the WLS
    // plate as a way to virtually extend the WLS plate. In this context, a realistic simulation
    // of the X-ARAPUCA with this "chunk" of X-ARAPUCA would entail setting plate_dx_ to L/dimples_no
    // (where L is the actual length of the WLS plate within a supercell and dimple_no is the number
    // of dimples, and plate_dz_ to the M/2, where M is the width of the WLS plate within a supercell.
    // Indeed, using the symmetry of one supercell we can collapse it all onto one chunk of supercell.
    // The reflections of this limited geometry generate the whole supercell. In this interpretation,
    // this chunk is the span of the whole supercell (except for the border chunks, which would
    // introduce some border effects)

    HamamatsuS133606050VE sipm;
    G4double sipm_thickn = sipm.GetThickness();

    // The reflective enclosure thickness is that of the board
    G4double enclosure_thickn = board_thickn_;

    // Compute the dimensions of the enclosure
    G4double enclosure_length = enclosure_thickn +plate_dx_ +gap_ +enclosure_thickn; //X
    if(with_real_sipm_){
        enclosure_length += sipm_thickn;
    }

    G4double enclosure_height   = board_height_+(2.*enclosure_thickn); //Y
    G4double enclosure_width    = plate_dz_ + 2.*enclosure_thickn; //Z
    // The reflective enclosure tightly fits the lateral faces of the WLS plate.  


    const G4String enclosure_name = "REF_ENCLOSURE";
    G4Box* enclosure_minuend = new G4Box("AUX_OUTTER", enclosure_length/2., enclosure_height/2., enclosure_width/2.);
    G4Box* enclosure_subtrahend = new G4Box("AUX_INNER",    enclosure_length/2. -board_thickn_/2., // The reflective board serves a cover, 
                                                                                                   // so here we just have to subtract the board_thickn one time
                                                            enclosure_height/2. -board_thickn_,
                                                            enclosure_width/2.  -board_thickn_);

    G4SubtractionSolid* enclosure_solid = new G4SubtractionSolid(   enclosure_name, enclosure_minuend, 
                                                                    enclosure_subtrahend, nullptr, 
                                                                    G4ThreeVector(board_thickn_/2., 0., 0.));

    G4LogicalVolume* enclosure_logic = new G4LogicalVolume(enclosure_solid, materials::FR4(), enclosure_name);
    
    /*
    // Set its color for visualization purposes
    G4VisAttributes board_col = nexus::WhiteAlpha();
    board_col.SetForceSolid(true);
    board_logic->SetVisAttributes(board_col);*/

    //Now create the reflective optical surface
    const G4String refsurf_name = "REF_SURFACE";
    G4OpticalSurface* refsurf_opsurf = 
    new G4OpticalSurface(refsurf_name, unified, ground, dielectric_metal, 1);
    
    refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
    new G4LogicalSkinSurface("ENCLOSURE_REF_SURFACE", enclosure_logic, refsurf_opsurf);   

    G4double x_pos = gap_/2.;
    if(with_real_sipm_){
        x_pos += sipm_thickn/2.;
    }

    new G4PVPlacement(  nullptr, G4ThreeVector(x_pos, 0., 0.),
                        enclosure_name, enclosure_logic,
                        mother_physical, false, 0, true);
    return;
  }

  G4bool DimpleTest::geometry_is_ill_formed() const
  {
    if(plate_dy_>= board_height_) return true;
    return false;
  }

  G4ThreeVector DimpleTest::GenerateVertex(const G4String&) const
  {
    G4double generation_vertex_height = ((plate_dy_/2.)+(board_height_/2.))/2.;
    return G4ThreeVector(0., generation_vertex_height, 0.);
  }
}