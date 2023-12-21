#include "APEX.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "WLSPlate.h"
#include "SiPMMPPC.h"
#include "HamamatsuS133606050VE.h"
#include "HamamatsuS133605075HQR.h"
#include "FbkNuvHdCryoTT.h"
#include "BroadcomAFBRS4N44P044M.h"
#include "PerfectSiPMMPPC.h"
#include "SiPMBoard.h"
#include "RandomUtils.h"
#include "Visibilities.h"

#include <algorithm>
#include <random>
#include <cmath>
#include <G4GenericMessenger.hh>
#include <G4UserLimits.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4SubtractionSolid.hh>
#include <G4ThreeVector.hh>
#include <G4MultiUnion.hh>
#include <G4VisAttributes.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>

#include <CLHEP/Units/SystemOfUnits.h>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(APEX, GeometryBase)

namespace nexus{

  APEX::APEX():
  GeometryBase(), 
  surrounding_media_                    ("lar"                        ),
  MLS_thickn_                           (0.010  *mm                   ),
  MLS_rindex_                           (1.68                         ),
  coating_thickn_                       (3.226  *um                   ),  // Based on arxiv.org/abs/1912.09191 and TDR vol.IX, section 5.8.3.1,
                                                                          // the pTP film thickness is such that there's 400 micrograms of pTP
                                                                          // deposited over each square centimeter of DF.
                                                                          // This, together with the pTP density, (1.24g/cm3, found in
                                                                          // en.wikipedia.org/wiki/Terphenyl), gives a pTP film thickness of 
                                                                          // 3.226 micrometers
  coating_rindex_                       (1.65                         ),  
  remove_coating_                       (false                        ),
  remove_MLS_                           (false                        ),
  plate_length_                         (450.   *mm                   ),   ///X // SiPMs are placed long the short side, which is 450 mm long
  plate_thickn_                         (6.0    *mm                   ),   ///Y // The values for these dimensions were measured from an step
  plate_width_                          (495.    *mm                  ),   ///Z // file which we received from F. Cavanna
  WLSp_rindex_                          (1.502                        ),
  secondary_wls_attlength_              (1.     *m                    ),
  cromophore_concentration_             (40.                          ),

  reflective_foil_thickn_               (0.065  *mm                   ),   /// Got foil thickness from isoltronic.ch/assets/of-m-vikuiti-esr-app-guide.pdf
  SiPM_code_                            (1                            ),
  num_phsensors_                        (30                           ),   /// This is seemingly the APEX baseline
  board_position_code_                  (1                            ),
  gap_                                  (0.5    *mm                   ),
  ref_phsensors_supports_               (true                         ), 
  with_dimples_                         (false                         ),
  dimple_type_                          ("cylindrical"                ),
  flat_dimple_width_                    (6.1    *mm                   ),
  flat_dimple_depth_                    (2.     *mm                   ),
  curvy_dimple_radius_                  (1.5    *mm                   ),
  generation_region_                    ("random"                     ),
  gen_x_                                (0.     *cm                   ),
  gen_z_                                (0.     *cm                   ),
  gen_diameter_                         (1.*cm                        ),
  path_to_inwards_dichroic_data_        (""                           ),
  path_to_outwards_dichroic_data_       (""                           ),
  world_extra_thickn_                   (100.   *cm                   )
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/APEX/",
				"Control commands of geometry APEX.");

    G4GenericMessenger::Command& sm_cmd =
      msg_->DeclareProperty("surrounding_media", surrounding_media_,
			    "Which media to place the APEX in");

    G4GenericMessenger::Command& mlst_cmd =
      msg_->DeclareProperty("MLS_thickn", MLS_thickn_,
			    "Thickness of the DF multilayer structure (MLS) which is deposited on top of the WLS plate.");
    mlst_cmd.SetUnitCategory("Length");
    mlst_cmd.SetParameterName("MLS_thickn", false);
    mlst_cmd.SetRange("MLS_thickn>0.");

    G4GenericMessenger::Command& mlsr_cmd =
      msg_->DeclareProperty("MLS_rindex", MLS_rindex_,
			    "Effective refractive index of the multi-layer structure.");
    mlsr_cmd.SetParameterName("MLS_thickn", false);
    mlsr_cmd.SetRange("MLS_thickn>=1.");

    G4GenericMessenger::Command& ptpct_cmd =
      msg_->DeclareProperty("coating_thickn", coating_thickn_,
			    "Thickness of the coating layer that is deposited over the MLS.");
    ptpct_cmd.SetUnitCategory("Length");
    ptpct_cmd.SetParameterName("coating_thickn", false);
    ptpct_cmd.SetRange("coating_thickn>0.");

    G4GenericMessenger::Command& ptpcr_cmd =
      msg_->DeclareProperty("coating_rindex", coating_rindex_,
			    "Refractive index of the coating layer that is deposited over the MLS.");
    ptpcr_cmd.SetParameterName("coating_rindex", false);
    ptpcr_cmd.SetRange("coating_rindex>1.");

    G4GenericMessenger::Command& rc_cmd =
      msg_->DeclareProperty("remove_coating", remove_coating_,
			    "Whether to remove the coating layer that is deposited over the MLS.");

    G4GenericMessenger::Command& rmls_cmd =
      msg_->DeclareProperty("remove_MLS", remove_MLS_,
			    "Whether to remove the DF (the MLS) together with the coating layer that is deposited on top of it.");

    G4GenericMessenger::Command& pl_cmd =
      msg_->DeclareProperty("plate_length", plate_length_,
			    "Length of the WLS plate.");
    pl_cmd.SetUnitCategory("Length");
    pl_cmd.SetParameterName("plate_length", false);
    pl_cmd.SetRange("plate_length>0.");

    G4GenericMessenger::Command& pw_cmd =
      msg_->DeclareProperty("plate_width", plate_width_,
			    "Width of the WLS plate.");
    pw_cmd.SetUnitCategory("Length");
    pw_cmd.SetParameterName("plate_width", false);
    pw_cmd.SetRange("plate_width>0.");

    G4GenericMessenger::Command& pt_cmd =
      msg_->DeclareProperty("plate_thickn", plate_thickn_,
			    "Thickness of the WLS plate.");
    pt_cmd.SetUnitCategory("Length");
    pt_cmd.SetParameterName("plate_thickn", false);
    pt_cmd.SetRange("plate_thickn>0.");

    G4GenericMessenger::Command& wlspr_cmd =
      msg_->DeclareProperty("WLSp_rindex", WLSp_rindex_,
			    "Refractive index of the wavelength shifting plate.");
    wlspr_cmd.SetParameterName("WLSp_rindex", false);
    wlspr_cmd.SetRange("WLSp_rindex>=1.");

    G4GenericMessenger::Command& swlsal_cmd =
      msg_->DeclareProperty("secondary_wls_attlength", secondary_wls_attlength_,
			    "Attenuation length of the secondary WLShifter (the WLS plate), in case EJ286 is used.");
    swlsal_cmd.SetUnitCategory("Length");
    swlsal_cmd.SetParameterName("secondary_wls_attlength", false);
    swlsal_cmd.SetRange("secondary_wls_attlength>0.");

    G4GenericMessenger::Command& crco_cmd =
      msg_->DeclareProperty("cromophore_concentration", cromophore_concentration_,
			    "Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter (the WLS plate), in case G2P_FB118 is used.");
    crco_cmd.SetParameterName("cromophore_concentration", false);
    crco_cmd.SetRange("cromophore_concentration>0.");

    G4GenericMessenger::Command& rft_cmd =
      msg_->DeclareProperty("reflective_foil_thickn", reflective_foil_thickn_,
			    "Reflective foil thickness.");
    rft_cmd.SetUnitCategory("Length");
    rft_cmd.SetParameterName("reflective_foil_thickn", false);
    rft_cmd.SetRange("reflective_foil_thickn>0.");

    G4GenericMessenger::Command& sc_cmd =
      msg_->DeclareProperty("SiPM_code", SiPM_code_,
			    "Integer signalling which SiPM to construct.");
    sc_cmd.SetParameterName("SiPM_code", false);
    sc_cmd.SetRange("SiPM_code>=1");

    G4GenericMessenger::Command& np_cmd =
      msg_->DeclareProperty("num_phsensors", num_phsensors_,
			    "Number of photosensors per board.");
    np_cmd.SetParameterName("num_phsensors", false);
    np_cmd.SetRange("num_phsensors>=0");

    G4GenericMessenger::Command& bpc_cmd =
      msg_->DeclareProperty("board_position_code", board_position_code_,
			    "Integer signalling where to place the SiPM board.");
    bpc_cmd.SetParameterName("board_position_code", false);
    bpc_cmd.SetRange("board_position_code>=1");

    G4GenericMessenger::Command& g_cmd =
      msg_->DeclareProperty("gap", gap_,
			    "Gap between the photosensors and the WLS plate. A negative gap can help modelate the immersion of the SiPMs into the dimples. Be careful not to collide the SiPMs into the plate.");
    g_cmd.SetUnitCategory("Length");
    //g_cmd.SetParameterName("gap", false);
    //g_cmd.SetRange("gap>0.");    // These are commented so that gap_ can help modelate the immersion of the SiPMs into the flat dimple

    G4GenericMessenger::Command& rps_cmd =
      msg_->DeclareProperty("ref_phsensors_supports", ref_phsensors_supports_,
			    "Whether photosensors supports are reflective.");

    G4GenericMessenger::Command& wd_cmd =
      msg_->DeclareProperty("with_dimples", with_dimples_,
			    "Whether the plate has carved dimples on it.");

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

    G4GenericMessenger::Command& gr_cmd =
      msg_->DeclareProperty("generation_region", generation_region_,
			    "Where to place the generation vertex.");

    G4GenericMessenger::Command& gx_cmd =
      msg_->DeclareProperty("gen_x", gen_x_,
			    "Average X-coordinate of the generation vertex if generation_region_=='custom'.");
    gx_cmd.SetUnitCategory("Length");
          
    G4GenericMessenger::Command& gz_cmd =
      msg_->DeclareProperty("gen_z", gen_z_,
			    "Average Z-coordinate of the generation vertex if generation_region_=='custom'.");
    gz_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command& gd_cmd =
      msg_->DeclareProperty("gen_diameter", gen_diameter_,
			    "Diameter of the circle where the generation vertex could be randomly sampled if generation_region_=='custom' is True.");
    gd_cmd.SetUnitCategory("Length");
    gd_cmd.SetParameterName("gen_diameter", false);
    gd_cmd.SetRange("gen_diameter>0.");

    G4GenericMessenger::Command& ptidd_cmd =
      msg_->DeclareProperty("path_to_inwards_dichroic_data", path_to_inwards_dichroic_data_,
			    "Absolute path to the dichroic data file that is to be sampled for the light trying to enter the WLS plate.");

    G4GenericMessenger::Command& ptodd_cmd =
      msg_->DeclareProperty("path_to_outwards_dichroic_data", path_to_outwards_dichroic_data_,
			    "Absolute path to the dichroic data file that is to be sampled for the light trying to escape the WLS plate.");


    // When testing WLS plates with opticalprops::noAbsLength_, it is possible that a photon
    // gets trapped within the plate (below the critical angle) into an infinite-bouncing-loop.
    // For the case of an EJ286 plate with DUNE supercells dimensions, with this absorption 
    // length and immersed into LAr, some analysis showed that particles that did not fall 
    // into this infinite loop had track lengths of less than one hundred meters.
    ul_ = new G4UserLimits();
    ul_->SetUserMaxTrackLength(100*m);
  }

  APEX::~APEX()
  {
      if(ul_) delete ul_;
  }


  void APEX::Construct()
  {

    // Compute internal attributes  // COME BACK HERE AND COMPUTE THESE !!!!
    // overall_length_ = ... ;  
    // overall_thickn_ = ... ;
    // overall_width_  = ... ;

    board_length_ = plate_length_;

    if(GeometryIsIllFormed()){
      G4Exception("[APEX]", "Construct()", FatalException,
      "The given dimensions do not describe a feasible APEX.");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    // The biggest volume is a vacuum box with the dimensions of the APEX device
    // plus 2*world_extra_thickn_, for each dimension. This volume is the effective 
    // world volume FOR THE NEXUS USER, but it is not the world volume of the overall
    // Geant4 application (afterwards nexus takes the user's biggest volume 
    // and places it inside another world volume whose dimensions are enough so as to 
    // fit the whole span of the biggest volume you implemented). The nexus user 
    // seems not to have access to the physical placement of the biggest volume (in our 
    // case, the vacuum box) which is implemented in line 80 of source/base/DetectorConstruction.cc, 
    // AFTER calling GeometryBase::Construct() (i.e. your geometry must be constructed 
    // before placing your biggest volume in nexus world volume,so you cannot possibly 
    // have access to the physical placement of your biggest volume when you Construct() it.).
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    // VACUUM CAPSULE
    const G4String world_name = "VACUUM_CAPSULE";

    G4Box* world_solid =
        new G4Box(world_name,
                (overall_length_/2.)+world_extra_thickn_,
                (overall_thickn_/2.)+world_extra_thickn_,
                (overall_width_/2.) +world_extra_thickn_);
                
    G4Material* vacuum =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

    G4LogicalVolume* world_logic = 
        new G4LogicalVolume(world_solid, vacuum, world_name, 0, 0, 0, true);
    world_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    this->SetLogicalVolume(world_logic);

    // surrounding_media_ box that contains all other volumes.

    G4String sm_name;
    G4MaterialPropertiesTable* mpt_ptr;
    if(surrounding_media_=="gar"){
      sm_name = "G4_Ar";
      mpt_ptr = opticalprops::GAr(10000/MeV);
    }
    else if(surrounding_media_=="air"){
      sm_name = "G4_AIR";
      mpt_ptr = opticalprops::Air();
    }
    else{
      sm_name = "G4_lAr";
      mpt_ptr = opticalprops::LAr();
    }

    const G4String sm_box_name = sm_name+"_BOX";

    G4Box* sm_box_solid =
      new G4Box(sm_box_name,
                (overall_length_+world_extra_thickn_)/2.,
                (overall_thickn_+world_extra_thickn_)/2.,
                (overall_width_ +world_extra_thickn_)/2.);

    G4Material* sm_material = G4NistManager::Instance()->FindOrBuildMaterial(sm_name);
    sm_material->SetMaterialPropertiesTable(mpt_ptr);

    G4LogicalVolume* sm_box_logic =
      new G4LogicalVolume(sm_box_solid, sm_material, sm_box_name);
    sm_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VPhysicalVolume* mother_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(new G4RotationMatrix(), 
                          G4ThreeVector(0., 0., 0.), 
                          sm_box_logic, 
                          sm_box_name, 
                          world_logic, 
                          false, 0, true));

    ConstructWLSPlate(mother_physical);
    ConstructSiPMSAndBoard(mother_physical);
    ConstructReflectiveFoil(mother_physical);         
    if(!remove_MLS_) ConstructDichroicFilter(mother_physical);


    return;
  }

  void APEX::ConstructWLSPlate(G4VPhysicalVolume* mother_physical) const
  { 
    WLSPlate* plate = new WLSPlate  (plate_length_, plate_thickn_, plate_width_, 
                                    opticalprops::G2P_FB118(cromophore_concentration_, WLSp_rindex_, true), false,
                                    //opticalprops::EJ286(secondary_wls_attlength_), false,
                                    with_dimples_, dimple_type_, num_phsensors_, 
                                    false, flat_dimple_width_, 
                                    flat_dimple_depth_, curvy_dimple_radius_);
    plate->Construct();
    G4LogicalVolume* plate_logic = plate->GetLogicalVolume();
    plate_logic->SetUserLimits(ul_);

    
    G4VisAttributes wlsp_col = nexus::LightBlueAlpha();
    wlsp_col.SetForceSolid(true);
    plate_logic->SetVisAttributes(wlsp_col);
    

    if (!plate_logic) {
      G4Exception("[APEX]", "ConstructWLSPlate()",
                  FatalException, "Null pointer to logical volume.");
    }

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), plate_logic->GetName(), 
                      plate_logic, mother_physical, false, 0, true);
    
    return;
  }

  void APEX::ConstructSiPMSAndBoard(G4VPhysicalVolume* mother_physical) const
  {
    // Construct first the SiPMs
    SiPMMPPC * sipm = nullptr;
    if(SiPM_code_==1){
      sipm = dynamic_cast<HamamatsuS133606050VE*>(sipm);
      sipm = new HamamatsuS133606050VE();
    }
    else if(SiPM_code_==2){
      sipm = dynamic_cast<HamamatsuS133605075HQR*>(sipm);
      sipm = new HamamatsuS133605075HQR();
    }
    else if(SiPM_code_==3){
      sipm = dynamic_cast<FbkNuvHdCryoTT*>(sipm);
      sipm = new FbkNuvHdCryoTT();
    }
    else if(SiPM_code_==4){
      sipm = dynamic_cast<BroadcomAFBRS4N44P044M*>(sipm);
      sipm = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm = dynamic_cast<PerfectSiPMMPPC*>(sipm);
      sipm = new PerfectSiPMMPPC();
    }
    sipm->SetReflectiveSupports(ref_phsensors_supports_);  
    sipm->Construct();
    G4double sipm_thickn = sipm->GetThickness();
    G4LogicalVolume* sipm_logic_vol = sipm->GetLogicalVolume();

    if (!sipm_logic_vol) {
      G4Exception("[APEX]", "ConstructSiPMSAndBoard()",
                  FatalException, "Null pointer to logical volume.");
    }

    G4RotationMatrix* sipm_rot = new G4RotationMatrix();
    G4ThreeVector base_pos;

    if(board_position_code_==1) // Board in the middle of a large face
    {
      sipm_rot->rotateX(0.0*deg);
      base_pos.set( (-1.*board_length_/2.) + (0.5*board_length_/num_phsensors_),
                    -1.*(plate_thickn_/2.)-1.*(sipm_thickn/2.)-gap_,    // Note that what's placed in the global origin of 
                                                                        // coordinates is the plate, not the reflective foil. 
                    0.);
    }
    else
    {
      sipm_rot->rotateX(-90.*deg);
      base_pos.set( (-1.*board_length_/2.) + (0.5*board_length_/num_phsensors_),
                    0.,
                    -1.*(plate_width_/2.)-1.*(sipm_thickn/2.)-gap_);
    }

    G4int phsensor_id = 0;
    for (G4int i=0; i<num_phsensors_; ++i) {
      new G4PVPlacement(sipm_rot, base_pos+G4ThreeVector(i*board_length_/num_phsensors_, 0., 0.),
                        sipm->GetModel(), sipm_logic_vol,
                        mother_physical, true, phsensor_id, true);

      phsensor_id += 1;
    }


    // Then construct the board (Vikuiti-coated FR4 piece)
    G4String board_name = "BOARD";

    G4double sipm_height = sipm->GetTransverseDim();
    G4double board_thickn = 1.*mm;

    G4Box* board_solid =
        new G4Box(  board_name, board_length_/2., 
                                sipm_height/2.,     // Board height matches that of the SiPMs
                                board_thickn/2.);

    G4LogicalVolume* board_logic = 
        new G4LogicalVolume(board_solid, materials::FR4(), board_name);

    G4VisAttributes board_col = nexus::White();
    //board_col.SetForceSolid(true);
    board_logic->SetVisAttributes(board_col);

    //VIKUITI coating for the board
    const G4String bc_name = "BOARD_COATING";
    G4OpticalSurface* board_coating = 
      new G4OpticalSurface(bc_name, unified, ground, dielectric_metal, 1);
    
    board_coating->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());
    new G4LogicalSkinSurface(bc_name, board_logic, board_coating); 

    G4RotationMatrix* board_rot = new G4RotationMatrix();
    G4ThreeVector board_pos;

    if(board_position_code_==1) // Board in the middle of a large face
    {
      board_rot->rotateX(90.0*deg);
      board_pos.set(0.,
                    -1.*(plate_thickn_/2.)-gap_
                    -sipm_thickn-1.*(board_thickn/2.),  // Note that what's placed in the global origin of 
                                                        // coordinates is the plate, not the reflective foil. 
                    0.);
    }
    else
    {
      board_rot->rotateX(0.0*deg);
      board_pos.set(0.,
                    0.,
                    -1.*(plate_width_/2.)-gap_
                    -sipm_thickn-1.*(board_thickn/2.));
    }

    //Place it
    new G4PVPlacement(board_rot, board_pos,
                      "COATED_BOARD", board_logic,
                      mother_physical,
                      false, 0, true);

    return;
  }

  void APEX::ConstructReflectiveFoil(G4VPhysicalVolume* mother_physical) const
  {
    const G4String ref_foil_name = "REF_FOIL";

    // The reflective foil covers every face of the plate but one
    // Get its volume as a subtraction solid from two boxes

    G4Box* aux_outer_box = new G4Box(  "AUX_OUTER_BOX", 
                                        (plate_length_  + (2.*reflective_foil_thickn_))/2., 
                                        (plate_thickn_+reflective_foil_thickn_)/2., 
                                        (plate_width_  + (2.*reflective_foil_thickn_))/2.);

    // Extra thickness to prevent boolean subtraction of solids with matching surfaces
    // See geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html#solids-made-by-boolean-operations

    G4double tolerance = 1.*mm; // To prevent matching surfaces in the boolean subtraction
    G4Box* aux_inner_box =  new G4Box(  "AUX_INNER_BOX", 
                                        plate_length_/2., 
                                        (plate_thickn_/2.)+tolerance, 
                                        plate_width_/2.);

    G4SubtractionSolid* ref_foil_solid =    new G4SubtractionSolid( ref_foil_name, 
                                                                    aux_outer_box, aux_inner_box, 
                                                                    nullptr, G4ThreeVector(0., (reflective_foil_thickn_/2.)+tolerance, 0.));
    SiPMMPPC* sipm_ptr = nullptr;
    if(SiPM_code_==1){
      sipm_ptr = new HamamatsuS133606050VE();
    }
    else if(SiPM_code_==2){
      sipm_ptr = new HamamatsuS133605075HQR();
    }
    else if(SiPM_code_==3){
      sipm_ptr = new FbkNuvHdCryoTT();
    }
    else if(SiPM_code_==4){
      sipm_ptr = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm_ptr = new PerfectSiPMMPPC();
    }

    G4double sipm_transverse_dim = sipm_ptr->GetTransverseDim();
    G4double sipm_thickness = sipm_ptr->GetThickness();

    G4double thickness_of_dummy_sipm = reflective_foil_thickn_+(0.1*mm);
    G4Box* dummy_sipm =  new G4Box( "DUMMY_SIPM", 
                                    sipm_transverse_dim/2., 
                                    thickness_of_dummy_sipm/2., // Setting here the reflective-foil thickness plus some tolerance so that:
                                                                //  1)  if board_position_code_==1, the carved hole is a pass-through hole
                                                                //  2)  if board_position_code_==2, we prevent matching surfaces in the boolean subtraction
                                                                //      In this second case, the value of the tolerance actually matters. It must be big 
                                                                //      enough so as to prevent matching surfaces, but small enough so as to not carve to 
                                                                //      match the horizontal portion of the reflective foil.
                                    sipm_transverse_dim/2.);

    G4MultiUnion* reflective_foil_holes = new G4MultiUnion("REF_FOIL_HOLES");

    G4double pos;
    G4Transform3D* transform_ptr = nullptr;
    G4RotationMatrix* rot = new G4RotationMatrix();
    if(board_position_code_!=1){
      rot->rotateX(90.0*deg);
    }

    for(G4int i=0; i<num_phsensors_; i++){

      pos = (-1.*board_length_/2.) + ((0.5 + i)*board_length_/num_phsensors_);    
      transform_ptr = new G4Transform3D(*rot, G4ThreeVector(pos, 0., 0.));
      reflective_foil_holes->AddNode(*dummy_sipm, *transform_ptr);    

    }

    reflective_foil_holes->Voxelize();

    G4ThreeVector vec = G4ThreeVector(0., 
                                      -1.*plate_thickn_/2., // Minus half the thickness of AUX_OUTER_BOX 
                                                            // plus half the reflective-foil thickness
                                      0.);
    if(board_position_code_!=1){

      SiPMBoard board;
      board.SetSiPMCode(SiPM_code_);

      vec = G4ThreeVector(0.,
                          0.,
                          -1.*(plate_length_/2.)-1.*(reflective_foil_thickn_/2.));  // Minus half the length of the plate
                                                                                    // minus half the reflective-foil thickness
    }

    ref_foil_solid = new G4SubtractionSolid(ref_foil_name, 
                                            ref_foil_solid, reflective_foil_holes, 
                                            nullptr, vec);
    G4LogicalVolume* ref_case_logic = 
      new G4LogicalVolume(ref_foil_solid, materials::FR4(), ref_foil_name);

    // Set its color for visualization purposes
    G4VisAttributes ref_case_col = nexus::WhiteAlpha();
    //ref_case_col.SetForceSolid(true);
    ref_case_logic->SetVisAttributes(ref_case_col);

    //Now create the reflectivie optical surface
    const G4String ref_surf_name = "REF_SURFACE";
    G4OpticalSurface* refsurf_opsurf = 
      new G4OpticalSurface(ref_surf_name, unified, ground, dielectric_metal, 1);
    
    // From geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#optical-photon-processes
    // The dielectric_metal->ground configuration of the unified model works as:
    // "Only reflection or absorption; No refraction: Reflection probability set by
    // reflectivity. If reflected, one of the four specular spike, backscatter,
    // lambertian or specular lobe reflection with respect to a FacetNormal takes
    // place according to the assigned probabilities."
    // So, make sure you have set the reflectivity and the probabilities for each
    // type of reflection. On the other hand, for this configuration, it does not matter
    // if you have set the transmission, since that option is already banned from the
    // configuration model.

    refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());
    new G4LogicalSkinSurface(ref_surf_name, ref_case_logic, refsurf_opsurf);   
    
    new G4PVPlacement(  nullptr, G4ThreeVector(0., -1.*reflective_foil_thickn_/2., 0.),
                        ref_foil_name, ref_case_logic, 
                        mother_physical,
                        false, 0, true);
    return;
  }

  void APEX::ConstructDichroicFilter(G4VPhysicalVolume* mother_physical) const
  {
      
    // -------------------------------- DICHROIC FILTER MODEL --------------------------------
    //
    // The DF model is implemented in the following manner:
    //
    //
    // _______________________________________________________________________________
    //
    //                                  PTP coating
    // _______________________________________________________________________________
    //
    //                              First half of the MLS
    //                                   (n = n_x)
    //                        
    // --------------------------- G4LogicalBorderSurface ----------------------------
    //
    //                              Second half of the MLS
    //                                   (n = n_x)
    // _______________________________________________________________________________
    //
    //                                   WLS Plate
    // _______________________________________________________________________________
    //
    //
    // When modelling the X-Arapuca (in XArapuca.cc), with actual DFs accounting for a 
    // substrate, we substracted the contribution from
    //  
    //  1) The Fresnel transmission (FT) from the lab r.index to the DF substrate and
    //  2) The FT from the MLS effective refractive index to the lab r.index,
    //
    //  from the measured transmission curve. For more information on why we introduced 
    //  this correction, check the documentation in XArapuca::ConstructDichroicAssemblies(), 
    //  in XArapuca.cc.
    //
    //  There are two options now:
    //
    //  1)  The first one is setting n_x = n_eff, where n_eff is the effective refractive 
    //      index of the MLS. In this case, G4 simulates the FT of the n_pTP->n_eff and 
    //      n_eff->n_WLSp interfaces:
    //
    //        On the n_pTP->n_eff interface: There are reasons to think both, that this is 
    //        OK and not OK. pTP-coated filters are not measured in the laboratory, so
    //        letting G4 simulate the FT in this interface might be the best approach. On 
    //        the other hand, one could think that it is not realistic to let G4 simulate 
    //        the FT in this interface, since the first layer (and every layer) of the MLS 
    //        is very thin and inteference phenomena with the subsequent layer may occurr.
    //    
    //        On the n_eff->n_WLSp interface: Letting G4 simulate this one is analogous to
    //        what we did in the XArapuca case, since there, we substracted the FT
    //        from the MLS to the lab r. index, to then simulate the FT from the MLS volume
    //        (with refractive index equal to that of the substrate, in order not to lose
    //        the snell information) to the LAr r. index. This was already an inaccuracy. 
    //        Now, we will be simulating the FT from the MLS to the WLS plate.
    //
    //  2)  The second one entails setting n_x = n_wlsp, where n_wlsp is the refractive 
    //      index of the wavelength shifting plate. This alternative makes sense if we
    //      stick to the following understanding:
    //      
    //        The intrinsic transmitance curve (ITC), which we compute by dividing the 
    //        measured transmission curve (TC) by the FT of the n_lab->n_DFsubs and by
    //        the FT of the n_MLS->n_lab, contains the TC information of the MLS alone.
    //
    //      Since APEX contains no DF substrate, but just its MLS, it is reasonable to 
    //      think that the DF implementation should allow the photons to just 'feel' the
    //      the ITC in the way inwards and outwards the WLS plate. This implementation, 
    //      i.e. n_x=n_wlsp, does so except for the fact that it simulates the FT from
    //      n_ptp->n_wlsp.
    //
    //        For photons that are emitted by pTP, they feel the FT of n_ptp->n_wlsp and
    //        the ITC.
    //
    //        For photons that are emitted by the WLSp, they only feel the ITC and, in
    //        case they are transmitted by the ITC, then they feel the n_ptp->n_wlsp FT.
    //        
    //      Again, there're reasons to think that letting G4 simulate the n_ptp<->n_wlsp
    //      is both correct and incorrect. They key point is that G4 won't ever correctly
    //      simulate the interface between some media and the first layer of a MLS, where 
    //      the layers are wavelength-order-of-magnitude thin. That's because G4 cannot 
    //      simulate wave interference. In this context, I think it is good enough (and
    //      maybe the best we can achieve with G4) to account for some reflectance in the 
    //      PTP->MLS interface (which is realized by the n_pTP->n_wlsp FT in this 
    //      implementation), which will (physically) happen in APEX since there's a change 
    //      of refractive index from PTP to the MLS, while still sticking to our 
    //      understanding that, when reaching the DF side of the WLS plate, from within the
    //      WLS plate, the photon should not 'feel' internal reflection, but just the ITC.
    //      That's why I am going to stick to this alternative, i.e. n_x=n_wlsp
    //
    //  Summary:  I am going to set n_x=n_wlsp, where n_wlsp is the refractive index of the
    //            wavelength shifting plate.
    //      
    //
    // ---------------------------------------------------------------------------------------

    // This function is called by APEX::Construct() only if !remove_MLS_

    G4Box* MLS_half_solid = new G4Box("AUX", plate_length_/2., MLS_thickn_/4., plate_width_/2.);
    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    mat->SetMaterialPropertiesTable(opticalprops::TunableRIMat(WLSp_rindex_));  // Change this to MLS_rindex_ if you
                                                                                // want to go for the first alternative
                                                                                // of DF implementation explained above
    G4LogicalVolume* MLS_half_logic = new G4LogicalVolume(MLS_half_solid, mat, "MLS_HALF");
            
    G4VisAttributes MLS_col = nexus::BloodRedAlpha();
    MLS_col.SetForceSolid(true);
    MLS_half_logic->SetVisAttributes(MLS_col);

    // Place the MLS
    G4VPhysicalVolume* MLS_first_half = dynamic_cast<G4VPhysicalVolume*>(   // This is the outermost one
        new G4PVPlacement(nullptr, G4ThreeVector(0., plate_thickn_/2.
                                                    +MLS_thickn_/2.         // Note that the thickness of MLS_half_solid is MLS_thickn_/2
                                                    +MLS_thickn_/4., 0.), 
                          "FIRST_MLS_HALF", MLS_half_logic, mother_physical, true, 0, true));

    G4VPhysicalVolume* MLS_second_half = dynamic_cast<G4VPhysicalVolume*>(  // This is the internal one
        new G4PVPlacement(nullptr, G4ThreeVector(0., plate_thickn_/2.
                                                    +MLS_thickn_/4., 0.), 
                          "SECOND_MLS_HALF", MLS_half_logic, mother_physical, true, 1, true));

    // Check that there's dichroic information for ingoing (wrt APEX) photons
    if(path_to_inwards_dichroic_data_==""){
        G4Exception("[APEX]", "ConstructDichroicFilter()",
                    FatalException, "The path to the inwards dichroic data file was not set.");
    }

    // Check that there's dichroic information for outgoing photons
    if(path_to_outwards_dichroic_data_==""){
        G4Exception("[APEX]", "ConstructDichroicFilter()",
                    FatalException, "The path to the outwards dichroic data file was not set.");
    }

    // Construct the ingoing optical surface
    setenv("G4DICHROICDATA", path_to_inwards_dichroic_data_, 1);
    G4OpticalSurface* df_inwards_opsurf =                 // G4OpticalSurface constructor loads the
        new G4OpticalSurface( "DICHROIC_INWARDS_OPSURF",  // dichroic information from the file which
                              dichroic,                   // is currently pointed to by the environment
                              polished,                   // variable G4DICHROICDATA
                              dielectric_dichroic);

    // Construct the outgoung optical surface
    setenv("G4DICHROICDATA", path_to_outwards_dichroic_data_, 1);   // Note that, if you did not compile the modified version of G4 code, 
                                                                    // then different G4 dichroic data cannot be loaded. Instead, the first 
                                                                    // one (i.e. the one I am setting from path_to_inwards_dichroic_data_), 
                                                                    // is the one that will apply for every dichroic boundary in the simulation.
    G4OpticalSurface* df_outwards_opsurf =   
        new G4OpticalSurface( "DICHROIC_OUTWARDS_OPSURF", 
                              dichroic, 
                              polished, 
                              dielectric_dichroic);

    // Endow the MLS_first_half->MLS_second_half surface with the ingoing optical surface
    new G4LogicalBorderSurface( "MLS1->MLS2", 
                                MLS_first_half, 
                                MLS_second_half, 
                                df_inwards_opsurf);

    // Endow the MLS_second_half->MLS_first_half surface with the outgoing optical surface
    new G4LogicalBorderSurface( "MLS2->MLS1", 
                                MLS_second_half, 
                                MLS_first_half, 
                                df_outwards_opsurf);      
    // pTP coating
    if(!remove_coating_)
    {
        G4Box* coating_solid = new G4Box( "COATING", 
                                          plate_length_/2., coating_thickn_/2., plate_width_/2.);

        G4Material* coating_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
        coating_mat->SetMaterialPropertiesTable(opticalprops::PTP(coating_rindex_));
        G4LogicalVolume* coating_logic = 
                            new G4LogicalVolume(coating_solid, coating_mat, "COATING");   

        G4VisAttributes coating_col = nexus::TitaniumGreyAlpha();
        coating_col.SetForceSolid(true);
        coating_logic->SetVisAttributes(coating_col);

        // Place the coating
        G4VPhysicalVolume* coating_physical = dynamic_cast<G4VPhysicalVolume*>(
            new G4PVPlacement(  nullptr, G4ThreeVector(0.,  plate_thickn_/2.
                                                            +MLS_thickn_
                                                            +coating_thickn_/2., 0.), 
                                "COATING", coating_logic, mother_physical, false, 0, true));

        // Make the LAR-coating interface rough, so that photons cannot be trapped within the coating
        G4OpticalSurface* coating_rough_surf =
                new G4OpticalSurface("COATING_ROUGH_SURFACE", glisur, ground, dielectric_dielectric, .01);
                // 0.01 is the polish value for glisur model that was measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
                // This is the best reference we have, since both PTP and TPB are the result of an evaporation+deposition process
        new G4LogicalBorderSurface( "SURROUNDINGS->COATING", 
                                    mother_physical, 
                                    coating_physical, 
                                    coating_rough_surf);
        new G4LogicalBorderSurface( "COATING->SURROUNDINGS", 
                                    coating_physical, 
                                    mother_physical, 
                                    coating_rough_surf);
        // We will also add roughness for the coating->MLS interface, but only with such ordering. The alternative case takes place
        // when the photon travels from the MLS to the coating. The MLS is supposed to be polished, so the photon may not see a rough 
        // surface.
        new G4LogicalBorderSurface( "COATING->MLS", 
                                    coating_physical, 
                                    MLS_first_half, 
                                    coating_rough_surf);
    }
    return;
  }

  void APEX::ConstructBoard(G4VPhysicalVolume* mother_physical) const
  {
    SiPMBoard board;
    board.SetBaseID(0);
    board.SetBoardLength(board_length_);
    board.SetSiPMCode(SiPM_code_);
    board.SetNumPhsensors(num_phsensors_);
    board.SetReflectiveSupports(ref_phsensors_supports_);
    board.Construct();
    G4LogicalVolume* board_logic_vol = board.GetLogicalVolume();

    G4RotationMatrix* rot = new G4RotationMatrix();
    G4ThreeVector pos;

    if(board_position_code_==1) // Board in the middle of a large face
    {
      rot->rotateX(-90.0*deg);
      pos.set(0.,
              -1.*(plate_thickn_/2.)-1.*(board.GetOverallThickness()/2.)-gap_,  // Note that what's placed in the global origin of 
                                                                                // coordinates is the plate, not the reflective foil. 
              0.);
    }
    else
    {
      rot->rotateY(+180.0*deg);
      pos.set(0., 
              0., 
              -1.*(plate_length_/2.)-1.*(board.GetOverallThickness()/2.));
    }
  
    new G4PVPlacement(rot, pos, "SIPMS_BOARD", board_logic_vol, 
                      mother_physical, false, 0, false);
    // SiPMBoard logical volume is an encasing volume which may collide into other volumes
    // No need to set pSurfCheck for that volume (dimples). As we are setting pSurfCheck=false
    // (so that no harmless-overlap warning pops up in a with-dimples configuration), you have to
    // be extra careful to examine when there's actually a problematic overlap of this volume with
    // another one (since Geant4 won't warn you).

    return;
  }

  G4ThreeVector APEX::GenerateVertex(const G4String&) const{

    G4double tolerance = 0.1*mm;    // Small distance over the dichroic filter from
                                    // which photons are launched. Also, the width of 
                                    // the outer border projected over the DF from 
                                    // which photons won't be launched (Just see the 
                                    // implementation in x_pos and z_pos below to 
                                    // understand its meaning)
    G4double x_pos, z_pos;
    G4double y_pos = plate_thickn_/2. +MLS_thickn_ +coating_thickn_ +tolerance;

    if(generation_region_=="custom"){
      G4double random_radius =  UniformRandomInRange(gen_diameter_/2., 0.);
      G4double random_angle =   UniformRandomInRange(twopi, 0.); 
      x_pos = gen_x_ +(random_radius*sin(random_angle));
      z_pos = gen_z_ +(random_radius*cos(random_angle));
    }
    else{ // Default behaviour is that of generation_region_=="random"
      x_pos = UniformRandomInRange( plate_length_/2.,
                                    -1.*plate_length_/2.);
      z_pos = UniformRandomInRange( plate_width_/2.,
                                    -1.*plate_width_/2.);
    }
    return G4ThreeVector(x_pos, y_pos, z_pos);
  }

  G4bool APEX::GeometryIsIllFormed()                ///< The only check to make is that the sipm thickness should be bigger or equal to the 
                                                    ///< reflective foil thickness. Otherwise, the SiPM surface may not make it to the WLS 
                                                    ///< plate surface depending on whether the board height is bigger or smaller than the sipm 
                                                    ///< height (collision of the SiPM board into the reflective foil may happen). For the rest 
                                                    ///< of it, if the given parameters comply with the range set to their 
                                                    ///< G4GenericMessenger::Command, then the geometry is always feasible. One could think of
                                                    ///< an exception regarding the number of photosensors that a SiPMBoard can allocate, but
                                                    ///< that check is already performed by SiPMBoard::GeometryIsIllFormed(). 
  {

    SiPMMPPC* sipm_ptr = nullptr;
    if(SiPM_code_==1){
      sipm_ptr = new HamamatsuS133606050VE();
    }
    else if(SiPM_code_==2){
      sipm_ptr = new HamamatsuS133605075HQR();
    }
    else if(SiPM_code_==3){
      sipm_ptr = new FbkNuvHdCryoTT();
    }
    else if(SiPM_code_==4){
      sipm_ptr = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm_ptr = new PerfectSiPMMPPC();
    }

    if(sipm_ptr->GetThickness()<reflective_foil_thickn_){
      return true;
    }

    return false;
  }

} //End namespace nexus
