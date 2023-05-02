#include "LAttMeas.h"
#include "PmtR7378A.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "Visibilities.h"

#include <algorithm>
#include <random>
#include <cmath>
#include <G4GenericMessenger.hh>
#include <G4UserLimits.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4ThreeVector.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4Tubs.hh>
#include <G4VisAttributes.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Orb.hh>

#include <CLHEP/Units/SystemOfUnits.h>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(LAttMeas, GeometryBase)

//Next thing to do: write geometry_is_ill_formed()

namespace nexus{

  LAttMeas::LAttMeas():
  GeometryBase(), 
  ///Get internal reflector cavity dimensions from arxiv.org/abs/1912.09191
  config_code_              (1),
  with_box_                 (true),
  black_box_internal_dx_    (70.*cm),     
  black_box_internal_dy_    (65.*cm),
  black_box_internal_dz_    (40.*cm),
  plate_x_pos_              (35.*cm),                      ///< Cartesian coordinates of the geometric
  plate_y_pos_              (10.*cm),                      ///< center of the WLS plate with respect to
  plate_z_pos_              (0.*cm),
  plate_dx_                 (50.*cm),
  plate_dy_                 (3.9*mm),
  plate_dz_                 (12.5*cm),
  wls_attlength_            (1.5*m),
  gap_                      (5.*mm),
  with_holder_              (false),
  holder_x_pos_             (45.*cm),
  holder_y_pos_             (10.*cm),
  holder_z_pos_             (0.*cm),
  holder_dx_                (7.6*cm),
  holder_dy_                (4.*cm),
  holder_dz_                (3.8*cm),
  pmt_x_pos_                (62.15*cm),
  pmt_y_pos_                (10.*cm),
  pmt_z_pos_                (0.*cm),
  pmt_case_depth_           (7.*cm),
  pmt_case_diameter_        (4.5*cm),
  pmt_sensarea_diameter_    (2.7*cm),
  with_adapter_             (false),
  adapter_x_pos_            (32.6*cm),
  adapter_y_pos_            (10.*cm),
  adapter_z_pos_            (0.*cm),
  adp_diameter_             (6.*cm),
  adp_thickn_               (5.*mm),
  adp_hole_height_          (26.4*mm),
  adp_slot_width_           (4.*mm),
  adp_slot_thickn_          (3.*mm),
  gen_vertex_x_             (0.*cm),
  gen_vertex_y_             (15.*cm),
  gen_vertex_z_             (15.*cm),
  world_length_             (1000.*m)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/LAttMeas/",
				"Control commands of geometry LAttMeas.");

    G4GenericMessenger::Command& cc_cmd =
      msg_->DeclareProperty("config_code", config_code_,
			    "Configuration code.");
    cc_cmd.SetParameterName("config_code", false);
    cc_cmd.SetRange("config_code>=1"); 

    G4GenericMessenger::Command& wb_cmd =
      msg_->DeclareProperty("with_box", with_box_,
			    "Whether to construct a black box surrounding the plate or not.");

    G4GenericMessenger::Command& bbidx_cmd =
      msg_->DeclareProperty("black_box_internal_dx", black_box_internal_dx_,
			    "Internal length of the black box.");
    bbidx_cmd.SetUnitCategory("Length");
    bbidx_cmd.SetParameterName("black_box_internal_dx", false);
    bbidx_cmd.SetRange("black_box_internal_dx>0.");

    G4GenericMessenger::Command& bbidy_cmd =
      msg_->DeclareProperty("black_box_internal_dy", black_box_internal_dy_,
			    "Internal height of the black box.");
    bbidy_cmd.SetUnitCategory("Length");
    bbidy_cmd.SetParameterName("black_box_internal_dy", false);
    bbidy_cmd.SetRange("black_box_internal_dy>0.");

    G4GenericMessenger::Command& bbidz_cmd =
      msg_->DeclareProperty("black_box_internal_dz", black_box_internal_dz_,
			    "Internal width of the black box.");
    bbidz_cmd.SetUnitCategory("Length");
    bbidz_cmd.SetParameterName("black_box_internal_dz", false);
    bbidz_cmd.SetRange("black_box_internal_dz>0.");
    
   G4GenericMessenger::Command& pxp_cmd =
      msg_->DeclareProperty("plate_x_pos", plate_x_pos_,
			    "X-Position of the WLS plate.");
    pxp_cmd.SetUnitCategory("Length");
    pxp_cmd.SetParameterName("plate_x_pos", false);

   G4GenericMessenger::Command& pyp_cmd =
      msg_->DeclareProperty("plate_y_pos", plate_y_pos_,
			    "Y-Position of the WLS plate.");
    pyp_cmd.SetUnitCategory("Length");
    pyp_cmd.SetParameterName("plate_y_pos", false);

    G4GenericMessenger::Command& pzp_cmd =
      msg_->DeclareProperty("plate_z_pos", plate_z_pos_,
			    "Z-Position of the WLS plate.");
    pzp_cmd.SetUnitCategory("Length");
    pzp_cmd.SetParameterName("plate_z_pos", false);

    G4GenericMessenger::Command& pdx_cmd =
      msg_->DeclareProperty("plate_dx", plate_dx_,
			    "Length of the WLS plate.");
    pdx_cmd.SetUnitCategory("Length");
    pdx_cmd.SetParameterName("plate_dx", false);
    pdx_cmd.SetRange("plate_dx>0.");

    G4GenericMessenger::Command& pdy_cmd =
      msg_->DeclareProperty("plate_dy", plate_dy_,
			    "Width of the WLS plate.");
    pdy_cmd.SetUnitCategory("Length");
    pdy_cmd.SetParameterName("plate_dy", false);
    pdy_cmd.SetRange("plate_dy>0.");

    G4GenericMessenger::Command& pdz_cmd =
      msg_->DeclareProperty("plate_dz", plate_dz_,
			    "Thickness of the WLS plate.");
    pdz_cmd.SetUnitCategory("Length");
    pdz_cmd.SetParameterName("plate_dz", false);
    pdz_cmd.SetRange("plate_dz>0.");

    G4GenericMessenger::Command& wlsal_cmd =
      msg_->DeclareProperty("wls_attlength", wls_attlength_,
			    "Attenuation length of the WLS plate.");
    wlsal_cmd.SetUnitCategory("Length");
    wlsal_cmd.SetParameterName("wls_attlength", false);
    wlsal_cmd.SetRange("wls_attlength>0.");

    G4GenericMessenger::Command& g_cmd =
      msg_->DeclareProperty("gap", gap_,
			    "Gap between the WLS plate and the adapter/PMT.");
    g_cmd.SetUnitCategory("Length");
    g_cmd.SetParameterName("gap", false);
    g_cmd.SetRange("gap>=0.");

    G4GenericMessenger::Command& wh_cmd =
      msg_->DeclareProperty("with_holder", with_holder_,
			    "Whether to construct the holder or not.");

    G4GenericMessenger::Command& hxp_cmd =
      msg_->DeclareProperty("holder_x_pos", holder_x_pos_,
          "X-Position of the holder plate.");
    hxp_cmd.SetUnitCategory("Length");
    hxp_cmd.SetParameterName("holder_x_pos", false);

   G4GenericMessenger::Command& hyp_cmd =
      msg_->DeclareProperty("holder_y_pos", holder_y_pos_,
          "Y-Position of the holder plate.");
    hyp_cmd.SetUnitCategory("Length");
    hyp_cmd.SetParameterName("holder_y_pos", false);

    G4GenericMessenger::Command& hzp_cmd =
      msg_->DeclareProperty("holder_z_pos", holder_z_pos_,
          "Z-Position of the holder plate.");
    hzp_cmd.SetUnitCategory("Length");
    hzp_cmd.SetParameterName("holder_z_pos", false);

    G4GenericMessenger::Command& hdx_cmd =
      msg_->DeclareProperty("holder_dx", holder_dx_,
          "Length of the holder plate.");
    hdx_cmd.SetUnitCategory("Length");
    hdx_cmd.SetParameterName("holder_dx", false);
    hdx_cmd.SetRange("holder_dx>0.");

    G4GenericMessenger::Command& hdy_cmd =
      msg_->DeclareProperty("holder_dy", holder_dy_,
          "Width of the holder plate.");
    hdy_cmd.SetUnitCategory("Length");
    hdy_cmd.SetParameterName("holder_dy", false);
    hdy_cmd.SetRange("holder_dy>0.");

    G4GenericMessenger::Command& hdz_cmd =
      msg_->DeclareProperty("holder_dz", holder_dz_,
          "Thickness of the holder plate.");
    hdz_cmd.SetUnitCategory("Length");
    hdz_cmd.SetParameterName("holder_dz", false);
    hdz_cmd.SetRange("holder_dz>0.");

    G4GenericMessenger::Command& pmtxp_cmd =
      msg_->DeclareProperty("pmt_x_pos", pmt_x_pos_,
          "X-Position of the PMT.");
    pmtxp_cmd.SetUnitCategory("Length");
    pmtxp_cmd.SetParameterName("pmt_x_pos", false);

   G4GenericMessenger::Command& pmtyp_cmd =
      msg_->DeclareProperty("pmt_y_pos", pmt_y_pos_,
          "Y-Position of the PMT.");
    pmtyp_cmd.SetUnitCategory("Length");
    pmtyp_cmd.SetParameterName("pmt_y_pos", false);

    G4GenericMessenger::Command& pmtzp_cmd =
      msg_->DeclareProperty("pmt_z_pos", pmt_z_pos_,
          "Z-Position of the PMT.");
    pmtzp_cmd.SetUnitCategory("Length");
    pmtzp_cmd.SetParameterName("pmt_z_pos", false);

    G4GenericMessenger::Command& pmtcde_cmd =
      msg_->DeclareProperty("pmt_case_depth", pmt_case_depth_,
			    "Depth of the PMT case.");
    pmtcde_cmd.SetUnitCategory("Length");
    pmtcde_cmd.SetParameterName("pmt_case_depth", false);
    pmtcde_cmd.SetRange("pmt_case_depth>0.");

    G4GenericMessenger::Command& pmtcdi_cmd =
      msg_->DeclareProperty("pmt_case_diameter", pmt_case_diameter_,
			    "Diameter of the PMT case.");
    pmtcdi_cmd.SetUnitCategory("Length");
    pmtcdi_cmd.SetParameterName("pmt_case_diameter", false);
    pmtcdi_cmd.SetRange("pmt_case_diameter>0.");

    G4GenericMessenger::Command& pmtsd_cmd =
      msg_->DeclareProperty("pmt_sensarea_diameter", pmt_sensarea_diameter_,
			    "Diameter of the PMT sensitive area (window).");
    pmtsd_cmd.SetUnitCategory("Length");
    pmtsd_cmd.SetParameterName("pmt_sensarea_diameter", false);
    pmtsd_cmd.SetRange("pmt_sensarea_diameter>0.");

    G4GenericMessenger::Command& wa_cmd =
      msg_->DeclareProperty("with_adapter", with_adapter_,
			    "Whether to place an adapter in between the WLS plate and the PMT.");

    G4GenericMessenger::Command& axp_cmd =
      msg_->DeclareProperty("adapter_x_pos", adapter_x_pos_,
          "X-Position of the PMT adapter.");
    axp_cmd.SetUnitCategory("Length");
    axp_cmd.SetParameterName("adapter_x_pos", false);

   G4GenericMessenger::Command& ayp_cmd =
      msg_->DeclareProperty("adapter_y_pos", adapter_y_pos_,
          "Y-Position of the PMT adapter.");
    ayp_cmd.SetUnitCategory("Length");
    ayp_cmd.SetParameterName("adapter_y_pos", false);

    G4GenericMessenger::Command& azp_cmd =
      msg_->DeclareProperty("adapter_z_pos", adapter_z_pos_,
          "Z-Position of the PMT adapter.");
    azp_cmd.SetUnitCategory("Length");
    azp_cmd.SetParameterName("adapter_z_pos", false);

    G4GenericMessenger::Command& ad_cmd =
      msg_->DeclareProperty("adp_diameter", adp_diameter_,
			    "Adapter diameter.");
    ad_cmd.SetUnitCategory("Length");
    ad_cmd.SetParameterName("adp_diameter", false);
    ad_cmd.SetRange("adp_diameter>0.");

    G4GenericMessenger::Command& at_cmd =
      msg_->DeclareProperty("adp_thickn", adp_thickn_,
			    "Adapter thickness.");
    at_cmd.SetUnitCategory("Length");
    at_cmd.SetParameterName("adp_thickn", false);
    at_cmd.SetRange("adp_thickn>0.");

    G4GenericMessenger::Command& ahh_cmd =
      msg_->DeclareProperty("adp_hole_height", adp_hole_height_,
			    "Height of the adapter hole.");
    ahh_cmd.SetUnitCategory("Length");
    ahh_cmd.SetParameterName("adp_hole_height", false);
    ahh_cmd.SetRange("adp_hole_height>0.");

    G4GenericMessenger::Command& asw_cmd =
      msg_->DeclareProperty("adp_slot_width", adp_slot_width_,
			    "Width of the adapter slot.");
    asw_cmd.SetUnitCategory("Length");
    asw_cmd.SetParameterName("adp_slot_width", false);
    asw_cmd.SetRange("adp_slot_width>0.");

    G4GenericMessenger::Command& ast_cmd =
      msg_->DeclareProperty("adp_slot_thickn", adp_slot_thickn_,
			    "Thickness of the adapter right at the slot.");
    ast_cmd.SetUnitCategory("Length");
    ast_cmd.SetParameterName("adp_slot_thickn", false);
    ast_cmd.SetRange("adp_slot_thickn>0.");

    G4GenericMessenger::Command& gvx_cmd =
      msg_->DeclareProperty("gen_vertex_x", gen_vertex_x_,
			    "X-position of the light generation vertex.");
    gvx_cmd.SetUnitCategory("Length");
    gvx_cmd.SetParameterName("gen_vertex_x", false);

    G4GenericMessenger::Command& gvy_cmd =
      msg_->DeclareProperty("gen_vertex_y", gen_vertex_y_,
			    "Y-position of the light generation vertex.");
    gvy_cmd.SetUnitCategory("Length");
    gvy_cmd.SetParameterName("gen_vertex_y", false);

    G4GenericMessenger::Command& gvz_cmd =
      msg_->DeclareProperty("gen_vertex_z", gen_vertex_z_,
			    "Z-position of the light generation vertex.");
    gvz_cmd.SetUnitCategory("Length");
    gvz_cmd.SetParameterName("gen_vertex_z", false);

    // When testing WLS plates with opticalprops::noAbsLength_, it is possible that a photon
    // gets trapped within the plate (below the critical angle) into an infinite-bouncing-loop.
    // For the case of an EJ286 plate with DUNE supercells dimensions, with this absorption 
    // length and immersed into LAr, some analysis revealed that particles that did not fall 
    // into this infinite loop had track lengths of less than one hundred meters.
    ul_ = new G4UserLimits();
    ul_->SetUserMaxTrackLength(100.*m);
  }

  LAttMeas::~LAttMeas()
  {
      if(ul_) delete ul_;
  }

  void LAttMeas::Construct()
  {

    // Compute internal attributes ----------
    gen_vertex_y_ = plate_y_pos_; // This one matches the height dimension
                                  // This one is the same for every configuration
    if(config_code_==1){
      // In this case, gen_vertex_x_ is taken from the configuration macro
      gen_vertex_z_ = plate_z_pos_ - plate_dz_/2. -2.*mm;
    }
    else if(config_code_==2){
      gen_vertex_x_ = plate_x_pos_ - plate_dx_/2. - 1.*cm;
      gen_vertex_z_ = plate_z_pos_;
    }
    else{ // config_code_==3 is default
      gen_vertex_z_ = plate_z_pos_;

      G4double aux = plate_dx_;
      plate_dx_=plate_dy_;
      plate_dy_=aux;

      gen_vertex_x_ = plate_x_pos_ - plate_dx_/2. - 1.*cm;
    }

    PmtR7378A* pmt = new PmtR7378A();
    pmt->Construct();
    G4double pmt_length = pmt->Length();
    delete pmt;

    pmt_x_pos_ = plate_x_pos_+(plate_dx_/2.)+(pmt_length/2.);
    if(with_adapter_) pmt_x_pos_ += adp_slot_thickn_;
    pmt_x_pos_ += gap_;

    pmt_y_pos_ = plate_y_pos_;
    pmt_z_pos_ = plate_z_pos_;

    adapter_x_pos_ = pmt_x_pos_ - (pmt->Length()/2.) - adp_thickn_/2.;
    adapter_y_pos_ = pmt_y_pos_;
    adapter_z_pos_ = pmt_z_pos_;

    // --------------------------------------

    /*
    if(geometry_is_ill_formed()){
      G4Exception("[LAttMeas]", "Construct()", FatalException,
      "The given dimensions do not describe a feasible setup.");
    }
    */

    // The biggest volume is a vacuum box with dimensions 
    // world_length_ x world_length_ x world_length_, which 
    // surrounds the setup. 
    
    // Vacuum box
    const G4String world_name = "WORLD";

    G4Box* world_solid =
        new G4Box(world_name,
                world_length_/2.,
                world_length_/2.,
                world_length_/2.);
                
    G4Material* vacuum =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

    G4LogicalVolume* world_logic = 
        new G4LogicalVolume(world_solid, vacuum, world_name, 0, 0, 0, true);
    world_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    this->SetLogicalVolume(world_logic);

    // Air atmosphere
    const G4String atmosphere_name = "AIR";
    
    G4Box* atmosphere_solid =
      new G4Box(atmosphere_name,
                world_length_/2.,
                world_length_/2.,
                world_length_/2.);

    G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Air());

    G4LogicalVolume* atmosphere_logic =
      new G4LogicalVolume(atmosphere_solid, air, atmosphere_name);
    atmosphere_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VPhysicalVolume* atmosphere_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement (new G4RotationMatrix(), G4ThreeVector(0., 0., 0.), atmosphere_logic, 
                        atmosphere_name, world_logic, false, 0, true));

    G4VPhysicalVolume* mother_ptr = nullptr;
    if(!with_box_)  mother_ptr = atmosphere_physical;
    else            mother_ptr = ConstructBlackBox(atmosphere_physical);

    ConstructWLSPlate(mother_ptr);
    if(with_adapter_) ConstructPlateAdapter(mother_ptr);
    ConstructPMT(mother_ptr);       
    if(with_holder_) ConstructPlateHolders(mother_ptr);
    return;

  }

  G4VPhysicalVolume* LAttMeas::ConstructBlackBox(G4VPhysicalVolume* mother_physical) const
  { 
     
    // Black box (filled with air)
    const G4String black_box_name = "BLACK_BOX";

    G4Box* black_box_solid =
      new G4Box(black_box_name,
                black_box_internal_dx_/2.,
                black_box_internal_dy_/2.,
                black_box_internal_dz_/2.);

    G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Air());

    G4LogicalVolume* black_box_logic =
      new G4LogicalVolume(black_box_solid, air, black_box_name);
    //black_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    return dynamic_cast<G4VPhysicalVolume*>(
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), black_box_name, 
               black_box_logic, mother_physical, false, 0, true));
    
  }

  void LAttMeas::ConstructWLSPlate(G4VPhysicalVolume* mother_physical) const
  { 
    
    const G4String plate_name = "PLATE";

    G4Box* plate_solid =
        new G4Box(plate_name,
                plate_dx_/2.,
                plate_dy_/2.,
                plate_dz_/2.);

    G4Material* plastic_scintillator = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    plastic_scintillator->SetMaterialPropertiesTable(opticalprops::G2P_FB118(wls_attlength_));

    G4LogicalVolume* plate_logic =
      new G4LogicalVolume(plate_solid, plastic_scintillator, plate_name);

    new G4PVPlacement(nullptr, G4ThreeVector(plate_x_pos_, plate_y_pos_, plate_z_pos_), plate_name, 
               plate_logic, mother_physical, false, 0, true);

    return;

    /*
    plate_logic->SetUserLimits(ul_);
    G4VisAttributes wlsp_col = nexus::BlueAlpha();
    wlsp_col.SetForceSolid(true);
    plate_logic->SetVisAttributes(wlsp_col);
    */
    
  }

  void LAttMeas::ConstructPMT(G4VPhysicalVolume* mother_physical) const
  {

    // Place the PMT
    const G4String pmt_name = "PMT";

    PmtR7378A* pmt = new PmtR7378A();
    pmt->Construct();
    G4LogicalVolume* pmt_logic = pmt->GetLogicalVolume();

    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateY(+90.*deg);

    new G4PVPlacement(  rot, 
                        G4ThreeVector(pmt_x_pos_, pmt_y_pos_, pmt_z_pos_), 
                        pmt_name, 
                        pmt_logic, 
                        mother_physical, 
                        false, 
                        0, 
                        true);

    return;

  }

  void LAttMeas::ConstructPlateHolders(G4VPhysicalVolume* mother_physical) const
  {

    G4double tolerance = 1*mm;

      const G4String auxo_name = "AUX_OUTTER";

      G4Box* encasing_solid =
          new G4Box(auxo_name,
              holder_dx_ / 2.,
              holder_dy_ / 2.,
              holder_dz_ / 2.);
     

      const G4String auxi_name = "AUX_INNER";

      G4Box* hole_solid =
          new G4Box(auxi_name,
              holder_dx_ / 2. + tolerance,
              holder_dy_ / 4.+ tolerance,
              plate_dz_ / 2.);

     const G4ThreeVector transVector = G4ThreeVector(0., holder_dy_/4.+ tolerance, 0.);

     G4SubtractionSolid* support_solid = 
        new G4SubtractionSolid ("PLATE_HOLDER", encasing_solid, hole_solid, new G4RotationMatrix{}, transVector);

      G4LogicalVolume* support_logic =
        new G4LogicalVolume(support_solid, 
                            materials::FR4(), 
                            "PLATE_HOLDER");

        new G4PVPlacement(nullptr,
                      G4ThreeVector( holder_x_pos_, 
                                     holder_y_pos_, 
                                     holder_z_pos_), 
                      "PLATE_HOLDER",
                      support_logic,
                      mother_physical, 
                      false, 
                      0, 
                      true);

  }

  void LAttMeas::ConstructPlateAdapter(G4VPhysicalVolume* mother_physical) const
  {

    G4double tolerance = 1*mm;

    const G4String upperT_name = "UPPER_T";
    G4Box* upperT_solid =
        new G4Box(upperT_name,
            adp_diameter_ /2. +tolerance,                         // tolerance should surpass each side of the tub
            (adp_thickn_ - adp_slot_thickn_)/2. +(tolerance/2.),  // tolerance should surpass the upper side of the slit
            adp_slot_width_/2.);

    const G4String lowerT_name = "LOWER_T";
    G4Box* lowerT_solid =
        new G4Box(lowerT_name,
            adp_hole_height_/2.,
            adp_slot_thickn_/2. +tolerance,       // tolerance into the upper_T part and tolerance 
                                                  // to surpass the lower side of the slit
            adp_slot_width_/2.);                  

    G4double relative_y_pos = -1.*((adp_thickn_ - adp_slot_thickn_)/2. +(tolerance/2.))   // minus half of the height of the upper T
                              -1.*(adp_slot_thickn_/2. +tolerance)                        // minus half of the height of the lower T
                              +tolerance;                                                 // plus one tolerance (overlap)

    G4UnionSolid* T_solid = 
      new G4UnionSolid( "AUXILIAR VOLUME", 
                        upperT_solid, 
                        lowerT_solid, 
                        new G4RotationMatrix{}, 
                        G4ThreeVector(0., relative_y_pos, 0.));

    G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateX(-90.*deg);

    G4double relative_z_pos =   -1.*adp_thickn_/2.
                                +adp_slot_thickn_
                                +((adp_thickn_ - adp_slot_thickn_)/2. +(tolerance/2.));

    const G4String adapter_name = "ADAPTER";
    G4double innerRadius = 0.*cm;
    G4double startAngle = 0.*deg;
    G4double spanningAngle = 360.*deg;

    G4Tubs* base_solid =
      new G4Tubs( adapter_name, 
                  innerRadius, 
                  adp_diameter_/2., 
                  adp_thickn_/2, 
                  startAngle, 
                  spanningAngle);

    G4SubtractionSolid* adapter_solid = 
      new G4SubtractionSolid( "ADAPTER", 
                              base_solid, 
                              T_solid, 
                              rot, 
                              G4ThreeVector(0., 0., relative_z_pos));

    G4LogicalVolume* adapter_logic =
      new G4LogicalVolume(adapter_solid, 
                          materials::FR4(), 
                          "PMT_ADAPTER");

    G4VisAttributes ref_case_col = nexus::WhiteAlpha();
    ref_case_col.SetForceSolid(true);
    adapter_logic->SetVisAttributes(ref_case_col);

    G4RotationMatrix* rot2 = new G4RotationMatrix();
    rot2->rotateX(+90.*deg);
    rot2->rotateY(+90.*deg);

    new G4PVPlacement(rot2,
                  G4ThreeVector( adapter_x_pos_, 
                                 adapter_y_pos_, 
                                 adapter_z_pos_), 
                  "PMT_ADAPTER",
                  adapter_logic,
                  mother_physical, 
                  false, 
                  0, 
                  true);
                  
  }



  G4bool LAttMeas::geometry_is_ill_formed()
  {
    return false;
  }

  G4ThreeVector LAttMeas::GenerateVertex(const G4String&) const{
    return G4ThreeVector(gen_vertex_x_, gen_vertex_y_, gen_vertex_z_);
  }

} //End namespace nexus
