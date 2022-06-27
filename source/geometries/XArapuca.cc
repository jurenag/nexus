#include "XArapuca.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "WLSPlate.h"
#include "HamamatsuS133606050VE.h"
#include "SiPMBoard.h"
#include "MomentumSD.h"
#include "RandomUtils.h"
#include "Visibilities.h"

#include <G4GenericMessenger.hh>
#include <G4UserLimits.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4ThreeVector.hh>
#include <G4VisAttributes.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4SDManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(XArapuca, GeometryBase)

namespace nexus{

  XArapuca::XArapuca():
  GeometryBase(), 
  ///Get internal reflector cavity dimensions from arxiv.org/abs/1912.09191
  internal_length_          (488.*mm),    ///X
  internal_width_           (100.*mm),     ///Z
  internal_thickn_          (8*mm),     ///Y
  plate_length_             (487.*mm),    ///X
  plate_width_              (93.*mm),     ///Z
  plate_thickn_             (3.5*mm),     ///Y
  case_thickn_              (1.*mm),      ///Get foil thickness from isoltronic.ch/assets/of-m-vikuiti-esr-app-guide.pdf
  df_thickn_                (1.*mm),
  gap_                      (0.5*mm),
  num_phsensors_            (24),
  ref_phsensors_supports_   (true), 
  with_boards_              (false),
  double_sided_             (true),
  collectors_are_reflective_(false),
  random_generation_vertex_ (true),
  remove_filters_           (false),
  path_to_dichroic_data_(""),
  world_extra_thickn_       (100.*cm)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/XArapuca/",
				"Control commands of geometry XArapuca.");
    
    G4GenericMessenger::Command& il_cmd =
      msg_->DeclareProperty("internal_length", internal_length_,
			    "Internal length of the X-Arapuca cavity.");
    il_cmd.SetUnitCategory("Length");
    il_cmd.SetParameterName("internal_length", false);
    il_cmd.SetRange("internal_length>0.");

    G4GenericMessenger::Command& iw_cmd =
      msg_->DeclareProperty("internal_width", internal_width_,
			    "Internal width of the X-Arapuca cavity.");
    iw_cmd.SetUnitCategory("Length");
    iw_cmd.SetParameterName("internal_width", false);
    iw_cmd.SetRange("internal_width>0.");

    G4GenericMessenger::Command& it_cmd =
      msg_->DeclareProperty("internal_thickn", internal_thickn_,
			    "Internal thickness of the X-Arapuca cavity.");
    it_cmd.SetUnitCategory("Length");
    it_cmd.SetParameterName("internal_thickness", false);
    it_cmd.SetRange("internal_thickness>0.");

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

    G4GenericMessenger::Command& ct_cmd =
      msg_->DeclareProperty("case_thickn", case_thickn_,
			    "Thickness of the XArapuca case.");
    ct_cmd.SetUnitCategory("Length");
    ct_cmd.SetParameterName("case_thickn", false);
    ct_cmd.SetRange("case_thickn>0.");

    G4GenericMessenger::Command& dft_cmd =
      msg_->DeclareProperty("dichroic_filter_thickn", df_thickn_,
			    "Thickness of the dichroic filter.");
    dft_cmd.SetUnitCategory("Length");
    dft_cmd.SetParameterName("df_thickn", false);
    dft_cmd.SetRange("df_thickn>0.");

    G4GenericMessenger::Command& g_cmd =
      msg_->DeclareProperty("gap", gap_,
			    "Gap between the WLS plate and the photosensors.");
    g_cmd.SetUnitCategory("Length");
    g_cmd.SetParameterName("gap", false);
    g_cmd.SetRange("gap>0.");

    G4GenericMessenger::Command& np_cmd =
      msg_->DeclareProperty("num_phsensors", num_phsensors_,
			    "Number of photosensors per long side.");
    np_cmd.SetParameterName("num_phsensors", false);
    np_cmd.SetRange("num_phsensors>=0"); 

    G4GenericMessenger::Command& rps_cmd =
      msg_->DeclareProperty("ref_phsensors_supports", ref_phsensors_supports_,
			    "Whether the photosensors supports are VIKUITI-coated.");

    G4GenericMessenger::Command& ds_cmd =
      msg_->DeclareProperty("double_sided", double_sided_,
			    "Whether the dichroic filter is set on both sides of the WLS plate.");

    G4GenericMessenger::Command& cr_cmd =
      msg_->DeclareProperty("collectors_are_reflective", collectors_are_reflective_,
			    "Whether the test collectors are reflective.");

    G4GenericMessenger::Command& rgv_cmd =
      msg_->DeclareProperty("random_generation_vertex", random_generation_vertex_,
			    "Whether the generation vertex is sampled randomly or not.");

    G4GenericMessenger::Command& rf_cmd =
      msg_->DeclareProperty("remove_filters", remove_filters_,
			    "Whether to remove the filters or not.");

    G4GenericMessenger::Command& ptdd_cmd =
      msg_->DeclareProperty("path_to_dichroic_data", path_to_dichroic_data_,
			    "Absolute path to the file containing the transmission data of the dichroic filter.");

    // When testing WLS plates with opticalprops::noAbsLength_, it is possible that a photon
    // gets trapped within the plate (below the critical angle) into an infinite-bouncing-loop.
    // For the case of an EJ286 plate with DUNE supercells dimensions, with this absorption 
    // length and immersed into LAr, some analysis revealed that particles that did not fall 
    // into this infinite loop had track lengths of less than one hundred meters.
    ul_ = new G4UserLimits();
    ul_->SetUserMaxTrackLength(100*m);
  }

  XArapuca::~XArapuca()
  {
      if(ul_) delete ul_;
  }


  void XArapuca::Construct()
  {

    if(geometry_is_ill_formed()){
      G4Exception("[XArapuca]", "Construct()", FatalException,
      "The given dimensions do not describe a feasible X-Arapuca.");
    }

    // The biggest volume is a vacuum box with the dimensions of the X-ARAPUCA device
    // plus 2*world_extra_thickn_, for each dimension. This volume is the effective 
    // world volume FOR THE NEXUS USER, but it is not the world volume of the overall
    // Geant4 application (afterwards nexus takes the user's biggest volume 
    // and places it inside another world volume whose dimensions are enough so as to 
    // fit the whole span of the biggest volume you implemented). Placed within this
    // vacuum box, there's a LAr box with dimensions of the X-ARAPUCA device plus 
    // world_extra_thickn_, along each dimension. The previous version of this class
    // simply implemented a LAr box as the biggest volume, without the need to encapsulate
    // it within a bigger volume (apart from that of nexus). However, as of the implementation
    // of the dichroic filter, the use of G4LogicalBorderSurface class is necessary. 
    // Objects of such class must be given two ordered physical volumes, so that a particle 
    // flying from the first one to the second one, interacts with such an optical surface. 
    // Within this context, it is convenient that we specify the LAr volume which encapsulates 
    // the X-ARAPUCA as the first physical volume, and the filter physical volume as the second 
    // volume, so that the filter properties are applied in both cases: 1) when a photon tries 
    // to enter the dichroic filter volume from outside the X-ARAPUCA, and 2) when a photon 
    // tries to enter the dichroic filter volume from inside the X-ARAPUCA (indeed, take into 
    // account that X-ARAPUCA different parts are placed within this LAr box, thus giving room 
    // to LAr gaps within the X-ARAPUCA cavity that belong to the LAr box mother volume). Since 
    // the nexus user seems not to have access to the physical placement of the biggest volume 
    // (in our case, the vacuum box) which is implemented in line 80 of 
    // source/base/DetectorConstruction.cc, AFTER calling GeometryBase::Construct() (i.e. your
    // geometry must be constructed before placing your biggest volume in nexus world volume,
    // so you cannot possibly have access to the physical placement of your biggest volume
    // when you Construct() it.), one way around this issue is adding this intermediate LAr box.   
    // LAr ///////////////////////////////////////////////////////////
    // Liquid argon box that contains all other volumes.
    
    // VACUUM CAPSULE
    const G4String world_name = "VACUUM_CAPSULE";

    G4Box* world_solid =
        new G4Box(world_name,
                (internal_length_ +case_thickn_)/2. +world_extra_thickn_,
                (internal_thickn_ +case_thickn_)/2. +world_extra_thickn_,
                (internal_width_ +case_thickn_)/2. +world_extra_thickn_);
                
    G4Material* vacuum =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

    G4LogicalVolume* world_logic = 
        new G4LogicalVolume(world_solid, vacuum, world_name, 0, 0, 0, true);
    world_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    this->SetLogicalVolume(world_logic);

    // LAr BOX
    const G4String LAr_box_name = "LAr";

    G4Box* LAr_box_solid =
      new G4Box(LAr_box_name,
                (internal_length_ +case_thickn_ +world_extra_thickn_)/2.,
                (internal_thickn_ +case_thickn_ +world_extra_thickn_)/2.,
                (internal_width_ +case_thickn_ +world_extra_thickn_)/2.);

    G4Material* LAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    LAr->SetMaterialPropertiesTable(opticalprops::paulucci_LAr());

    G4LogicalVolume* LAr_box_logic =
      new G4LogicalVolume(LAr_box_solid, LAr, LAr_box_name);
    LAr_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VPhysicalVolume* mother_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement (new G4RotationMatrix(), G4ThreeVector(0., 0., 0.), LAr_box_logic, 
                        LAr_box_name, world_logic, false, 0, true));

    ConstructWLSPlate(mother_physical);
    ConstructReflectiveBox(mother_physical);
    if(with_boards_) ConstructBoards(mother_physical);
    else ConstructPhotosensors(mother_physical);
    //ConstructCollectors(mother_physical); 
    if(!remove_filters_) ConstructDichroicFilters(mother_physical);

    return;
  }


  void XArapuca::ConstructWLSPlate(G4VPhysicalVolume* mother_physical) const
  { 

    WLSPlate* plate = new WLSPlate(plate_length_, plate_thickn_, plate_width_, opticalprops::EJ286(), false);
    //WLSPlate* plate = new WLSPlate(plate_length_, plate_thickn_, plate_width_, false);
    plate->Construct();
    G4LogicalVolume* plate_logic = plate->GetLogicalVolume();
    plate_logic->SetUserLimits(ul_);

    G4VisAttributes wlsp_col = nexus::BlueAlpha();
    wlsp_col.SetForceSolid(true);
    plate_logic->SetVisAttributes(wlsp_col);

    if (!plate_logic) {
      G4Exception("[XArapuca]", "ConstructWLSPlate()",
                  FatalException, "Null pointer to logical volume.");
    }

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), plate_logic->GetName(), 
                    plate_logic, mother_physical, false, 0, true);
    
    return;
  }
  

  void XArapuca::ConstructPhotosensors(G4VPhysicalVolume* mother_physical) const
  {
    HamamatsuS133606050VE sipm;
    sipm.SetReflectiveSupports(ref_phsensors_supports_);
    sipm.Construct();
    G4double sipm_thickn = sipm.GetThickness();
    G4LogicalVolume* sipm_logic_vol = sipm.GetLogicalVolume();

    G4VisAttributes sipm_col = nexus::Red();
    sipm_col.SetForceSolid(true);
    sipm_logic_vol->SetVisAttributes(sipm_col);

    if (!sipm_logic_vol) {
      G4Exception("[XArapuca]", "ConstructPhotosensors()",
                  FatalException, "Null pointer to logical volume.");
    }

    G4int phsensor_id = 0;

    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateX(-90.*deg);

    for (G4int i=0; i<num_phsensors_; ++i) {
      phsensor_id = i;
      G4ThreeVector pos(-plate_length_/2. + (0.5 + i) * plate_length_/num_phsensors_,
                        0.,
                        -plate_width_/2. - sipm_thickn/2. -gap_);

      new G4PVPlacement(rot, pos,
                        "S133606050VE_MPPC", sipm_logic_vol,
                        mother_physical, false, phsensor_id, true);
    }

    G4RotationMatrix* rot2 = new G4RotationMatrix();
    rot2->rotateX(90.*deg);

    for (G4int i=0; i<num_phsensors_; ++i) {
      phsensor_id = num_phsensors_ +i;
      G4ThreeVector pos(-plate_length_/2. + (0.5 + i) * plate_length_/num_phsensors_,
                        0.,
                        +plate_width_/2. +sipm_thickn/2. +gap_);

      new G4PVPlacement(rot2, pos,
                        "S133606050VE_MPPC", sipm_logic_vol,
                        mother_physical, false, phsensor_id, true);
    }
    
    return;
  }

  void XArapuca::ConstructBoards(G4VPhysicalVolume* mother_physical) const
  {

    SiPMBoard board1;
    board1.SetBaseID(0);
    board1.SetNumPhsensors(num_phsensors_);
    board1.SetReflectiveSupports(ref_phsensors_supports_);
    board1.Construct();
    G4LogicalVolume* board1_logic_vol = board1.GetLogicalVolume();

    G4double z_pos = (plate_width_/2.) +(board1.GetOverallThickness()/2.) +gap_;
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., z_pos),
                    "SIPMS_BOARD", board1_logic_vol, mother_physical, true, 0., true);

    SiPMBoard board2;
    board2.SetBaseID(24);
    board2.SetNumPhsensors(num_phsensors_);
    board2.SetReflectiveSupports(ref_phsensors_supports_);
    board2.Construct();
    G4LogicalVolume* board2_logic_vol = board2.GetLogicalVolume();

    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateX(+180.*deg);

    new G4PVPlacement(rot, G4ThreeVector(0., 0., -1.*z_pos),
                    "SIPMS_BOARD", board2_logic_vol, mother_physical, true, 0., true);

    return;

  }

  void XArapuca::ConstructReflectiveBox(G4VPhysicalVolume* mother_physical) const
  {
    // REFLECTIVE FOILS ///////////////////////////////////////////////////////
    
    //Based on arxiv.org/abs/1912.09191 and TDR vol.IX, the reflective foils
    //fully cover the internal cavity of the X-Arapuca.

    const G4String foil_name = "REF_FOIL";

    //Get reflective foil volume as subtraction solid from two boxes
    G4Box* aux_outter_box = 
      new G4Box("AUX_OUTTER_BOX",  
      internal_length_/2. +case_thickn_,
      internal_thickn_/2. +case_thickn_, 
      internal_width_/2.  +case_thickn_);  

    //Extra thickness to prevent boolean subtraction of solids with matching surfaces
    //See geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html#solids-made-by-boolean-operations
    G4double tolerance = 1.*cm;

/*
    G4SubtractionSolid* ref_foils_solid = nullptr;
    if(!double_sided_)
    {
        G4Box * aux_inner_box = 
        new G4Box("AUX_INNER_BOX",  
        internal_length_/2., 
        (internal_thickn_ +case_thickn_ +tolerance)/2.,
        internal_width_/2.);

        ref_foils_solid = 
        new G4SubtractionSolid("REFLECTIVE_FOIL", aux_outter_box, aux_inner_box, 
        new G4RotationMatrix(), G4ThreeVector(0., (case_thickn_+tolerance)/2., 0.));  

    }
    else
    {
        G4Box * aux_inner_box = 
        new G4Box("AUX_INNER_BOX",  
        internal_length_/2., 
        (internal_thickn_/2.) +case_thickn_ +(tolerance/2.),
        internal_width_/2.);

      ref_foils_solid = 
      new G4SubtractionSolid("REFLECTIVE_FOIL", aux_outter_box, aux_inner_box, 
      new G4RotationMatrix(), G4ThreeVector(0., 0., 0.));  
    }
    */

    G4Box* aux_aux_inner_box_1 = 
    new G4Box("AUX_INNER_BOX",  
    internal_length_/2., 
    internal_thickn_/2.,
    internal_width_/2.);

    G4SubtractionSolid* ref_foils_solid = nullptr;
    if(!double_sided_)
    {

        G4Box* aux_aux_inner_box_2 = 
        new G4Box("AUX_INNER_BOX",
        plate_length_/2.,
        (internal_thickn_+case_thickn_+tolerance)/2.,
        plate_width_/2.);

        G4UnionSolid* aux_inner_box = 
        new G4UnionSolid("AUX_INNER_BOX", aux_aux_inner_box_1, aux_aux_inner_box_2, 
        new G4RotationMatrix(), G4ThreeVector(0., (case_thickn_+tolerance)/2., 0.));

        ref_foils_solid = 
        new G4SubtractionSolid("REFLECTIVE_FOIL", aux_outter_box, aux_inner_box, 
        new G4RotationMatrix(), G4ThreeVector(0., 0., 0.));  

    }
    else
    {

        G4Box* aux_aux_inner_box_2 = 
        new G4Box("AUX_INNER_BOX",
        plate_length_/2.,
        (internal_thickn_+tolerance)/2. + case_thickn_,
        plate_width_/2.);

        G4UnionSolid* aux_inner_box = 
        new G4UnionSolid("AUX_INNER_BOX", aux_aux_inner_box_1, aux_aux_inner_box_2, 
        new G4RotationMatrix(), G4ThreeVector(0., 0., 0.));

        ref_foils_solid = 
        new G4SubtractionSolid("REFLECTIVE_FOIL", aux_outter_box, aux_inner_box, 
        new G4RotationMatrix(), G4ThreeVector(0., 0., 0.));  
    }

    G4LogicalVolume* ref_foils_logic = 
      new G4LogicalVolume(ref_foils_solid, materials::FR4(), foil_name);

    G4VisAttributes ref_foils_col = nexus::WhiteAlpha();
    ref_foils_col.SetForceSolid(true);
    ref_foils_logic->SetVisAttributes(ref_foils_col);

    //Now create the surface
    const G4String refsurf_name = "REF_SURFACE";
    G4OpticalSurface* refsurf_opsurf = 
      new G4OpticalSurface(refsurf_name, unified, polishedfrontpainted, dielectric_dielectric, 1);
    
    refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::VIKUITI());
    new G4LogicalSkinSurface("REF_FOIL_SURFACE", ref_foils_logic, refsurf_opsurf);   
    
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.),
          foil_name, ref_foils_logic, 
          mother_physical,
          false, 0, true);

    return;
  }

  void XArapuca::ConstructCollectors(G4VPhysicalVolume* mother_physical) const
  {

    // Box which is made out of the same material of the world
    // that encloses the XArapuca cavity, i.e. replaces the dichroic filter.

    const G4String collector_name = "COLLECTOR";

    G4Box* collector_solid = new G4Box(collector_name, plate_length_/2., 
                                    case_thickn_/2., plate_width_/2.);

    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    //  The collector will be coated with opticalprops::PerfectPhotonCollector()
    //  No need to set the optical properties of the collector volume
    G4LogicalVolume* collector_logic = 
                        new G4LogicalVolume(collector_solid, lAr, collector_name);

    /*
    G4VisAttributes collector_col = nexus::LightGreenAlpha();
    collector_col.SetForceSolid(true);
    collector_logic->SetVisAttributes(collector_col);
    */

    //  Detection coating
    G4OpticalSurface* collector_opsurf =
      new G4OpticalSurface("COLLECTOR_OPSURF", unified, polished, dielectric_metal);
    if(!collectors_are_reflective_){
        collector_opsurf->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());
    }
    else{
        collector_opsurf->SetMaterialPropertiesTable(opticalprops::PerfectPhotonReflector());
    }
    new G4LogicalSkinSurface("COLLECTOR_OPSURF", collector_logic, collector_opsurf);

    /*
    //  Sensitive Detector
    MomentumSD* msd = new MomentumSD("MomentumSD");
    G4SDManager* SDMPointer = G4SDManager::GetSDMpointer();
    if(!SDMPointer){
      G4Exception("[XArapuca]", "ConstructCollectors()", FatalException,
      "Could not retrieve the Sensitive Detector Manager pointer.");
    }
    SDMPointer->AddNewDetector(msd);
    collector_logic->SetSensitiveDetector(msd);
    */
    
    //  Place the collector(s)
    new G4PVPlacement(nullptr, G4ThreeVector(0., internal_thickn_/2. + case_thickn_/2., 0.), 
                        collector_name, collector_logic, mother_physical, false, 0, true);
    if(double_sided_)
    {
        new G4PVPlacement(nullptr, G4ThreeVector(0., -internal_thickn_/2. -case_thickn_/2., 0.), 
                        collector_name, collector_logic, mother_physical, false, 0, true);
    }

    return;

  }

  
  void XArapuca::ConstructDichroicFilters(G4VPhysicalVolume* mother_physical) const
  {
      
    // Dichroic filters. They have the right dimensions so as to enclose the 
    // X-ARAPUCA internal cavity.

    const G4String df_name = "DICHROIC_FILTER";

    G4Box* df_solid = new G4Box(df_name, plate_length_/2., 
                                    df_thickn_/2., plate_width_/2.);

    G4Material* df_substrate = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    // The filter will be implemented as a G4LogicalBorderSurface from the surrounding LAr physical volume
    // towards the filter physical volume. For photons that are transmited from outside the X-ARAPUCA towards
    // the dichroic filter volume, they shall not be further reflected or refracted. To avoid this, the 
    // dichroic filter substrate optical properties must match those of the surrounding LAr. This was checked 
    // in an alternative setup.
    df_substrate->SetMaterialPropertiesTable(opticalprops::paulucci_LAr());
    G4LogicalVolume* df_logic = 
                        new G4LogicalVolume(df_solid, df_substrate, df_name);

    G4VisAttributes df_col = nexus::LightGreenAlpha();
    df_col.SetForceSolid(true);
    df_logic->SetVisAttributes(df_col);
    
    //  Place the filters
    G4VPhysicalVolume* first_df_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, G4ThreeVector(0., internal_thickn_/2. +df_thickn_/2., 0.), 
                        df_name, df_logic, mother_physical, true, 0, true));

    // Add dichroic specifications
    if(path_to_dichroic_data_==""){
        G4Exception("[XArapuca]", "ConstructDichroicFilters()",
                    FatalException, "The path to the dichroic data file was not set.");
    }

    setenv("G4DICHROICDATA", path_to_dichroic_data_, 1);
    G4OpticalSurface* df_opsurf =   
        new G4OpticalSurface("DICHROIC_OPSURF", dichroic, polished, dielectric_dichroic);

    new G4LogicalBorderSurface("LAr->DICHROIC1", mother_physical, first_df_physical, df_opsurf);
    
    if(double_sided_){
        G4VPhysicalVolume* second_df_physical = dynamic_cast<G4VPhysicalVolume*>(
            new G4PVPlacement(nullptr, G4ThreeVector(0., -internal_thickn_/2. -df_thickn_/2., 0.), 
                            df_name, df_logic, mother_physical, true, 1, true));
        new G4LogicalBorderSurface("LAr->DICHROIC2", mother_physical, second_df_physical, df_opsurf);
    }

    return;
  }

  G4bool XArapuca::geometry_is_ill_formed()
  {
    // What we need to make sure here is that there's room enough within the XArapuca internal cavity
    // to allocate everything that we intending to put in it. To do so, we may find the span of the 
    // internal geometry along each axis, and then compare it to the dimensions of the internal cavity.

    G4double internal_geom_length_span, internal_geom_width_span, internal_geom_thickn_span;

    if(with_boards_){
        SiPMBoard board;
        internal_geom_length_span = std::max(plate_length_, board.GetBoardLength());
        internal_geom_width_span = plate_width_+(2.*gap_)+(2.*board.GetOverallThickness());
        internal_geom_thickn_span = std::max(plate_thickn_, board.GetOverallHeight());
    }
    else{
        HamamatsuS133606050VE sipm;
        internal_geom_length_span = plate_length_;  //If you check XArapuca::ConstructPhotosensors() you 
                                                    //can check that num_phsensors_ sensors are uniformly 
                                                    //distributed along the length of the plate. Therefore,
                                                    //in this case it is necessary to make one additional
                                                    //check: the plate must be long enough so as to fit
                                                    //num_phsensors_ sensors.

        if(num_phsensors_*sipm.GetWidth()>plate_length_) { return true; }

        internal_geom_width_span = plate_width_+(2.*gap_)+(2.*sipm.GetThickness());
        internal_geom_thickn_span = std::max(plate_thickn_, sipm.GetHeight());
    }

    G4bool check1 = internal_length_ >= internal_geom_length_span;
    G4bool check2 = internal_width_  >= internal_geom_width_span;
    G4bool check3 = internal_thickn_ >= internal_geom_thickn_span;

    return !(check1*check2*check3);
  }

  G4ThreeVector XArapuca::GenerateVertex(const G4String&) const{

    G4double tolerance = 0.1*mm;    // Short distance over the dichroic filter (DF)
                                    // from which photons are launched. Also, the 
                                    // width of the outter border projected over the
                                    // DF from which photons won't be launched
                                    // (Just see the implementation in x_pos and
                                    // z_pos below to understand its meaning)

    G4double y_pos = (internal_thickn_/2.) +df_thickn_ +tolerance;
    if(!random_generation_vertex_){
        return G4ThreeVector(0., y_pos, 0.);
    }
    else{
        G4double x_pos = UniformRandomInRange(-1.*plate_length_/2. +tolerance, 
                                                plate_length_/2. -tolerance);
        G4double z_pos = UniformRandomInRange(-1.*plate_width_/2. +tolerance, 
                                                plate_width_/2. -tolerance);
        return G4ThreeVector(x_pos, y_pos, z_pos);
    }
  }

} //End namespace nexus
