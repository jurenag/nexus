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
#include <G4ExtrudedSolid.hh>
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
  shape_code_                           (0                            ),    
  detach_DF_                            (false                        ),
  wlsp_DF_gap_                          (1.13   *mm                   ),  // Assuming a 8.26 mm tall SiPM (that's the height of the BroadcomAFBRS4N44P044M
                                                                          // SiPMs) and a 6 mm thick WLS plate, 1.13 mm is the required gap for the DF to 
                                                                          // fully enclose the LAr gap that's left between the WLS plate, the DF and the SiPMs.
  DF_substrate_thickn_                  (1.000  *mm                   ),
  DF_substrate_mpt_                     (opticalprops::FusedSilica()  ),  // This one is still hardcoded. Be careful to choose the one you want to
  //DF_substrate_mpt_                     (opticalprops::SCHOTT_B270()  ),// actually simulate, in accordance to the transmission curves you are setting.
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
  secondary_wls_attlength_              (-1.     *m                   ),
  cromophore_concentration_             (40.                          ),
  cryogenic_temperature_                (false                        ),
  reflective_foil_thickn_               (0.065  *mm                   ),   /// Got foil thickness from isoltronic.ch/assets/of-m-vikuiti-esr-app-guide.pdf
  remove_reflective_foil_               (false                        ),
  SiPM_code_                            (1                            ),
  num_phsensors_                        (30                           ),   /// This is seemingly the APEX baseline
  board_position_code_                  (1                            ),
  align_lower_edges_of_plate_and_SiPMs_ (false                        ),
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

    G4GenericMessenger::Command& shc_cmd =
      msg_->DeclareProperty("shape_code", shape_code_,
			    "The shape of the built APEX depends on this parameter (0 - rectangular, 1 - triangular).");
    shc_cmd.SetParameterName("shape_code", false);
    shc_cmd.SetRange("shape_code>=0");
    shc_cmd.SetRange("shape_code<=1");

    G4GenericMessenger::Command& ddf_cmd =
      msg_->DeclareProperty("detach_DF", detach_DF_,
			    "Whether to detach the DF from the WLS plate.");    

    G4GenericMessenger::Command& wlspdfg_cmd =
      msg_->DeclareProperty("wlsp_DF_gap", wlsp_DF_gap_,
			    "Thickness of the gap which is left between the DF and the WLS plate.");
    wlspdfg_cmd.SetUnitCategory("Length");
    wlspdfg_cmd.SetParameterName("wlsp_DF_gap", false);
    wlspdfg_cmd.SetRange("wlsp_DF_gap>0.");

    G4GenericMessenger::Command& dfst_cmd =
      msg_->DeclareProperty("DF_substrate_thickn", DF_substrate_thickn_,
			    "Thickness of the DF substrate.");
    dfst_cmd.SetUnitCategory("Length");
    dfst_cmd.SetParameterName("DF_substrate_thickn", false);
    dfst_cmd.SetRange("DF_substrate_thickn>0.");

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
			    "Constant (wavelength indepedent) attenuation length of the secondary WLShifter.");
    swlsal_cmd.SetUnitCategory("Length");
    // Allow negative values for this parameter, so that it can be used as a flag
    // to signal that G2P_FB118() should use its own attenuation length spectrum.

    G4GenericMessenger::Command& crco_cmd =
      msg_->DeclareProperty("cromophore_concentration", cromophore_concentration_,
			    "Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter (the WLS plate), in case G2P_FB118 is used.");
    crco_cmd.SetParameterName("cromophore_concentration", false);
    crco_cmd.SetRange("cromophore_concentration>0.");

    G4GenericMessenger::Command& crte_cmd =
      msg_->DeclareProperty("cryogenic_temperature", cryogenic_temperature_,
			    "Whether the secondary WLShifter is at cryogenic temperature or not. It only makes a difference if G2P_FB118 is used.");

    G4GenericMessenger::Command& rft_cmd =
      msg_->DeclareProperty("reflective_foil_thickn", reflective_foil_thickn_,
			    "Reflective foil thickness.");
    rft_cmd.SetUnitCategory("Length");
    rft_cmd.SetParameterName("reflective_foil_thickn", false);
    rft_cmd.SetRange("reflective_foil_thickn>0.");

    G4GenericMessenger::Command& rrf_cmd =
      msg_->DeclareProperty("remove_reflective_foil", remove_reflective_foil_,
			    "If true, the reflective foil is not constructed.");

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

    G4GenericMessenger::Command& aleopas_cmd =
      msg_->DeclareProperty("align_lower_edges_of_plate_and_SiPMs", align_lower_edges_of_plate_and_SiPMs_,
			    "It only makes a difference if board_position_code_ is set to 2 or 3. If set to true, the lower edges of the WLS plate and the SiPMs are aligned.");

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

    // Compute internal attributes
    // Open issue: The overall dimensions of the APEX device as
    // a function of its attributes are not computed here yet.
    // They are set here to 10 meters as a quick workaround.
    overall_length_ = 10.*m ;  
    overall_thickn_ = 10.*m ;
    overall_width_  = 10.*m ;

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
    // before placing your biggest volume in nexus world volume, so you cannot possibly 
    // have access to the physical placement of your biggest volume when you Construct() it.).
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    // VACUUM CAPSULE
    const G4String world_name = "VACUUM_CAPSULE";

    G4Box* world_solid =
        new G4Box(
          world_name,
          (overall_length_/2.)+world_extra_thickn_,
          (overall_thickn_/2.)+world_extra_thickn_,
          (overall_width_/2.) +world_extra_thickn_
        );
                
    G4Material* vacuum =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

    G4LogicalVolume* world_logic = 
        new G4LogicalVolume(
          world_solid,
          vacuum,
          world_name,
          0,
          0,
          0,
          true
        );

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
      new G4Box(
        sm_box_name,
        (overall_length_+world_extra_thickn_)/2.,
        (overall_thickn_+world_extra_thickn_)/2.,
        (overall_width_ +world_extra_thickn_)/2.
      );

    G4Material* sm_material = G4NistManager::Instance()->FindOrBuildMaterial(sm_name);
    sm_material->SetMaterialPropertiesTable(mpt_ptr);

    G4LogicalVolume* sm_box_logic =
      new G4LogicalVolume(
        sm_box_solid,
        sm_material,
        sm_box_name
      );
    sm_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VPhysicalVolume* mother_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
          new G4PVPlacement(
            new G4RotationMatrix(), 
            G4ThreeVector(0., 0., 0.), 
            sm_box_logic, 
            sm_box_name, 
            world_logic, 
            false,
            0,
            true
          )
        );

    ConstructWLSPlate(mother_physical);
    ConstructSiPMSAndBoard(mother_physical);
    if(!remove_reflective_foil_) ConstructReflectiveFoil(mother_physical);
    if(!remove_MLS_)
    {
      if(!detach_DF_)
      {
        ConstructAttachedDichroicFilter(mother_physical);
      }
      else
      {
        ConstructDetachedDichroicFilter(mother_physical);
      }
    }

    return;
  }

  void APEX::ConstructWLSPlate(G4VPhysicalVolume* mother_physical) const
  { 
    G4bool dimples_at_z_minus = false;
    G4bool dimples_at_z_plus = false;

    if(shape_code_==0 && with_dimples_)
    {
      // Note that board_position_code_ is limited to >=1 
      // via a G4GenericMessenger::Command
      if(board_position_code_>=2)
      {
        dimples_at_z_minus = true;

        if(board_position_code_>=3)
        {
          dimples_at_z_plus = true;
        }
      }
    }

    // The ternary operators in the first arguments of the WLSPlate
    // constructor fix the fact that the WLSPlate class aligns the 
    // base of the triangular plate (if shape_code_ equals 1) with
    // the Z-axis, while the APEX class aligns the base of the
    // triangular plate with the X-axis. Take into account that the
    // documentation in APEX.h says that, for shape_code_==1,
    // plate_length_ should give the length of the triangular plate base.
    WLSPlate* plate = new WLSPlate(
      shape_code_==0 ? plate_length_ : plate_width_,
      plate_thickn_,
      shape_code_==0 ? plate_width_ : plate_length_,
      opticalprops::G2P_FB118(
        cromophore_concentration_,
        secondary_wls_attlength_,
        // WLSp_rindex_, No longer used here, this is related to an open issue (search for occurrences of 'Open issue')
        cryogenic_temperature_,
        true
      ),
      //opticalprops::EJ286(secondary_wls_attlength_),
      shape_code_,
      false,
      false,                // dimples_at_x_plus_
      false,                // dimples_at_x_minus_
      dimples_at_z_plus,    // dimples_at_z_plus_
      dimples_at_z_minus,   // dimples_at_z_minus_
      dimple_type_, 
      num_phsensors_, 
      flat_dimple_width_, 
      flat_dimple_depth_, 
      curvy_dimple_radius_
    );

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

    G4RotationMatrix* rot = new G4RotationMatrix();
    if(shape_code_==1)
    { 
      // The rotation about the X axis compensates the fact that the extrusion of the
      // triangular polygon performed by the WLSPlate class is aligned with the Z axis.
      // For more information check the WLSPlate::ConstructWLSPlate() method
      // documentation. The rotation about the Y axis is needed to align the
      // the triangular-plate base with the X axis.
      rot->rotateY(90.*deg);
      rot->rotateX(-90.*deg);
    }

    new G4PVPlacement(
      rot,
      // For the triangular plate case, the WLSPlate class places the plate in the
      // centroid of the (isosceles) triangular polygon. Correcting its position to
      // make its center coincide with half the triangle height, is more convenient
      // to reuse the rectangular-plate code. Particularly, the SiPMs-board code
      // and the vertex generation code.
      shape_code_==0 ? G4ThreeVector(0., 0., 0.) : G4ThreeVector(0., 0., -plate_width_/6.),
      plate_logic->GetName(), 
      plate_logic,
      mother_physical,
      false,
      0,
      true
    );
    
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

    // Board in the middle of a large face, only applicable to the shape_code_==0 case
    if(board_position_code_==1 && shape_code_==0)
    {
      sipm_rot->rotateX(0.0*deg);
      base_pos.set(
        (-1.*board_length_/2.) + (0.5*board_length_/num_phsensors_),
        -1.*(plate_thickn_/2.)-1.*(sipm_thickn/2.)-gap_,    // Note that what's placed in the global origin of 
        0.                                                  // coordinates is the plate, not the reflective foil. 
      );

      G4int phsensor_id = 0;
      for (G4int i=0; i<num_phsensors_; ++i)
      {
        new G4PVPlacement(
          sipm_rot,
          base_pos+G4ThreeVector(i*board_length_/num_phsensors_, 0., 0.),
          sipm->GetModel(),
          sipm_logic_vol,
          mother_physical,
          true,
          phsensor_id,
          true
        );
        phsensor_id += 1;
      }
    }
    else  // board_position_code_ is 2 or 3
    {
      G4double sipms_y_pos = 0.;
      if(align_lower_edges_of_plate_and_SiPMs_) 
      {
        // Assuming that the WLS plate is placed at the origin of coordinates
        sipms_y_pos = (-1.*plate_thickn_/2.)+(sipm->GetTransverseDim()/2.);
      }

      sipm_rot->rotateX(-90.*deg);
      base_pos.set(
        (-1.*board_length_/2.) + (0.5*board_length_/num_phsensors_),
        sipms_y_pos,
        -1.*(plate_width_/2.)-1.*(sipm_thickn/2.)-gap_
      );

      G4int phsensor_id = 0;
      for (G4int i=0; i<num_phsensors_; ++i) {
        new G4PVPlacement(
          sipm_rot,
          base_pos+G4ThreeVector(i*board_length_/num_phsensors_, 0., 0.),
          sipm->GetModel(),
          sipm_logic_vol,
          mother_physical,
          true,
          phsensor_id,
          true
        );
        phsensor_id += 1;
      }

      if(board_position_code_>=3 && shape_code_==0)
      {
        G4RotationMatrix* sipm_rot_2 = new G4RotationMatrix();
        G4ThreeVector base_pos_2;

        sipm_rot_2->rotateX(+90.*deg);
        base_pos_2.set(
          (-1.*board_length_/2.) + (0.5*board_length_/num_phsensors_),
          sipms_y_pos,
          (plate_width_/2.)+(sipm_thickn/2.)+gap_
        );

        phsensor_id = 0;
        for (G4int i=0; i<num_phsensors_; ++i) {
          new G4PVPlacement(
            sipm_rot_2,
            base_pos_2+G4ThreeVector(i*board_length_/num_phsensors_, 0., 0.),
            sipm->GetModel(),
            sipm_logic_vol,
            mother_physical,
            true,
            phsensor_id,
            true
          );
          phsensor_id += 1;
        }
      }
    }

    // Then construct the board (Vikuiti-coated FR4 piece)
    G4String board_name = "BOARD";

    G4double sipm_height = sipm->GetTransverseDim();
    G4double board_thickn = 1.*mm;

    G4Box* board_solid =
        new G4Box(
          board_name,
          board_length_/2., 
          sipm_height/2.,     // Board height matches that of the SiPMs
          board_thickn/2.
        );

    G4LogicalVolume* board_logic = 
        new G4LogicalVolume(
          board_solid,
          materials::FR4(),
          board_name
        );

    G4VisAttributes board_col = nexus::White();
    //board_col.SetForceSolid(true);
    board_logic->SetVisAttributes(board_col);

    //VIKUITI coating for the board
    const G4String bc_name = "BOARD_COATING";
    G4OpticalSurface* board_coating = 
      new G4OpticalSurface(
        bc_name,
        unified,
        ground,
        dielectric_metal,
        1
      );
    
    board_coating->SetMaterialPropertiesTable(opticalprops::Vikuiti());
    new G4LogicalSkinSurface(
      bc_name,
      board_logic,
      board_coating
    ); 

    G4RotationMatrix* board_rot = new G4RotationMatrix();
    G4ThreeVector board_pos;

    // Board in the middle of a large face, only applicable to the shape_code_==0 case
    if(board_position_code_==1 && shape_code_==0)
    {
      board_rot->rotateX(90.0*deg);
      board_pos.set(
        0.,
        -1.*(plate_thickn_/2.)-gap_
        -sipm_thickn-1.*(board_thickn/2.),  // Note that what's placed in the global origin of 
                                            // coordinates is the plate, not the reflective foil. 
        0.
      );

      //Place it
      new G4PVPlacement(
        board_rot,
        board_pos,
        "COATED_BOARD",
        board_logic,
        mother_physical,
        false,
        0,
        true
      );
    }
    else  // board_position_code_ is 2 or 3
    {
      G4double board_y_pos = 0.;
      if(align_lower_edges_of_plate_and_SiPMs_) 
      {
        // Assuming that the WLS plate is placed at the origin of coordinates
        board_y_pos = (-1.*plate_thickn_/2.)+(sipm->GetTransverseDim()/2.);
      }

      board_rot->rotateX(0.0*deg);
      board_pos.set(
        0.,
        board_y_pos,
        -1.*(plate_width_/2.)-gap_
        -sipm_thickn-1.*(board_thickn/2.)
      );
      //Place it
      new G4PVPlacement(
        board_rot,
        board_pos,
        "COATED_BOARD",
        board_logic,
        mother_physical,
        false,
        0,
        true
      );

      if(board_position_code_>=3 && shape_code_==0)
      {
        G4RotationMatrix* board_rot_2 = new G4RotationMatrix();
        G4ThreeVector board_pos_2;

        board_rot_2->rotateX(0.0*deg);
        board_pos_2.set(
          0.,
          board_y_pos,
          (plate_width_/2.)+gap_
          +sipm_thickn+(board_thickn/2.)
        );
        //Place it
        new G4PVPlacement(
          board_rot_2,
          board_pos_2,
          "COATED_BOARD",
          board_logic,
          mother_physical,
          false,
          0,
          true
        );
      }
    }
    return;
  }

  void APEX::ConstructReflectiveFoil(G4VPhysicalVolume* mother_physical) const
  {
    // I probably introduced some technical debt here when extending this
    // method to support the triangular-APEX case. This method should be
    // revised and refactored. 

    const G4String ref_foil_name = "REF_FOIL";

    // The reflective foil covers every face of the plate but one
    // Get its volume as a subtraction solid from two boxes
    G4double tolerance = 1.*mm; // To prevent matching surfaces in the boolean subtraction
    G4double outer_box_half_thickn = (plate_thickn_+reflective_foil_thickn_)/2.;

    // Extra thickness to prevent boolean subtraction of solids with matching surfaces
    // See geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html#solids-made-by-boolean-operations
    G4double inner_box_half_thickn = (plate_thickn_/2.)+tolerance;

    G4VSolid* aux_outer_box = nullptr;
    G4VSolid* aux_inner_box = nullptr;

    if(shape_code_==0){
      aux_outer_box = dynamic_cast<G4VSolid*>(
        new G4Box(
          "AUX_OUTER_BOX", 
          (plate_length_  + (2.*reflective_foil_thickn_))/2., 
          outer_box_half_thickn,
          (plate_width_  + (2.*reflective_foil_thickn_))/2.
        )
      );

      aux_inner_box = dynamic_cast<G4VSolid*>(
        new G4Box(
          "AUX_INNER_BOX",
          plate_length_/2.,
          inner_box_half_thickn,
          plate_width_/2.
        )
      );
    }
    else if(shape_code_==1)
    {
      std::vector<G4TwoVector> outer_prism_base;

      // Point A of sketch of the APEX::GetTriangularPlatePrismBase() method, but
      // displaced in X and Y axes by one or two times the thickness of the
      // reflective foil. Displacing all of them at the same time by one (or two)
      // times the thickness of the reflective foil gives visualization problems,
      // and we would not get the exact specified thickness of the reflective foil,
      // anyway. I.e. implementing such foil thickness exactly for the two edges
      // which are not aligned with any axes, would involve some geometric calculations.
      // It is important, though, that the Y-points are only displaced by one
      // thickness of the reflective foil, since that's the side on which the
      // windows for the SiPMs are carved (and such carvings are computed, further
      // in the code, assuming that the reflective foil is exactly
      // reflective_foil_thickn_ thick).
      outer_prism_base.push_back(
        G4TwoVector(
          -1.*plate_length_/2. -(2.*reflective_foil_thickn_),
          -1.*plate_width_/2. -reflective_foil_thickn_
        )
      );
      
      // Point B
      outer_prism_base.push_back(
        G4TwoVector(
          0.,
          plate_width_/2. +(2.*reflective_foil_thickn_)
        )
      );

      // Point C
      outer_prism_base.push_back(
        G4TwoVector(
          plate_length_/2. +(2.*reflective_foil_thickn_),
          -1.*plate_width_/2. -reflective_foil_thickn_
        )
      );

      aux_outer_box = dynamic_cast<G4VSolid*>(
        new G4ExtrudedSolid(
          "AUX_OUTER_BOX",
          outer_prism_base,
          outer_box_half_thickn
        )
      );

      aux_inner_box = dynamic_cast<G4VSolid*>(
        new G4ExtrudedSolid(
          "AUX_INNER_BOX",
          this->GetTriangularPlatePrismBase(),
          inner_box_half_thickn
        )
      );
    }
    else
    {
      G4Exception("[APEX]", "ConstructReflectiveFoil()",
                  FatalException, "The given shape code is not recognized.");
    }

    G4SubtractionSolid* ref_foil_solid = new G4SubtractionSolid(
      ref_foil_name,
      aux_outer_box,
      aux_inner_box,
      nullptr,
      G4ThreeVector(
        0.,
        // The thickness of the rectangular-APEX solids are already aligned
        // with the Y axis, but the triangular-APEX solids are aligned with
        // with the extrusion direction (i.e. the Z axis).
        shape_code_==0 ? (reflective_foil_thickn_/2.)+tolerance : 0.,
        shape_code_==0 ? 0. : -(reflective_foil_thickn_/2.)-tolerance
      )
    );

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
    G4Box* dummy_sipm =  new G4Box(
      "DUMMY_SIPM", 
      sipm_transverse_dim/2., 
      thickness_of_dummy_sipm/2., // Setting here the reflective-foil thickness plus some tolerance so that:
                                  //  1)  if board_position_code_==1, the carved hole is a pass-through hole
                                  //  2)  if board_position_code_>=2, we prevent matching surfaces in the boolean subtraction
                                  //      In this second case, the value of the tolerance actually matters. It must be big 
                                  //      enough so as to prevent matching surfaces, but small enough so as to not carve too 
                                  //      much the horizontal portion of the reflective foil.
      sipm_transverse_dim/2.
    );

    G4MultiUnion* reflective_foil_holes = new G4MultiUnion("REF_FOIL_HOLES");

    G4double pos;
    G4Transform3D* transform_ptr = nullptr;
    G4RotationMatrix* rot = new G4RotationMatrix();

    // The only case where the rotation is not needed is for
    // the rectangular APEX when the board_position_code_ is 1.
    if(board_position_code_!=1 || shape_code_==1)
    {
      rot->rotateX(90.0*deg);
    }

    for(G4int i=0; i<num_phsensors_; i++){

      pos = (-1.*board_length_/2.) + ((0.5 + i)*board_length_/num_phsensors_);    
      transform_ptr = new G4Transform3D(*rot, G4ThreeVector(pos, 0., 0.));
      reflective_foil_holes->AddNode(*dummy_sipm, *transform_ptr);    

    }

    reflective_foil_holes->Voxelize();

    G4ThreeVector vec = G4ThreeVector(
      0., 
      -1.*plate_thickn_/2., // Minus half the thickness of AUX_OUTER_BOX 
                            // plus half the reflective-foil thickness
      0.
    );

    G4double sipms_y_pos = 0.;
    if(align_lower_edges_of_plate_and_SiPMs_) 
    {
      // I was expecting an error in this case, due to matching surfaces in the boolean subtraction 
      // of ref_foil_solid minus reflective_foil_holes. Particularly, since the dummy_sipm which
      // creates the holes is slightly thicker than the reflective foil (there is an explanation why
      // this is done in the comments of the dummy_sipm definition), I was expecting that for the
      // the case where align_lower_edges_of_plate_and_SiPMs_ is true, the lower face of the 
      // dummy_sipm geometry would partially match the upper face of the bottom of the reflective foil.
      // However, I saw no warning while Geant4 builds the geomeotry, and no error while running for
      // a million photons.
      sipms_y_pos = (-1.*plate_thickn_/2.)+(sipm_transverse_dim/2.);
    }

    G4RotationMatrix* rot2 = new G4RotationMatrix();

    if(board_position_code_!=1 && shape_code_==0)
    {
      vec = G4ThreeVector(
        0.,
        sipms_y_pos,
        -1.*(plate_width_/2.)-1.*(reflective_foil_thickn_/2.)     // Minus half the width of the plate
      );                                                          // minus half the reflective-foil thickness
    }
    else if(shape_code_==1)
    {
      vec = G4ThreeVector(
        0.,
        -1.*(plate_width_/2.)-1.*(reflective_foil_thickn_/2.),    // Same coordinates as for the rectangular case,
        -sipms_y_pos                                              // but inverting Y and Z since the reflective foil
      );                                                          // is not rotated until its placement.

      rot2->rotateX(-90.*deg);
    }

    ref_foil_solid = new G4SubtractionSolid(
      ref_foil_name, 
      ref_foil_solid,
      reflective_foil_holes, 
      rot2,
      vec
    );

    // If board_position_code_ is 3 for the rectangular APEX
    // case, then also carve the holes for a second strip of SiPMs
    if(shape_code_==0 && board_position_code_>=3){
      
      G4ThreeVector vec_2 = G4ThreeVector(
        0.,
        sipms_y_pos,
        (plate_width_/2.)+(reflective_foil_thickn_/2.)    // Minus half the width of the plate
      );                                                  // minus half the reflective-foil thickness

      ref_foil_solid = new G4SubtractionSolid(
        ref_foil_name, 
        ref_foil_solid,
        reflective_foil_holes, 
        nullptr,
        vec_2
      );
    }
    
    G4LogicalVolume* ref_case_logic = 
      new G4LogicalVolume(
        ref_foil_solid,
        materials::FR4(),
        ref_foil_name
      );

    // Set its color for visualization purposes
    G4VisAttributes ref_case_col = nexus::WhiteAlpha();
    //ref_case_col.SetForceSolid(true);
    ref_case_logic->SetVisAttributes(ref_case_col);

    //Now create the reflectivie optical surface
    const G4String ref_surf_name = "REF_SURFACE";
    G4OpticalSurface* refsurf_opsurf = 
      new G4OpticalSurface(
        ref_surf_name,
        unified,
        ground,
        dielectric_metal,
        1
      );
    
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

    refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
    new G4LogicalSkinSurface(
      ref_surf_name,
      ref_case_logic,
      refsurf_opsurf
    );

    G4RotationMatrix* rot3 = new G4RotationMatrix();
    if(shape_code_==1)
    { 
      // In case a triangular plate (and so, a triangular relfective foil) is used,
      // rotate it about the X axis so that its thickness is aligned with the Z axis.
      // Note that this volume was the result of an extrusion, which happens along
      // the Z axis, by definition of G4ExtrudedSolid.
      rot3->rotateX(-90.*deg);
    }  
    
    new G4PVPlacement(
      rot3,
      G4ThreeVector(
        0.,
        -1.*reflective_foil_thickn_/2.,
        0.
      ),
      ref_foil_name,
      ref_case_logic,
      mother_physical,
      false,
      0,
      true
    );

    return;
  }

  void APEX::ConstructAttachedDichroicFilter(G4VPhysicalVolume* mother_physical) const
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

    // This function is called by APEX::Construct() only if !remove_MLS_ and !detach_DF_
    
    G4VSolid* MLS_half_solid = nullptr;
    if(shape_code_==0)
    {
      MLS_half_solid = dynamic_cast<G4VSolid*>(
        new G4Box(
          "AUX",
          plate_length_/2.,
          MLS_thickn_/4.,
          plate_width_/2.
        )
      );
    }
    else if(shape_code_==1)
    {
      // Extrude the triangular
      MLS_half_solid = dynamic_cast<G4VSolid*>(
        new G4ExtrudedSolid(
          "AUX",
          this->GetTriangularPlatePrismBase(),
          MLS_thickn_/4.
        )
      );
    }
    else
    {
      G4Exception("[APEX]", "ConstructAttachedDichroicFilter()",
                  FatalException, "The given shape code is not recognized.");
    }

    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    // Open issue: Some time ago, the refractive index of the opticalprops::G2P_FB118() G4MaterialPropertiesTable pointer
    // was set via a G4double parameter, i.e. its refractive index was constant and fixed to the given G4double input.
    // The refractive index of the MLS material (which should be set to that of the wavelength shifting plate according
    // to the second alternative of DF implementation explained above) was set to the WLSp_rindex_ attribute, which was,
    // at the same time, set to the refractive index of the opticalprops::G2P_FB118() G4MaterialPropertiesTable. However,
    // at some point (before 13/03/2025), the refractive index of the opticalprops::G2P_FB118() G4MaterialPropertiesTable
    // was set to some hardcoded array which depends on the wavelength (which was done for the sake of realism), while the
    // MLS rindex remained set to the constant value given to the WLSp_rindex_ array (which is no longer the refractive
    // index of the WLS plate). Although both refractive indices are similar, they are not exactly the same. Thus, to
    // stick to a more realistic implementation, we should try to retrieve this hardcoded wavelength-dependent refractive
    // index from opticalprops::G2P_FB118(), and give it to opticalprops::TunableRIMat() here.

    mat->SetMaterialPropertiesTable(opticalprops::TunableRIMat(WLSp_rindex_));  // Change this to MLS_rindex_ if you
                                                                                // want to go for the first alternative
                                                                                // of DF implementation explained above

                                                                                // This opticalprops::TunableRIMat() has 
                                                                                // no defined absorption length, just a
                                                                                // defined refractive index. The direct
                                                                                // comparison with respect to the 
                                                                                // detach_DF_==false case is, then, only 
                                                                                // possible if the DF_substrate_mpt_
                                                                                // attributes points to an MPT whose 
                                                                                // absorption length is also undefined
                                                                                // or practically infinte. P.e. 
                                                                                // opticalprops::FusedSilica() works
                                                                                // for our case, since its absorption
                                                                                // length below 6.5 eV (i.e. above 190
                                                                                // nm) is set to opticalprops::noAbsLength_.

    G4LogicalVolume* MLS_half_logic = new G4LogicalVolume(
      MLS_half_solid,
      mat,
      "MLS_HALF"
    );
            
    G4VisAttributes MLS_col = nexus::BloodRedAlpha();
    //MLS_col.SetForceSolid(true);
    MLS_half_logic->SetVisAttributes(MLS_col);

    G4RotationMatrix* rot = new G4RotationMatrix();
    if(shape_code_==1)
    { 
      // In case a triangular plate (and so, a triangular MLS) is used, rotate it
      // about the X axis so that its thickness is aligned with the Z axis. Note
      // that this volume was the result of an extrusion, which happens along the
      // Z axis, by definition of G4ExtrudedSolid.
      rot->rotateX(-90.*deg);
    }

    // Place the MLS
    G4VPhysicalVolume* MLS_first_half = dynamic_cast<G4VPhysicalVolume*>(   // This is the outermost one
        new G4PVPlacement(
          rot,
          G4ThreeVector(
            0.,
            plate_thickn_/2.
            +MLS_thickn_/2.         // Note that the thickness of MLS_half_solid is MLS_thickn_/2
            +MLS_thickn_/4.,
            0.
          ), 
          "FIRST_MLS_HALF",
          MLS_half_logic,
          mother_physical,
          true,
          0,
          true
        )
      );

    G4VPhysicalVolume* MLS_second_half = dynamic_cast<G4VPhysicalVolume*>(  // This is the internal one
        new G4PVPlacement(
          rot,
          G4ThreeVector(
            0.,
            plate_thickn_/2.
            +MLS_thickn_/4.,
            0.
          ), 
          "SECOND_MLS_HALF",
          MLS_half_logic,
          mother_physical,
          true,
          1,
          true
        )
      );

    // Check that there's dichroic information for ingoing (wrt APEX) photons
    if(path_to_inwards_dichroic_data_==""){
        G4Exception("[APEX]", "ConstructAttachedDichroicFilter()",
                    FatalException, "The path to the inwards dichroic data file was not set.");
    }

    // Check that there's dichroic information for outgoing photons
    if(path_to_outwards_dichroic_data_==""){
        G4Exception("[APEX]", "ConstructAttachedDichroicFilter()",
                    FatalException, "The path to the outwards dichroic data file was not set.");
    }

    // Construct the ingoing optical surface
    setenv("G4DICHROICDATA", path_to_inwards_dichroic_data_, 1);
    G4OpticalSurface* df_inwards_opsurf =
        new G4OpticalSurface(           // G4OpticalSurface constructor loads the
          "DICHROIC_INWARDS_OPSURF",    // dichroic information from the file which
          dichroic,                     // is currently pointed to by the environment
          polished,                     // variable G4DICHROICDATA
          dielectric_dichroic
        );

    // Construct the outgoung optical surface
    setenv("G4DICHROICDATA", path_to_outwards_dichroic_data_, 1);   // Note that, if you did not compile the modified version of G4 code, 
                                                                    // then different G4 dichroic data cannot be loaded. Instead, the first 
                                                                    // one (i.e. the one I am setting from path_to_inwards_dichroic_data_), 
                                                                    // is the one that will apply for every dichroic boundary in the simulation.
    G4OpticalSurface* df_outwards_opsurf =   
        new G4OpticalSurface(
          "DICHROIC_OUTWARDS_OPSURF", 
          dichroic, 
          polished, 
          dielectric_dichroic
        );

    // Endow the MLS_first_half->MLS_second_half surface with the ingoing optical surface
    new G4LogicalBorderSurface(
      "MLS1->MLS2",
      MLS_first_half,
      MLS_second_half,
      df_inwards_opsurf
    );

    // Endow the MLS_second_half->MLS_first_half surface with the outgoing optical surface
    new G4LogicalBorderSurface(
      "MLS2->MLS1",
      MLS_second_half,
      MLS_first_half,
      df_outwards_opsurf
    );

    // pTP coating
    if(!remove_coating_)
    {
        G4VSolid* coating_solid = nullptr;

        if(shape_code_==0)
        {
          coating_solid = dynamic_cast<G4VSolid*>(
            new G4Box(
              "COATING",
              plate_length_/2.,
              coating_thickn_/2.,
              plate_width_/2.
            )
          );
        }
        else // It has been checked before (within this same method) that shape_code_ is either 0 or 1
        {
          // Extrude the same polygon (prism base) that we used for the MLS solids
          coating_solid = dynamic_cast<G4VSolid*>(
            new G4ExtrudedSolid(
              "COATING",
              this->GetTriangularPlatePrismBase(),
              coating_thickn_/2.
            )
          );
        }

        G4Material* coating_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
        coating_mat->SetMaterialPropertiesTable(opticalprops::PTP(coating_rindex_));
        G4LogicalVolume* coating_logic = new G4LogicalVolume(
          coating_solid,
          coating_mat,
          "COATING"
        );

        G4VisAttributes coating_col = nexus::TitaniumGreyAlpha();
        coating_col.SetForceSolid(true);
        coating_logic->SetVisAttributes(coating_col);

        // Place the coating
        G4VPhysicalVolume* coating_physical = dynamic_cast<G4VPhysicalVolume*>(
          new G4PVPlacement(
            // Using the same rotation matrix as the one used for the MLS
            rot,
            G4ThreeVector(
              0.,
              plate_thickn_/2.
              +MLS_thickn_
              +coating_thickn_/2.,
              0.
            ),
            "COATING",
            coating_logic,
            mother_physical,
            false,
            0,
            true
          )
        );

        // Make the LAR-coating interface rough, so that photons cannot be trapped within the coating
        G4OpticalSurface* coating_rough_surf = new G4OpticalSurface(
          "COATING_ROUGH_SURFACE",
          glisur,
          ground,
          dielectric_dielectric,
          .01                     // 0.01 is the polish value for glisur model that was
        );                        // measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
                                  // This is the best reference we have, since both PTP and
                                  // TPB are the result of an evaporation+deposition process
                
        new G4LogicalBorderSurface(
          "SURROUNDINGS->COATING",
          mother_physical,
          coating_physical,
          coating_rough_surf
        );

        new G4LogicalBorderSurface(
          "COATING->SURROUNDINGS",
          coating_physical,
          mother_physical,
          coating_rough_surf
        );
        // We will also add roughness for the coating->MLS interface, but only with such ordering. The alternative case takes place
        // when the photon travels from the MLS to the coating. The MLS is supposed to be polished, so the photon may not see a rough 
        // surface.
        new G4LogicalBorderSurface(
          "COATING->MLS",
          coating_physical,
          MLS_first_half,
          coating_rough_surf
        );
    }
    return;
  }

  void APEX::ConstructDetachedDichroicFilter(G4VPhysicalVolume* mother_physical) const
  {
  
    // ------------------------------------ DICHROIC FILTER MODEL ------------------------------------
    //
    // The DF model implemented here is the same as the one implemented in 
    // XArapuca::ConstructDichroicAssemblies . 
    //
    //                                         LAr (Cryostat)
    //
    //  --------------------------------------------------------------------------------------------
    //                                            coating
    //  --------------------------------------------------------------------------------------------
    //                                                               
    //                                       substrate (n=n_subs)
    //                                                               
    //  -----------------------------------G4LogicalBorderSurfaces----------------------------------
    //                                                               
    //                                   MLS (with substrate rindex) (n=n_subs)
    //  ____________________________________________________________________________________________
    //                                                
    //                                 LAr (X-ARAPUCA internal cavity)
    //
    //  ____________________________________________________________________________________________
    //
    //                                            WLS plate
    //
    //  ____________________________________________________________________________________________
    //
    //
    // This function is called by APEX::Construct() only if !remove_MLS_ and detach_DF_
    

    // DF substrate
    G4VSolid* DF_substrate_solid = nullptr;
    if(shape_code_==0)
    {
      DF_substrate_solid = dynamic_cast<G4VSolid*>(
        new G4Box(
          "DICHROIC_FILTER_SUBSTRATE", 
          plate_length_/2.,
          DF_substrate_thickn_/2.,
          plate_width_/2.
        )
      );
    }
    else if(shape_code_==1)
    {
      // Extrude the triangular
      DF_substrate_solid = dynamic_cast<G4VSolid*>(
        new G4ExtrudedSolid(
          "DICHROIC_FILTER_SUBSTRATE",
          this->GetTriangularPlatePrismBase(),
          DF_substrate_thickn_/2.
        )
      );
    }
    else
    {
      G4Exception("[APEX]", "ConstructAttachedDichroicFilter()",
                  FatalException, "The given shape code is not recognized.");
    }

    G4Material* DF_substrate_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    DF_substrate_mat->SetMaterialPropertiesTable(DF_substrate_mpt_);

    G4LogicalVolume* DF_substrate_logic = new G4LogicalVolume(
      DF_substrate_solid,
      DF_substrate_mat,
      "DICHROIC_FILTER_SUBSTRATE"
    );

    G4RotationMatrix* rot = new G4RotationMatrix();
    if(shape_code_==1)
    { 
      // In case a triangular plate (and so, a triangular DF) is used, rotate it
      // about the X axis so that its thickness is aligned with the Z axis. Note
      // that this volume was the result of an extrusion, which happens along the
      // Z axis, by definition of G4ExtrudedSolid.
      rot->rotateX(-90.*deg);
    }

    G4VPhysicalVolume* DF_substrate_physical = dynamic_cast<G4VPhysicalVolume*>(
      new G4PVPlacement(
        rot,
        G4ThreeVector(
          0.,
          plate_thickn_/2.
          +wlsp_DF_gap_
          +MLS_thickn_
          +DF_substrate_thickn_/2.,
          0.
        ),
        "DICHROIC_FILTER_SUBSTRATE",
        DF_substrate_logic,
        mother_physical,
        false,
        0,
        true
      )
    );

    // DF MLS
    G4VSolid* MLS_solid = nullptr;
    if(shape_code_==0)
    {
      MLS_solid = dynamic_cast<G4VSolid*>(
        new G4Box(
          "MLS",
          plate_length_/2.,
          MLS_thickn_/2.,
          plate_width_/2.
        )
      );
    }
    else  // It has been checked before (within this same method) that shape_code_ is either 0 or 1
    {
      // Extrude the triangular
      MLS_solid = dynamic_cast<G4VSolid*>(
        new G4ExtrudedSolid(
          "MLS",
          this->GetTriangularPlatePrismBase(),
          MLS_thickn_/2.
        )
      );
    }

    G4LogicalVolume* MLS_logic = new G4LogicalVolume(
      MLS_solid, 
      DF_substrate_mat,   // Yes, in the detached DF model, the
      "MLS"               // MLS MPT is that of the DF substrate
    );
    
    G4VisAttributes MLS_col = nexus::BloodRedAlpha();
    //MLS_col.SetForceSolid(true);
    MLS_logic->SetVisAttributes(MLS_col);

    G4VPhysicalVolume* MLS_physical = dynamic_cast<G4VPhysicalVolume*>(
      new G4PVPlacement(
        rot,
        G4ThreeVector(
          0.,
          plate_thickn_/2.
          +wlsp_DF_gap_
          +MLS_thickn_/2.,
          0.
        ), 
        "MLS",
        MLS_logic,
        mother_physical,
        false,
        0,
        true
      )
    );

    // Check that there's dichroic information for ingoing (wrt APEX) photons
    if(path_to_inwards_dichroic_data_==""){
        G4Exception("[APEX]", "ConstructDetachedDichroicFilter()",
                    FatalException, "The path to the inwards dichroic data file was not set.");
    }

    // Check that there's dichroic information for outgoing photons
    if(path_to_outwards_dichroic_data_==""){
        G4Exception("[APEX]", "ConstructDetachedDichroicFilter()",
                    FatalException, "The path to the outwards dichroic data file was not set.");
    }

    // Construct the ingoing optical surface
    setenv("G4DICHROICDATA", path_to_inwards_dichroic_data_, 1);
    G4OpticalSurface* df_inwards_opsurf =
        new G4OpticalSurface(             // G4OpticalSurface constructor loads the
          "DICHROIC_INWARDS_OPSURF",      // dichroic information from the file which
          dichroic,                       // is currently pointed to by the environment
          polished,                       // variable G4DICHROICDATA
          dielectric_dichroic
        );

    // Construct the outgoung optical surface
    setenv("G4DICHROICDATA", path_to_outwards_dichroic_data_, 1);   // Note that, if you did not compile the modified version of G4 code, 
                                                                    // then different G4 dichroic data cannot be loaded. Instead, the first 
                                                                    // one (i.e. the one I am setting from path_to_inwards_dichroic_data_), 
                                                                    // is the one that will apply for every dichroic boundary in the simulation.
    G4OpticalSurface* df_outwards_opsurf =   
        new G4OpticalSurface(
          "DICHROIC_OUTWARDS_OPSURF",
          dichroic,
          polished,
          dielectric_dichroic
        );

    // Endow the DF_substrate_physical->MLS_physical surface with the ingoing optical surface
    new G4LogicalBorderSurface(
      "DF SUBSTRATE->DF MLS", 
      DF_substrate_physical, 
      MLS_physical, 
      df_inwards_opsurf
    );

    // Endow the MLS_physical->DF_substrate_physical surface with the outgoing optical surface
    new G4LogicalBorderSurface(
      "DF MLS->DF SUBSTRATE",
      MLS_physical,
      DF_substrate_physical,
      df_outwards_opsurf
    );

    // pTP coating
    if(!remove_coating_)
    {
        G4VSolid* coating_solid = nullptr;
        if(shape_code_==0)
        {
          coating_solid = dynamic_cast<G4VSolid*>(
            new G4Box(
              "COATING",
              plate_length_/2.,
              coating_thickn_/2.,
              plate_width_/2.
            )
          );
        }
        else  // It has been checked before (within this same method) that shape_code_ is either 0 or 1
        {
          // Extrude the triangular
          coating_solid = dynamic_cast<G4VSolid*>(
            new G4ExtrudedSolid(
              "MLS",
              this->GetTriangularPlatePrismBase(),
              coating_thickn_/2.
            )
          );
        }

        G4Material* coating_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
        coating_mat->SetMaterialPropertiesTable(opticalprops::PTP(coating_rindex_));
        G4LogicalVolume* coating_logic = new G4LogicalVolume(
          coating_solid,
          coating_mat,
          "COATING"
        );

        G4VisAttributes coating_col = nexus::TitaniumGreyAlpha();
        coating_col.SetForceSolid(true);
        coating_logic->SetVisAttributes(coating_col);

        // Place the coating
        G4VPhysicalVolume* coating_physical = dynamic_cast<G4VPhysicalVolume*>(
          new G4PVPlacement(
            rot,
            G4ThreeVector(
              0.,
              plate_thickn_/2.
              +wlsp_DF_gap_
              +MLS_thickn_
              +DF_substrate_thickn_
              +coating_thickn_/2.,
              0.
            ), 
            "COATING",
            coating_logic,
            mother_physical,
            false,
            0,
            true
          )
        );

        // Make the LAR-coating interface rough, so that photons cannot be trapped within the coating
        G4OpticalSurface* coating_rough_surf =
                new G4OpticalSurface(
                  "COATING_ROUGH_SURFACE",
                  glisur,
                  ground,
                  dielectric_dielectric,
                  .01
                );
                // 0.01 is the polish value for glisur model that was measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
                // This is the best reference we have, since both PTP and TPB are the result of an evaporation+deposition process
        new G4LogicalBorderSurface(
          "SURROUNDINGS->COATING",
          mother_physical,
          coating_physical,
          coating_rough_surf
        );
        new G4LogicalBorderSurface(
          "COATING->SURROUNDINGS", 
          coating_physical,
          mother_physical,
          coating_rough_surf
        );
        // We will also add roughness for the coating->DF substrate, but only with such ordering. The alternative case takes place
        // when the photon travels from the DF substrate to the coating. The DF substrate is supposed to be polished, so the photon may not see a rough 
        // surface.
        new G4LogicalBorderSurface(
          "COATING->DF SUBSTRATE",
          coating_physical,
          DF_substrate_physical,
          coating_rough_surf
        );
    }
    return;
  }

  // This function is no longer supported. It is kept here for future reference.
  void APEX::ConstructBoard(G4VPhysicalVolume* mother_physical) const ///< Deprecated
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

    // Board in the middle of a large face
    if(board_position_code_==1)
    {
      rot->rotateX(-90.0*deg);
      pos.set(
        0.,
        -1.*(plate_thickn_/2.)-1.*(board.GetOverallThickness()/2.)-gap_,  // Note that what's placed in the global origin of 
                                                                          // coordinates is the plate, not the reflective foil. 
        0.
      );
    }
    else
    {
      rot->rotateY(+180.0*deg);
      pos.set(
        0.,
        0.,
        -1.*(plate_length_/2.)-1.*(board.GetOverallThickness()/2.)
      );
    }
  
    new G4PVPlacement(
      rot,
      pos,
      "SIPMS_BOARD",
      board_logic_vol,
      mother_physical,
      false,
      0,
      false
    );
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

    if(detach_DF_)
    {
      y_pos += wlsp_DF_gap_ + DF_substrate_thickn_;
    }

    if(generation_region_=="custom"){
      G4double random_radius =  UniformRandomInRange(gen_diameter_/2., 0.);
      G4double random_angle =   UniformRandomInRange(twopi, 0.); 
      x_pos = gen_x_ +(random_radius*sin(random_angle));
      z_pos = gen_z_ +(random_radius*cos(random_angle));
    }
    else{ // Default behaviour is that of generation_region_=="random"

      if(shape_code_==0)
      {
        x_pos = UniformRandomInRange(
          plate_length_/2.,
          -1.*plate_length_/2.
        );
        z_pos = UniformRandomInRange(
          plate_width_/2.,
          -1.*plate_width_/2.
        );
      }
      else if(shape_code_==1)
      {
        // The parameterization in the body of this conditional block takes into
        // account that the triangle point which is placed at the coordinates
        // system origin is the point which
        //
        //  1) belongs to the symmetry axis of the isosceles triangle, and
        //  2) is placed at half the triangle height (plate_width_) from the base of the triangle
        //
        // Note that this parameterization is different to the one used in the
        // WLSPlate::GenerateVertex() method (at the time of writing).
        //
        // The situation that we got now is something like this:
        //
        //                            +z
        //                            /\
        //                            |
        //                            |
        //                            x   
        //                   z_1(x) / | \ z_2(x)      
        //                            |         
        //     -x <--------------/----+----\------------------> +x
        //                            |         
        //                   /________|________\
        //                            |
        //                            |
        //                            |
        //                            v
        //                            -z
        //
        // Imagine that we parameterize the left (resp. right) side of the triangle with
        // the function z_1(x) (resp. z_2(x)). This dependency can be inverted to find
        // x_1(z) (resp. x_2(z)). Then, after having generated a random z_pos, in the
        // (-plate_width_/2., plate_width_/2.) range, we can calculate x_1(z_pos), which
        // is basically the (negative) half-width of the isosceles triangle at a height
        // given by the sampled z_pos. Now, to not bias the generation of the vertex
        // towards the upper peak of the triangle, we need to accept the generated z_pos
        // with a probability proportional to the ratio of the computed margin to half
        // of the triangle base. Since the margin we computed is negative, we just invert
        // its sign and compute its ratio to (plate_length_/2.). Note that the rejection
        // probability approaches 0, when the generated z_pos landed very close to the
        // triangle base, and it is almost 1 when the generated z_pos is very close to
        // the triangle peak. This is as it should be, since the area near the base of
        // the triangle is much bigger than the area near the peak, and a random generation
        // should give an uniform superficial density of photon hits throughout the
        // triangle surface.

        G4double negative_margin;
        G4bool accepted_z_pos = false;

        while(!accepted_z_pos)
        {
          z_pos = UniformRandomInRange(
            (plate_width_/2.)-tolerance,
            (-1.*plate_width_/2.)+tolerance
          );

          negative_margin = ((z_pos*plate_length_)/(2.*plate_width_))-(plate_length_/4);

          if(UniformRandomInRange(1., 0.) < (-1.*negative_margin)/(plate_length_/2.))
          {
            accepted_z_pos = true;
          }
        }

        // Under the assumption that the triangle height (plate_width_) is
        // bigger than the tolerance (t), i.e. plate_width_>t, and that the
        // tolerance is smaller than 1.0, then you can prove that, for
        //
        //  ((z_pos*plate_length_)/(2.*plate_width_))-(plate_length_/4)+(k*tolerance) < 0     (1)
        //
        // to hold in the whole range of
        //
        //    z_pos \in [(-plate_width_/2.)+tolerance, (plate_width_/2.)-tolerance],
        //
        // (which is the range of random generation of z_pos), it is needed
        // that k < plate_length_/(2*plate_width_). The reason why we need
        // the inequality (1) above is that UniformRandomInRange(x, y) works
        // for x>y (otherwise we are inverting the range, which makes no sense).
        // That's why we are introducing the factor k in the range limits of
        // the x_pos random generation.

        G4double k = 0.5 * plate_length_/(2.*plate_width_); // Smaller than plate_length_/(2*plate_width_)
        negative_margin += k*tolerance;

        x_pos = UniformRandomInRange(
          -1.*negative_margin,
          negative_margin
        );
      }
      else
      {
        G4Exception("[APEX]", "GenerateVertex()",
                    FatalException, "The given shape code is not recognized.");
      }
    }
    return G4ThreeVector(x_pos, y_pos, z_pos);
  }

  std::vector<G4TwoVector> APEX::GetTriangularPlatePrismBase() const{

    std::vector<G4TwoVector> prism_base;
    // plate_length_ (resp. plate_width_) is the base (resp. height) of
    // the triangle. The 2D polygon which we are creating in the XY plane
    // is the following:
    //
    //                            +y
    //                            /\
    //                            |
    //                            |B
    //                            x
    //                          / | \
    //                            |
    //     -x <--------------/----+----\------------------> +x
    //                            |
    //                   /________|________\
    //                 A          |          C
    //                            |
    //                            |
    //                            v
    //                            -y
    //
    // where the horizontal line (parallel to the x axis) which goes from
    // point A to point C is, the base of the isosceles triangle. It is
    // important to note two things:
    //
    //    1)  The triangle is centered about the Y-axis, while point A (or
    //        C) is placed at a vertical distance of plate_width_/2 from
    //        the origin of coordinates.
    //
    //    2)  After extrusion in the Z direction, we will need to rotate
    //        the triangle 90º degrees about the X-axis, so that the
    //        thickness dimension of the plate is set along the Y-axis,
    //        as for the rectangular plate case.

    // Point A of the sketch above
    prism_base.push_back(G4TwoVector(-1.*plate_length_/2., -1.*plate_width_/2.));
    // Point B of the sketch above
    prism_base.push_back(G4TwoVector(0., plate_width_/2.));
    // Point C of the sketch above
    prism_base.push_back(G4TwoVector(plate_length_/2., -1.*plate_width_/2.));

    return prism_base;
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
