#include "XArapucaDS.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "WLSPlate.h"
#include "HamamatsuS133606050VE.h"
#include "HamamatsuS133605075HQR.h"
#include "FbkNuvHdCryoTT.h"
#include "PerfectSiPMMPPC.h"
#include "ScalableHamamatsuS133606050VE.h"
#include "SiPMBoard.h"
#include "RandomUtils.h"
#include "Visibilities.h"

#include <algorithm>
#include <random>
#include <cmath>
#include <G4GenericMessenger.hh>
#include <G4UserLimits.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Orb.hh>
#include <G4Transform3D.hh>
#include <G4MultiUnion.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
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

REGISTER_CLASS(XArapucaDS, GeometryBase)

namespace nexus{

  XArapucaDS::XArapucaDS():
  GeometryBase(), 
  ///Get internal reflector cavity dimensions from arxiv.org/abs/1912.09191
  config_code_                          (1                            ),
  surrounding_media_                    ("lar"                        ),
  internal_length_                      (488.   *mm                   ),   ///X
  internal_width_                       (100.   *mm                   ),   ///Z
  internal_thickn_                      (8      *mm                   ),   ///Y
  CSA_thickn_                           (1.     *mm                   ),
  CS_thickn_                            (1.000  *mm                   ),
  CS_mpt_                               (opticalprops::FusedSilica()  ),
  //CS_mpt_                     (opticalprops::SCHOTT_B270()  ),
  CS_pos_wrt_CSA_pos_                   (0.                           ),
  substrates_are_coated_                (true                         ),
  coating_thickn_                       (3.226  *um                   ),  // Based on arxiv.org/abs/1912.09191 and TDR vol.IX, section 5.8.3.1,
                                                                          // the pTP film thickness is such that there's 400 micrograms of pTP
                                                                          // deposited over each square centimeter of DF.
                                                                          // This, together with the pTP density, (1.24g/cm3, found in
                                                                          // en.wikipedia.org/wiki/Terphenyl), gives a pTP film thickness of 
                                                                          // 3.226 micrometers
  outter_frame_width_along_wlsplength_  (10.    *mm                   ),
  outter_frame_width_along_wlspwidth_   (10.    *mm                   ),
  inner_frames_width_along_wlsplength_  (6.     *mm                   ),
  inner_frames_width_along_wlspwidth_   (2.     *mm                   ),
  cs_no_along_wlsplength_               (2                            ),
  cs_no_along_wlspwidth_                (3                            ),
  CSA_frame_is_reflective_              (false                        ),
  CSA_frame_is_specular_                (true                         ),
  remove_CSs_                           (false                        ),  
  remove_CSA_frame_                     (false                        ),
  SS_cromophore_concentration_          (40.                          ),
  TS_cromophore_concentration_          (40.                          ),
  case_thickn_                          (1.     *mm                   ),   ///Get foil thickness from isoltronic.ch/assets/of-m-vikuiti-esr-app-guide.pdf
  SiPM_code_                            (1                            ),
  num_phsensors_                        (24                           ),
  gap_                                  (0.5    *mm                   ),
  ref_phsensors_supports_               (true                         ), 
  double_sided_                         (true                         ),
  collectors_are_reflective_            (false                        ),
  generation_region_                    ("random"                     ),
  gen_x_                                (0.     *cm                   ),
  gen_z_                                (0.     *cm                   ),
  gen_diameter_                         (1.*cm                        ),
  world_extra_thickn_                   (100.   *cm                   ),
  plate_length_                         (487.   *mm                   ),   ///X
  plate_width_                          (93.    *mm                   ),   ///Z
  plate_thickn_                         (3.5    *mm                   ),   ///Y
  ring_thickness_                       (20     *mm                   ),
  plate_ring_gap_                       (1      *mm                   ),
  only_sipms_along_long_sides_          (true                         ),
  with_boards_                          (true                         ),
  PS_config_code_                       (1                            ),
  with_dimples_                         (true                         ),
  dimple_type_                          ("cylindrical"                ),
  flat_dimple_width_                    (6.1    *mm                   ),
  flat_dimple_depth_                    (2.     *mm                   ),
  curvy_dimple_radius_                  (1.5    *mm                   )
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/XArapucaDS/",
				"Control commands of geometry XArapucaDS.");

    G4GenericMessenger::Command& cc_cmd =
      msg_->DeclareProperty("config_code", config_code_,
			    "Configuration code.");
    cc_cmd.SetParameterName("config_code", false);
    cc_cmd.SetRange("config_code>=1"); 

    G4GenericMessenger::Command& sm_cmd =
      msg_->DeclareProperty("surrounding_media", surrounding_media_,
			    "Which media to place the XArapuca DS in");
    
    G4GenericMessenger::Command& il_cmd =
      msg_->DeclareProperty("internal_length", internal_length_,
			    "Internal length of the X-Arapuca DS cavity.");
    il_cmd.SetUnitCategory("Length");
    il_cmd.SetParameterName("internal_length", false);
    il_cmd.SetRange("internal_length>0.");

    G4GenericMessenger::Command& iw_cmd =
      msg_->DeclareProperty("internal_width", internal_width_,
			    "Internal width of the X-Arapuca DS cavity.");
    iw_cmd.SetUnitCategory("Length");
    iw_cmd.SetParameterName("internal_width", false);
    iw_cmd.SetRange("internal_width>0.");

    G4GenericMessenger::Command& it_cmd =
      msg_->DeclareProperty("internal_thickn", internal_thickn_,
			    "Internal thickness of the X-Arapuca DS cavity.");
    it_cmd.SetUnitCategory("Length");
    it_cmd.SetParameterName("internal_thickness", false);
    it_cmd.SetRange("internal_thickness>0.");

    G4GenericMessenger::Command& csat_cmd =
      msg_->DeclareProperty("CSA_thickn", CSA_thickn_,
			    "Stands for Coated Substrates Assembly. It is the frame thickness.");
    csat_cmd.SetUnitCategory("Length");
    csat_cmd.SetParameterName("CSA_thickn", false);
    csat_cmd.SetRange("CSA_thickn>0.");

    G4GenericMessenger::Command& cst_cmd =
      msg_->DeclareProperty("CS_thickn", CS_thickn_,
			    "Thickness of the coated substrates.");
    cst_cmd.SetUnitCategory("Length");
    cst_cmd.SetParameterName("CS_thickn", false);
    cst_cmd.SetRange("CS_thickn>0.");

    G4GenericMessenger::Command& cpwcp_cmd =
      msg_->DeclareProperty("CS_pos_wrt_CSA_pos", CS_pos_wrt_CSA_pos_,
			    "Position (height) of the coated substrates with respect to the CSA position (height). This parameter can take values from 0 to 1.");
    cpwcp_cmd.SetParameterName("CS_pos_wrt_CSA_pos", false);
    cpwcp_cmd.SetRange("CS_pos_wrt_CSA_pos>=0.");
    cpwcp_cmd.SetRange("CS_pos_wrt_CSA_pos<=1.");

    G4GenericMessenger::Command& sac_cmd =
      msg_->DeclareProperty("substrates_are_coated", substrates_are_coated_,
			    "Whether the substrates are coated with PTP or not.");

    G4GenericMessenger::Command& cot_cmd =
      msg_->DeclareProperty("coating_thickn", coating_thickn_,
			    "Thickness of the coating layer that is deposited over the substrates.");
    cot_cmd.SetUnitCategory("Length");
    cot_cmd.SetParameterName("coating_thickn", false);
    cot_cmd.SetRange("coating_thickn>0.");

    G4GenericMessenger::Command& ofwawl_cmd =
      msg_->DeclareProperty("outter_frame_width_along_wlsplength", outter_frame_width_along_wlsplength_,
			    "Outter frame dimension along the length of the X-ARAPUCA DS.");
    ofwawl_cmd.SetUnitCategory("Length");
    ofwawl_cmd.SetParameterName("outter_frame_width_along_wlsplength", false);
    ofwawl_cmd.SetRange("outter_frame_width_along_wlsplength>0.");

    G4GenericMessenger::Command& ofwaww_cmd =
      msg_->DeclareProperty("outter_frame_width_along_wlspwidth", outter_frame_width_along_wlspwidth_,
			    "Outter frame dimension along the width of the X-ARAPUCA DS.");
    ofwaww_cmd.SetUnitCategory("Length");
    ofwaww_cmd.SetParameterName("outter_frame_width_along_wlspwidth", false);
    ofwaww_cmd.SetRange("outter_frame_width_along_wlspwidth>0.");

    G4GenericMessenger::Command& ifwawl_cmd =
      msg_->DeclareProperty("inner_frames_width_along_wlsplength", inner_frames_width_along_wlsplength_,
			    "Inner frames dimension along the length of the X-ARAPUCA DS.");
    ifwawl_cmd.SetUnitCategory("Length");
    ifwawl_cmd.SetParameterName("inner_frames_width_along_wlsplength", false);
    ifwawl_cmd.SetRange("inner_frames_width_along_wlsplength>0.");

    G4GenericMessenger::Command& ifwaww_cmd =
      msg_->DeclareProperty("inner_frames_width_along_wlspwidth", inner_frames_width_along_wlspwidth_,
			    "Inner frames dimension along the width of the X-ARAPUCA DS.");
    ifwaww_cmd.SetUnitCategory("Length");
    ifwaww_cmd.SetParameterName("inner_frames_width_along_wlspwidth", false);
    ifwaww_cmd.SetRange("inner_frames_width_along_wlspwidth>0.");

    G4GenericMessenger::Command& csnawl_cmd =
      msg_->DeclareProperty("cs_no_along_wlsplength", cs_no_along_wlsplength_,
			    "Number of coated substrates along the length of the X-ARAPUCA DS.");
    csnawl_cmd.SetParameterName("cs_no_along_wlsplength", false);
    csnawl_cmd.SetRange("cs_no_along_wlsplength>=0"); 

    G4GenericMessenger::Command& csnaww_cmd =
      msg_->DeclareProperty("cs_no_along_wlspwidth", cs_no_along_wlspwidth_,
			    "Number of coated substrates along the width of the X-ARAPUCA DS.");
    csnaww_cmd.SetParameterName("cs_no_along_wlspwidth", false);
    csnaww_cmd.SetRange("cs_no_along_wlspwidth>=0"); 

    G4GenericMessenger::Command& csafir_cmd =
      msg_->DeclareProperty("CSA_frame_is_reflective", CSA_frame_is_reflective_,
			    "Whether the FR4 CSA frame is vikuiti-coated or not.");

    G4GenericMessenger::Command& csafis_cmd =
      msg_->DeclareProperty("CSA_frame_is_specular", CSA_frame_is_specular_,
			    "Whether the vikuiti coating of the CSA frame is specular-spikely reflective or diffusively reflective. Only makes a difference if CSA_frame_is_reflective_==True.");

    G4GenericMessenger::Command& rcss_cmd =
      msg_->DeclareProperty("remove_CSs", remove_CSs_,
			    "Whether to remove the coated substrates or not.");

    G4GenericMessenger::Command& rcsaf_cmd =
      msg_->DeclareProperty("remove_CSA_frame", remove_CSA_frame_,
			    "Whether to remove the CSA frame.");

    G4GenericMessenger::Command& sscrco_cmd =
      msg_->DeclareProperty("SS_cromophore_concentration", SS_cromophore_concentration_,
			    "Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter, in case G2P_FB118 is used.");
    sscrco_cmd.SetParameterName("SS_cromophore_concentration", false);
    sscrco_cmd.SetRange("SS_cromophore_concentration>0.");

    G4GenericMessenger::Command& tscrco_cmd =
      msg_->DeclareProperty("TS_cromophore_concentration", TS_cromophore_concentration_,
			    "Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the tertiary WLShifter, in case FakeG2P_FB118 is used.");
    tscrco_cmd.SetParameterName("TS_cromophore_concentration", false);
    tscrco_cmd.SetRange("TS_cromophore_concentration>0.");

    G4GenericMessenger::Command& ct_cmd =
      msg_->DeclareProperty("case_thickn", case_thickn_,
			    "Reflective foils thickness.");
    ct_cmd.SetUnitCategory("Length");
    ct_cmd.SetParameterName("case_thickn", false);
    ct_cmd.SetRange("case_thickn>0.");

    G4GenericMessenger::Command& sc_cmd =
      msg_->DeclareProperty("SiPM_code", SiPM_code_,
			    "Integer signalling which SiPM to construct.");
    sc_cmd.SetParameterName("SiPM_code", false);
    sc_cmd.SetRange("SiPM_code>=1"); 

    G4GenericMessenger::Command& np_cmd =
      msg_->DeclareProperty("num_phsensors", num_phsensors_,
			    "If config_code_==1, this is the number of SiPMs per long side.");
    np_cmd.SetParameterName("num_phsensors", false);
    np_cmd.SetRange("num_phsensors>=0"); 

    G4GenericMessenger::Command& g_cmd =
      msg_->DeclareProperty("gap", gap_,
			    "Gap between photosensors and WLS ring.");
    g_cmd.SetUnitCategory("Length");
    //g_cmd.SetParameterName("gap", false);
    //g_cmd.SetRange("gap>0.");    // These are commented so that gap_ can help modelate the immersion of the SiPMs into the flat dimple

    G4GenericMessenger::Command& rps_cmd =
      msg_->DeclareProperty("ref_phsensors_supports", ref_phsensors_supports_,
			    "Whether the photosensors supports are VIKUITI-coated.");

    G4GenericMessenger::Command& ds_cmd =
      msg_->DeclareProperty("double_sided", double_sided_,
			    "Whether there are coated substrates on both sides of the WLS plate.");

    G4GenericMessenger::Command& cr_cmd =
      msg_->DeclareProperty("collectors_are_reflective", collectors_are_reflective_,
			    "Whether the collectors that replace the coated substrates are reflective or not.");

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
			    "Diameter of the circle where the GV could be randomly sampled. It is only used if generation_region_=='custom' is True.");
    gd_cmd.SetUnitCategory("Length");
    gd_cmd.SetParameterName("gen_diameter", false);
    gd_cmd.SetRange("gen_diameter>0.");

    G4GenericMessenger::Command& pl_cmd =
      msg_->DeclareProperty("plate_length", plate_length_,
			    "Length of the secondary WLShifter plate dimensions.");
    pl_cmd.SetUnitCategory("Length");
    pl_cmd.SetParameterName("plate_length", false);
    pl_cmd.SetRange("plate_length>0.");

    G4GenericMessenger::Command& pw_cmd =
      msg_->DeclareProperty("plate_width", plate_width_,
			    "Width of the secondary WLShifter plate dimensions.");
    pw_cmd.SetUnitCategory("Length");
    pw_cmd.SetParameterName("plate_width", false);
    pw_cmd.SetRange("plate_width>0.");

    G4GenericMessenger::Command& pt_cmd =
      msg_->DeclareProperty("plate_thickn", plate_thickn_,
			    "Thickness of the secondary WLShifter plate dimensions.");
    pt_cmd.SetUnitCategory("Length");
    pt_cmd.SetParameterName("plate_thickn", false);
    pt_cmd.SetRange("plate_thickn>0.");

    G4GenericMessenger::Command& rt_cmd =
      msg_->DeclareProperty("ring_thickness", ring_thickness_,
			    "Tertiary WLShifter ring thickness.");
    rt_cmd.SetUnitCategory("Length");
    rt_cmd.SetParameterName("ring_thickness", false);
    rt_cmd.SetRange("ring_thickness>0.");

    G4GenericMessenger::Command& prg_cmd =
      msg_->DeclareProperty("plate_ring_gap", plate_ring_gap_,
			    "Gap in between the WLS plate and the WLS ring.");
    prg_cmd.SetUnitCategory("Length");
    prg_cmd.SetParameterName("plate_ring_gap", false);
    prg_cmd.SetRange("plate_ring_gap>0.");

    G4GenericMessenger::Command& osals_cmd =
      msg_->DeclareProperty("only_sipms_along_long_sides", only_sipms_along_long_sides_,
			    "Whether the photo sensors are installed only along the long sides of the X-ARAPUCA DS (i.e. the 'length' dimension), or along every side.");

    G4GenericMessenger::Command& wb_cmd =
      msg_->DeclareProperty("with_boards", with_boards_,
			    "Whether the photosensors are mounted on a board.");

    G4GenericMessenger::Command& pscc_cmd =
      msg_->DeclareProperty("PS_config_code", PS_config_code_,
			    "Photo sensors configuration code.");
    pscc_cmd.SetParameterName("PS_config_code", false);
    pscc_cmd.SetRange("PS_config_code>=1"); 

    G4GenericMessenger::Command& wd_cmd =
      msg_->DeclareProperty("with_dimples", with_dimples_,
			    "Whether the ring has carved dimples on it.");

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

    // When testing WLS plates with opticalprops::noAbsLength_, it is possible that a photon
    // gets trapped within the plate (below the critical angle) into an infinite-bouncing-loop.
    // For the case of an EJ286 plate with DUNE supercells dimensions, with this absorption 
    // length and immersed into LAr, some analysis revealed that particles that did not fall 
    // into this infinite loop had track lengths of less than one hundred meters.
    ul_ = new G4UserLimits();
    ul_->SetUserMaxTrackLength(100*m);
  }

  XArapucaDS::~XArapucaDS()
  {
      if(ul_) delete ul_;
  }

  void XArapucaDS::Construct()
  {

    // Compute internal attributes
    overall_length_ = internal_length_  + 2*case_thickn_;
    overall_thickn_ = internal_thickn_  + 2*CSA_thickn_;    // geometry_is_ill_formed() will make sure that CS_thickn_<=CSA_thickn_ 
    overall_width_  = internal_width_   + 2*case_thickn_;
    CSA_length_ = overall_length_;                          // For now, these two match. 
    CSA_width_  = overall_width_;                           // For now, these two match.  
    CS_length_  = (CSA_length_  -(2.*outter_frame_width_along_wlsplength_)  -((cs_no_along_wlsplength_-1.)*inner_frames_width_along_wlsplength_))   /(1.*cs_no_along_wlsplength_);   
    CS_width_   = (CSA_width_   -(2.*outter_frame_width_along_wlspwidth_)   -((cs_no_along_wlspwidth_-1.)*inner_frames_width_along_wlspwidth_))     /(1.*cs_no_along_wlspwidth_);

    G4cout << "CS_length_=" << CS_length_ << G4endl;
    G4cout << "CS_width_=" << CS_width_ << G4endl;

    if(geometry_is_ill_formed()){
      G4Exception("[XArapucaDS]", "Construct()", FatalException,
      "The given dimensions do not describe a feasible X-Arapuca.");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // The biggest volume is a vacuum box with the dimensions of the X-ARAPUCA DS device
    // plus 2*world_extra_thickn_, for each dimension. This volume is the effective 
    // world volume FOR THE NEXUS USER, but it is not the world volume of the overall
    // Geant4 application (afterwards nexus takes the user's biggest volume 
    // and places it inside another world volume whose dimensions are enough so as to 
    // fit the whole span of the biggest volume you implemented). Placed within this
    // vacuum box, there's a box, made out of surrounding_media_, with dimensions of 
    // the X-ARAPUCA DS device plus world_extra_thickn_, along each dimension. The previous 
    // version of this class simply implemented a surrounding_media_ box as the biggest 
    // volume, without the need to encapsulate it within a bigger volume (apart from 
    // that of nexus). However, as of the implementation of the dichroic filter, the 
    // use of G4LogicalBorderSurface class is necessary. Objects of such class must 
    // be given two ordered physical volumes, so that a particle flying from the first 
    // one to the second one, interacts with such an optical surface. Since the nexus 
    // user seems not to have access to the physical placement of the biggest volume 
    // (in our case, the vacuum box) which is implemented in line 80 of 
    // source/base/DetectorConstruction.cc, AFTER calling GeometryBase::Construct() 
    // (i.e. your geometry must be constructed before placing your biggest volume in 
    // nexus world volume,so you cannot possibly have access to the physical placement 
    // of your biggest volume when you Construct() it.), one way around this issue is 
    // adding this intermediate surrounding_media_ box.   
    /////////////////////////////////////////////////////////////////////////////////////
    
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
                (internal_length_ +case_thickn_ +world_extra_thickn_)/2.,
                (internal_thickn_ +case_thickn_ +world_extra_thickn_)/2.,
                (internal_width_ +case_thickn_ +world_extra_thickn_)/2.);

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

    // ConstructReflectiveCase(mother_physical);
    // //ConstructCollectors(mother_physical); 
    // if(!remove_CSs_ || !remove_CSA_frame_) ConstructSubstratesAssemblies(mother_physical);
    if(config_code_==1){
        ConstructWLSPlate(mother_physical);
        ConstructWLSRing(mother_physical);
        if(with_boards_) ConstructBoards(mother_physical);
        else ConstructPhotosensors(mother_physical);
    }
    else{
      G4Exception("[XArapucaDS]", "XArapucaDS()",
                  FatalException, "Not allowed configuration code.");
    }

    return;
  }

  void XArapucaDS::ConstructWLSPlate(G4VPhysicalVolume* mother_physical) const
  { 

    WLSPlate* plate = new WLSPlate( plate_length_, 
                                    plate_thickn_, 
                                    plate_width_, 
                                    opticalprops::G2P_FB118(SS_cromophore_concentration_, 1.502, true), 
                                    false,
                                    only_sipms_along_long_sides_ ? false : true,
                                    only_sipms_along_long_sides_ ? false : true,
                                    true,
                                    true,
                                    "", 
                                    0,          // For XArapuca DS, dimples are not carved
                                    0.,         // on the secondary WLSplate. If there are
                                    0.,         // dimples, they are carved on the WLSring.
                                    0.);    

    plate->Construct();
    G4LogicalVolume* plate_logic = plate->GetLogicalVolume();
    plate_logic->SetUserLimits(ul_);

    /*
    G4VisAttributes wlsp_col = nexus::BlueAlpha();
    wlsp_col.SetForceSolid(true);
    plate_logic->SetVisAttributes(wlsp_col);
    */

    if (!plate_logic) {
      G4Exception("[XArapucaDS]", "ConstructWLSPlate()",
                  FatalException, "Null pointer to logical volume.");
    }

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), plate_logic->GetName(), 
                    plate_logic, mother_physical, false, 0, true);
    
    return;
  }  

  void XArapucaDS::ConstructWLSRing(G4VPhysicalVolume* mother_physical) const
  { 

    const G4String ring_name = "WLS_RING";

    G4Box* aux_outter_box = new G4Box("AUX_OUTTER_BOX", 
                                      plate_length_/2. +plate_ring_gap_+ring_thickness_, 
                                      plate_thickn_/2., 
                                      plate_width_/2. +plate_ring_gap_+ring_thickness_);  
    G4double tolerance = 1.*cm;
    G4Box* aux_inner_box =  new G4Box("AUX_INNER_BOX", 
                                      plate_length_/2. +plate_ring_gap_, 
                                      (plate_thickn_/2. +tolerance)/2., 
                                      plate_width_/2. +plate_ring_gap_);

    G4SubtractionSolid* ring_solid =    new G4SubtractionSolid( ring_name, aux_outter_box, aux_inner_box, 
                                                                nullptr, G4ThreeVector(0., 0., 0.));

    G4VSolid* geometry_solid = nullptr; 
    G4int how_many_dimples = num_phsensors_;  
    if(with_dimples_ && how_many_dimples>=1)
    {
        G4double tolerance = 0.5*mm; // To avoid boolean subtraction of matching surfaces
        G4VSolid* carving_solid = nullptr;
        if(dimple_type_=="flat"){
            // Yes, using 2*flat_dimple_depth_ as the whole dimension of the carving along the z-axis is ok 
            // (the carvings are subtracted from the very edge of the plate)
            carving_solid = dynamic_cast<G4VSolid*>(new G4Box("AUX", flat_dimple_width_/2., plate_thickn_/2. +tolerance, flat_dimple_depth_));
        }
        else if(dimple_type_=="spherical"){
            carving_solid = dynamic_cast<G4VSolid*>(new G4Orb("AUX", curvy_dimple_radius_));
        }
        else{ //Default is cylindrical dimples
            carving_solid = dynamic_cast<G4VSolid*>(
                                        new G4Tubs("AUX", 0., curvy_dimple_radius_, plate_thickn_/2. +tolerance, 0., twopi));
        }

        G4RotationMatrix* rot = new G4RotationMatrix();
        if(dimple_type_=="cylindrical"){
            rot->rotateX(+90.*deg);
        }

        G4Transform3D* transform_ptr = nullptr;
        G4MultiUnion* carvings_multiunion_solid = new G4MultiUnion("CARVINGS");
        G4double aux_x = plate_length_  + (2.*plate_ring_gap_) +(2.*ring_thickness_);
        G4double aux_z = plate_width_   + (2.*plate_ring_gap_) +(2.*ring_thickness_); 
        for(G4int i=0; i<how_many_dimples; i++){
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector((-1.*aux_x/2.)+(1.*(0.5+i)*aux_x/(1.*how_many_dimples)), 0., +aux_z/2.));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector((-1.*aux_x/2.)+(1.*(0.5+i)*aux_x/(1.*how_many_dimples)), 0., -1.*aux_z/2.));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
        }

        if(dimple_type_=="flat"){
            rot->rotateY(+90.*deg);
        }

        if(!only_sipms_along_long_sides_){
            for(G4int i=0; i<how_many_dimples; i++){
                transform_ptr = new G4Transform3D(*rot, G4ThreeVector(+aux_x/2., 0., (-1.*aux_z/2.)+(1.*(0.5+i)*aux_z/(1.*how_many_dimples))));
                carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
                transform_ptr = new G4Transform3D(*rot, G4ThreeVector(-1.*aux_x/2., 0., (-1.*aux_z/2.)+(1.*(0.5+i)*aux_z/(1.*how_many_dimples))));
                carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
            }
        }
        carvings_multiunion_solid->Voxelize();

        G4SubtractionSolid* dimpled_ring_solid = new G4SubtractionSolid(ring_name, 
                                                                        ring_solid, carvings_multiunion_solid);
        geometry_solid = dynamic_cast<G4VSolid*>(dimpled_ring_solid);
    }
    else{
        geometry_solid = dynamic_cast<G4VSolid*>(ring_solid);
    }

    G4Material* pmma = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS");
    pmma->SetMaterialPropertiesTable(opticalprops::FakeG2P_FB118(TS_cromophore_concentration_, 1.502, true));
    G4LogicalVolume* ring_logic = new G4LogicalVolume(geometry_solid, pmma, ring_name, 0, 0, 0, true);

    new G4PVPlacement(  nullptr, G4ThreeVector(0., 0., 0.),
                        ring_name, ring_logic,
                        mother_physical, false, 0, true); 

    return;
  }  

  void XArapucaDS::ConstructPhotosensors(G4VPhysicalVolume* mother_physical) const
  {
    G4VisAttributes sipm_col = nexus::Red();
    sipm_col.SetForceSolid(true);

    if(PS_config_code_==1){
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
        else{
          sipm_ptr = new PerfectSiPMMPPC();
        }

        sipm_ptr->SetReflectiveSupports(ref_phsensors_supports_);
        sipm_ptr->Construct();
        G4double sipm_thickn = sipm_ptr->GetThickness();
        G4LogicalVolume* sipm_logic_vol = sipm_ptr->GetLogicalVolume();

        if (!sipm_logic_vol) {
        G4Exception("[XArapucaDS]", "ConstructPhotosensors()",
                    FatalException, "Null pointer to logical volume.");
        }

        sipm_logic_vol->SetVisAttributes(sipm_col);

        G4int phsensor_id = 0;

        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateX(-90.*deg);
        G4double aux = plate_length_ + 2.*plate_ring_gap_ +2.*ring_thickness_;
        for (G4int i=0; i<num_phsensors_; ++i) {
            G4ThreeVector pos(  -1.*aux/2. + (0.5 + i) * aux/num_phsensors_,
                                0.,
                                -plate_width_/2. -ring_thickness_ -sipm_thickn/2. -gap_);

            new G4PVPlacement(  rot, pos,
                                sipm_ptr->GetModel(), sipm_logic_vol,
                                mother_physical, false, phsensor_id, true);
            phsensor_id += 1;
        }

        G4RotationMatrix* rot2 = new G4RotationMatrix();
        rot2->rotateX(90.*deg);
        for (G4int i=0; i<num_phsensors_; ++i) {
            G4ThreeVector pos(  -1.*aux/2. + (0.5 + i) * aux/num_phsensors_,
                                0.,
                                +plate_width_/2. +ring_thickness_ +sipm_thickn/2. +gap_);

            new G4PVPlacement(rot2, pos,
                                sipm_ptr->GetModel(), sipm_logic_vol,
                                mother_physical, false, phsensor_id, true);
            phsensor_id += 1;
        }

        if(!only_sipms_along_long_sides_){
            G4RotationMatrix* rot3 = new G4RotationMatrix();
            rot3->rotateZ(90.*deg);
            aux = plate_width_ + 2.*plate_ring_gap_ +2.*ring_thickness_;
            for (G4int i=0; i<num_phsensors_; ++i) {
                G4ThreeVector pos(  -plate_length_/2. -ring_thickness_ -sipm_thickn/2. -gap_,
                                    0.,
                                    -1.*aux/2. + (0.5 + i) * aux/num_phsensors_);

                new G4PVPlacement(  rot3, pos,
                                    sipm_ptr->GetModel(), sipm_logic_vol,
                                    mother_physical, false, phsensor_id, true);
                phsensor_id += 1;
            }

            G4RotationMatrix* rot4 = new G4RotationMatrix();
            rot4->rotateZ(-90.*deg);
            for (G4int i=0; i<num_phsensors_; ++i) {
                G4ThreeVector pos(  +plate_length_/2. +ring_thickness_ +sipm_thickn/2. +gap_,
                                    0.,
                                    -aux/2. + (0.5 + i) * aux/num_phsensors_);

                new G4PVPlacement(  rot4, pos,
                                    sipm_ptr->GetModel(), sipm_logic_vol,
                                    mother_physical, false, phsensor_id, true);
                phsensor_id += 1;
            }
        }

    }
    else if(PS_config_code_==2){

        ScalableHamamatsuS133606050VE sipm_along_length(plate_length_, plate_thickn_);
        sipm_along_length.SetReflectiveSupports(ref_phsensors_supports_);
        sipm_along_length.Construct();
        G4double sipm_thickn = sipm_along_length.GetThickness();
        G4LogicalVolume* sipm_along_length_logic_vol = sipm_along_length.GetLogicalVolume();
        if(!sipm_along_length_logic_vol){
        G4Exception("[XArapucaDS]", "ConstructPhotosensors()",
                    FatalException, "Null pointer to logical volume.");
        }
        sipm_along_length_logic_vol->SetVisAttributes(sipm_col);

        G4RotationMatrix* rot0 = new G4RotationMatrix();
        rot0->rotateX(-90.*deg);
        G4ThreeVector pos0( 0.,
                            0.,
                            -plate_width_/2. -ring_thickness_ -sipm_thickn/2. -gap_);
        new G4PVPlacement(rot0, pos0,
                            "SS133606050VE_MPPC", sipm_along_length_logic_vol,
                            mother_physical, false, 0, true);

        G4RotationMatrix* rot1 = new G4RotationMatrix();
        rot1->rotateX(90.*deg);
        G4ThreeVector pos1( 0.,
                            0.,
                            +plate_width_/2. +ring_thickness_ +sipm_thickn/2. +gap_);
        new G4PVPlacement(rot1, pos1,
                            "SS133606050VE_MPPC", sipm_along_length_logic_vol,
                            mother_physical, false, 1, true);

        if(!only_sipms_along_long_sides_){

            ScalableHamamatsuS133606050VE sipm_along_width(plate_width_, plate_thickn_);
            sipm_along_width.SetReflectiveSupports(ref_phsensors_supports_);
            sipm_along_width.Construct();
            G4double sipm_thickn = sipm_along_width.GetThickness();
            G4LogicalVolume* sipm_along_width_logic_vol = sipm_along_width.GetLogicalVolume();
            if (!sipm_along_width_logic_vol) {
            G4Exception("[XArapucaDS]", "ConstructPhotosensors()",
                        FatalException, "Null pointer to logical volume.");
            }
            sipm_along_width_logic_vol->SetVisAttributes(sipm_col);

            G4RotationMatrix* rot2 = new G4RotationMatrix();
            rot2->rotateY(-90.*deg);
            rot2->rotateX(-90.*deg);
            G4ThreeVector pos2(  -plate_length_/2. -ring_thickness_ -sipm_thickn/2. -gap_,
                                0.,
                                0.);
            new G4PVPlacement(rot2, pos2,
                                "SS133606050VE_MPPC", sipm_along_width_logic_vol,
                                mother_physical, false, 2, true);

            G4RotationMatrix* rot3 = new G4RotationMatrix();
            rot3->rotateY(-90.*deg);
            rot3->rotateX(90.*deg);
            G4ThreeVector pos3(  +plate_length_/2. +ring_thickness_ +sipm_thickn/2. +gap_,
                                0.,
                                0.);
            new G4PVPlacement(rot3, pos3,
                                "SS133606050VE_MPPC", sipm_along_width_logic_vol,
                                mother_physical, false, 3, true);
        }

    }
    else{
        G4Exception(    "[XArapucaDS]", "ConstructPhotosensors()",
                        FatalException, "Not allowed configuration code.");

    }

    return;
  }

  void XArapucaDS::ConstructBoards(G4VPhysicalVolume* mother_physical) const
  {

    if(config_code_==1){

        G4double aux = 2.*plate_ring_gap_ +2.*ring_thickness_; 

        SiPMBoard board1;
        board1.SetBaseID(0);
        board1.SetBoardLength(plate_length_+aux);
        board1.SetSiPMCode(SiPM_code_);
        board1.SetNumPhsensors(num_phsensors_);
        board1.SetReflectiveSupports(ref_phsensors_supports_);
        board1.Construct();
        G4LogicalVolume* board1_logic_vol = board1.GetLogicalVolume();

        G4double z_pos = (plate_width_/2.) +plate_ring_gap_ +ring_thickness_ +(board1.GetOverallThickness()/2.) +gap_;
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., z_pos),
                        "SIPMS_BOARD", board1_logic_vol, mother_physical, false, 0, false);
        // SiPMBoard logical volume is an encasing volume which may collide into other volumes
        // No need to set pSurfCheck for that volume (dimples). As we are setting pSurfCheck=false
        // (so that no harmless overlap warning pops up in a with-dimples configuration), you have to
        // be extra careful to examine when there's actually a problematic overlap of this volume with
        // another one (since Geant4 won't warn you).

        SiPMBoard board2;
        board2.SetBaseID(num_phsensors_);
        board2.SetBoardLength(plate_length_+aux);
        board2.SetSiPMCode(SiPM_code_);
        board2.SetNumPhsensors(num_phsensors_);
        board2.SetReflectiveSupports(ref_phsensors_supports_);
        board2.Construct();
        G4LogicalVolume* board2_logic_vol = board2.GetLogicalVolume();

        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateY(+180.*deg);

        z_pos = z_pos = (plate_width_/2.) +plate_ring_gap_ +ring_thickness_ +(board1.GetOverallThickness()/2.) +gap_;
        new G4PVPlacement(rot, G4ThreeVector(0., 0., -1.*z_pos),
                        "SIPMS_BOARD", board2_logic_vol, mother_physical, true, 0, false);

        if(!only_sipms_along_long_sides_){
            SiPMBoard board3;
            board3.SetBaseID(2*num_phsensors_);
            board3.SetBoardLength(plate_width_+aux);
            board3.SetSiPMCode(SiPM_code_);
            board3.SetNumPhsensors(num_phsensors_);
            board3.SetReflectiveSupports(ref_phsensors_supports_);
            board3.Construct();
            G4LogicalVolume* board3_logic_vol = board3.GetLogicalVolume();

            G4RotationMatrix* rot2 = new G4RotationMatrix();
            rot2->rotateY(-90.*deg);

            G4double x_pos = (plate_length_/2.) +plate_ring_gap_ +ring_thickness_ +(board3.GetOverallThickness()/2.) +gap_;
            new G4PVPlacement(rot2, G4ThreeVector(x_pos, 0., 0.),
                            "SIPMS_BOARD", board3_logic_vol, mother_physical, true, 0, false);

            SiPMBoard board4;
            board4.SetBaseID(3*num_phsensors_);
            board4.SetBoardLength(plate_width_+aux);
            board4.SetSiPMCode(SiPM_code_);
            board4.SetNumPhsensors(num_phsensors_);
            board4.SetReflectiveSupports(ref_phsensors_supports_);
            board4.Construct();
            G4LogicalVolume* board4_logic_vol = board4.GetLogicalVolume();

            G4RotationMatrix* rot3 = new G4RotationMatrix();
            rot3->rotateY(+90.*deg);

            x_pos = (plate_length_/2.) +plate_ring_gap_ +ring_thickness_ +(board3.GetOverallThickness()/2.) +gap_;
            new G4PVPlacement(rot3, G4ThreeVector(-1.*x_pos, 0., 0.),
                            "SIPMS_BOARD", board4_logic_vol, mother_physical, true, 0, false);
        }
    }
    else{
        G4Exception("[XArapucaDS]", "ConstructBoards()",
                    FatalException, "Not allowed configuration code.");
    }

    return;

  }

  void XArapucaDS::ConstructReflectiveCase(G4VPhysicalVolume* mother_physical) const
  {
    
    // Based on arxiv.org/abs/1912.09191 and TDR vol.IX, the reflective foils
    // fully cover the internal cavity of the X-Arapuca.

    const G4String ref_case_name = "REF_CASE";

    // Get reflective case volume as subtraction solid from two boxes
    G4Box* aux_outter_box = new G4Box("AUX_OUTTER_BOX", 
                            overall_length_/2., internal_thickn_/2., overall_width_/2.);  

    // Extra thickness to prevent boolean subtraction of solids with matching surfaces
    // See geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html#solids-made-by-boolean-operations
    G4double tolerance = 1.*cm;

    G4Box* aux_inner_box =  new G4Box("AUX_INNER_BOX", 
                            internal_length_/2., (internal_thickn_+tolerance)/2., internal_width_/2.);
    G4SubtractionSolid* ref_case_solid =    new G4SubtractionSolid(ref_case_name, aux_outter_box, aux_inner_box, 
                                            nullptr, G4ThreeVector(0., 0., 0.));

    // The reflective case is 'open' on both sides. If the X-ARAPUCA DS is not 
    // double sided, then Construct() will take care of placing a reflective cover.

    G4LogicalVolume* ref_case_logic = 
      new G4LogicalVolume(ref_case_solid, materials::FR4(), ref_case_name);

    // Set its color for visualization purposes
    G4VisAttributes ref_case_col = nexus::WhiteAlpha();
    //ref_case_col.SetForceSolid(true);
    ref_case_logic->SetVisAttributes(ref_case_col);

    //Now create the reflective optical surface
    const G4String refsurf_name = "REF_SURFACE";
    G4OpticalSurface* refsurf_opsurf = 
      new G4OpticalSurface(refsurf_name, unified, ground, dielectric_metal, 1);
    
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
    new G4LogicalSkinSurface("REF_CASE_SURFACE", ref_case_logic, refsurf_opsurf);   
    
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.),
          ref_case_name, ref_case_logic, 
          mother_physical,
          false, 0, true);

    return;

  }
  
  void XArapucaDS::ConstructSubstratesAssemblies(G4VPhysicalVolume* mother_physical) const
  {
      
    // ------------------------------------ SUBSTRATES MODEL ------------------------------------
    //
    // Coating-substrates assembly (CSA) is made up of one-piece frame and a bunch of coated 
    // substrates which are inserted in the frame. The volume that faces the internal cavity of 
    // the X-Arapuca DS is meant to be the substrate (a piece of glass with good transmitance for 
    // the coating light). The volume that faces the cryostat is the coating volume, with the 
    // coating optical properties, such as WLS absorption length, refractive index and the 
    // following rough surfaces: LAr->coating, coating->LAr and coating->substrate. substrate 
    // is a volume with the substrate optical properties, such as absorption length and refractive 
    // index. Since both the coating and the substrate have well-defined refractive indices, this 
    // model is able to simulate the light behaviour in the coating-substrate interface, which has 
    // not been measured in the laboratory.  The resulting model is:
    //
    //                                            LAr (Cryostat)
    //
    //                       __________                               __________
    //                                |                               |
    //                                |                               |
    //                                |                               |
    //                                |                               |
    //                                |                               |
    //                                |                               |
    //                                |                               |
    //                                |-------------------------------|
    //                                |            coating            |
    //                                |-------------------------------|
    //                                |                               |
    //                        frame   |            substrate          | frame
    //                       _________|_______________________________|_________
    //                                                
    //                                 LAr (X-ARAPUCA internal cavity)
    //
    // For a common configuration, the coating is PTP with refractive index equal to 1.65, and
    // the substrate is fused silica (R. index approximately equal to 1.47 for PTP light).
    //
    // -----------------------------------------------------------------------------------------------

    G4Box* cover_solid = new G4Box("AUX", CSA_length_/2., CSA_thickn_/2., CSA_width_/2.);
    G4double tolerance = 5.*mm; // Tolerance to prevent matching surfaces when subtracting carvings from frame
    
    G4Box* carving_solid                    =  new G4Box("FRAME_CARVING",               CS_length_/2., CSA_thickn_/2. +tolerance,   CS_width_/2.); 
    G4Box* substrate_solid                  =  new G4Box("SUBSTRATE",   CS_length_/2., CS_thickn_/2.,               CS_width_/2.); 
    G4Box* coating_solid                    =  nullptr;
    if(substrates_are_coated_) coating_solid=  new G4Box("DF_COATING", CS_length_/2., coating_thickn_/2., CS_width_/2.); 

    G4MultiUnion* frame_carvings_multiunion_solid       = new G4MultiUnion("FRAME_CARVINGS");
    G4MultiUnion* substrates_multiunion_solid           = new G4MultiUnion("SUBSTRATES");
    G4MultiUnion* coatings_multiunion_solid             = nullptr;
    if(substrates_are_coated_) coatings_multiunion_solid= new G4MultiUnion("COATINGS"); 

    G4double x_pos, z_pos;
    G4Transform3D* transform_ptr = nullptr;
    G4RotationMatrix* rot = new G4RotationMatrix();

    for(G4int i=0; i<cs_no_along_wlsplength_; i++){
        for(G4int j=0; j<cs_no_along_wlspwidth_; j++){
            x_pos = (-1.*(CSA_length_/2.))  +outter_frame_width_along_wlsplength_   +(i*inner_frames_width_along_wlsplength_)   +((i+0.5)*CS_length_);
            z_pos = (-1.*(CSA_width_/2.))   +outter_frame_width_along_wlspwidth_    +(j*inner_frames_width_along_wlspwidth_)    +((j+0.5)*CS_width_);
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector(x_pos, 0., z_pos));

            frame_carvings_multiunion_solid                     ->AddNode(*carving_solid,       *transform_ptr);
            substrates_multiunion_solid                         ->AddNode(*substrate_solid,     *transform_ptr);
            if(substrates_are_coated_) coatings_multiunion_solid->AddNode(*coating_solid,       *transform_ptr);
        }
    }
    
    frame_carvings_multiunion_solid                     ->Voxelize();
    substrates_multiunion_solid                         ->Voxelize();
    if(substrates_are_coated_) coatings_multiunion_solid->Voxelize();

    // FRAME //
    if(!remove_CSA_frame_){
        const G4String frame_name = "COATED-SUBSTRATES_ASSEMBLY_FRAME";
        G4SubtractionSolid* frame_solid = new G4SubtractionSolid(frame_name, cover_solid,
                                                                frame_carvings_multiunion_solid, 
                                                                rot, G4ThreeVector(0., 0., 0.));

        // At some point there was a segfault bug with the logical frame. I am suspicious that it has to do with
        // the fact that for the case where the internal cavity thickness matches the height of the SiPM board, there are
        // two logical skin surfaces matching in the space: The one of the SiPM board and the one of the frame.
        // To prevent this, make internal_thickn_ slightly bigger than the height of the SiPM board.

        G4LogicalVolume* frame_logic = new G4LogicalVolume(frame_solid, materials::FR4(), frame_name);

        // Set its color for visualization purposes
        G4VisAttributes frame_col = nexus::DarkGreen();
        //frame_col.SetForceSolid(true);
        frame_logic->SetVisAttributes(frame_col);
        //frame_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

        if(CSA_frame_is_reflective_){
            const G4String refcoat_name = "REF_COATING";
            G4OpticalSurface* refcoat_opsurf = 
            new G4OpticalSurface(refcoat_name, unified, ground, dielectric_metal, 1);
            if(CSA_frame_is_specular_)  refcoat_opsurf->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());
            else                        refcoat_opsurf->SetMaterialPropertiesTable(opticalprops::diffusiveVIKUITI());
            new G4LogicalSkinSurface(refcoat_name, frame_logic, refcoat_opsurf);
        }
        
        G4PVPlacement* frame_physical = new G4PVPlacement(nullptr, G4ThreeVector(0., (internal_thickn_+CSA_thickn_)/2. +0.*cm, 0.),
                                        frame_name, frame_logic, mother_physical, true, 0, true);

        if(double_sided_){
            new G4PVPlacement(nullptr, G4ThreeVector(0., -1.*(internal_thickn_+CSA_thickn_)/2., 0.), 
                                frame_name, frame_logic, mother_physical, true, 1, true);
        }
        else{
            // If the X-ARAPUCA DS is not double-sided, then the dummy plane is reflective
            const G4String refsurf_name = "REF_SURFACE";
            G4OpticalSurface* refsurf_opsurf = 
            new G4OpticalSurface(refsurf_name, unified, ground, dielectric_metal, 1);
            refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());

            const G4String cover_name = "REFLECTIVE_COVER";
            G4LogicalVolume* cover_logic = new G4LogicalVolume(cover_solid, materials::FR4(), cover_name);
            cover_logic->SetVisAttributes(frame_col);
            new G4LogicalSkinSurface(refsurf_name, cover_logic, refsurf_opsurf);   
            new G4PVPlacement(nullptr, G4ThreeVector(0., -1.*(internal_thickn_+CSA_thickn_)/2., 0.), 
                                cover_name, cover_logic, mother_physical, true, 1, true);
        }
    }

    // DICHROIC FILTERS //
    if(!remove_CSs_){
        G4Material* substrate_mat    = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
        substrate_mat->SetMaterialPropertiesTable(CS_mpt_);
        G4LogicalVolume* substrates_multiunion_logic = 
                            new G4LogicalVolume(substrates_multiunion_solid, substrate_mat, "SUBSTRATES");

        // Place the filters
        G4VPhysicalVolume* first_substrates_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
            new G4PVPlacement(nullptr, G4ThreeVector(0., (internal_thickn_/2.) 
                                                        +(CS_thickn_/2.) 
                                                        +(CS_pos_wrt_CSA_pos_*(CSA_thickn_-CS_thickn_)), 0.), 
                            "SUBSTRATES_1", substrates_multiunion_logic, mother_physical, true, 0, true));
            
        G4VPhysicalVolume* second_substrates_multiunion_physical = nullptr;
        if(double_sided_){
            second_substrates_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                new G4PVPlacement(nullptr, G4ThreeVector(0.,    -1.*(
                                                                (internal_thickn_/2.) 
                                                                +(CS_thickn_/2.) 
                                                                +(CS_pos_wrt_CSA_pos_*(CSA_thickn_-CS_thickn_))), 0.), 
                                "SUBSTRATES_2", substrates_multiunion_logic, mother_physical, true, 1, true));
        }

        // COATINGS //
        if(substrates_are_coated_){
            G4Material* coatings_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
            coatings_mat->SetMaterialPropertiesTable(opticalprops::PTP());
            G4LogicalVolume* coatings_multiunion_logic = 
                                new G4LogicalVolume(coatings_multiunion_solid, coatings_mat, "COATINGS");   

            G4VisAttributes coatings_col = nexus::TitaniumGreyAlpha();
            coatings_col.SetForceSolid(true);
            coatings_multiunion_logic->SetVisAttributes(coatings_col);

            // Place the coatings
            G4VPhysicalVolume* first_coatings_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                new G4PVPlacement(  nullptr, G4ThreeVector(0.,  (internal_thickn_/2.)
                                                                +CS_thickn_
                                                                +(CS_pos_wrt_CSA_pos_*(CSA_thickn_-CS_thickn_)) 
                                                                +(coating_thickn_/2.), 0.), 
                                    "COATINGS_1", coatings_multiunion_logic, mother_physical, true, 0, true));

            // Make the LAR-coating interface rough, so that photons cannot be trapped within the coating
            G4OpticalSurface* coating_rough_surf =
                    new G4OpticalSurface("COATING_ROUGH_SURFACE", glisur, ground, dielectric_dielectric, .01);
                    // 0.01 is the polish value for glisur model that was measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
                    // This is the best reference we have, since both PTP and TPB are the result of an evaporation+deposition process
            new G4LogicalBorderSurface( "SURROUNDINGS->COATING1", 
                                        mother_physical, 
                                        first_coatings_multiunion_physical, 
                                        coating_rough_surf);
            new G4LogicalBorderSurface( "COATING1->SURROUNDINGS", 
                                        first_coatings_multiunion_physical, 
                                        mother_physical, 
                                        coating_rough_surf);
            // We will also add roughness for the coating->substrate interface, but only with such ordering. The alternative case takes 
            // place when the photon travels from the substrate to the coating. The glass is supposed to be polished, so the photon may 
            // not see a rough surface.
            new G4LogicalBorderSurface( "COATING1->SUBSTRATE1", 
                                        first_coatings_multiunion_physical, 
                                        first_substrates_multiunion_physical, 
                                        coating_rough_surf);

            if(double_sided_){
                G4VPhysicalVolume* second_coatings_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                    new G4PVPlacement(  nullptr, G4ThreeVector(0.,  -1.*((internal_thickn_/2.)
                                                                    +CS_thickn_
                                                                    +(CS_pos_wrt_CSA_pos_*(CSA_thickn_-CS_thickn_)) 
                                                                    +(coating_thickn_/2.)), 0.), 
                                        "COATINGS_2", coatings_multiunion_logic, mother_physical, true, 1, true));
                new G4LogicalBorderSurface( "SURROUNDINGS->COATING2", 
                                            mother_physical, 
                                            second_coatings_multiunion_physical, 
                                            coating_rough_surf);
                new G4LogicalBorderSurface( "COATING2->SURROUNDINGS", 
                                            second_coatings_multiunion_physical, 
                                            mother_physical, 
                                            coating_rough_surf);
                new G4LogicalBorderSurface( "COATING2->SUBSTRATE2", 
                                            second_coatings_multiunion_physical, 
                                            second_substrates_multiunion_physical, 
                                            coating_rough_surf);
            }
        }
    }
    return;
  }

  // COME BACK HERE: check the dimensions of the big photosensor in case PS_config_code_==2 is true
  // COME BACK HERE: REWRITE THE CONSTRUCTION OF THE BOARDS/PHOTOSENSORS

  G4bool XArapucaDS::geometry_is_ill_formed()
  {

    // Some checks are independent of the configuration code. Perform those tests here:
    if(CS_thickn_>CSA_thickn_){
        G4Exception("[XArapucaDS]", "geometry_is_ill_formed()", FatalException,
        "The coated-substrates thickness cannot be bigger than the frame thickness.");
    }

    if(CS_pos_wrt_CSA_pos_<0. || CS_pos_wrt_CSA_pos_>1.){
        G4Exception("[XArapucaDS]", "geometry_is_ill_formed()", FatalException,
        "The dichroic filter position with respect to the CSA position must belong to the [0, 1] interval.");
    }

    // What we need to make sure here is that there's room enough within the XArapuca DS internal cavity
    // to allocate everything that we intending to put in it. To do so, we may find the span of the 
    // internal geometry along each axis, and then compare it to the dimensions of the internal cavity.

    if(internal_length_<internal_width_){
        G4Exception("[XArapucaDS]", "geometry_is_ill_formed()", FatalException,
        "For consistency with the code, internal_length_ must be greater than internal_width_.");
    }
    
    G4double internal_geom_length_span, internal_geom_width_span, internal_geom_thickn_span;
    if(config_code_==1){
        if(plate_length_<plate_width_){
            G4Exception("[XArapucaDS]", "geometry_is_ill_formed()", FatalException,
            "For consistency with the code, plate_length_ must be greater than plate_width_.");
        }

        if(with_boards_){
            // If with_boards_==true, then PS_config_code_ is ignored, the geometry is loaded
            // as if PS_config_code_==1, so there's no need to distinguish cases here
            SiPMBoard board;
            internal_geom_length_span   = plate_length_ +2.*plate_ring_gap_ +2.*ring_thickness_;
            if(!only_sipms_along_long_sides_){
                internal_geom_length_span   += (2.*gap_)+(2.*board.GetOverallThickness());    
            }
            internal_geom_width_span    = plate_width_ +2.*plate_ring_gap_ +2.*ring_thickness_+(2.*gap_)+(2.*board.GetOverallThickness());
            internal_geom_thickn_span   = std::max(plate_thickn_, board.GetOverallHeight());
        }
        else{
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
            else{
              sipm_ptr = new PerfectSiPMMPPC();
            }
            if(PS_config_code_==1){
                G4double aux = plate_length_ +2.*plate_ring_gap_ +2.*ring_thickness_;
                if(num_phsensors_*sipm_ptr->GetTransverseDim()>aux) { return true; }
                aux = plate_width_ +2.*plate_ring_gap_ +2.*ring_thickness_;
                if(!only_sipms_along_long_sides_){
                    if(num_phsensors_*sipm_ptr->GetTransverseDim()>aux) { return true; }
                }
            }
            // Else if PS_config_code_==2, the scalable SiPM is build with the exact plate 
            // dimensions, so there's no need to care about such case here

            internal_geom_length_span = plate_length_ +2.*plate_ring_gap_   // Either the HamamatsuS133606050VE and the ScalableHamamatsuS133606050VE
                                        +2.*ring_thickness_ +(2.*gap_)      // have the same thickness, so there's no need to distinguish between both
                                        +(2.*sipm_ptr->GetThickness());     // according to PS_config_code_ value
                                                                                                                                            
            internal_geom_width_span =  plate_width_+2.*plate_ring_gap_
                                        +2.*ring_thickness_+(2.*gap_)
                                        +(2.*sipm_ptr->GetThickness());

            internal_geom_thickn_span = std::max(plate_thickn_, sipm_ptr->GetTransverseDim());
            delete sipm_ptr;
        }
    }
    else{
        G4Exception("[XArapucaDS]", "geometry_is_ill_formed()",
                    FatalException, "Not allowed configuration code.");
    }
    
    G4bool check1 = internal_length_ >= internal_geom_length_span;
    G4bool check2 = internal_thickn_ >= internal_geom_thickn_span;
    G4bool check3 = internal_width_  >= internal_geom_width_span;

    // Uncomment the following chunk for testing purposes
    if (!check1) {
        G4cout << "WARNING: First condition is not met" << G4endl;
        G4cout << "internal_length_=" << internal_length_ << ", internal_geom_length_span=" << internal_geom_length_span << G4endl;
    }
    if (!check2) {
        G4cout << "WARNING: Second condition is not met" << G4endl;
        G4cout << "internal_thickn_=" << internal_thickn_ << ", internal_geom_thickn_span=" << internal_geom_thickn_span << G4endl;
    }
    if (!check3) {
        G4cout << "WARNING: Third condition is not met" << G4endl;
        G4cout << "internal_width_=" << internal_width_ << ", internal_geom_width_span=" << internal_geom_width_span << G4endl;
    }

    return !(check1*check2*check3);
  }

  G4ThreeVector XArapucaDS::GenerateVertex(const G4String&) const{

    G4double tolerance = 0.1*mm;    // Small distance over the CSA from which
                                    // photons are launched. Also, the width 
                                    // of the outter border projected over the 
                                    // CS from which photons won't be launched 
                                    // (Just see the implementation in x_pos and 
                                    // z_pos below to understand its meaning)

    G4double x_pos, z_pos;
    G4double y_pos = (internal_thickn_/2.) +CSA_thickn_ +tolerance;

    if(generation_region_=="dichroic"){
      std::random_device rd;  //Will be used to obtain a seed for the random number engine
      std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
      std::uniform_int_distribution<> dist(0, cs_no_along_wlsplength_*cs_no_along_wlspwidth_ -1);
      G4int filter_no = dist(gen);

      G4int column_no = filter_no%cs_no_along_wlsplength_;                // \in[0, df_no_along_wlsplength -1]  
      G4int row_no = (int) std::floor(filter_no/cs_no_along_wlsplength_); // \in[0, df_no_along_wlspwidth -1]

      G4double selected_filter_x_center =  (-1.*(CSA_length_/2.)) +outter_frame_width_along_wlsplength_   +(column_no*inner_frames_width_along_wlsplength_)   +((column_no+0.5)*CS_length_);
      G4double selected_filter_z_center =  (-1.*(CSA_width_/2.))  +outter_frame_width_along_wlspwidth_    +(row_no*inner_frames_width_along_wlspwidth_)       +((row_no+0.5)*CS_width_);

      x_pos = UniformRandomInRange(  selected_filter_x_center -(CS_length_/2.) +tolerance, 
                                              selected_filter_x_center +(CS_length_/2.) -tolerance    );
      z_pos = UniformRandomInRange(  selected_filter_z_center -(CS_width_/2.) +tolerance, 
                                              selected_filter_z_center +(CS_width_/2.) -tolerance     );
    }
    else if(generation_region_=="custom"){
      G4double random_radius =  UniformRandomInRange(gen_diameter_/2., 0.);
      G4double random_angle =   UniformRandomInRange(twopi, 0.); 
      x_pos = gen_x_ +(random_radius*sin(random_angle));
      z_pos = gen_z_ +(random_radius*cos(random_angle));
    }
    else{ // Default behaviour is that of generation_region_=="random"
      x_pos = UniformRandomInRange(  -1.*CSA_length_/2.,
                                      CSA_length_/2.      );
      z_pos = UniformRandomInRange(  -1.*CSA_width_/2.,
                                      CSA_width_/2.       );
    }
    return G4ThreeVector(x_pos, y_pos, z_pos);
  }

} //End namespace nexus
