#include "XArapucaRods.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
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
#include <G4Para.hh>
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

#include <CLHEP/Units/SystemOfUnits.h>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(XArapucaRods, GeometryBase)

namespace nexus{

  XArapucaRods::XArapucaRods():
  GeometryBase(), 
  ///Get internal reflector cavity dimensions from arxiv.org/abs/1912.09191
  surrounding_media_                    ("lar"                        ),
  DFA_thickn_                           (1.     *mm                   ),
  DF_thickn_                            (1.010  *mm                   ),
  DF_substrate_thickn_                  (1.000  *mm                   ),
  DF_substrate_mpt_                     (opticalprops::FusedSilica()  ),
  //DF_substrate_mpt_                     (opticalprops::SCHOTT_B270()  ),      // EY! YOU CHANGED THIS!
  DF_pos_wrt_DFA_pos_                   (0.                           ),
  DF_are_coated_                        (true                         ),
  coating_thickn_                       (3.226  *um                   ),  // Based on arxiv.org/abs/1912.09191 and TDR vol.IX, section 5.8.3.1,
                                                                          // the pTP film thickness is such that there's 400 micrograms of pTP
                                                                          // deposited over each square centimeter of DF.
                                                                          // This, together with the pTP density, (1.24g/cm3, found in
                                                                          // en.wikipedia.org/wiki/Terphenyl), gives a pTP film thickness of 
                                                                          // 3.226 micrometers
  coating_rindex_                       (1.65                         ),  
  outter_frame_width_along_wlsplength_  (10.    *mm                   ),
  outter_frame_width_along_wlspwidth_   (10.    *mm                   ),
  inner_frames_width_along_wlsplength_  (6.     *mm                   ),
  inner_frames_width_along_wlspwidth_   (2.     *mm                   ),
  df_no_along_wlsplength_               (4                            ),
  df_no_along_wlspwidth_                (4                            ),
  DFA_config_code_                      (1                            ),
  DFA_frame_is_reflective_              (false                        ),
  DFA_frame_is_specular_                (false                        ),
  remove_DFs_                           (false                        ),  
  remove_DFA_frame_                     (false                        ),
  secondary_wls_attlength_              (1.     *m                    ),
  cromophore_concentration_             (16.                          ),
  case_thickn_                          (1.     *mm                   ),   ///Get foil thickness from isoltronic.ch/assets/of-m-vikuiti-esr-app-guide.pdf
  SiPM_code_                            (1                            ),
  num_phsensors_                        (24                           ),
  rod_sipm_gap_                         (0.5    *mm                   ),
  ref_phsensors_supports_               (true                         ), 
  double_sided_                         (true                         ),
  generation_region_                    ("random"                     ),
  gen_x_                                (0.     *cm                   ),
  gen_z_                                (0.     *cm                   ),
  gen_diameter_                         (1.     *cm                   ),
  path_to_inwards_dichroic_data_        (""                           ),
  path_to_outwards_dichroic_data_       (""                           ),
  world_extra_thickn_                   (100.   *cm                   ),
  tunneling_probability_                (0.0                          ),
  add_blocks_between_sipms_             (false                        ),
  rod_length_                           (500.   *mm                   ),
  rod_height_                           (6.     *mm                   ),
  rod_width_                            (6.     *mm                   ),
  rods_no_                              (80                           ),
  add_reflective_separators_            (true                         ),
  rods_lateral_gap_                     (2.*mm                        ),
  along_long_side_                      (true                         ),
  cut_rods_                             (false                        ),
  cut_angle_                            (45.*deg                      ),
  cut_thickness_                        (1.*mm                        ),
  place_separator_at_the_cut_           (false                        )
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/XArapucaRods/",
				"Control commands of geometry XArapucaRods.");

    G4GenericMessenger::Command& sm_cmd =
      msg_->DeclareProperty("surrounding_media", surrounding_media_,
			    "Which media to place the XArapucaRods in.");

    G4GenericMessenger::Command& dfat_cmd =
      msg_->DeclareProperty("DFA_thickn", DFA_thickn_,
			    "Thickness of the dichroic filters assembly.");
    dfat_cmd.SetUnitCategory("Length");
    dfat_cmd.SetParameterName("DFA_thickn", false);
    dfat_cmd.SetRange("DFA_thickn>0.");

    G4GenericMessenger::Command& dft_cmd =
      msg_->DeclareProperty("DF_thickn", DF_thickn_,
			    "Overall thickness of the dichroic filters (MLS+substrate). Must be smaller than the frame thickness.");
    dft_cmd.SetUnitCategory("Length");
    dft_cmd.SetParameterName("DF_thickn", false);
    dft_cmd.SetRange("DF_thickn>0.");

    G4GenericMessenger::Command& dfst_cmd =
      msg_->DeclareProperty("DF_substrate_thickn", DF_substrate_thickn_,
			    "Thickness of the dichroic filters substrate.");
    dfst_cmd.SetUnitCategory("Length");
    dfst_cmd.SetParameterName("DF_substrate_thickn", false);
    dfst_cmd.SetRange("DF_substrate_thickn>0.");

    G4GenericMessenger::Command& dpwdp_cmd =
      msg_->DeclareProperty("DF_pos_wrt_DFA_pos", DF_pos_wrt_DFA_pos_,
			    "Position (height) of the dichroic filters with respect to the DFA position (height). This parameter can take values from 0 to 1.");
    dpwdp_cmd.SetParameterName("DF_pos_wrt_DFA_pos", false);
    dpwdp_cmd.SetRange("DF_pos_wrt_DFA_pos>=0.");
    dpwdp_cmd.SetRange("DF_pos_wrt_DFA_pos<=1.");

    G4GenericMessenger::Command& dfac_cmd =
      msg_->DeclareProperty("DF_are_coated", DF_are_coated_,
			    "Whether the dichroic filter is set on both sides of the secondary-WLS material.");

    G4GenericMessenger::Command& ptpct_cmd =
      msg_->DeclareProperty("coating_thickn", coating_thickn_,
			    "Thickness of the coating layer that is deposited over the dichroic filters.");
    ptpct_cmd.SetUnitCategory("Length");
    ptpct_cmd.SetParameterName("coating_thickn", false);
    ptpct_cmd.SetRange("coating_thickn>0.");

    G4GenericMessenger::Command& ptpcr_cmd =
      msg_->DeclareProperty("coating_rindex", coating_rindex_,
			    "Refractive index of the coating layer.");
    ptpcr_cmd.SetParameterName("coating_rindex", false);
    ptpcr_cmd.SetRange("coating_rindex>1.");

    G4GenericMessenger::Command& ofwawl_cmd =
      msg_->DeclareProperty("outter_frame_width_along_wlsplength", outter_frame_width_along_wlsplength_,
			    "Width of the outter DFA frame along the secondary-WLS material length.");
    ofwawl_cmd.SetUnitCategory("Length");
    ofwawl_cmd.SetParameterName("outter_frame_width_along_wlsplength", false);
    ofwawl_cmd.SetRange("outter_frame_width_along_wlsplength>0.");

    G4GenericMessenger::Command& ofwaww_cmd =
      msg_->DeclareProperty("outter_frame_width_along_wlspwidth", outter_frame_width_along_wlspwidth_,
			    "Width of the outter DFA frame along the WLS secondary-WLS material width.");
    ofwaww_cmd.SetUnitCategory("Length");
    ofwaww_cmd.SetParameterName("outter_frame_width_along_wlspwidth", false);
    ofwaww_cmd.SetRange("outter_frame_width_along_wlspwidth>0.");

    G4GenericMessenger::Command& ifwawl_cmd =
      msg_->DeclareProperty("inner_frames_width_along_wlsplength", inner_frames_width_along_wlsplength_,
			    "Width of the inner DFA frames along the WLS secondary-WLS material length.");
    ifwawl_cmd.SetUnitCategory("Length");
    ifwawl_cmd.SetParameterName("inner_frames_width_along_wlsplength", false);
    ifwawl_cmd.SetRange("inner_frames_width_along_wlsplength>0.");

    G4GenericMessenger::Command& ifwaww_cmd =
      msg_->DeclareProperty("inner_frames_width_along_wlspwidth", inner_frames_width_along_wlspwidth_,
			    "Width of the inner DFA frames along the WLS secondary-WLS material width.");
    ifwaww_cmd.SetUnitCategory("Length");
    ifwaww_cmd.SetParameterName("inner_frames_width_along_wlspwidth", false);
    ifwaww_cmd.SetRange("inner_frames_width_along_wlspwidth>0.");

    G4GenericMessenger::Command& dfnawl_cmd =
      msg_->DeclareProperty("df_no_along_wlsplength", df_no_along_wlsplength_,
			    "Number of dichroic filters along the secondary-WLS material length.");
    dfnawl_cmd.SetParameterName("df_no_along_wlsplength", false);
    dfnawl_cmd.SetRange("df_no_along_wlsplength>=0"); 

    G4GenericMessenger::Command& dfnaww_cmd =
      msg_->DeclareProperty("df_no_along_wlspwidth", df_no_along_wlspwidth_,
			    "Number of dichroic filters along the secondary-WLS material width.");
    dfnaww_cmd.SetParameterName("df_no_along_wlspwidth", false);
    dfnaww_cmd.SetRange("df_no_along_wlspwidth>=0"); 

    G4GenericMessenger::Command& dfacc_cmd =
      msg_->DeclareProperty("DFA_config_code", DFA_config_code_,
			    "DFA configuration code.");
    dfacc_cmd.SetParameterName("DFA_config_code", false);
    dfacc_cmd.SetRange("DFA_config_code>=1"); 

    G4GenericMessenger::Command& dfafir_cmd =
      msg_->DeclareProperty("DFA_frame_is_reflective", DFA_frame_is_reflective_,
			    "Whether the FR4 DFA frame is vikuiti-coated or not.");

    G4GenericMessenger::Command& dfafis_cmd =
      msg_->DeclareProperty("DFA_frame_is_specular", DFA_frame_is_specular_,
			    "Whether the vikuiti coating of the DFA frame is specular-spikely reflective or diffusively reflective. Only makes a difference if DFA_frame_is_reflective_==True.");

    G4GenericMessenger::Command& rdfs_cmd =
      msg_->DeclareProperty("remove_DFs", remove_DFs_,
			    "Whether to remove the dichroic filters or not.");

    G4GenericMessenger::Command& rdfaf_cmd =
      msg_->DeclareProperty("remove_DFA_frame", remove_DFA_frame_,
			    "Whether to remove the dichroic filters assembly or not.");

    G4GenericMessenger::Command& swlsal_cmd =
      msg_->DeclareProperty("secondary_wls_attlength", secondary_wls_attlength_,
			    "Attenuation length of the secondary WLShifter, in case EJ286 is used.");
    swlsal_cmd.SetUnitCategory("Length");
    swlsal_cmd.SetParameterName("secondary_wls_attlength", false);
    swlsal_cmd.SetRange("secondary_wls_attlength>0.");

    G4GenericMessenger::Command& crco_cmd =
      msg_->DeclareProperty("cromophore_concentration", cromophore_concentration_,
			    "Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter, in case G2P_FB118 is used.");
    crco_cmd.SetParameterName("cromophore_concentration", false);
    crco_cmd.SetRange("cromophore_concentration>0.");

    G4GenericMessenger::Command& ct_cmd =
      msg_->DeclareProperty("case_thickn", case_thickn_,
			    "Thickness of the XArapucaRods case.");
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
			    "Number of photosensors per long side.");
    np_cmd.SetParameterName("num_phsensors", false);
    np_cmd.SetRange("num_phsensors>=0"); 

    G4GenericMessenger::Command& rsg_cmd =
      msg_->DeclareProperty("rod_sipm_gap", rod_sipm_gap_,
			    "Gap between photosensors and the secondary-WLS material.");
    rsg_cmd.SetUnitCategory("Length");
    //rsg_cmd.SetParameterName("rod_sipm_gap", false);
    //rsg_cmd.SetRange("rod_sipm_gap>0.");    // These are commented so that rod_sipm_gap_ can help modelate the immersion of the SiPMs into the flat dimple

    G4GenericMessenger::Command& rps_cmd =
      msg_->DeclareProperty("ref_phsensors_supports", ref_phsensors_supports_,
			    "Whether the photosensors supports are VIKUITI-coated.");

    G4GenericMessenger::Command& ds_cmd =
      msg_->DeclareProperty("double_sided", double_sided_,
			    "Whether the dichroic filter is set on both sides of the secondary-WLS material.");

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
			    "Diameter of the circle where the GV could be randomly sampled if generation_region_=='custom' is True.");
    gd_cmd.SetUnitCategory("Length");
    gd_cmd.SetParameterName("gen_diameter", false);
    gd_cmd.SetRange("gen_diameter>0.");

    G4GenericMessenger::Command& ptidd_cmd =
      msg_->DeclareProperty("path_to_inwards_dichroic_data", path_to_inwards_dichroic_data_,
			    "Absolute path to dichroic data file that is to be sampled for the light trying to enter the XAR cavity.");

    G4GenericMessenger::Command& ptodd_cmd =
      msg_->DeclareProperty("path_to_outwards_dichroic_data", path_to_outwards_dichroic_data_,
			    "Absolute path to dichroic data file that is to be sampled for the light trying to escape the XAR cavity.");

    G4GenericMessenger::Command& tp_cmd =
      msg_->DeclareProperty("tunneling_probability", tunneling_probability_,
			    "Probability that a photon tunnels through the surface of the secondary-WLS material (inwards or outwards).");
    tp_cmd.SetParameterName("tunneling_probability", false);
    tp_cmd.SetRange("tunneling_probability>=0.");
    tp_cmd.SetRange("tunneling_probability<=1.");

    G4GenericMessenger::Command& abbs_cmd =
      msg_->DeclareProperty("add_blocks_between_sipms", add_blocks_between_sipms_,
			    "Whether to add reflective blocks filling the gaps in between the SiPMs.");

    G4GenericMessenger::Command& rl_cmd =
      msg_->DeclareProperty("rod_length", rod_length_,
			    "Length of each rod (X).");
    rl_cmd.SetUnitCategory("Length");
    rl_cmd.SetParameterName("rod_length", false);
    rl_cmd.SetRange("rod_length>0.");

    G4GenericMessenger::Command& rh_cmd =
      msg_->DeclareProperty("rod_height", rod_height_,
			    "Height of each rod (Y).");
    rh_cmd.SetUnitCategory("Length");
    rh_cmd.SetParameterName("rod_height", false);
    rh_cmd.SetRange("rod_height>0.");

    G4GenericMessenger::Command& rw_cmd =
      msg_->DeclareProperty("rod_width", rod_width_,
			    "Width of each rod (Z).");
    rw_cmd.SetUnitCategory("Length");
    rw_cmd.SetParameterName("rod_width", false);
    rw_cmd.SetRange("rod_width>0.");

    G4GenericMessenger::Command& rn_cmd =
      msg_->DeclareProperty("rods_no", rods_no_,
			    "Number of rods that are installed within the XArapucaRods.");
    rn_cmd.SetParameterName("rods_no", false);
    rn_cmd.SetRange("rods_no>=1");

    G4GenericMessenger::Command& ars_cmd =
      msg_->DeclareProperty("add_reflective_separators", add_reflective_separators_,
			    "Whether to add reflective blocks filling the gaps in between the SiPMs.");

    G4GenericMessenger::Command& rlg_cmd =
      msg_->DeclareProperty("rods_lateral_gap", rods_lateral_gap_,
			    "Lateral gap between any two adjacent rods.");
    rlg_cmd.SetUnitCategory("Length");
    rlg_cmd.SetParameterName("rods_lateral_gap", false);
    rlg_cmd.SetRange("rods_lateral_gap>0.");

    G4GenericMessenger::Command& als_cmd =
      msg_->DeclareProperty("along_long_side", along_long_side_,
			    "Whether the rods are parallel to the long side of the XArapucaRods or not.");

    G4GenericMessenger::Command& cro_cmd =
      msg_->DeclareProperty("cut_rods", cut_rods_,
			    "Whether to split up each rod into two pieces.");

    G4GenericMessenger::Command& ca_cmd =
      msg_->DeclareProperty("cut_angle", cut_angle_,
			    "Angle of the cut, with respect to the Z axis.");
    ca_cmd.SetParameterName("cut_angle", false);
    ca_cmd.SetRange("cut_angle>=0.");
    ca_cmd.SetRange("cut_angle<=1.57079632679");  // pi/2 ~ 1.57079632679

    G4GenericMessenger::Command& cut_cmd =
      msg_->DeclareProperty("cut_thickness", cut_thickness_,
			    "Thickness of each cut/crack that is carved from the rod.");
    cut_cmd.SetUnitCategory("Length");
    cut_cmd.SetParameterName("cut_thickness", false);
    cut_cmd.SetRange("cut_thickness>0.");

    G4GenericMessenger::Command& psatc_cmd =
      msg_->DeclareProperty("place_separator_at_the_cut", place_separator_at_the_cut_,
			    "Whether to insert a reflective separator in the rod cut.");

    // When testing WLS plates with opticalprops::noAbsLength_, it is possible that a photon
    // gets trapped within the plate (below the critical angle) into an infinite-bouncing-loop.
    // For the case of an EJ286 plate with DUNE supercells dimensions, with this absorption 
    // length and immersed into LAr, some analysis revealed that particles that did not fall 
    // into this infinite loop had track lengths of less than one hundred meters.
    ul_ = new G4UserLimits();
    ul_->SetUserMaxTrackLength(100*m);
  }


  XArapucaRods::~XArapucaRods()
  {
      if(ul_) delete ul_;
  }


  void XArapucaRods::Construct()
  {

    SiPMBoard board;
    board.SetSiPMCode(SiPM_code_);

    // Compute internal attributes

    internal_length_  = rod_length_ +(2.*rod_sipm_gap_)+(2.*board.GetOverallThickness());
    internal_width_   = (rod_width_*rods_no_)   +(rods_lateral_gap_*(rods_no_-1))   +(rods_lateral_gap_*2);
                        // Width of the rods    // Width of the gap between rods    // Width of the extremal gaps
    internal_thickn_  = std::max(rod_height_, board.GetOverallHeight());
    overall_length_ = internal_length_  + 2*case_thickn_;
    overall_thickn_ = internal_thickn_  + 2*DFA_thickn_;    // geometry_is_ill_formed() will make sure that DF_thickn_<=DFA_thickn_
    overall_width_  = internal_width_   + 2*case_thickn_;
    DFA_length_ = overall_length_;                          // For now, these two match. 
    DFA_width_  = overall_width_;                           // For now, these two match.  
    DF_length_  = (DFA_length_  -(2.*outter_frame_width_along_wlsplength_)  -((df_no_along_wlsplength_-1.)*inner_frames_width_along_wlsplength_))   /(1.*df_no_along_wlsplength_);   
    DF_width_   = (DFA_width_   -(2.*outter_frame_width_along_wlspwidth_)   -((df_no_along_wlspwidth_-1.)*inner_frames_width_along_wlspwidth_))     /(1.*df_no_along_wlspwidth_);
    residual_displacement_ = 1e-6*cut_thickness_;

    G4cout << "DF_length_=" << DF_length_ << G4endl;
    G4cout << "DF_width_=" << DF_width_ << G4endl;

    if(geometry_is_ill_formed()){
      G4Exception("[XArapucaRods]", "Construct()", FatalException,
      "The given dimensions do not describe a feasible X-Arapuca.");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // The biggest volume is a vacuum box with the dimensions of the XArapucaRods device
    // plus 2*world_extra_thickn_, for each dimension. This volume is the effective 
    // world volume FOR THE NEXUS USER, but it is not the world volume of the overall
    // Geant4 application (afterwards nexus takes the user's biggest volume 
    // and places it inside another world volume whose dimensions are enough so as to 
    // fit the whole span of the biggest volume you implemented). Placed within this
    // vacuum box, there's a box, made out of surrounding_media_, with dimensions of 
    // the XArapucaRods device plus world_extra_thickn_, along each dimension. The previous 
    // version of this class simply implemented a surrounding_media_ box as the biggest 
    // volume, without the need to encapsulate it within a bigger volume (apart from 
    // that of nexus). However, as of the implementation of the dichroic filter, the 
    // use of G4LogicalBorderSurface class is necessary. Objects of such class must 
    // be given two ordered physical volumes, so that a particle flying from the first 
    // one to the second one, interacts with such an optical surface. Within this 
    // context, it is convenient that we specify the surrounding_media_ volume which 
    // encapsulates the XArapucaRods as the first physical volume, and the filter physical 
    // volume as the second volume, so that the filter properties are applied in both 
    // cases: 1) when a photon tries to enter the dichroic filter volume from outside 
    // the XArapucaRods, and 2) when a photon tries to enter the dichroic filter volume 
    // from inside the XArapucaRods (indeed, take into account that XArapucaRods different 
    // parts are placed within this surrounding_media_ box, thus giving room to 
    // surrounding_media_ gaps within the XArapucaRods cavity that belong to the 
    // surrounding_media_ box mother volume). Since the nexus user seems not to have 
    // access to the physical placement of the biggest volume (in our case, the vacuum 
    // box) which is implemented in line 80 of source/base/DetectorConstruction.cc, 
    // AFTER calling GeometryBase::Construct() (i.e. your geometry must be constructed 
    // before placing your biggest volume in nexus world volume,so you cannot possibly 
    // have access to the physical placement of your biggest volume when you Construct() 
    // it.), one way around this issue is adding this intermediate surrounding_media_ box.   
    /////////////////////////////////////////////////////////////////////////////////////
    
    // VACUUM CAPSULE
    const G4String world_name = "VACUUM_CAPSULE";

    G4Box* world_solid =
        new G4Box(world_name,
                (internal_length_)/2. +case_thickn_+world_extra_thickn_,
                (internal_thickn_)/2. +DFA_thickn_ +world_extra_thickn_,
                (internal_width_ )/2. +case_thickn_+world_extra_thickn_);
                
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
                (internal_length_)/2. +case_thickn_+(world_extra_thickn_/2.),
                (internal_thickn_)/2. +DFA_thickn_ +(world_extra_thickn_/2.),
                (internal_width_ )/2. +case_thickn_+(world_extra_thickn_/2.));

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


    ConstructReflectiveCase(mother_physical);
    if(!remove_DFs_ || !remove_DFA_frame_) ConstructDichroicAssemblies(mother_physical);
    ConstructRods(mother_physical);
    ConstructBoards(mother_physical);
    if(add_reflective_separators_) ConstructReflectiveSeparators(mother_physical);

    return;
  }

  void XArapucaRods::ConstructRods(G4VPhysicalVolume* mother_physical) const
  {
    const G4String rod_name = "ROD";
    G4VSolid* geometry_solid = nullptr;
    G4Box* whole_rod_solid = new G4Box(rod_name,  rod_length_/2., 
                                            rod_height_/2., 
                                            rod_width_/2.); 
    if(cut_rods_){
      G4Para* subtrahend_solid = new G4Para("SUBTRAHEND",
                                            (rod_length_/2.)+(cut_thickness_/2.),
                                            2.*(rod_height_/2.),  // Setting here twice the dimensions to prevent
                                            2.*(rod_width_/2.),   // matching-surfaces in boolean subtraction
                                            0.0,
                                            cut_angle_,                   // Parallelepiped angle with respect to 
                                                                          // the WLS-rod side which is dz_ long
                                            0.0);

      G4SubtractionSolid* half_rod_solid = new G4SubtractionSolid(rod_name, 
                                                                  whole_rod_solid, 
                                                                  subtrahend_solid, 
                                                                  nullptr, 
                                                                  G4ThreeVector(rod_length_/2., 0., 0.)); 

      // Placing the geometric center of the parallelepiped (subtrahend) right onto the 
      // edge of the WLS rod, so that we are left with one half of the rod minus half 
      // of the cut. The other half of the cut is carved from the other half of the rod.

      G4MultiUnion* multiunion_geometry_solid = new G4MultiUnion(rod_name);

      G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateY(180.*deg);

      multiunion_geometry_solid->AddNode( half_rod_solid, 
                                          G4Transform3D(  G4RotationMatrix{}, 
                                                          G4ThreeVector(0., 0., 0.)));  // Note that the geometric center of half_rod_solid 
                                                                                        // is that of the original (uncut) WLS rod, which is 
                                                                                        // geometrically outside the half_rod_solid volume.
      multiunion_geometry_solid->AddNode( half_rod_solid, 
                                          G4Transform3D(  *rot,
                                                          G4ThreeVector(residual_displacement_,   // G4Transform3D is a typedef of HepGeom::Transform3D, whose .cc
                                                                        0.,                       // source code is here:
                                                                        0.)));                    // apc.u-paris.fr/~franco/g4doxy4.10/html/_transform3_d_8cc_source.html#l00142
                                                                                                  // I had to add a residual displacement, because otherwise, apparently the 
                                                                                                  // affine transformation has no inverse (i.e. its inverse transformation has 
                                                                                                  // zero determinant) and the error at line 157 of Transform3D.cc pops up. 
                                                                                                  // I.e. an affine transformation which consists just of a 180º rotation about 
                                                                                                  // the Y axis apparently has no inverse.
      multiunion_geometry_solid->Voxelize();                                                                                  
      geometry_solid = dynamic_cast<G4VSolid*>(multiunion_geometry_solid);
    }
    else{
      geometry_solid = dynamic_cast<G4VSolid*>(whole_rod_solid);
    }

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    //pvt->SetMaterialPropertiesTable(opticalprops::EJ286(secondary_wls_attlength_));
    pvt->SetMaterialPropertiesTable(opticalprops::G2P_FB118(cromophore_concentration_));

    G4LogicalVolume* rod_logic = new G4LogicalVolume( geometry_solid, 
                                                      pvt, 
                                                      rod_name);
    rod_logic->SetUserLimits(ul_);

    G4double running_z_pos =  (-1.*internal_width_/2.)
                              +rods_lateral_gap_
                              +(rod_width_/2.);

    G4double z_step = rods_lateral_gap_+rod_width_;
    G4int k = 0;
    G4ThreeVector pos(0., 0., 0.);
    G4RotationMatrix* rot = new G4RotationMatrix();    
    rot->rotateY(0.0*deg);
    for(G4int i=0; i<rods_no_; i++){
        pos.setZ(running_z_pos);
        new G4PVPlacement(  rot, pos, rod_name, rod_logic, 
                            mother_physical, true, k, true);
        running_z_pos += z_step;
        k++;
    }
    return;

    // !!!! HAVE YOU USED residual_displacement_ TO ACCOMODATE gap_, internal_length_ AND SO ON? !!!! 
    // YOU NEED TO IMPLEMENT THE TUNNELING PROBABILITY YET HERE !!!!
  }


  void XArapucaRods::ConstructBoards(G4VPhysicalVolume* mother_physical) const
  {

    G4double sipm_protrusion = 0.1*mm;  // According to the add_blocks_between_sipms_ attribute documentation 
                                        // in XArapucaRods.h. Check such documentation for more information.
    SiPMBoard board1;
    board1.SetBaseID(0);
    board1.SetSiPMCode(SiPM_code_);
    board1.SetNumPhsensors(num_phsensors_);
    board1.SetReflectiveSupports(ref_phsensors_supports_);
    board1.SetAddBlocks(add_blocks_between_sipms_);
    board1.SetSiPMProtrusion(sipm_protrusion);

    SiPMBoard board2;
    board2.SetBaseID(num_phsensors_);
    board2.SetSiPMCode(SiPM_code_);
    board2.SetNumPhsensors(num_phsensors_);
    board2.SetReflectiveSupports(ref_phsensors_supports_);
    board2.SetAddBlocks(add_blocks_between_sipms_);
    board2.SetSiPMProtrusion(sipm_protrusion);

    board1.SetBoardLength(internal_width_-(rods_lateral_gap_*2)); // Equivalent to ((rod_width_*rods_no_)+(rods_lateral_gap_*(rods_no_-1)))
    board1.Construct();
    G4LogicalVolume* board1_logic_vol = board1.GetLogicalVolume();
        
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateY(-90.*deg);

    G4double x_pos = (rod_length_/2.) + (board1.GetOverallThickness()/2.) +rod_sipm_gap_;
    G4ThreeVector pos(x_pos, 0., 0.);
    new G4PVPlacement(rot, 
                      pos, 
                      "SIPMS_BOARD", 
                      board1_logic_vol, 
                      mother_physical, 
                      true, 0, true);

    board2.SetBoardLength(internal_width_-(rods_lateral_gap_*2)); // Equivalent to ((rod_width_*rods_no_)+(rods_lateral_gap_*(rods_no_-1)))
    board2.Construct();
    G4LogicalVolume* board2_logic_vol = board2.GetLogicalVolume();

    G4RotationMatrix* rot_2 = new G4RotationMatrix();
    rot_2->rotateY(+90.*deg);

    x_pos = (rod_length_/2.) + (board2.GetOverallThickness()/2.) +rod_sipm_gap_;
    G4ThreeVector pos_2(-1.*x_pos, 0., 0.);
    new G4PVPlacement(rot_2, 
                      pos_2, 
                      "SIPMS_BOARD", 
                      board2_logic_vol, 
                      mother_physical, 
                      true, 0, true);
    return;
  }


  void XArapucaRods::ConstructReflectiveCase(G4VPhysicalVolume* mother_physical) const
  {
    
    // Based on arxiv.org/abs/1912.09191 and TDR vol.IX, the reflective foils
    // fully cover the internal cavity of the X-Arapuca.

    const G4String ref_case_name = "REF_CASE";

    // Get reflective case volume as subtraction solid from two boxes
    G4Box* aux_outter_box = new G4Box("AUX_OUTTER_BOX", 
                                      overall_length_/2., 
                                      internal_thickn_/2., 
                                      overall_width_/2.);  

    // Extra thickness to prevent boolean subtraction of solids with matching surfaces
    // See geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html#solids-made-by-boolean-operations
    G4double tolerance = 1.*cm;

    G4Box* aux_inner_box =  new G4Box("AUX_INNER_BOX", 
                                      internal_length_/2., 
                                      (internal_thickn_+tolerance)/2., 
                                      internal_width_/2.);

    G4SubtractionSolid* ref_case_solid =    new G4SubtractionSolid( ref_case_name, 
                                                                    aux_outter_box, 
                                                                    aux_inner_box, 
                                                                    nullptr, 
                                                                    G4ThreeVector(0., 0., 0.));

    // The reflective case is 'open' on both sides. If the X-ARAPUCA is not 
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

    refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
    new G4LogicalSkinSurface("REF_CASE_SURFACE", ref_case_logic, refsurf_opsurf);   
    
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.),
          ref_case_name, ref_case_logic, 
          mother_physical,
          false, 0, true);

    return;

  }
  
 
  void XArapucaRods::ConstructDichroicAssemblies(G4VPhysicalVolume* mother_physical) const
  {
      
    // ------------------------------------ DICHROIC FILTER MODEL ------------------------------------
    //
    // Dichroic-filters assembly is made up of one-piece frame and a bunch of dichroic filters 
    // inserted in the frame. The DF (without WLS coating) is implemented as two stacked volumes. 
    // The volume that faces the cryostat is meant to be the DF substrate (a piece of glass 
    // with good transmitance for the coating light). The volume that faces the intenal cavity 
    // of the X-Arapuca is meant to be the multilayer structure (MLS). The DF transmitance curve 
    // is implemented as a properly-placed G4LogicalBorderSurface's. When the DF are coated, a 
    // third volume is stacked on top of the filter (on top of the substrate volume).
    // Here, coating is a volume with the coating optical properties, such as WLS absorption length,
    // refractive index and the following rough surfaces: LAr->coating, coating->LAr and 
    // coating->substrate. substrate is a volume with the DF substrate optical properties, such 
    // as absorption length and refractive index. Since both the coating and the substrate have 
    // well-defined refractive indices, this model is able to simulate the light behaviour in the 
    // coating-substrate interface, which has not been measured in the laboratory. The DF transmission
    // curve (TC) is implemented in the substrate-MLS interface. The resulting volume is:
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
    //                                |                               |
    //                                |----G4LogicalBorderSurfaces----|
    //                                |                               |
    //                                |  MLS (with substrate rindex)  |
    //                       _________|_______________________________|_________
    //                                                
    //                                 LAr (X-ARAPUCA internal cavity)
    //
    // For a common configuration, the coating is PTP with refractive index equal to 1.65, and
    // the substrate is fused silica (R. index approximately equal to 1.47 for PTP light). The
    // transmitance curve of the DF must be "seen" not only by photons that try to enter the 
    // X-Arapuca coming from the cryostat, but also by those that are already trapped within the 
    // X-Arapuca. Such transmitance curve must be seen either in the substrate-MLS interface, 
    // or in the MLS-LAr interface.   
    //
    // For photons trying to enter the X-Arapuca for the first time, we want them to undergo a 
    // DF->LAr refraction for the sake of realism (See **). So, either we set the transmission curve 
    // in the substrate-MLS interface and let the DF-LAr refraction happen in the MLS-LAr interface, 
    // or we set the transmission curve in the MLS-LAr interface and let the DF-LAr refraction happen 
    // in the substrate-MLS interface. A priori, both options are valid, but if we go for the second one, 
    // where the DF->LAr refraction happens in the substrate-MLS interface, we need to set the MLS 
    // G4MaterialPropertiesTable as if it was LAr. For the sake of keeping the DF source code 
    // independent of the outter media, which, in this case is LAr but may change for other-detectors 
    // simulations, I am going to implement the second option, i.e.:
    //
    //                          1) the photons see the transmitance curve in the substrate-MLS interface
    //                          2) the DF->LAr refraction happens in the MLS-LAr interface
    //
    // Now, up to Geant4 functioning, in the interface where such transmission curve is implemented, 
    // no refraction is simulated. The photon is either reflected or transmitted with no deviation 
    // with respect to its original trajectory. This must be kept in mind to set the refractive 
    // indices of the volumes so that the Snell's invariant (*) is not destroyed, i.e. refractive 
    // indices of the two volumes whose interface is endowed with the TC, must be the same.
    // In our case, this means that the MLS volume must be endowed with the substrate refractive index.
    //
    // There's another thing to take into account. The PTP-substrate and the MLS-LAr interfaces 
    // might transmit or reflect photons up to Fresnel probabilities based on the AOI and the 
    // refractive indices of the surrounding media. However, all of the reflection-or-transmission 
    // information should be already taken into account in the transmission curve. Indeed, when 
    // we measure a filter in the laboratory, the transmitance (or reflectance) that this "external"
    // interfaces add are convoluted with the DF TC itself. Endowing the substrate-MLS interface
    // with the measured TC is problematic because, not only would we be accounting for the 
    // external interfaces twice (one of them in the simulation, because Geant4 will simulate
    // the transmitance in those surfaces using Fresnel equations, and the other one in the
    // TC we measured, which already take those into account) but also, in the simulation we may
    // not want to simulate air->DF->air external interfaces, but maybe PTP->DF->LAr. A way to 
    // solve this is to de-convolute the TCs we measure in the lab. I.e., assume that the curves
    // we measure, MTC=MTC(wavelength, AOI)=MTC(wl, AOI), is the product of the Fresnel transmitance
    // from the lab media to the DF, FT1=FT1(AOI), the intrinsic transmission curve of the 
    // dichroic filter (that of the MLS), TC=TC(wl, AOI) and the Fresnel transmitance from the DF
    // to the lab media, FT2=FT2(AOI):
    //
    //  MTC(wl, AOI)=FT1(AOI)*TC(wl, AOI)*FT2(AOI)
    //
    // And so, we can compute the intrinsic TC 
    //
    // TC(wl, AOI)=MTC(wl, AOI)/(FT1(AOI)*FT2(AOI))
    //
    // And endow the substrate-MLS interface with such transmission curve. Overall, with the
    // previously specified DF-model and the intrinsic TC, Geant4 will be able to realistically 
    // simulate the coated DF immersed in LAr, because photons interacting with it will (rightfully)
    // undergo the coating->DF and the DF->LAr interfaces.
    // 
    //
    // * For a light beam traversing a parallel-faces stacking, the product n*sin(theta) is a constant
    // which is preserved for every volume, i.e. n_i*sin(theta_i)=K \forall i. If, at some interface,
    // the light is allowed to traverse from one volume to another one (with different refractive)
    // with no refraction, the snell's invariant will be destroyed. Since this is not physically
    // realistic, we should prevent this from happening. A way to do this so is to set the same 
    // refractive index for the volumes whose interface is endowed with the TC.
    //
    // ** In the lab, when we measure the transmitance curve (TC) for a DF at a certain AOI, the DF is
    // not coated. This means that the light impinges into the DF from air, and after interacting with
    // it, the light emerges back to air. Since the filter is a stacking of face-parallel layers, 
    // the transmitted light beam is parallel to the incident light beam. However, for a coated DF, 
    // the incident light beam is (hopefully) WLSed by the coating, which, in turn, emits the light
    // isotropic-ly. Thus, the light that impinges into the DF not only has a random angle within the
    // coating, but also such angle won't match the angle of the photon in case it is transmitted
    // through the DF into the surrounding media, say LAr. Such final angle should be computed applying
    // snell's law towards LAr, from the coating volume or any other volume where the photon still
    // preserves the snell's invariant ( {n, \theta}, for face-parallel geometries) of the coating-emission.
    // -----------------------------------------------------------------------------------------------

    G4double df_MLS_thickn = DF_thickn_-DF_substrate_thickn_;

    G4Box* cover_solid = new G4Box("AUX", DFA_length_/2., DFA_thickn_/2., DFA_width_/2.);
    G4double tolerance = 5.*mm; // Tolerance to prevent matching surfaces when subtracting carvings from frame
    
    G4Box* carving_solid                =  new G4Box("FRAME_CARVING",               DF_length_/2., DFA_thickn_/2. +tolerance,   DF_width_/2.); 
    G4Box* df_substrate_solid           =  new G4Box("DICHROIC_FILTER_SUBSTRATE",   DF_length_/2., DF_substrate_thickn_/2.,     DF_width_/2.); 
    G4Box* df_MLS_solid                 =  new G4Box("DICHROIC_FILTER_MLS",         DF_length_/2., df_MLS_thickn/2.,            DF_width_/2.); 
    G4Box* coating_solid                =  nullptr;
    if(DF_are_coated_) coating_solid    =  new G4Box("DF_COATING", DF_length_/2., coating_thickn_/2., DF_width_/2.); 

    G4MultiUnion* frame_carvings_multiunion_solid   = new G4MultiUnion("FRAME_CARVINGS");
    G4MultiUnion* df_substrates_multiunion_solid    = new G4MultiUnion("DICHROIC_FILTER_SUBSTRATES");
    G4MultiUnion* df_MLSs_multiunion_solid          = new G4MultiUnion("DICHROIC_FILTER_MLSs");
    G4MultiUnion* coatings_multiunion_solid         = nullptr;
    if(DF_are_coated_) coatings_multiunion_solid    = new G4MultiUnion("DF_COATINGS"); 

    G4double x_pos, z_pos;
    G4Transform3D* transform_ptr = nullptr;
    G4RotationMatrix* rot = new G4RotationMatrix();

    for(G4int i=0; i<df_no_along_wlsplength_; i++){
        for(G4int j=0; j<df_no_along_wlspwidth_; j++){
            x_pos = (-1.*(DFA_length_/2.))  +outter_frame_width_along_wlsplength_   +(i*inner_frames_width_along_wlsplength_)   +((i+0.5)*DF_length_);
            z_pos = (-1.*(DFA_width_/2.))   +outter_frame_width_along_wlspwidth_    +(j*inner_frames_width_along_wlspwidth_)    +((j+0.5)*DF_width_);
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector(x_pos, 0., z_pos));

            frame_carvings_multiunion_solid                 ->AddNode(*carving_solid,       *transform_ptr);
            df_substrates_multiunion_solid                  ->AddNode(*df_substrate_solid,  *transform_ptr);
            df_MLSs_multiunion_solid                        ->AddNode(*df_MLS_solid,        *transform_ptr);
            if(DF_are_coated_) coatings_multiunion_solid    ->AddNode(*coating_solid,       *transform_ptr);
        }
    }
    
    frame_carvings_multiunion_solid                 ->Voxelize();
    df_substrates_multiunion_solid                  ->Voxelize();
    df_MLSs_multiunion_solid                        ->Voxelize();
    if(DF_are_coated_) coatings_multiunion_solid    ->Voxelize();

    // FRAME //
    if(!remove_DFA_frame_){
        const G4String frame_name = "DICHROIC-FILTERS_ASSEMBLY_FRAME";
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

        if(DFA_frame_is_reflective_){
            const G4String refcoat_name = "REF_COATING";
            G4OpticalSurface* refcoat_opsurf = 
            new G4OpticalSurface(refcoat_name, unified, ground, dielectric_metal, 1);
            if(DFA_frame_is_specular_)  refcoat_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());  // Fix this !!!!
            else                        refcoat_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
            new G4LogicalSkinSurface(refcoat_name, frame_logic, refcoat_opsurf);
        }
        
        G4PVPlacement* frame_physical = new G4PVPlacement(nullptr, G4ThreeVector(0., (internal_thickn_+DFA_thickn_)/2. +0.*cm, 0.),
                                        frame_name, frame_logic, mother_physical, true, 0, true);

        if(double_sided_){
            new G4PVPlacement(nullptr, G4ThreeVector(0., -1.*(internal_thickn_+DFA_thickn_)/2., 0.), 
                                frame_name, frame_logic, mother_physical, true, 1, true);
        }
        else{
            // If the X-ARAPUCA is not double-sided, then the dummy plane is reflective
            const G4String refsurf_name = "REF_SURFACE";
            G4OpticalSurface* refsurf_opsurf = 
            new G4OpticalSurface(refsurf_name, unified, ground, dielectric_metal, 1);
            refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());

            const G4String cover_name = "REFLECTIVE_COVER";
            G4LogicalVolume* cover_logic = new G4LogicalVolume(cover_solid, materials::FR4(), cover_name);
            cover_logic->SetVisAttributes(frame_col);
            new G4LogicalSkinSurface(refsurf_name, cover_logic, refsurf_opsurf);   
            new G4PVPlacement(nullptr, G4ThreeVector(0., -1.*(internal_thickn_+DFA_thickn_)/2., 0.), 
                                cover_name, cover_logic, mother_physical, true, 1, true);
        }
    }

    // DICHROIC FILTERS //
    if(!remove_DFs_){
        G4Material* DF_substrate_mat    = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
        DF_substrate_mat->SetMaterialPropertiesTable(DF_substrate_mpt_);
        G4LogicalVolume* df_substrates_multiunion_logic = 
                            new G4LogicalVolume(df_substrates_multiunion_solid, DF_substrate_mat, "DICHROIC_FILTER_SUBSTRATES");
        G4LogicalVolume* df_MLSs_multiunion_logic = 
                            new G4LogicalVolume(df_MLSs_multiunion_solid,       DF_substrate_mat, "DICHROIC_FILTER_MLSs");
            
        // Visually transparent substrate, and visually red MLs
        G4VisAttributes df_MLSs_col = nexus::BloodRedAlpha();
        df_MLSs_col.SetForceSolid(true);
        df_MLSs_multiunion_logic->SetVisAttributes(df_MLSs_col);

        // Place the filters
        G4VPhysicalVolume* first_df_substrates_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
            new G4PVPlacement(nullptr, G4ThreeVector(0., (internal_thickn_/2.) 
                                                        +df_MLS_thickn 
                                                        +(DF_substrate_thickn_/2.) 
                                                        +(DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_)), 0.), 
                            "DICHROIC_FILTER_SUBSTRATES_1", df_substrates_multiunion_logic, mother_physical, true, 0, true));

        G4VPhysicalVolume* first_df_MLSs_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
            new G4PVPlacement(nullptr, G4ThreeVector(0.,  (internal_thickn_/2.)
                                                          +(df_MLS_thickn/2.)
                                                          +(DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_)), 0.),
                              "DICHROIC_FILTER_MLSs_1", df_MLSs_multiunion_logic, mother_physical, true, 0, true));
        
        // Check that there's dichroic information for ingoing (wrt XA) photons
        if(path_to_inwards_dichroic_data_==""){
            G4Exception("[XArapucaRods]", "ConstructDichroicAssemblies()",
                        FatalException, "The path to the inwards dichroic data file was not set.");
        }

        // Check that there's dichroic information for outgoing photons
        if(path_to_outwards_dichroic_data_==""){
            G4Exception("[XArapucaRods]", "ConstructDichroicAssemblies()",
                        FatalException, "The path to the outwards dichroic data file was not set.");
        }

        // Construct the ingoing optical surface
        setenv("G4DICHROICDATA", path_to_inwards_dichroic_data_, 1);
        G4OpticalSurface* dfs_inwards_opsurf =                // G4OpticalSurface constructor loads the
            new G4OpticalSurface( "DICHROIC_INWARDS_OPSURF",  // dichroic information from the file which
                                  dichroic,                   // is currently pointed to by the environment
                                  polished,                   // variable G4DICHROICDATA
                                  dielectric_dichroic);

        // Construct the outgoung optical surface
        setenv("G4DICHROICDATA", path_to_outwards_dichroic_data_, 1);
        G4OpticalSurface* dfs_outwards_opsurf =   
            new G4OpticalSurface( "DICHROIC_OUTWARDS_OPSURF", 
                                  dichroic, 
                                  polished, 
                                  dielectric_dichroic);

        // Endow the substrate->MLS surface with the ingoing optical surface
        new G4LogicalBorderSurface( "Substrate1->MLS1", 
                                    first_df_substrates_multiunion_physical, 
                                    first_df_MLSs_multiunion_physical, 
                                    dfs_inwards_opsurf);

        // Endow the MLS->substrate surface with the outgoing optical surface
        new G4LogicalBorderSurface( "MLS1->Substrate1", 
                                    first_df_MLSs_multiunion_physical, 
                                    first_df_substrates_multiunion_physical, 
                                    dfs_outwards_opsurf);
            
        G4VPhysicalVolume* second_df_substrates_multiunion_physical = nullptr;
        G4VPhysicalVolume* second_df_MLSs_multiunion_physical = nullptr;
        if(double_sided_){
            second_df_substrates_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                new G4PVPlacement(nullptr, G4ThreeVector(0.,    -1.*(
                                                                (internal_thickn_/2.) 
                                                                +df_MLS_thickn 
                                                                +(DF_substrate_thickn_/2.) 
                                                                +(DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_))), 0.), 
                                "DICHROIC_FILTER_SUBSTRATES_2", df_substrates_multiunion_logic, mother_physical, true, 1, true));

            second_df_MLSs_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                new G4PVPlacement(nullptr, G4ThreeVector(0.,  -1.*(
                                                              (internal_thickn_/2.)
                                                              +(df_MLS_thickn/2.)
                                                              +(DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_))), 0.),
                                  "DICHROIC_FILTER_MLSs_2", df_MLSs_multiunion_logic, mother_physical, true, 0, true));

            new G4LogicalBorderSurface( "Substrate2->MLS2", 
                                        second_df_substrates_multiunion_physical, 
                                        second_df_MLSs_multiunion_physical, 
                                        dfs_inwards_opsurf);

            new G4LogicalBorderSurface( "MLS2->Substrate2", 
                                        second_df_MLSs_multiunion_physical, 
                                        second_df_substrates_multiunion_physical, 
                                        dfs_outwards_opsurf);
        }

        // COATINGS //
        if(DF_are_coated_){
            G4Material* coatings_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
            coatings_mat->SetMaterialPropertiesTable(opticalprops::PTP(coating_rindex_));
            G4LogicalVolume* coatings_multiunion_logic = 
                                new G4LogicalVolume(coatings_multiunion_solid, coatings_mat, "DF_COATINGS");   

            G4VisAttributes coatings_col = nexus::TitaniumGreyAlpha();
            //coatings_col.SetForceSolid(true);
            coatings_multiunion_logic->SetVisAttributes(coatings_col);

            // Place the coatings
            G4VPhysicalVolume* first_coatings_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                new G4PVPlacement(  nullptr, G4ThreeVector(0.,  (internal_thickn_/2.)
                                                                +DF_thickn_
                                                                +(DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_)) 
                                                                +(coating_thickn_/2.), 0.), 
                                    "DICHROIC_FILTER_COATINGS_1", coatings_multiunion_logic, mother_physical, true, 0, true));

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
            // We will also add roughness for the coating->substrate interface, but only with such ordering. The alternative case takes place
            // when the photon travels from the DF glass substrate to the coating. The glass is supposed to be polished, so the photon may not
            // see a rough surface.
            new G4LogicalBorderSurface( "COATING1->SUBSTRATE1", 
                                        first_coatings_multiunion_physical, 
                                        first_df_substrates_multiunion_physical, 
                                        coating_rough_surf);

            if(double_sided_){
                G4VPhysicalVolume* second_coatings_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                    new G4PVPlacement(  nullptr, G4ThreeVector(0.,  -1.*((internal_thickn_/2.)
                                                                    +DF_thickn_
                                                                    +(DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_)) 
                                                                    +(coating_thickn_/2.)), 0.), 
                                        "DICHROIC_FILTER_COATINGS_2", coatings_multiunion_logic, mother_physical, true, 1, true));
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
                                            second_df_substrates_multiunion_physical, 
                                            coating_rough_surf);
            }
        }
    }
    return;
  }


  void XArapucaRods::ConstructReflectiveSeparators(G4VPhysicalVolume* mother_physical) const
  {

    const G4String separator_name = "SEPARATOR";

    SiPMMPPC * sipm_ptr = ChooseSiPM(SiPM_code_);
    G4double separator_length = rod_length_ +(2.*sipm_ptr->GetThickness());

    G4Box* separator_solid = new G4Box(separator_name,  separator_length/2., 
                                                        rod_height_/2., 
                                                        rods_lateral_gap_/4.); 

    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");

    G4LogicalVolume* separator_logic = new G4LogicalVolume( separator_solid, 
                                                            mat, 
                                                            separator_name);
    //Now create the reflective optical surface
    const G4String refsurf_name = "REF_SURFACE";
    G4OpticalSurface* refsurf_opsurf = 
      new G4OpticalSurface(refsurf_name, unified, ground, dielectric_metal, 1);

    refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
    new G4LogicalSkinSurface("REF_SEPARATOR_SURFACE", separator_logic, refsurf_opsurf);   

    G4double running_z_pos =  (-1.*internal_width_/2.)
                              +rods_lateral_gap_        // Separators are placed
                              +rod_width_               // in between every two
                              +(rods_lateral_gap_/2.);  // adjacent rods
    G4double z_step = rods_lateral_gap_+rod_width_;
    G4int k = 0;
    G4ThreeVector pos(0., 0., 0.);
    G4RotationMatrix* rot = new G4RotationMatrix();    
    rot->rotateY(0.0*deg);
    for(G4int i=0; i<(rods_no_-1); i++){    // There are rods_no_-1 separators
        pos.setZ(running_z_pos);
        new G4PVPlacement(  rot, pos, separator_name, separator_logic, 
                            mother_physical, true, k, true);
        running_z_pos += z_step;
        k++;
    }

    if(cut_rods_ && place_separator_at_the_cut_){ // In this case, XArapucaRods::ConstructReflectiveSeparators 
                                                  // is also responsible for placing the cut separators

      const G4String cut_separator_name = "CUT_SEPARATOR";
      G4Para* cut_separator_solid = new G4Para( cut_separator_name,
                                                (cut_thickness_/2.)/2.,
                                                rod_height_/2., 
                                                (rod_width_/2.) +(rods_lateral_gap_/4.),  // So that the cut-separators fit tightly
                                                                                          // between two rods-separators
                                                0.0,
                                                cut_angle_,                   // Parallelepiped angle with respect to the
                                                                              // WLS-rod side which is rod_width_ long
                                                0.0);
      G4LogicalVolume* cut_separator_logic = 
        new G4LogicalVolume(cut_separator_solid, materials::FR4(), cut_separator_name);

      G4OpticalSurface* opsurf = 
        new G4OpticalSurface("CUT_SEPARATOR_REF_SURFACE", unified, ground, dielectric_metal, 1);
    
      opsurf->SetMaterialPropertiesTable(opticalprops::Vikuiti());
      new G4LogicalSkinSurface("CUT_SEPARATOR_REF_SURFACE", cut_separator_logic, opsurf);

      running_z_pos = (-1.*internal_width_/2.)
                      +rods_lateral_gap_
                      +(rod_width_/2.);

      z_step = rods_lateral_gap_+rod_width_;
      G4int l = 0;
      G4ThreeVector pos_2(0., 0., 0.);
      for(G4int i=0; i<rods_no_; i++){
          pos_2.setZ(running_z_pos);
          new G4PVPlacement(  nullptr, pos_2, cut_separator_name, cut_separator_logic, 
                              mother_physical, true, l, true);
          running_z_pos += z_step;
          l++;
      }
    }

    return;
  }

  G4bool XArapucaRods::geometry_is_ill_formed()
  {

    if(DF_substrate_thickn_>=DF_thickn_){
        G4Exception("[XArapucaRods]", "geometry_is_ill_formed()", FatalException,
        "The dichroic filters substrate thickness cannot be bigger or equal to the overall dichroic filter thickness.");
    }

    if(DF_thickn_>DFA_thickn_){
        G4Exception("[XArapucaRods]", "geometry_is_ill_formed()", FatalException,
        "The dichroic filters thickness cannot be bigger than the frame thickness.");
    }

    if(internal_length_<internal_width_){
        G4Exception("[XArapucaRods]", "geometry_is_ill_formed()", FatalException,
        "For consistency with the code, internal_length_ must be greater than internal_width_.");
    }

    return false;
  }


  G4ThreeVector XArapucaRods::GenerateVertex(const G4String&) const{

    G4double tolerance = 0.1*mm;    // Small distance over the dichroic filter assembly 
                                    // (DFA) from which photons are launched. Also, the 
                                    // width of the outter border projected over the DF 
                                    // from which photons won't be launched (Just see the 
                                    // implementation in x_pos and z_pos below to understand 
                                    // its meaning)

    G4double x_pos, z_pos;
    G4double y_pos = (internal_thickn_/2.) +DFA_thickn_ +tolerance;

    if(generation_region_=="dichroic"){
      std::random_device rd;  //Will be used to obtain a seed for the random number engine
      std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
      std::uniform_int_distribution<> dist(0, df_no_along_wlsplength_*df_no_along_wlspwidth_ -1);
      G4int filter_no = dist(gen);

      G4int column_no = filter_no%df_no_along_wlsplength_;                // \in[0, df_no_along_wlsplength -1]  
      G4int row_no = (int) std::floor(filter_no/df_no_along_wlsplength_); // \in[0, df_no_along_wlspwidth -1]

      G4double selected_filter_x_center =  (-1.*(DFA_length_/2.)) +outter_frame_width_along_wlsplength_   +(column_no*inner_frames_width_along_wlsplength_)   +((column_no+0.5)*DF_length_);
      G4double selected_filter_z_center =  (-1.*(DFA_width_/2.))  +outter_frame_width_along_wlspwidth_    +(row_no*inner_frames_width_along_wlspwidth_)       +((row_no+0.5)*DF_width_);

      x_pos = UniformRandomInRange(  selected_filter_x_center -(DF_length_/2.) +tolerance, 
                                              selected_filter_x_center +(DF_length_/2.) -tolerance    );
      z_pos = UniformRandomInRange(  selected_filter_z_center -(DF_width_/2.) +tolerance, 
                                              selected_filter_z_center +(DF_width_/2.) -tolerance     );
    }
    else if(generation_region_=="custom"){
      G4double random_radius =  UniformRandomInRange(gen_diameter_/2., 0.);
      G4double random_angle =   UniformRandomInRange(twopi, 0.); 
      x_pos = gen_x_ +(random_radius*sin(random_angle));
      z_pos = gen_z_ +(random_radius*cos(random_angle));
    }
    else{ // Default behaviour is that of generation_region_=="random"
      x_pos = UniformRandomInRange( DFA_length_/2.,
                                    -1.*DFA_length_/2.);
      z_pos = UniformRandomInRange( DFA_width_/2.,
                                    -1.*DFA_width_/2.);
    }
    return G4ThreeVector(x_pos, y_pos, z_pos);
  }


  SiPMMPPC * XArapucaRods::ChooseSiPM(G4int input) const{
    SiPMMPPC* sipm_ptr = nullptr;
    if(input==1){
      sipm_ptr = new HamamatsuS133606050VE();
    }
    else if(input==2){
      sipm_ptr = new HamamatsuS133605075HQR();
    }
    else if(input==3){
      sipm_ptr = new FbkNuvHdCryoTT();
    }
    else if(input==4){
      sipm_ptr = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm_ptr = new PerfectSiPMMPPC();
    }
    return sipm_ptr;
  }
} //End namespace nexus
