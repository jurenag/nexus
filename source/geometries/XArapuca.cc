#include "XArapuca.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "WLSPlate.h"
#include "HamamatsuS133606050VE.h"
#include "ScalableHamamatsuS133606050VE.h"
#include "SiPMBoard.h"
#include "MomentumSD.h"
#include "RandomUtils.h"
#include "Visibilities.h"

#include <algorithm>
#include <random>
#include <cmath>
#include <G4GenericMessenger.hh>
#include <G4UserLimits.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
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

REGISTER_CLASS(XArapuca, GeometryBase)

namespace nexus{

  XArapuca::XArapuca():
  GeometryBase(), 
  ///Get internal reflector cavity dimensions from arxiv.org/abs/1912.09191
  config_code_                          (1          ),
  internal_length_                      (488.   *mm ),   ///X
  internal_width_                       (100.   *mm ),   ///Z
  internal_thickn_                      (8      *mm ),   ///Y
  DFA_thickn_                           (1.     *mm ),
  DF_thickn_                            (1.     *mm ),
  DF_pos_wrt_DFA_pos_                   (0.         ),
  DF_are_coated_                        (true       ),
  coating_thickn_                       (3.226  *um ),  // Based on arxiv.org/abs/1912.09191 and TDR vol.IX, section 5.8.3.1,
                                                        // the pTP film thickness is such that there's 400 micrograms of pTP
                                                        // deposited over each square centimeter of DF.
                                                        // This, together with the pTP density, (1.24g/cm3, found in
                                                        // en.wikipedia.org/wiki/Terphenyl), gives a pTP film thickness of 
                                                        // 3.226 micrometers
  outter_frame_width_along_wlsplength_  (10.    *mm ),
  outter_frame_width_along_wlspwidth_   (10.    *mm ),
  inner_frames_width_along_wlsplength_  (6.     *mm ),
  inner_frames_width_along_wlspwidth_   (2.     *mm ),
  df_no_along_wlsplength_               (2          ),
  df_no_along_wlspwidth_                (3          ),
  DFA_frame_is_reflective_              (false      ),
  DFA_frame_is_specular_                (true       ),
  remove_DFA_                           (false      ),
  remove_DFs_                           (false      ),    
  case_thickn_                          (1.     *mm ),   ///Get foil thickness from isoltronic.ch/assets/of-m-vikuiti-esr-app-guide.pdf
  PS_config_code_                       (1          ),
  num_phsensors_                        (24         ),
  gap_                                  (0.5    *mm ),
  ref_phsensors_supports_               (true       ), 
  double_sided_                         (true       ),
  collectors_are_reflective_            (false      ),
  generation_vertex_over_df_            (true       ),
  path_to_dichroic_data_                (""         ),
  world_extra_thickn_                   (100.   *cm ),
  plate_length_                         (487.   *mm ),   ///X
  plate_width_                          (93.    *mm ),   ///Z
  plate_thickn_                         (3.5    *mm ),   ///Y
  only_sipms_along_long_sides_          (true       ),
  with_boards_                          (true       ),
  with_dimples_                         (true       ),
  dimple_type_                          ("cylindrical"),
  flat_dimple_width_                    (6.1    *mm ),
  flat_dimple_depth_                    (2.     *mm ),
  curvy_dimple_radius_                  (1.5    *mm ),
  fibers_no_                            (6          ),
  fiber_planes_no_                      (1          ),    
  fiber_radius_                         (5.     *mm ),
  fiber_length_                         (400.   *mm ),
  along_long_side_                      (true       )
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/XArapuca/",
				"Control commands of geometry XArapuca.");

    G4GenericMessenger::Command& cc_cmd =
      msg_->DeclareProperty("config_code", config_code_,
			    "Configuration code.");
    cc_cmd.SetParameterName("config_code", false);
    cc_cmd.SetRange("config_code>=1"); 
    
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

    G4GenericMessenger::Command& dfat_cmd =
      msg_->DeclareProperty("DFA_thickn", DFA_thickn_,
			    "Thickness of the dichroic filters assembly.");
    dfat_cmd.SetUnitCategory("Length");
    dfat_cmd.SetParameterName("DFA_thickn", false);
    dfat_cmd.SetRange("DFA_thickn>0.");

    G4GenericMessenger::Command& dft_cmd =
      msg_->DeclareProperty("DF_thickn", DF_thickn_,
			    "Thickness of the dichroic filters.");
    dft_cmd.SetUnitCategory("Length");
    dft_cmd.SetParameterName("DF_thickn", false);
    dft_cmd.SetRange("DF_thickn>0.");

    G4GenericMessenger::Command& dpwdp_cmd =
      msg_->DeclareProperty("DF_pos_wrt_DFA_pos", DF_pos_wrt_DFA_pos_,
			    "Position (height) of the dichroic filters with respect to the DFA position (height). This parameter can take values from 0 to 1.");
    dpwdp_cmd.SetParameterName("DF_pos_wrt_DFA_pos", false);
    dpwdp_cmd.SetRange("DF_pos_wrt_DFA_pos>=0.");
    dpwdp_cmd.SetRange("DF_pos_wrt_DFA_pos<=1.");

    G4GenericMessenger::Command& dfac_cmd =
      msg_->DeclareProperty("DF_are_coated", DF_are_coated_,
			    "Whether the dichroic filter is set on both sides of the WLS plate.");

    G4GenericMessenger::Command& ptpct_cmd =
      msg_->DeclareProperty("coating_thickn", coating_thickn_,
			    "Thickness of the coating layer that is deposited over the dichroic filters.");
    ptpct_cmd.SetUnitCategory("Length");
    ptpct_cmd.SetParameterName("coating_thickn", false);
    ptpct_cmd.SetRange("coating_thickn>0.");

    G4GenericMessenger::Command& ofwawl_cmd =
      msg_->DeclareProperty("outter_frame_width_along_wlsplength", outter_frame_width_along_wlsplength_,
			    "Width of the outter DFA frame along the WLS plate length.");
    ofwawl_cmd.SetUnitCategory("Length");
    ofwawl_cmd.SetParameterName("outter_frame_width_along_wlsplength", false);
    ofwawl_cmd.SetRange("outter_frame_width_along_wlsplength>0.");

    G4GenericMessenger::Command& ofwaww_cmd =
      msg_->DeclareProperty("outter_frame_width_along_wlspwidth", outter_frame_width_along_wlspwidth_,
			    "Width of the outter DFA frame along the WLS plate width.");
    ofwaww_cmd.SetUnitCategory("Length");
    ofwaww_cmd.SetParameterName("outter_frame_width_along_wlspwidth", false);
    ofwaww_cmd.SetRange("outter_frame_width_along_wlspwidth>0.");

    G4GenericMessenger::Command& ifwawl_cmd =
      msg_->DeclareProperty("inner_frames_width_along_wlsplength", inner_frames_width_along_wlsplength_,
			    "Width of the inner DFA frames along the WLS plate length.");
    ifwawl_cmd.SetUnitCategory("Length");
    ifwawl_cmd.SetParameterName("inner_frames_width_along_wlsplength", false);
    ifwawl_cmd.SetRange("inner_frames_width_along_wlsplength>0.");

    G4GenericMessenger::Command& ifwaww_cmd =
      msg_->DeclareProperty("inner_frames_width_along_wlspwidth", inner_frames_width_along_wlspwidth_,
			    "Width of the inner DFA frames along the WLS plate width.");
    ifwaww_cmd.SetUnitCategory("Length");
    ifwaww_cmd.SetParameterName("inner_frames_width_along_wlspwidth", false);
    ifwaww_cmd.SetRange("inner_frames_width_along_wlspwidth>0.");

    G4GenericMessenger::Command& dfnawl_cmd =
      msg_->DeclareProperty("df_no_along_wlsplength", df_no_along_wlsplength_,
			    "Number of dichroic filters along the WLS plate length.");
    dfnawl_cmd.SetParameterName("df_no_along_wlsplength", false);
    dfnawl_cmd.SetRange("df_no_along_wlsplength>=0"); 

    G4GenericMessenger::Command& dfnaww_cmd =
      msg_->DeclareProperty("df_no_along_wlspwidth", df_no_along_wlspwidth_,
			    "Number of dichroic filters along the WLS plate width.");
    dfnaww_cmd.SetParameterName("df_no_along_wlspwidth", false);
    dfnaww_cmd.SetRange("df_no_along_wlspwidth>=0"); 

    G4GenericMessenger::Command& ct_cmd =
      msg_->DeclareProperty("case_thickn", case_thickn_,
			    "Thickness of the XArapuca case.");
    ct_cmd.SetUnitCategory("Length");
    ct_cmd.SetParameterName("case_thickn", false);
    ct_cmd.SetRange("case_thickn>0.");

    G4GenericMessenger::Command& pscc_cmd =
      msg_->DeclareProperty("PS_config_code", PS_config_code_,
			    "Photo sensors configuration code.");
    pscc_cmd.SetParameterName("PS_config_code", false);
    pscc_cmd.SetRange("PS_config_code>=1"); 

    G4GenericMessenger::Command& np_cmd =
      msg_->DeclareProperty("num_phsensors", num_phsensors_,
			    "Number of photosensors per long side.");
    np_cmd.SetParameterName("num_phsensors", false);
    np_cmd.SetRange("num_phsensors>=0"); 

    G4GenericMessenger::Command& g_cmd =
      msg_->DeclareProperty("gap", gap_,
			    "Gap between the WLS plate and the photosensors.");
    g_cmd.SetUnitCategory("Length");
    //g_cmd.SetParameterName("gap", false);
    //g_cmd.SetRange("gap>0.");    // These are commented so that gap_ can help modelate the immersion of the SiPMs into the flat dimple

    G4GenericMessenger::Command& rps_cmd =
      msg_->DeclareProperty("ref_phsensors_supports", ref_phsensors_supports_,
			    "Whether the photosensors supports are VIKUITI-coated.");

    G4GenericMessenger::Command& ds_cmd =
      msg_->DeclareProperty("double_sided", double_sided_,
			    "Whether the dichroic filter is set on both sides of the WLS plate.");

    G4GenericMessenger::Command& cr_cmd =
      msg_->DeclareProperty("collectors_are_reflective", collectors_are_reflective_,
			    "Whether the test collectors are reflective.");

    G4GenericMessenger::Command& gvod_cmd =
      msg_->DeclareProperty("generation_vertex_over_df", generation_vertex_over_df_,
			    "Whether the generation vertex is randomly sampled over any dichroic filter. If not, it is randomly sampled over the whohle DFA, including the frame itself.");

    G4GenericMessenger::Command& dfafir_cmd =
      msg_->DeclareProperty("DFA_frame_is_reflective", DFA_frame_is_reflective_,
			    "Whether the FR4 DFA frame is vikuiti-coated or not.");

    G4GenericMessenger::Command& dfafis_cmd =
      msg_->DeclareProperty("DFA_frame_is_specular", DFA_frame_is_specular_,
			    "Whether the vikuiti coating of the DFA frame is specular-spikely reflective or diffusively reflective. Only makes a difference if DFA_frame_is_reflective_==True.");

    G4GenericMessenger::Command& rdfa_cmd =
      msg_->DeclareProperty("remove_DFA", remove_DFA_,
			    "Whether to remove the dichroic filters assembly or not.");

    G4GenericMessenger::Command& rdfs_cmd =
      msg_->DeclareProperty("remove_DFs", remove_DFs_,
			    "Whether to remove the dichroic filters or not.");

    G4GenericMessenger::Command& ptdd_cmd =
      msg_->DeclareProperty("path_to_dichroic_data", path_to_dichroic_data_,
			    "Absolute path to the file containing the transmission data of the dichroic filter.");

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

    G4GenericMessenger::Command& osals_cmd =
      msg_->DeclareProperty("only_sipms_along_long_sides", only_sipms_along_long_sides_,
			    "Whether SiPM boards are installed only along the long sides of the X-ARAPUCA, or along every side.");

    G4GenericMessenger::Command& wb_cmd =
      msg_->DeclareProperty("with_boards", with_boards_,
			    "Whether the SiPMs are mounted on boards or floating.");

    G4GenericMessenger::Command& wd_cmd =
      msg_->DeclareProperty("with_dimples", with_dimples_,
			    "Whether the plate has carved dimples.");

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

    G4GenericMessenger::Command& fn_cmd =
      msg_->DeclareProperty("fibers_number", fibers_no_,
			    "Number of fibers within the X-ARAPUCA.");
    fn_cmd.SetParameterName("fibers_no", false);
    fn_cmd.SetRange("fibers_no>=0"); 

    G4GenericMessenger::Command& fp_cmd =
      msg_->DeclareProperty("fiber_planes", fiber_planes_no_,
			    "Number of fiber planes within the X-ARAPUCA.");
    fp_cmd.SetParameterName("fiber_planes", false);
    fp_cmd.SetRange("fiber_planes>=0"); 

    G4GenericMessenger::Command& fr_cmd =
      msg_->DeclareProperty("fiber_radius", fiber_radius_,
			    "Radius of the fibers.");
    fr_cmd.SetUnitCategory("Length");
    fr_cmd.SetParameterName("fiber_radius", false);
    fr_cmd.SetRange("fiber_radius>0.");

    G4GenericMessenger::Command& fl_cmd =
      msg_->DeclareProperty("fiber_length", fiber_length_,
			    "Length of the fibers.");
    fl_cmd.SetUnitCategory("Length");
    fl_cmd.SetParameterName("fiber_length", false);
    fl_cmd.SetRange("fiber_length>0.");

    G4GenericMessenger::Command& als_cmd =
      msg_->DeclareProperty("along_long_side", along_long_side_,
			    "Whether the fibers are parallel to the long side of the X-ARAPUCA.");

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

    // Compute internal attributes
    overall_length_ = internal_length_  + 2*case_thickn_;
    overall_thickn_ = internal_thickn_  + 2*DFA_thickn_;    // geometry_is_ill_formed() will make sure that DF_thickn_<=DFA_thickn_
    overall_width_  = internal_width_   + 2*case_thickn_;
    DFA_length_ = overall_length_;                          // For now, these two match. 
    DFA_width_  = overall_width_;                           // For now, these two match.  
    DF_length_  = (DFA_length_  -(2.*outter_frame_width_along_wlsplength_)  -((df_no_along_wlsplength_-1.)*inner_frames_width_along_wlsplength_))   /(1.*df_no_along_wlsplength_);   
    DF_width_   = (DFA_width_   -(2.*outter_frame_width_along_wlspwidth_)   -((df_no_along_wlspwidth_-1.)*inner_frames_width_along_wlspwidth_))     /(1.*df_no_along_wlspwidth_);

    G4cout << "DF_length_=" << DF_length_ << G4endl;
    G4cout << "DF_width_=" << DF_width_ << G4endl;

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
    LAr->SetMaterialPropertiesTable(opticalprops::LAr());

    G4LogicalVolume* LAr_box_logic =
      new G4LogicalVolume(LAr_box_solid, LAr, LAr_box_name);
    LAr_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VPhysicalVolume* mother_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement (new G4RotationMatrix(), G4ThreeVector(0., 0., 0.), LAr_box_logic, 
                        LAr_box_name, world_logic, false, 0, true));

    ConstructReflectiveCase(mother_physical);
    //ConstructCollectors(mother_physical); 
    if(!remove_DFA_) ConstructDichroicAssemblies(mother_physical);
    if(config_code_==1){
        ConstructWLSPlate(mother_physical);
        if(with_boards_) ConstructBoards(mother_physical);
        else ConstructPhotosensors(mother_physical);
    }
    else if(config_code_==2){
        ConstructFibers(mother_physical);
        ConstructBoards(mother_physical);
    }
    else{
      G4Exception("[XArapuca]", "XArapuca()",
                  FatalException, "Not allowed configuration code.");
    }

    return;
  }


  void XArapuca::ConstructWLSPlate(G4VPhysicalVolume* mother_physical) const
  { 

    WLSPlate* plate = new WLSPlate  (plate_length_, plate_thickn_, plate_width_, 
                                    opticalprops::EJ286(), false,
                                    with_dimples_, dimple_type_, num_phsensors_, 
                                    !only_sipms_along_long_sides_, flat_dimple_width_, 
                                    flat_dimple_depth_, curvy_dimple_radius_);


    //WLSPlate* plate = new WLSPlate(plate_length_, plate_thickn_, plate_width_, false);
    plate->Construct();
    G4LogicalVolume* plate_logic = plate->GetLogicalVolume();
    plate_logic->SetUserLimits(ul_);

    /*
    G4VisAttributes wlsp_col = nexus::BlueAlpha();
    wlsp_col.SetForceSolid(true);
    plate_logic->SetVisAttributes(wlsp_col);
    */

    if (!plate_logic) {
      G4Exception("[XArapuca]", "ConstructWLSPlate()",
                  FatalException, "Null pointer to logical volume.");
    }

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), plate_logic->GetName(), 
                    plate_logic, mother_physical, false, 0, true);
    
    return;
  }

  void XArapuca::ConstructFibers(G4VPhysicalVolume* mother_physical) const
  {

    const G4String fiber_name = "FIBER";
    G4Tubs* fiber_solid = new G4Tubs(fiber_name, 0., fiber_radius_, fiber_length_/2., 0., twopi);
    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    pvt->SetMaterialPropertiesTable(opticalprops::EJ286());
    G4LogicalVolume* fiber_logic = new G4LogicalVolume(fiber_solid, pvt, fiber_name);
    fiber_logic->SetUserLimits(ul_);

    G4double spanning_length_perp_to_fiber = internal_width_;
    if(!along_long_side_){
        spanning_length_perp_to_fiber = internal_length_;
    }

    G4double aux_h = spanning_length_perp_to_fiber/(1.*(fibers_no_/fiber_planes_no_));
    G4double aux_v = internal_thickn_/(1.*fiber_planes_no_);

    G4int k = 0;
    G4ThreeVector pos(0., 0., 0.);
    if(along_long_side_){
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateY(+90.*deg);
        
        for(G4int i=0; i<fiber_planes_no_; i++){
            // geometry_is_ill_formed() ensures the divisibility of fibers_no_ by fiber_planes_no_
            for(G4int j=0; j<(fibers_no_/fiber_planes_no_); j++){
                pos.setY((-1.*internal_thickn_/2)+((0.5+i)*aux_v));
                pos.setZ((-1.*internal_width_/2)+((0.5+j)*aux_h));
                new G4PVPlacement(rot, pos, fiber_name, fiber_logic, 
                                mother_physical, true, k, true);
                k++;
            }
        }

    }
    else{
        for(G4int i=0; i<fiber_planes_no_; i++){
            for(G4int j=0; j<(fibers_no_/fiber_planes_no_); j++){
                pos.setY((-1.*internal_thickn_/2)+((0.5+i)*aux_v));
                pos.setX((-1.*internal_width_/2)+((0.5+j)*aux_h));
                new G4PVPlacement(nullptr, pos, fiber_name, fiber_logic, 
                                mother_physical, true, k, true);
                k++;
            }
        }
    }
    return;
  }
  

  void XArapuca::ConstructPhotosensors(G4VPhysicalVolume* mother_physical) const
  {
    G4VisAttributes sipm_col = nexus::Red();
    sipm_col.SetForceSolid(true);

    if(PS_config_code_==1){

        HamamatsuS133606050VE sipm;
        sipm.SetReflectiveSupports(ref_phsensors_supports_);
        sipm.Construct();
        G4double sipm_thickn = sipm.GetThickness();
        G4LogicalVolume* sipm_logic_vol = sipm.GetLogicalVolume();

        if (!sipm_logic_vol) {
        G4Exception("[XArapuca]", "ConstructPhotosensors()",
                    FatalException, "Null pointer to logical volume.");
        }

        sipm_logic_vol->SetVisAttributes(sipm_col);

        G4int phsensor_id = 0;

        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateX(-90.*deg);
        for (G4int i=0; i<num_phsensors_; ++i) {
            G4ThreeVector pos(  -plate_length_/2. + (0.5 + i) * plate_length_/num_phsensors_,
                                0.,
                                -plate_width_/2. -sipm_thickn/2. -gap_);

            new G4PVPlacement(  rot, pos,
                                "S133606050VE_MPPC", sipm_logic_vol,
                                mother_physical, false, phsensor_id, true);
            phsensor_id += 1;
        }

        G4RotationMatrix* rot2 = new G4RotationMatrix();
        rot2->rotateX(90.*deg);
        for (G4int i=0; i<num_phsensors_; ++i) {
            G4ThreeVector pos(-plate_length_/2. + (0.5 + i) * plate_length_/num_phsensors_,
                                0.,
                                +plate_width_/2. +sipm_thickn/2. +gap_);

            new G4PVPlacement(rot2, pos,
                                "S133606050VE_MPPC", sipm_logic_vol,
                                mother_physical, false, phsensor_id, true);
            phsensor_id += 1;
        }

        if(!only_sipms_along_long_sides_){
            G4RotationMatrix* rot3 = new G4RotationMatrix();
            rot3->rotateZ(90.*deg);
            for (G4int i=0; i<num_phsensors_; ++i) {
                G4ThreeVector pos(  -plate_length_/2. -sipm_thickn/2. -gap_,
                                    0.,
                                    -plate_width_/2. + (0.5 + i) * plate_width_/num_phsensors_);

                new G4PVPlacement(  rot3, pos,
                                    "S133606050VE_MPPC", sipm_logic_vol,
                                    mother_physical, false, phsensor_id, true);
                phsensor_id += 1;
            }

            G4RotationMatrix* rot4 = new G4RotationMatrix();
            rot4->rotateZ(-90.*deg);
            for (G4int i=0; i<num_phsensors_; ++i) {
                G4ThreeVector pos(  +plate_length_/2. +sipm_thickn/2. +gap_,
                                    0.,
                                    -plate_width_/2. + (0.5 + i) * plate_width_/num_phsensors_);

                new G4PVPlacement(  rot4, pos,
                                    "S133606050VE_MPPC", sipm_logic_vol,
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
        G4Exception("[XArapuca]", "ConstructPhotosensors()",
                    FatalException, "Null pointer to logical volume.");
        }
        sipm_along_length_logic_vol->SetVisAttributes(sipm_col);

        G4RotationMatrix* rot0 = new G4RotationMatrix();
        rot0->rotateX(-90.*deg);
        G4ThreeVector pos0( 0.,
                            0.,
                            -plate_width_/2. -sipm_thickn/2. -gap_);
        new G4PVPlacement(rot0, pos0,
                            "SS133606050VE_MPPC", sipm_along_length_logic_vol,
                            mother_physical, false, 0, true);

        G4RotationMatrix* rot1 = new G4RotationMatrix();
        rot1->rotateX(90.*deg);
        G4ThreeVector pos1( 0.,
                            0.,
                            +plate_width_/2. +sipm_thickn/2. +gap_);
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
            G4Exception("[XArapuca]", "ConstructPhotosensors()",
                        FatalException, "Null pointer to logical volume.");
            }
            sipm_along_width_logic_vol->SetVisAttributes(sipm_col);

            G4RotationMatrix* rot2 = new G4RotationMatrix();
            rot2->rotateY(-90.*deg);
            rot2->rotateX(-90.*deg);
            G4ThreeVector pos2(  -plate_length_/2. -sipm_thickn/2. -gap_,
                                0.,
                                0.);
            new G4PVPlacement(rot2, pos2,
                                "SS133606050VE_MPPC", sipm_along_width_logic_vol,
                                mother_physical, false, 2, true);

            G4RotationMatrix* rot3 = new G4RotationMatrix();
            rot3->rotateY(-90.*deg);
            rot3->rotateX(90.*deg);
            G4ThreeVector pos3(  +plate_length_/2. +sipm_thickn/2. +gap_,
                                0.,
                                0.);
            new G4PVPlacement(rot3, pos3,
                                "SS133606050VE_MPPC", sipm_along_width_logic_vol,
                                mother_physical, false, 3, true);
        }

    }
    else{
        G4Exception(    "[XArapuca]", "ConstructPhotosensors()",
                        FatalException, "Not allowed configuration code.");

    }

    return;
  }

  void XArapuca::ConstructBoards(G4VPhysicalVolume* mother_physical) const
  {

    if(config_code_==1){
        SiPMBoard board1;
        board1.SetBaseID(0);
        board1.SetBoardLength(plate_length_);
        board1.SetNumPhsensors(num_phsensors_);
        board1.SetReflectiveSupports(ref_phsensors_supports_);
        board1.Construct();
        G4LogicalVolume* board1_logic_vol = board1.GetLogicalVolume();

        G4double z_pos = (plate_width_/2.) +(board1.GetOverallThickness()/2.) +gap_;
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., z_pos),
                        "SIPMS_BOARD", board1_logic_vol, mother_physical, false, 0, false);
        // SiPMBoard logical volume is an encasing volume which may collide into other volumes
        // No need to set pSurfCheck for that volume (dimples). As we are setting pSurfCheck=false
        // (so that no harmless overlap warning pops up in a with-dimples configuration), you have to
        // be extra careful to examine when there's actually a problematic overlap of this volume with
        // another one (since Geant4 won't warn you).

        SiPMBoard board2;
        board2.SetBaseID(num_phsensors_);
        board2.SetBoardLength(plate_length_);
        board2.SetNumPhsensors(num_phsensors_);
        board2.SetReflectiveSupports(ref_phsensors_supports_);
        board2.Construct();
        G4LogicalVolume* board2_logic_vol = board2.GetLogicalVolume();

        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateY(+180.*deg);

        z_pos = (plate_width_/2.) +(board2.GetOverallThickness()/2.) +gap_;
        new G4PVPlacement(rot, G4ThreeVector(0., 0., -1.*z_pos),
                        "SIPMS_BOARD", board2_logic_vol, mother_physical, true, 0, false);

        if(!only_sipms_along_long_sides_){
            SiPMBoard board3;
            board3.SetBaseID(2*num_phsensors_);
            board3.SetBoardLength(plate_width_);
            board3.SetNumPhsensors(num_phsensors_);
            board3.SetReflectiveSupports(ref_phsensors_supports_);
            board3.Construct();
            G4LogicalVolume* board3_logic_vol = board3.GetLogicalVolume();

            G4RotationMatrix* rot2 = new G4RotationMatrix();
            rot2->rotateY(-90.*deg);

            G4double x_pos = (plate_length_/2.) +(board3.GetOverallThickness()/2.) +gap_;
            new G4PVPlacement(rot2, G4ThreeVector(x_pos, 0., 0.),
                            "SIPMS_BOARD", board3_logic_vol, mother_physical, true, 0, false);

            SiPMBoard board4;
            board4.SetBaseID(3*num_phsensors_);
            board4.SetBoardLength(plate_width_);
            board4.SetNumPhsensors(num_phsensors_);
            board4.SetReflectiveSupports(ref_phsensors_supports_);
            board4.Construct();
            G4LogicalVolume* board4_logic_vol = board4.GetLogicalVolume();

            G4RotationMatrix* rot3 = new G4RotationMatrix();
            rot3->rotateY(+90.*deg);

            x_pos = (plate_length_/2.) +(board4.GetOverallThickness()/2.) +gap_;
            new G4PVPlacement(rot3, G4ThreeVector(-1.*x_pos, 0., 0.),
                            "SIPMS_BOARD", board4_logic_vol, mother_physical, true, 0, false);
        }
    }
    else if(config_code_==2){
        SiPMBoard board1;
        board1.SetBaseID(0);
        board1.SetNumPhsensors(num_phsensors_);
        board1.SetReflectiveSupports(ref_phsensors_supports_);

        SiPMBoard board2;
        board2.SetBaseID(num_phsensors_);
        board2.SetNumPhsensors(num_phsensors_);
        board2.SetReflectiveSupports(ref_phsensors_supports_);

        if(along_long_side_){
            board1.SetBoardLength(internal_width_);
            board1.Construct();
            G4LogicalVolume* board1_logic_vol = board1.GetLogicalVolume();
        
            G4RotationMatrix* rot = new G4RotationMatrix();
            rot->rotateY(-90.*deg);

            G4double x_pos = (fiber_length_/2.) + (board1.GetOverallThickness()/2.) +gap_;
            G4ThreeVector pos(x_pos, 0., 0.);
            new G4PVPlacement(rot, pos, "SIPMS_BOARD", board1_logic_vol, 
                            mother_physical, true, 0, true);

            board2.SetBoardLength(internal_width_);
            board2.Construct();
            G4LogicalVolume* board2_logic_vol = board2.GetLogicalVolume();

            G4RotationMatrix* rot_2 = new G4RotationMatrix();
            rot_2->rotateY(+90.*deg);

            x_pos = (fiber_length_/2.) + (board2.GetOverallThickness()/2.) +gap_;
            G4ThreeVector pos_2(-1.*x_pos, 0., 0.);
            new G4PVPlacement(rot_2, pos_2, "SIPMS_BOARD", board2_logic_vol, 
                            mother_physical, true, 0, true);
        }
        else{
            board1.SetBoardLength(internal_length_);
            board1.Construct();
            G4LogicalVolume* board1_logic_vol = board1.GetLogicalVolume();

            G4double z_pos = (fiber_length_/2.) + (board1.GetOverallThickness()/2.) +gap_;
            G4ThreeVector pos(0., 0., z_pos);
            new G4PVPlacement(nullptr, pos, "SIPMS_BOARD", board1_logic_vol, 
                            mother_physical, true, 0, true);

            board2.SetBoardLength(internal_length_);
            board2.Construct();
            G4LogicalVolume* board2_logic_vol = board2.GetLogicalVolume();

            G4RotationMatrix* rot = new G4RotationMatrix();
            rot->rotateY(+180.*deg);

            z_pos = (fiber_length_/2.) + (board2.GetOverallThickness()/2.) +gap_;
            G4ThreeVector pos_2(0., 0., -1.*z_pos);
            new G4PVPlacement(rot, pos_2, "SIPMS_BOARD", board2_logic_vol, 
                            mother_physical, true, 0, true);
        }
    }
    else{
        G4Exception("[XArapuca]", "ConstructBoards()",
                    FatalException, "Not allowed configuration code.");
    }

    return;

  }

  void XArapuca::ConstructReflectiveCase(G4VPhysicalVolume* mother_physical) const
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

    // The reflective case is 'open' on both sides. If the X-ARAPUCA is not double sided, then ConstructDichroicAssemblies will take care of placing
    // a reflective cover. 

    G4LogicalVolume* ref_case_logic = 
      new G4LogicalVolume(ref_case_solid, materials::FR4(), ref_case_name);

    // Set its color for visualization purposes
    G4VisAttributes ref_case_col = nexus::WhiteAlpha();
    ref_case_col.SetForceSolid(true);
    ref_case_logic->SetVisAttributes(ref_case_col);

    //Now create the reflectivie optical surface
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

  
  void XArapuca::ConstructDichroicAssemblies(G4VPhysicalVolume* mother_physical) const
  {
      
    // Dichroic-filters assembly. It is made up of one-piece frame, and the dichroic filters.
    G4Box* cover_solid = new G4Box("AUX", DFA_length_/2., DFA_thickn_/2., DFA_width_/2.);
    G4double tolerance = 5.*mm; // Tolerance to prevent matching surfaces when subtracting carvings from frame
    
    G4Box* carving_solid                =  new G4Box("FRAME_CARVING", DF_length_/2., DFA_thickn_/2. +tolerance, DF_width_/2.); 
    G4Box* df_solid                     =  new G4Box("DICHROIC_FILTER", DF_length_/2., DF_thickn_/2., DF_width_/2.); 
    G4Box* coating_solid                =  nullptr;
    if(DF_are_coated_) coating_solid    = new G4Box("DF_COATING", DF_length_/2., coating_thickn_/2., DF_width_/2.); 

    G4MultiUnion* filters_multiunion_solid          = new G4MultiUnion("DICHROIC_FILTERS");
    G4MultiUnion* frame_carvings_multiunion_solid   = new G4MultiUnion("FRAME_CARVINGS");
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
            filters_multiunion_solid->AddNode(*df_solid, *transform_ptr);
            frame_carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
            if(DF_are_coated_) coatings_multiunion_solid->AddNode(*coating_solid, *transform_ptr);
        }
    }
    filters_multiunion_solid->Voxelize();
    frame_carvings_multiunion_solid->Voxelize();
    if(DF_are_coated_) coatings_multiunion_solid->Voxelize();

    // FRAME //

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
    frame_col.SetForceSolid(true);
    frame_logic->SetVisAttributes(frame_col);

    if(DFA_frame_is_reflective_){
        const G4String refcoat_name = "REF_COATING";
        G4OpticalSurface* refcoat_opsurf = 
        new G4OpticalSurface(refcoat_name, unified, ground, dielectric_metal, 1);
        if(DFA_frame_is_specular_)  refcoat_opsurf->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());
        else                        refcoat_opsurf->SetMaterialPropertiesTable(opticalprops::diffusiveVIKUITI());
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
        refsurf_opsurf->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());

        const G4String cover_name = "REFLECTIVE_COVER";
        G4LogicalVolume* cover_logic = new G4LogicalVolume(cover_solid, materials::FR4(), cover_name);
        cover_logic->SetVisAttributes(frame_col);
        new G4LogicalSkinSurface(refsurf_name, cover_logic, refsurf_opsurf);   
        new G4PVPlacement(nullptr, G4ThreeVector(0., -1.*(internal_thickn_+DFA_thickn_)/2., 0.), 
                            cover_name, cover_logic, mother_physical, true, 1, true);
    }

    // DICHROIC FILTERS //
    if(!remove_DFs_){

            G4Material* dfs_substrate = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
            // The filter will be implemented as a G4LogicalBorderSurface from the surrounding LAr physical volume
            // towards the filter physical volume. For photons that are transmited from outside the X-ARAPUCA towards
            // the dichroic filter volume, they shall not be further reflected or refracted. To avoid this, the 
            // dichroic filter substrate optical properties must match those of the surrounding LAr. This was checked 
            // in an alternative setup.
            dfs_substrate->SetMaterialPropertiesTable(opticalprops::LAr());
            G4LogicalVolume* filters_multiunion_logic = 
                                new G4LogicalVolume(filters_multiunion_solid, dfs_substrate, "DICHROIC_FILTERS");

            G4VisAttributes dfs_col = nexus::BloodRedAlpha();
            dfs_col.SetForceSolid(true);
            filters_multiunion_logic->SetVisAttributes(dfs_col);

            // Place the filters
            G4VPhysicalVolume* first_df_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                new G4PVPlacement(nullptr, G4ThreeVector(0., (internal_thickn_+DF_thickn_)/2.+ (DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_)), 0.), 
                                "DICHROIC_FILTERS", filters_multiunion_logic, mother_physical, true, 0, true));

            // Add dichroic specifications
            if(path_to_dichroic_data_==""){
                G4Exception("[XArapuca]", "ConstructDichroicAssemblies()",
                            FatalException, "The path to the dichroic data file was not set.");
            }

            setenv("G4DICHROICDATA", path_to_dichroic_data_, 1);
            G4OpticalSurface* dfs_opsurf =   
                new G4OpticalSurface("DICHROIC_OPSURF", dichroic, polished, dielectric_dichroic);
            new G4LogicalBorderSurface("LAr->DICHROIC1", mother_physical, first_df_multiunion_physical, dfs_opsurf);
            
            if(double_sided_){
                G4VPhysicalVolume* second_df_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                    new G4PVPlacement(nullptr, G4ThreeVector(0., -1.*((internal_thickn_+DF_thickn_)/2.+ (DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_))), 0.), 
                                    "DICHROIC_FILTERS", filters_multiunion_logic, mother_physical, true, 1, true));
                new G4LogicalBorderSurface("LAr->DICHROIC2", mother_physical, second_df_multiunion_physical, dfs_opsurf);
            }

        // COATINGS //

        G4double ptp_df_gap = 1.*um;
        // As far as the G4LogicalBorderSurface that implements the dichroic transmitance curve is now written (mother_physical->df_multiunion_physical), there 
        // should be a gap in between the PTP and the DF. The volumes-evolution that a PTP-WLSed photon follows is:
        // (1) LAr world (mother_physical) -> 
        // (2) PTP (WLSed and luckily emitted towards the XA internal cavity) -> 
        // (3) LAr world (mother_physical) ->            
        // (4) DF -> 
        // (5) LAr world -> 
        // (6) WLS plate -> ...
        // The gap in (3) is needed because, in this way, the photon "sees" the DF transmitance curve in the transition (3)->(4)

        if(DF_are_coated_){
            G4Material* coatings_substrate = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
            coatings_substrate->SetMaterialPropertiesTable(opticalprops::LArPTPArtifact());
            G4LogicalVolume* coatings_multiunion_logic = 
                                new G4LogicalVolume(coatings_multiunion_solid, coatings_substrate, "DF_COATINGS");   

            G4VisAttributes coatings_col = nexus::TitaniumGreyAlpha();
            coatings_col.SetForceSolid(true);
            coatings_multiunion_logic->SetVisAttributes(coatings_col);

            // Place the coatings
            G4VPhysicalVolume* first_coatings_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                new G4PVPlacement(  nullptr, G4ThreeVector(0.,  (internal_thickn_+DF_thickn_)/2.+ (DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_)) 
                                                                +(DF_thickn_/2.) +(coating_thickn_/2.) +ptp_df_gap, 0.), 
                                    "DF_COATINGS", coatings_multiunion_logic, mother_physical, true, 0, true));

            if(double_sided_){
                G4VPhysicalVolume* second_coatings_multiunion_physical = dynamic_cast<G4VPhysicalVolume*>(
                    new G4PVPlacement(  nullptr, G4ThreeVector(0.,  -1.*((internal_thickn_+DF_thickn_)/2.+ (DF_pos_wrt_DFA_pos_*(DFA_thickn_-DF_thickn_)) 
                                                                    +(DF_thickn_/2.) +(coating_thickn_/2.) +ptp_df_gap), 0.), 
                                        "DF_COATINGS", coatings_multiunion_logic, mother_physical, true, 1, true));
            }
        }
    }


    return;
  }

  G4bool XArapuca::geometry_is_ill_formed()
  {
    // Some checks are independent of the configuration code. Perform those tests here:
    if(DF_thickn_>DFA_thickn_){
        G4Exception("[XArapuca]", "geometry_is_ill_formed()", FatalException,
        "The dichroic filters thickness cannot be bigger than the frame thickness.");
    }

    if(DF_pos_wrt_DFA_pos_<0. || DF_pos_wrt_DFA_pos_>1.){
        G4Exception("[XArapuca]", "geometry_is_ill_formed()", FatalException,
        "The dichroic filter position with respect to the DFA position must belong to the [0, 1] interval.");
    }

    // What we need to make sure here is that there's room enough within the XArapuca internal cavity
    // to allocate everything that we intending to put in it. To do so, we may find the span of the 
    // internal geometry along each axis, and then compare it to the dimensions of the internal cavity.

    if(internal_length_<internal_width_){
        G4Exception("[XArapuca]", "geometry_is_ill_formed()", FatalException,
        "For consistency with the code, internal_length_ must be greater than internal_width_.");
    }
    
    G4double internal_geom_length_span, internal_geom_width_span, internal_geom_thickn_span;
    if(config_code_==1){
        if(plate_length_<plate_width_){
            G4Exception("[XArapuca]", "geometry_is_ill_formed()", FatalException,
            "For consistency with the code, plate_length_ must be greater than plate_width_.");

        }

        if(with_boards_){
            // If with_boards_==true, then PS_config_code_ is ignored, the geometry is loaded
            // as if PS_config_code_==1, so there's no need to distinguish cases here
            SiPMBoard board;
            internal_geom_length_span   = plate_length_ +(2.*gap_)+(2.*board.GetOverallThickness());
            internal_geom_width_span    = plate_width_  +(2.*gap_)+(2.*board.GetOverallThickness());
            internal_geom_thickn_span   = std::max(plate_thickn_, board.GetOverallHeight());
        }
        else{
            HamamatsuS133606050VE sipm; 
            if(PS_config_code_==1){
                if(num_phsensors_*sipm.GetWidth()>plate_length_) { return true; }
                if(!only_sipms_along_long_sides_){
                    if(num_phsensors_*sipm.GetWidth()>plate_width_) { return true; }
                }

            }
            // Else if PS_config_code_==2, the scalable SiPM is build with the exact plate 
            // dimensions, so there's no need to care about such case here

            internal_geom_length_span = plate_length_+(2.*gap_)+(2.*sipm.GetThickness());   // Either the HamamatsuS133606050VE and the ScalableHamamatsuS133606050VE
                                                                                            // have the same thickness, so there's no need to distinguish between both
                                                                                            // according to PS_config_code_ value
            internal_geom_width_span = plate_width_+(2.*gap_)+(2.*sipm.GetThickness());
            internal_geom_thickn_span = std::max(plate_thickn_, sipm.GetHeight());
        }
    }
    else if(config_code_==2)
    {

        if(fibers_no_%fiber_planes_no_!=0){
            G4Exception("[XArapuca]", "geometry_is_ill_formed()",
                        FatalException, "The number of fibers to install must be divisible by the number of fiber planes.");
        }
        G4double fibers_per_plane = fibers_no_/fiber_planes_no_;

        SiPMBoard board;

        internal_geom_thickn_span = std::max(fiber_planes_no_*2*fiber_radius_, board.GetOverallHeight());
        internal_geom_length_span = fiber_length_+(2.*board.GetOverallThickness());
        internal_geom_width_span = fibers_per_plane*2*fiber_radius_;
        
        if(!along_long_side_){
            std::swap(internal_geom_length_span, internal_geom_width_span);
        }
    }
    else{
        G4Exception("[XArapuca]", "geometry_is_ill_formed()",
                    FatalException, "Not allowed configuration code.");
    }
    
    G4bool check1 = internal_length_ >= internal_geom_length_span;
    G4bool check2 = internal_thickn_ >= internal_geom_thickn_span;
    G4bool check3 = internal_width_  >= internal_geom_width_span;

    // Uncomment the following chunk for testing purposes
    /*
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
    */

    return !(check1*check2*check3);
  }

  G4ThreeVector XArapuca::GenerateVertex(const G4String&) const{

    G4double tolerance = 0.1*mm;    // Short distance over the dichroic filter (DF)
                                    // from which photons are launched. Also, the 
                                    // width of the outter border projected over the
                                    // DF from which photons won't be launched
                                    // (Just see the implementation in x_pos and
                                    // z_pos below to understand its meaning)

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dist(0, df_no_along_wlsplength_*df_no_along_wlspwidth_ -1);
    G4int filter_no = dist(gen);
    
    G4int column_no = filter_no%df_no_along_wlsplength_;                // \in[0, df_no_along_wlsplength -1]  
    G4int row_no = (int) std::floor(filter_no/df_no_along_wlsplength_); // \in[0, df_no_along_wlspwidth -1]

    G4double selected_filter_x_center =  (-1.*(DFA_length_/2.)) +outter_frame_width_along_wlsplength_   +(column_no*inner_frames_width_along_wlsplength_)   +((column_no+0.5)*DF_length_);
    G4double selected_filter_z_center =  (-1.*(DFA_width_/2.))  +outter_frame_width_along_wlspwidth_    +(row_no*inner_frames_width_along_wlspwidth_)       +((row_no+0.5)*DF_width_);
    G4double x_pos, z_pos;
    G4double y_pos = (internal_thickn_/2.) +DFA_thickn_ +tolerance;

    if(!generation_vertex_over_df_){
        x_pos = UniformRandomInRange(  -1.*DFA_length_/2.,
                                        DFA_length_/2.      );
        z_pos = UniformRandomInRange(  -1.*DFA_width_/2.,
                                        DFA_width_/2.       );
    }
    else{
        x_pos = UniformRandomInRange(  selected_filter_x_center -(DF_length_/2.) +tolerance, 
                                                selected_filter_x_center +(DF_length_/2.) -tolerance    );
        z_pos = UniformRandomInRange(  selected_filter_z_center -(DF_width_/2.) +tolerance, 
                                                selected_filter_z_center +(DF_width_/2.) -tolerance     );
    }
    return G4ThreeVector(x_pos, y_pos, z_pos);
  }

} //End namespace nexus
