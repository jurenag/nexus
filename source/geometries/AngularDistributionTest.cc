#include "AngularDistributionTest.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "Visibilities.h"
#include "RandomUtils.h"

#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4SubtractionSolid.hh>
#include <G4PVPlacement.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4GenericMessenger.hh>
#include <G4SystemOfUnits.hh>

using namespace nexus;

REGISTER_CLASS(AngularDistributionTest, GeometryBase)

namespace nexus{

  AngularDistributionTest::AngularDistributionTest():
  GeometryBase(),
  config_code_              (1                  ),
  filter_code_              (1                  ),
  world_radius_             (1.*m               ),
  world_height_             (1.*cm              ),
  sample_LDP_distance_      (0.5*cm             ),
  lgv_zpos_wrt_world_center_(2.5*mm             ),
  DF_thickn_                (1.*mm              ),
  ptp_thickn_               (1.*um              ),
  DF_x_size_                (10.*cm             ),
  DF_y_size_                (10.*cm             ),
  rib_height_               (5.*mm              ),
  shallowness_              (0.                 ),
  sm_mpt_                   (opticalprops::LAr())
  {

    msg_ = new G4GenericMessenger(this, "/Geometry/AngularDistributionTest/",
				"Control commands of geometry AngularDistributionTest.");

    G4GenericMessenger::Command& cc_cmd =
      msg_->DeclareProperty("config_code", config_code_,
			    "Configuration code.");
    cc_cmd.SetParameterName("config_code", false);
    cc_cmd.SetRange("config_code>=1"); 

    G4GenericMessenger::Command& fc_cmd =
      msg_->DeclareProperty("filter_code", filter_code_,
			    "Which filter to use.");
    fc_cmd.SetParameterName("filter_code", false);
    fc_cmd.SetRange("filter_code>=1"); 

    G4GenericMessenger::Command& wr_cmd =
      msg_->DeclareProperty("world_radius", world_radius_,
			    "World and LDP radius.");
    wr_cmd.SetUnitCategory("Length");
    wr_cmd.SetParameterName("world_radius", false);
    wr_cmd.SetRange("world_radius>0.");

    G4GenericMessenger::Command& wh_cmd =
      msg_->DeclareProperty("world_height", world_height_,
			    "World height.");
    wh_cmd.SetUnitCategory("Length");
    wh_cmd.SetParameterName("world_height", false);
    wh_cmd.SetRange("world_height>0.");

    G4GenericMessenger::Command& sldpd_cmd =
      msg_->DeclareProperty("sample_LDP_distance", sample_LDP_distance_,
			    "Distance from the LDP to the sample.");
    sldpd_cmd.SetUnitCategory("Length");
    sldpd_cmd.SetParameterName("sample_LDP_distance", false);
    sldpd_cmd.SetRange("sample_LDP_distance>0.");

    G4GenericMessenger::Command& lgvhdfs_cmd =
      msg_->DeclareProperty("lgv_zpos_wrt_world_center", lgv_zpos_wrt_world_center_,
			    "Light-generation-vertex z coordinate wrt the center of the world.");
    lgvhdfs_cmd.SetUnitCategory("Length");
    lgvhdfs_cmd.SetParameterName("lgv_zpos_wrt_world_center", false);
    lgvhdfs_cmd.SetRange("lgv_zpos_wrt_world_center>0.");

    G4GenericMessenger::Command& dft_cmd =
      msg_->DeclareProperty("DF_thickn", DF_thickn_,
			    "Dichroic filter thickness.");
    dft_cmd.SetUnitCategory("Length");
    dft_cmd.SetParameterName("DF_thickn", false);
    dft_cmd.SetRange("DF_thickn>0.");

    G4GenericMessenger::Command& ptpt_cmd =
      msg_->DeclareProperty("ptp_thickn", ptp_thickn_,
			    "PTP coating thickness.");
    ptpt_cmd.SetUnitCategory("Length");
    ptpt_cmd.SetParameterName("ptp_thickn", false);
    ptpt_cmd.SetRange("ptp_thickn>0.");

    G4GenericMessenger::Command& dfxs_cmd =
      msg_->DeclareProperty("DF_x_size", DF_x_size_,
			    "Dichroic Filter size along x dimension.");
    dfxs_cmd.SetUnitCategory("Length");
    dfxs_cmd.SetParameterName("DF_x_size", false);
    dfxs_cmd.SetRange("DF_x_size>0.");

    G4GenericMessenger::Command& dfys_cmd =
      msg_->DeclareProperty("DF_y_size", DF_y_size_,
			    "Dichroic Filter size along y dimension.");
    dfys_cmd.SetUnitCategory("Length");
    dfys_cmd.SetParameterName("DF_y_size", false);
    dfys_cmd.SetRange("DF_y_size>0.");

    G4GenericMessenger::Command& rh_cmd =
      msg_->DeclareProperty("rib_height", rib_height_,
			    "Height of the ribs that surround the dichroic filter.");
    rh_cmd.SetUnitCategory("Length");
    rh_cmd.SetParameterName("rib_height", false);
    rh_cmd.SetRange("rib_height>0.");

    G4GenericMessenger::Command& s_cmd =
      msg_->DeclareProperty("shallowness", shallowness_,
			    "Shallowness parameter of the DF within the ribs.");
    s_cmd.SetParameterName("shallowness", false);
    s_cmd.SetRange("shallowness>=0.");
    s_cmd.SetRange("shallowness<=1.");

    G4GenericMessenger::Command& rslgv_cmd =
      msg_->DeclareProperty("randomly_sample_lgv", randomly_sample_lgv_,
			    "Whether the generation vertex is randomly sampled over the coating filter surface.");    

  }

  AngularDistributionTest::~AngularDistributionTest()
  {
    if(msg_) delete msg_;
  }

  void AngularDistributionTest::Construct()
  {

    if(geometry_is_ill_formed()){
      G4Exception("AngularDistributionTest::Construct()", "0",
                  FatalException, "Geometry is ill formed.");
    }

    // WORLD ENCASING
    G4double world_tolerance = 1.*m;
    const G4String world_encasing_name = "WORLD_ENCASING";
    G4Tubs* world_encasing_solid =
      new G4Tubs(world_encasing_name, 0., world_radius_+world_tolerance, (world_height_+world_tolerance)/2., 0., 360.*degree);

    G4Material* vacuum =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4LogicalVolume* world_encasing_logic = 
        new G4LogicalVolume(world_encasing_solid, vacuum, world_encasing_name, 0, 0, 0, true);
        
    //world_encasing_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(world_encasing_logic);

    // WORLD (Fits tight within its encasing)
    G4String world_name = "WORLD";
    G4Tubs* world_solid = 
      new G4Tubs(world_name, 0., world_radius_, world_height_/2., 0., 360.*degree);

    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    lAr->SetMaterialPropertiesTable(sm_mpt_);

    G4LogicalVolume* world_logic = new G4LogicalVolume(world_solid, lAr, world_name);
    //world_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    // Place the world within its encasing
    G4VPhysicalVolume* world_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement (new G4RotationMatrix(), G4ThreeVector(0., 0., 0.), world_logic, 
                        world_name, world_encasing_logic, false, 0, true));

    ConstructLDP(world_physical);
    ConstructSample(world_physical); 
    return;
  }

  void AngularDistributionTest::ConstructSample(G4VPhysicalVolume* mother_physical)
  {
    if(config_code_==1) ConstructWorldTightDichroicFilter(mother_physical);
    else if(config_code_==2) ConstructFramedDichroicFilter(mother_physical);
    else if(config_code_==3) ConstructFramedPTPAlone(mother_physical);
    else{
      G4Exception("AngularDistributionTest::ConstructSample()", "0",
                  FatalException, "Not recognized configuration code.");
    }
    return;
  }

  void AngularDistributionTest::ConstructWorldTightDichroicFilter(G4VPhysicalVolume* mother_physical)
  {

    // When a manufacturer produces a filter and gives you the transmission curve,
    // that curve accounts for the whole filter, i.e. the dichroic deposition 
    // plus the dichroic substrate. So there's no need to separate these two volumes
    // here, i.e. you are implementing the whole dichroic filter by assigning the
    // transmission curve over the volume. However, there's an implementation 
    // peculiarity you have to take into account here. If you assign the dichroic
    // transmission curve in the transition LAr_volume->dichroic_volume, then
    // the photons will enter the dichroic_volume or get reflected back into the LAr
    // according to the transmission curve you provided. Now, if a reflection takes
    // place, then everything is fine. The angle of the reflected photon is the same
    // as the angle of the incident photon that is reflected. However, for a transmitted
    // photon, it is flying within the dichroic volume and it will eventually escape it.
    // When trying to escape it, the dichroic_volume needs to have a defined rindex, 
    // otherwise the photon will get killed in the boundary. Now, if the defined rindex
    // is different from that of the LAr, in case of refraction, the photon will get
    // deflected according to snell's law. Thus giving a real transmission curve
    // that is different from that implemented by the dichroic data. There are two
    // ways to solve this. Either place the second detector WITHIN the dichroic volume, 
    // or give the dichroic_volume the same optical properties as the surrounding volume.
    // Let us try the second one, which is the one that could be more easily implemented
    // in the eventual X-ARAPUCA geometry.

    const G4String df_name = "DICHROIC_SUBSTRATE";

    G4Tubs* df_solid =
      new G4Tubs(df_name, 0., world_radius_, DF_thickn_/2., 0., 360.*degree);

    G4Material* substrate_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
    substrate_mat->SetMaterialPropertiesTable(opticalprops::SCHOTT_B270());

    G4LogicalVolume* df_logic =
      new G4LogicalVolume(df_solid, substrate_mat, df_name);

    G4VisAttributes df_col = nexus::RedAlpha();
    df_col.SetForceSolid(true);
    df_logic->SetVisAttributes(df_col);

  
    G4VPhysicalVolume* df_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), df_name, df_logic,
        mother_physical, false, 0, true));

    // PTP COATING
    const G4String ptpd_name = "PTP_LAYER";
    G4Tubs* ptpd_solid =
      new G4Tubs(ptpd_name, 0., world_radius_, ptp_thickn_/2., 0., 360.*degree);

    G4Material* ptpd_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
    ptpd_mat->SetMaterialPropertiesTable(opticalprops::PTP());

    G4LogicalVolume* ptpd_logic =
      new G4LogicalVolume(ptpd_solid, ptpd_mat, ptpd_name, nullptr, nullptr, nullptr, true);

    G4VisAttributes ptpd_col = nexus::TitaniumGreyAlpha();
    ptpd_col.SetForceSolid(true);
    ptpd_logic->SetVisAttributes(ptpd_col);

    G4VPhysicalVolume* ptpd_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., (DF_thickn_/2.)+(ptp_thickn_/2.)), ptpd_name, ptpd_logic,
        mother_physical, false, 0, true));

    
    
    G4OpticalSurface* coating_rough_surf =
          new G4OpticalSurface("COATING_ROUGH_SURFACE", glisur, ground, dielectric_dielectric, .01);
            // 0.01 is the polish value for glisur model that was measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
            // This is the best reference we have, since both PTP and TPB are the result of an evaporation+deposition process
    new G4LogicalBorderSurface("SURROUNDINGS->COATING", mother_physical, ptpd_physical, coating_rough_surf);
    new G4LogicalBorderSurface("COATING->SURROUNDINGS", ptpd_physical, mother_physical, coating_rough_surf);
    new G4LogicalBorderSurface("COATING->DF", ptpd_physical, df_physical, coating_rough_surf);
    

    // Add dichroic specifications
    if(filter_code_==1){
      setenv("G4DICHROICDATA", "/home/jurenag/installations/nexus-source/data/full_transmission_test_dichroic_data", 1);
    }
    else if(filter_code_==2){
      setenv("G4DICHROICDATA", "/home/jurenag/cernbox/work_stuff/documentation_and_notes/materiales/filters_data/opto/ific_meas/g4_dichroic_data_from_schottb270", 1);
    }
    else{
      G4Exception("AngularDistributionTest::ConstructFramedDichroicFilter()", "0",
            FatalException, "Not recognized filter code.");
    }

    G4OpticalSurface* dichroic_opsurf =
        new G4OpticalSurface("DICHROIC_OPSURF", dichroic, polished, dielectric_dichroic);
    
    G4LogicalBorderSurface* first_lbs = new G4LogicalBorderSurface("LAr->DICHROIC1", df_physical, mother_physical, dichroic_opsurf);
    
    return;
  }

  void AngularDistributionTest::ConstructFramedDichroicFilter(G4VPhysicalVolume* mother_physical)
  {
    // When a manufacturer produces a filter and gives you the transmission curve,
    // that curve accounts for the whole filter, i.e. the dichroic deposition 
    // plus the dichroic substrate. So there's no need to separate these two volumes
    // here, i.e. you are implementing the whole dichroic filter by assigning the
    // transmission curve over the volume. However, there's an implementation 
    // peculiarity you have to take into account here. If you assign the dichroic
    // transmission curve in the transition LAr_volume->dichroic_volume, then
    // the photons will enter the dichroic_volume or get reflected back into the LAr
    // according to the transmission curve you provided. Now, if a reflection takes
    // place, then everything is fine. The angle of the reflected photon is the same
    // as the angle of the incident photon that is reflected. However, for a transmitted
    // photon, it is flying within the dichroic volume and it will eventually escape it.
    // When trying to escape it, the dichroic_volume needs to have a defined rindex, 
    // otherwise the photon will get killed in the boundary. Now, if the defined rindex
    // is different from that of the LAr, in case of refraction, the photon will get
    // deflected according to snell's law. Thus giving a real transmission curve
    // that is different from that implemented by the dichroic data. There are two
    // ways to solve this. Either place the second detector WITHIN the dichroic volume, 
    // or give the dichroic_volume the same optical properties as the surrounding volume.
    // Let us try the second one, which is the one that could be more easily implemented
    // in the eventual X-ARAPUCA geometry.


    // DICHROIC FILTER SURROUNDING-RIBS
    G4Box* aux_ribs_outer =
    new G4Box("AUX_RIBS_OUTER", (DF_x_size_/2.)+rib_height_, (DF_y_size_/2.)+rib_height_, rib_height_/2.);

    G4double geometric_tolerance = 1.*mm; // To avoid subtraction of solids with matching surfaces
    G4Box* aux_ribs_inner =
    new G4Box("AUX_RIBS_INNER", DF_x_size_/2., DF_y_size_/2., (rib_height_/2.)+geometric_tolerance);

    const G4String ribs_name = "RIBS";
    G4SubtractionSolid* ribs_solid = new G4SubtractionSolid(ribs_name, aux_ribs_outer, aux_ribs_inner, 
                                            nullptr, G4ThreeVector(0., 0., 0.));

    G4LogicalVolume* ribs_logic = new G4LogicalVolume(ribs_solid, materials::FR4(), ribs_name);

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), ribs_name, ribs_logic,
                      mother_physical, false, 0, true); // Last one is pSurfChk

    // DICHROIC FILTER
    const G4String df_name = "DICHROIC_FILTER";
    G4Box* df_solid =
    new G4Box(df_name, DF_x_size_/2., DF_y_size_/2., DF_thickn_/2.);

    G4Material* substrate_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
    substrate_mat->SetMaterialPropertiesTable(opticalprops::SCHOTT_B270());

    G4LogicalVolume* df_logic =
      new G4LogicalVolume(df_solid, substrate_mat, df_name);

    G4VisAttributes df_col = nexus::RedAlpha();
    df_col.SetForceSolid(true);
    df_logic->SetVisAttributes(df_col);

    G4VPhysicalVolume* df_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, 
                          G4ThreeVector(0., 0., ((rib_height_-DF_thickn_)/2.)*((2.*shallowness_)-1.)), 
                          df_name, df_logic,
                          mother_physical, false, 0, true));

    // Add dichroic specifications
    if(filter_code_==1){
      setenv("G4DICHROICDATA", "/home/jurenag/installations/nexus-source/data/full_transmission_test_dichroic_data", 1);
    }
    else if(filter_code_==2){
      setenv("G4DICHROICDATA", "/home/jurenag/cernbox/work_stuff/documentation_and_notes/materiales/filters_data/opto/ific_meas/g4_dichroic_data_from_schottb270", 1);
    }
    else{
      G4Exception("AngularDistributionTest::ConstructFramedDichroicFilter()", "0",
            FatalException, "Not recognized filter code.");
    }
    
    G4OpticalSurface* dichroic_opsurf =
        new G4OpticalSurface("DICHROIC_OPSURF", dichroic, polished, dielectric_dichroic);
    
    G4LogicalBorderSurface* first_lbs = new G4LogicalBorderSurface("DF->LAr", df_physical, mother_physical, dichroic_opsurf);

    // PTP COATING
    const G4String ptpd_name = "PTP_LAYER";
    G4Box* ptpd_solid =
      new G4Box(ptpd_name, DF_x_size_/2., DF_y_size_/2., ptp_thickn_/2.);

    G4Material* ptpd_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
    ptpd_mat->SetMaterialPropertiesTable(opticalprops::PTP());

    G4LogicalVolume* ptpd_logic =
      new G4LogicalVolume(ptpd_solid, ptpd_mat, ptpd_name, nullptr, nullptr, nullptr, true);

    G4VisAttributes ptpd_col = nexus::TitaniumGreyAlpha();
    ptpd_col.SetForceSolid(true);
    ptpd_logic->SetVisAttributes(ptpd_col);

    G4VPhysicalVolume* ptpd_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, 
                          G4ThreeVector(0., 0., ((rib_height_-DF_thickn_)/2.)*((2.*shallowness_)-1.)
                                                +(DF_thickn_/2.)
                                                +(ptp_thickn_/2.)), 
                          ptpd_name, 
                          ptpd_logic,
                          mother_physical, 
                          false, 0, true));

    G4OpticalSurface* coating_rough_surf =
          new G4OpticalSurface("COATING_ROUGH_SURFACE", glisur, ground, dielectric_dielectric, .01);
            // 0.01 is the polish value for glisur model that was measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
            // This is the best reference we have, since both PTP and TPB are the result of an evaporation+deposition process
    new G4LogicalBorderSurface("SURROUNDINGS->COATING", mother_physical, ptpd_physical, coating_rough_surf);
    new G4LogicalBorderSurface("COATING->SURROUNDINGS", ptpd_physical, mother_physical, coating_rough_surf);
    new G4LogicalBorderSurface("COATING->DF", ptpd_physical, df_physical, coating_rough_surf);
    
    return;

  }

    void AngularDistributionTest::ConstructFramedPTPAlone(G4VPhysicalVolume* mother_physical)
  {
    // When a manufacturer produces a filter and gives you the transmission curve,
    // that curve accounts for the whole filter, i.e. the dichroic deposition 
    // plus the dichroic substrate. So there's no need to separate these two volumes
    // here, i.e. you are implementing the whole dichroic filter by assigning the
    // transmission curve over the volume. However, there's an implementation 
    // peculiarity you have to take into account here. If you assign the dichroic
    // transmission curve in the transition LAr_volume->dichroic_volume, then
    // the photons will enter the dichroic_volume or get reflected back into the LAr
    // according to the transmission curve you provided. Now, if a reflection takes
    // place, then everything is fine. The angle of the reflected photon is the same
    // as the angle of the incident photon that is reflected. However, for a transmitted
    // photon, it is flying within the dichroic volume and it will eventually escape it.
    // When trying to escape it, the dichroic_volume needs to have a defined rindex, 
    // otherwise the photon will get killed in the boundary. Now, if the defined rindex
    // is different from that of the LAr, in case of refraction, the photon will get
    // deflected according to snell's law. Thus giving a real transmission curve
    // that is different from that implemented by the dichroic data. There are two
    // ways to solve this. Either place the second detector WITHIN the dichroic volume, 
    // or give the dichroic_volume the same optical properties as the surrounding volume.
    // Let us try the second one, which is the one that could be more easily implemented
    // in the eventual X-ARAPUCA geometry.

    G4double LDP_thickn = world_height_/1000.;

    // PTP COATING
    const G4String ptpd_name = "PTP_LAYER";
    G4Box* ptpd_solid =
      new G4Box(ptpd_name, DF_x_size_/2., DF_y_size_/2., ptp_thickn_/2.);

    G4Material* ptpd_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
    ptpd_mat->SetMaterialPropertiesTable(opticalprops::PTP());

    G4LogicalVolume* ptpd_logic =
      new G4LogicalVolume(ptpd_solid, ptpd_mat, ptpd_name, nullptr, nullptr, nullptr, true);

    G4VisAttributes ptpd_col = nexus::TitaniumGreyAlpha();
    ptpd_col.SetForceSolid(true);
    ptpd_logic->SetVisAttributes(ptpd_col);

    G4double ptp_zpos = ((rib_height_-DF_thickn_)/2.)*((2.*shallowness_)-1.)
                        +(DF_thickn_/2.)
                        +(ptp_thickn_/2.);

    G4VPhysicalVolume* ptpd_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, 
                          G4ThreeVector(0., 0., ptp_zpos), 
                          ptpd_name, 
                          ptpd_logic,
                          mother_physical, 
                          false, 0, true));

    G4OpticalSurface* coating_rough_surf =
          new G4OpticalSurface("COATING_ROUGH_SURFACE", glisur, ground, dielectric_dielectric, .01);
            // 0.01 is the polish value for glisur model that was measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
            // This is the best reference we have, since both PTP and TPB are the result of an evaporation+deposition process
    new G4LogicalBorderSurface("SURROUNDINGS->COATING", mother_physical, ptpd_physical, coating_rough_surf);
    new G4LogicalBorderSurface("COATING->SURROUNDINGS", ptpd_physical, mother_physical, coating_rough_surf);

    // FAKE AFTERDICHROIC VOLUME
    const G4String afterdichroic_name = "AFTER_DICHROIC";

    G4double afterdichroic_thickn = (world_height_/2.)
                                    -LDP_thickn
                                    +ptp_zpos
                                    -ptp_thickn_/2.;

    G4Tubs* afterdichroic_solid = 
      new G4Tubs(afterdichroic_name, 0., world_radius_, afterdichroic_thickn/2., 0., 360.*degree);

    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    lAr->SetMaterialPropertiesTable(sm_mpt_);

    G4LogicalVolume* afterdichroic_logic =
      new G4LogicalVolume(afterdichroic_solid, lAr, afterdichroic_name, nullptr, nullptr, nullptr, true);

    /*
    G4VisAttributes afterdichroic_col = nexus::Yellow();
    afterdichroic_col.SetForceSolid(true);
    afterdichroic_logic->SetVisAttributes(afterdichroic_col);
    */

    G4VPhysicalVolume* afterdichroic_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, 
                          G4ThreeVector(0., 0., -(world_height_/2.)
                                                +LDP_thickn
                                                +(afterdichroic_thickn/2.)), 
                          afterdichroic_name, 
                          afterdichroic_logic,
                          mother_physical, 
                          false, 0, true));

    // Add dichroic specifications
    if(filter_code_==1){
      setenv("G4DICHROICDATA", "/home/jurenag/installations/nexus-source/data/full_transmission_test_dichroic_data", 1);
    }
    else if(filter_code_==2){
      setenv("G4DICHROICDATA", "/home/jurenag/cernbox/work_stuff/documentation_and_notes/materiales/filters_data/opto/ific_meas/g4_dichroic_data_from_PTP", 1);
    }
    else{
      G4Exception("AngularDistributionTest::ConstructFramedPTPAlone()", "0",
            FatalException, "Not recognized filter code.");
    }
    
    G4OpticalSurface* dichroic_opsurf =
        new G4OpticalSurface("DICHROIC_OPSURF", dichroic, polished, dielectric_dichroic);
    
    G4LogicalBorderSurface* first_dichroic_lbs = new G4LogicalBorderSurface("PTP->AFTERDICHROIC", ptpd_physical, afterdichroic_physical, dichroic_opsurf);
    //G4LogicalBorderSurface* second_dichroic_lbs = new G4LogicalBorderSurface("PTP->AFTERDICHROIC", ptpd_physical, afterdichroic_physical, dichroic_opsurf);
    // The second one is actually useless in this geometry


    // DF FRAME ------------------------------------------------------------------------------
    G4double geometric_tolerance = 1.*mm; // To avoid subtraction of solids with matching surfaces
    G4Box* aux_frame_inner =
    new G4Box("AUX_FRAME_INNER", DF_x_size_/2., DF_y_size_/2., (rib_height_/2.)+geometric_tolerance);

    
    G4double upper_frame_thickn = 0.;
    // DF UPPER FRAME ------------------------------------------------------------------------
    if(shallowness_<1.){ // If shallowness_==1., then the upper frame is suppresed
      upper_frame_thickn = rib_height_/2. 
                                    -ptp_zpos
                                    +(ptp_thickn_/2.);

      G4Box* aux_upper_frame_outer =
      new G4Box("AUX_UPPER_FRAME_OUTER", (DF_x_size_/2.)+rib_height_, // This is the part of
                                  (DF_y_size_/2.)+rib_height_,        // frame that is outside
                                  upper_frame_thickn/2.);             // the after-dichroic volume

      const G4String upper_frame_name = "UPPER_FRAME";
      G4SubtractionSolid* upper_frame_solid = new G4SubtractionSolid( upper_frame_name, 
                                                                      aux_upper_frame_outer, 
                                                                      aux_frame_inner, 
                                                                      nullptr, 
                                                                      G4ThreeVector(0., 0., 0.));

      G4LogicalVolume* upper_frame_logic = new G4LogicalVolume( upper_frame_solid, 
                                                                materials::FR4(), 
                                                                upper_frame_name);
                                                                                          
      new G4PVPlacement(nullptr, G4ThreeVector( 0., 0., (rib_height_/2.)
                                                        -(upper_frame_thickn/2.)), 
                                                upper_frame_name, 
                                                upper_frame_logic,
                                                mother_physical, 
                                                false, 
                                                0, 
                                                true); // Last one is pSurfChk
    }
    // ---------------------------------------------------------------------------------------

    // DF LOWER FRAME ------------------------------------------------------------------------
    G4double lower_frame_thickn = rib_height_-upper_frame_thickn;
    G4Box* aux_lower_frame_outer =
    new G4Box("AUX_LOWER_FRAME_OUTER", (DF_x_size_/2.)+rib_height_, // This is the part of
                                (DF_y_size_/2.)+rib_height_,        // frame that is immersed
                                lower_frame_thickn/2.);             // in the after-dichroic volume

    const G4String lower_frame_name = "LOWER_FRAME";
    G4SubtractionSolid* lower_frame_solid = new G4SubtractionSolid( lower_frame_name, 
                                                                    aux_lower_frame_outer, 
                                                                    aux_frame_inner, 
                                                                    nullptr, 
                                                                    G4ThreeVector(0., 0., 0.));

    G4LogicalVolume* lower_frame_logic = new G4LogicalVolume( lower_frame_solid, 
                                                              materials::FR4(), 
                                                              lower_frame_name);
                                                                                         
    new G4PVPlacement(nullptr, G4ThreeVector( 0., 0., (afterdichroic_thickn/2)
                                                      -(lower_frame_thickn/2.)),  // Wrt to after-dichroic
                                                                                  // volume coordinates
                                                                                  // system!
                                              lower_frame_name, 
                                              lower_frame_logic,
                                              afterdichroic_physical, // Lower frame is placed within
                                              false,                  // the after-dichroic volume
                                              0, 
                                              true); // Last one is pSurfChk
    // ---------------------------------------------------------------------------------------

    return;

  }

  void AngularDistributionTest::ConstructLDP(G4VPhysicalVolume* mother_physical)
  {
    // Light Detecting Plane (Tightly fitted at the bottom of the cylinder-world)

    G4double LDP_thickn = world_height_/1000.;

    G4String LDP_name = "LDP";
    G4Tubs* LDP_solid = 
      new G4Tubs(LDP_name, 0., world_radius_, LDP_thickn/2., 0., 360.*degree);

    G4Material* LDP_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

    G4LogicalVolume* LDP_logic = new G4LogicalVolume(LDP_solid, LDP_mat, LDP_name);

    /*
    G4VisAttributes LDP_col = nexus::LightGreen();
    LDP_col.SetForceSolid(true);
    LDP_logic->SetVisAttributes(LDP_col);
    */

    G4OpticalSurface* LDP_surface = new G4OpticalSurface("LDP_coating", glisur, polished, dielectric_metal);
    LDP_surface->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());
    new G4LogicalSkinSurface("LDP_coating", LDP_logic, LDP_surface);
    
    new G4PVPlacement (new G4RotationMatrix(), G4ThreeVector(0., 0., -1.*(world_height_/2.)+(LDP_thickn/2.)), LDP_name, LDP_logic, 
                    mother_physical, false, 0, true);

    return;

  }

  G4bool AngularDistributionTest::geometry_is_ill_formed(){
    // By inner geometry I mean the geometry that is mean to fit inside the cylinder world
    G4double inner_geom_height_span;      // By height, I mean the direction of the world-cylinder height
    G4double inner_geom_transverse_span;  // By transverse, I mean the direction that is perpendicular to the world-cylinder height
    if(config_code_==1){
      inner_geom_height_span = DF_thickn_+ptp_thickn_+(world_height_/1000.);
      inner_geom_transverse_span = 2*world_radius_;
    }
    else if(config_code_==2){
      G4double aux = std::max(rib_height_, DF_thickn_+ptp_thickn_);
      inner_geom_height_span = aux+(world_height_/1000.); // The last term is the thickness of the LDP
      inner_geom_transverse_span = sqrt(  ((DF_x_size_+(2.*rib_height_))*(DF_x_size_+(2.*rib_height_)))
                                          +((DF_y_size_+(2.*rib_height_))*(DF_y_size_+(2.*rib_height_))));
    }
    else if(config_code_==3){
      G4double aux = std::max(rib_height_, ptp_thickn_);
      inner_geom_height_span = aux+(world_height_/1000.);
      inner_geom_transverse_span = sqrt(  ((DF_x_size_+(2.*rib_height_))*(DF_x_size_+(2.*rib_height_)))
                                          +((DF_y_size_+(2.*rib_height_))*(DF_y_size_+(2.*rib_height_))));
    }
    else{
      G4Exception("AngularDistributionTest::ConstructSample()", "0",
                  FatalException, "Not recognized configuration code.");
    }
    G4bool c1 = inner_geom_height_span <= world_height_;
    G4bool c2 = inner_geom_transverse_span <= 2*world_radius_;
    if(!(c1 && c2)) return true;
    else            return false;
  }

  G4ThreeVector AngularDistributionTest::GenerateVertex(const G4String&) const
  {
    if(config_code_==1){
      return G4ThreeVector(0., 0., lgv_zpos_wrt_world_center_);
    }
    else if(config_code_==2 || config_code_==3){
      if(randomly_sample_lgv_){
        G4double tolerance = 0.1*mm;
        G4double x_pos = UniformRandomInRange((DF_x_size_/2.)-tolerance,
                                              -(DF_x_size_/2.)+tolerance);
        G4double y_pos = UniformRandomInRange((DF_y_size_/2.)-tolerance,
                                              -(DF_y_size_/2.)+tolerance);
        return G4ThreeVector(x_pos, y_pos, lgv_zpos_wrt_world_center_);
      }
      else{
        return G4ThreeVector(0., 0., lgv_zpos_wrt_world_center_);
      }
    }
    else{
      G4Exception("AngularDistributionTest::ConstructSample()", "0",
                  FatalException, "Not recognized configuration code.");
    }
    return G4ThreeVector(0., 0., 0.);
  }

}