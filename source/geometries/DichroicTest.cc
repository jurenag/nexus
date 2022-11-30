#include "DichroicTest.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "LArSphere.h"
#include "SensorSD.h"
#include "Visibilities.h"

#include <G4Tubs.hh>
#include <G4Orb.hh>
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

REGISTER_CLASS(DichroicTest, GeometryBase)

namespace nexus{

  DichroicTest::DichroicTest():
  GeometryBase(),
  radius_(1.*m),
  thickn_(1.*mm),
  det_thickn_(5.*mm),
  det_depth_(2.*cm),
  with_ptp_deposition_(false),
  ptp_thickn_(3.*um),
  mpt_(opticalprops::SCHOTT_B270())
  {

    msg_ = new G4GenericMessenger(this, "/Geometry/DichroicTest/",
				"Control commands of geometry DichroicTest.");
    
    G4GenericMessenger::Command& dfr_cmd =
      msg_->DeclareProperty("dichroic_radius", radius_,
			    "Radius of the dichroic filter.");
    dfr_cmd.SetUnitCategory("Length");
    dfr_cmd.SetParameterName("dichroic_radius", false);
    dfr_cmd.SetRange("dichroic_radius>0.");

    G4GenericMessenger::Command& dft_cmd =
      msg_->DeclareProperty("dichroic_thickness", thickn_,
			    "Thickness of the dichroic filter.");
    dft_cmd.SetUnitCategory("Length");
    dft_cmd.SetParameterName("dichroic_thickness", false);
    dft_cmd.SetRange("dichroic_thickness>0.");

    G4GenericMessenger::Command& dett_cmd =
      msg_->DeclareProperty("detector_thickness", det_thickn_,
			    "Thickness of the detectors.");
    dett_cmd.SetUnitCategory("Length");
    dett_cmd.SetParameterName("det_thickness", false);
    dett_cmd.SetRange("det_thickness>0.");

    G4GenericMessenger::Command& detd_cmd =
      msg_->DeclareProperty("detector_depth", det_depth_,
			    "Depth of the detectors.");
    detd_cmd.SetUnitCategory("Length");
    detd_cmd.SetParameterName("det_depth", false);
    detd_cmd.SetRange("det_depth>0.");

    G4GenericMessenger::Command& wpd_cmd =
      msg_->DeclareProperty("with_ptp_deposition", with_ptp_deposition_,
			    "Whether the dichroic filter has a PTP coating or not.");

    G4GenericMessenger::Command& ptpt_cmd =
      msg_->DeclareProperty("ptp_thickn", ptp_thickn_,
			    "Thickness of PTP deposition.");
    dett_cmd.SetUnitCategory("Length");
    dett_cmd.SetParameterName("ptp_thickn", false);
    dett_cmd.SetRange("ptp_thickn>0.");

  }

  DichroicTest::~DichroicTest()
  {
    if(msg_) delete msg_;
  }

  void DichroicTest::Construct()
  {
    // VACUUM CAPSULE
    const G4String world_name = "VACUUM_CAPSULE";
    G4Orb* world_solid =
        new G4Orb(world_name, 2*radius_ +thickn_ +(10.*m));
    G4Material* vacuum =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4LogicalVolume* world_logic = 
        new G4LogicalVolume(world_solid, vacuum, world_name, 0, 0, 0, true);
    world_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(world_logic);

    // LAR SPHERE
    G4String sphere_name = "LAR_SPHERE";
    G4Orb* sphere_solid = new G4Orb(sphere_name, 2*radius_ +thickn_ +(5.*m));
    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    lAr->SetMaterialPropertiesTable(opticalprops::LAr());
    G4LogicalVolume* sphere_logic = new G4LogicalVolume(sphere_solid, lAr, sphere_name);
    sphere_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    // Place the LAr sphere within the vacuum world
    G4VPhysicalVolume* mother_physical = 
        dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement (new G4RotationMatrix(), G4ThreeVector(0., 0., 0.), sphere_logic, 
                        sphere_name, world_logic, false, 0, false));

    ConstructDetectors(mother_physical);
    ConstructDichroicFilter(mother_physical);
    
    return;

  }

  void DichroicTest::ConstructDichroicFilter(G4VPhysicalVolume* mother_physical)
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
    // in the eventual X-ARAPUCA geometry. If, in this second-solution case, you wanted
    // to implement the absorption length of the dichroic substrate but still implementing
    // the actual transmission curve, you could devise a G4MaterialPropertiesTable, 
    // something like opticalprops::FakeDichroicFilter, which has the same rindex as 
    // the surrounding LAr, but implements the actual absorption length of the 
    // dichroic filter substrate.

    const G4String df_name = "DICHROIC_SUBSTRATE";

    G4Tubs* df_solid =
      new G4Tubs(df_name, 0., radius_, thickn_/2., 0., 360.*degree);

    G4Material* substrate_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");

    if(mpt_){
      substrate_mat->SetMaterialPropertiesTable(mpt_);
    }
    else{
      G4Exception("[DichroicTest]", "ConstructDichroicFilter()", JustWarning,
      "The optical properties of the dichroic filter substrate were not set");
    }

    G4LogicalVolume* df_logic =
      new G4LogicalVolume(df_solid, substrate_mat, df_name);

    G4VisAttributes df_col = nexus::RedAlpha();
    df_col.SetForceSolid(true);
    df_logic->SetVisAttributes(df_col);

    G4VPhysicalVolume* df_physical = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -1.*thickn_/2.), df_name, df_logic,
        mother_physical, false, 0, false));


    G4VPhysicalVolume* ptpd_physical = nullptr;
    if(with_ptp_deposition_){
      const G4String ptpd_name = "PTP_LAYER";
      G4Tubs* ptpd_solid =
        new G4Tubs(ptpd_name, 0., radius_, ptp_thickn_/2., 0., 360.*degree);

      G4Material* ptpd_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TERPHENYL");
      ptpd_mat->SetMaterialPropertiesTable(opticalprops::PTP());

      G4LogicalVolume* ptpd_logic =
        new G4LogicalVolume(ptpd_solid, ptpd_mat, ptpd_name);

      G4VisAttributes ptpd_col = nexus::TitaniumGreyAlpha();
      ptpd_col.SetForceSolid(true);
      ptpd_logic->SetVisAttributes(ptpd_col);

      ptpd_physical = dynamic_cast<G4VPhysicalVolume*>(
          new G4PVPlacement(nullptr, G4ThreeVector(0., 0., +ptp_thickn_/2.), ptpd_name, ptpd_logic,
          mother_physical, false, 0, false));

      G4OpticalSurface* coating_rough_surf =
            new G4OpticalSurface("COATING_ROUGH_SURFACE", glisur, ground, dielectric_dielectric, .01);
            // 0.01 is the polish value for glisur model that was measured for TPB in doi.org/10.1140/epjc/s10052-018-5807-z
            // This is the best reference we have, since both PTP and TPB are the result of an evaporation+deposition process
      new G4LogicalBorderSurface("SURROUNDINGS->COATING", mother_physical, ptpd_physical, coating_rough_surf);
      new G4LogicalBorderSurface("COATING->SURROUNDINGS", ptpd_physical, mother_physical, coating_rough_surf);
      //new G4LogicalBorderSurface("COATING->DF", ptpd_physical, df_physical, coating_rough_surf);
    }

    // Add dichroic specifications
    setenv("G4DICHROICDATA", "data/full_transmission_test_dichroic_data", 1);
    G4OpticalSurface* dichroic_opsurf =
        new G4OpticalSurface("DICHROIC_OPSURF", dichroic, polished, dielectric_dichroic);
    
    G4LogicalBorderSurface* first_lbs = nullptr;   
    if(!with_ptp_deposition_) first_lbs = new G4LogicalBorderSurface("LAr->DICHROIC1", mother_physical, df_physical, dichroic_opsurf);
    else                      first_lbs = new G4LogicalBorderSurface("LAr->DICHROIC1", ptpd_physical, df_physical, dichroic_opsurf);
    
    return;
    

    /* CODE FOR ADDING A SECOND DF
    G4Material* substrate_mat_2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    if(mpt_){
      substrate_mat_2->SetMaterialPropertiesTable(mpt_);
    }
    else{
      G4Exception("[DichroicTest]", "ConstructDichroicFilter()", JustWarning,
      "The optical properties of the dichroic filter substrate were not set");
    }
    
    G4LogicalVolume* df_logic_2 =
      new G4LogicalVolume(df_solid, substrate_mat_2, df_name);
    
    G4VPhysicalVolume* df_physical_2 = dynamic_cast<G4VPhysicalVolume*>(
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.5*cm +thickn_/2.), df_name, df_logic_2,
        mother_physical, false, 0, false));
    
    //G4cout << getenv("G4DICHROICDATA") << G4endl;
    //unsetenv("G4DICHROICDATA");
    setenv("G4DICHROICDATA", "data/full_reflection_test_dichroic_data", 1);
    G4OpticalSurface* dichroic_opsurf_2 =
        new G4OpticalSurface("DICHROIC_OPSURF_2", dichroic, polished, dielectric_dichroic);
    G4LogicalBorderSurface* second_lbs = 
                new G4LogicalBorderSurface("LAr->DICHROIC2", mother_physical, df_physical_2, dichroic_opsurf_2);
    */

    // *Bug1: Digamos que si el primer fotón que interactúa con alguna de las G4LogicalBorderSurface
    // interactúa con la i-ésima G4LogicalBorderSurface. Entonces, todas las G4LogicalBorderSurface
    // adquieren las curvas de transmisión de dicha G4LogicalBorderSurface.
    
    return;

  }

  void DichroicTest::ConstructDetectors(G4VPhysicalVolume* mother_physical)
  {

    const G4String det_name = "DETECTOR";

    G4Tubs* aux_outter_det_solid =
      new G4Tubs("AUX", 0., radius_+det_thickn_, det_depth_/2., 0., 360.*degree);

    G4double tolerance = 1.*mm;

    G4Tubs* aux_inner_det_solid = 
      new G4Tubs("AUX", 0., radius_, (det_depth_-det_thickn_)/2. +tolerance/2., 0., 360.*degree);

    G4SubtractionSolid* det_solid = 
        new G4SubtractionSolid(det_name, aux_outter_det_solid, aux_inner_det_solid,
                            nullptr, G4ThreeVector(0., 0., det_thickn_/2. +tolerance/2.));

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    pvt->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());

    G4LogicalVolume* det_logic =
      new G4LogicalVolume(det_solid, pvt, det_name);

    G4VisAttributes det_col = nexus::WhiteAlpha();
    det_col.SetForceSolid(true);
    det_logic->SetVisAttributes(det_col);

    // Detection coating
    G4OpticalSurface* det_opsurf =
      new G4OpticalSurface("DETECTOR_OPSURF", unified, polished, dielectric_metal);
    det_opsurf->SetMaterialPropertiesTable(opticalprops::PerfectPhotonCollector());
    new G4LogicalSkinSurface("DETECTOR_OPSURF", det_logic, det_opsurf);

    // Sensitive detector
    SensorSD* sensdet = new SensorSD("/COLLECTOR");
    sensdet->SetDetectorVolumeDepth(0);
    sensdet->SetTimeBinning(1.*s);
    G4SDManager::GetSDMpointer()->AddNewDetector(sensdet);
    det_logic->SetSensitiveDetector(sensdet);

    // Placement of both detectors
    // Detectors and dichroic filter are placed so that final_z of a photon is <0 if the photon
    // has been reflected, and >0 if it has been transmited
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -det_depth_/2.), "R_DETECTOR", det_logic,
        mother_physical, false, 0, true);

    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateX(180.*degree);
    new G4PVPlacement(rot, G4ThreeVector(0., 0., +det_depth_/2), "T_DETECTOR", det_logic,
        mother_physical, false, 1, true);


    return;

  }

  G4ThreeVector DichroicTest::GenerateVertex(const G4String&) const
  {

    G4double distance_from_filter = 0.5*cm;
    if(with_ptp_deposition_) distance_from_filter += ptp_thickn_;
    return G4ThreeVector(0., 0., distance_from_filter);
  }

}