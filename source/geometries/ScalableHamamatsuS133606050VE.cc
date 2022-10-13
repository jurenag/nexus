#include "ScalableHamamatsuS133606050VE.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "SensorSD.h"

#include <G4Box.hh>
#include <G4SubtractionSolid.hh>
#include <G4RotationMatrix.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SystemOfUnits.hh>
#include <G4SDManager.hh>

namespace nexus {

  ScalableHamamatsuS133606050VE::ScalableHamamatsuS133606050VE(G4double height, G4double width):
    GeometryBase(),
    height_(height),
    width_(width),
    visibility_(true), 
    reflective_support_(false)
  {
  }

  ScalableHamamatsuS133606050VE::~ScalableHamamatsuS133606050VE()
  {
  }

  //Define static attributes
  const G4double ScalableHamamatsuS133606050VE::thickness_           = 1.35*mm;
  const G4double ScalableHamamatsuS133606050VE::window_thickness_    = 0.1*mm;
  const G4double ScalableHamamatsuS133606050VE::sensarea_thickness_  = 0.1*mm;
 
  //Define static methods
  G4double ScalableHamamatsuS133606050VE::GetThickness() { return thickness_; }
  

  void ScalableHamamatsuS133606050VE::Construct()
  {
    /// PHOTOSENSOR ENCASING

    G4String name = "SCALABLE_HAMAMATSU_S13360_6050VE_MPPC";

    G4Box* encasing_solid =
      new G4Box(name, height_/2., thickness_/2., width_/2.);

    G4LogicalVolume* encasing_logic =
      new G4LogicalVolume(encasing_solid, materials::FR4(), name);

    this->SetLogicalVolume(encasing_logic);

    /// FR4 SUPPORT

    name = "AUX";
    G4double aux = thickness_ -window_thickness_;
    G4Box* aux_solid =
      new G4Box("AUX", height_/2., aux/2., width_/2.);

    G4double tolerance = 0.1*mm;

    name = "MPPC_SENSAREA";
    G4Box* sensarea_solid =
      new G4Box(name, (height_/2.)-tolerance, sensarea_thickness_/2., (width_/2.)-tolerance);
      //Tolerance for preventing matching surfaces in subtraction solid

    name = "MPPC_SUPPORT";
    G4SubtractionSolid* support_solid = 
        new G4SubtractionSolid(name, aux_solid, sensarea_solid, 
                            new G4RotationMatrix{},
                            G4ThreeVector(0., (aux/2.) -(sensarea_thickness_/2.), 0.));

    G4LogicalVolume* support_logic =
      new G4LogicalVolume(support_solid, materials::FR4(), name);

    if(reflective_support_){
        //VIKUITI coating for the support
        const G4String sc_name = "SUPPORT_COATING";
        G4OpticalSurface* support_coating = 
        new G4OpticalSurface(sc_name, unified, ground, dielectric_metal, 1);

        support_coating->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());
        new G4LogicalSkinSurface(sc_name, support_logic, support_coating); 
    }

    new G4PVPlacement(nullptr, G4ThreeVector(0., (-1.*thickness_/2.) +(aux/2.), 0.),
                      support_logic, name, encasing_logic, false, 0, true);

    /// OPTICAL WINDOW

    name = "MPPC_WINDOW";

    G4Material* op_sil = materials::OpticalSilicone();
    op_sil->SetMaterialPropertiesTable(opticalprops::HamamatsuEpoxy());

    G4Box* window_solid =
      new G4Box(name, height_/2., window_thickness_/2., width_/2.);

    G4LogicalVolume* window_logic =
      new G4LogicalVolume(window_solid, op_sil, name);

    new G4PVPlacement(nullptr, G4ThreeVector(0., thickness_/2. -window_thickness_/2., 0.),
                      window_logic, name, encasing_logic, false, 0, true);

    /// PHOTOSENSITIVE AREA

    name = "MPPC_SENSAREA";
    G4LogicalVolume* sensarea_logic =
      new G4LogicalVolume(sensarea_solid,
                          G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"),
                          name);

    new G4PVPlacement(nullptr, G4ThreeVector(0., thickness_/2. -window_thickness_ -sensarea_thickness_/2., 0.),          
                      sensarea_logic, name, encasing_logic, false, 0, true);

    /// OPTICAL PROPERTIES

    name = "MPPC_OPSURF";

    G4double energy[]       = {opticalprops::optPhotMinE_, 1.384*eV, 1.4003*eV, 1.4255*eV, 1.4516*eV, 
                              1.4787*eV, 1.5068*eV, 1.536*eV, 1.5663*eV, 1.5979*eV, 1.6279*eV, 
                              1.6563*eV, 1.6856*eV, 1.716*eV, 1.7475*eV, 1.7802*eV, 1.8107*eV, 
                              1.8387*eV, 1.8676*eV, 1.8974*eV, 1.9281*eV, 1.9559*eV, 1.9803*eV, 
                              2.0053*eV, 2.031*eV, 2.0574*eV, 2.0844*eV, 2.1122*eV, 2.1407*eV, 
                              2.17*eV, 2.2001*eV, 2.2311*eV, 2.2629*eV, 2.2956*eV, 2.3293*eV, 
                              2.3641*eV, 2.4059*eV, 2.462*eV, 2.5341*eV, 2.6177*eV, 2.7071*eV, 
                              2.8028*eV, 2.8966*eV, 2.969*eV, 3.0256*eV, 3.0845*eV, 3.1353*eV, 
                              3.1771*eV, 3.2201*eV, 3.2644*eV, 3.3098*eV, 3.3565*eV, 3.3924*eV, 
                              3.4291*eV, 3.4792*eV, 3.5049*eV, 3.5308*eV, 3.5572*eV, 3.5841*eV, 
                              3.5976*eV, 3.6112*eV, 3.625*eV, 3.6389*eV, 3.6528*eV, 3.6669*eV, 
                              3.6907*eV, 3.7195*eV, 3.7488*eV, 3.7836*eV, 3.8138*eV, 3.8447*eV, 
                              opticalprops::optPhotMaxE_};

    G4double reflectivity[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    G4double efficiency[]   = {0.0, 0.0367, 0.0408, 0.0474, 0.0544, 0.0623, 0.07, 0.0783, 0.0871, 
                              0.0961, 0.1056, 0.1151, 0.1248, 0.1348, 0.145, 0.1554, 0.166, 0.1758, 
                              0.1855, 0.1961, 0.2073, 0.2197, 0.2321, 0.2424, 0.2519, 0.2606, 0.2709, 
                              0.2809, 0.2912, 0.3018, 0.3124, 0.3224, 0.3333, 0.3433, 0.3536, 0.3628, 
                              0.3738, 0.3849, 0.393, 0.399, 0.4017, 0.3951, 0.3856, 0.3758, 0.3661, 
                              0.3563, 0.3465, 0.3352, 0.3217, 0.3063, 0.2904, 0.2734, 0.2612, 0.2498, 
                              0.2319, 0.2201, 0.2103, 0.1993, 0.1859, 0.1741, 0.1656, 0.1546, 0.1461, 
                              0.1351, 0.127, 0.1107, 0.0928, 0.0752, 0.0586, 0.044, 0.0318, 0.0};

    G4double test_energy[]  = {opticalprops::optPhotMinE_, opticalprops::optPhotMaxE_};
    G4double test_efficiency[]    = {1., 1.};

    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();
    photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, 72);
    photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   72);

    G4OpticalSurface* photosensor_opsurf =
      new G4OpticalSurface(name, unified, polished, dielectric_metal);

    photosensor_opsurf->SetMaterialPropertiesTable(photosensor_mpt);
    new G4LogicalSkinSurface(name, sensarea_logic, photosensor_opsurf); 
   
    /// SENSITIVE DETECTOR
    G4SDManager* SDMptr = G4SDManager::GetSDMpointer();
    SensorSD* sensdet = dynamic_cast<SensorSD*>(SDMptr->FindSensitiveDetector("/PHOTOSENSOR/SHS133606050VE"));
    if(!sensdet){
        sensdet = new SensorSD("/PHOTOSENSOR/SHS133606050VE");
        sensdet->SetDetectorVolumeDepth(1);
        sensdet->SetTimeBinning(1.*s);
        G4SDManager::GetSDMpointer()->AddNewDetector(sensdet);
        sensarea_logic->SetSensitiveDetector(sensdet);
    }
    else{
        sensarea_logic->SetSensitiveDetector(sensdet);
    }

    return;
  }
}