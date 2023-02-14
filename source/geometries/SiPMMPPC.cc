#include "SiPMMPPC.h"

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

  SiPMMPPC::SiPMMPPC():
  GeometryBase(), 
  visibility_(true), 
  reflective_support_(false)
  {
  }

  SiPMMPPC::~SiPMMPPC()
  {
  }

  void SiPMMPPC::Construct()
  {
    /// PHOTOSENSOR ENCASING

    G4Box* encasing_solid =
      new G4Box(  this->GetModel(), 
                  (this->GetTransverseDim())/2., 
                  (this->GetThickness())/2., 
                  (this->GetTransverseDim())/2.);

    G4LogicalVolume* encasing_logic =
      new G4LogicalVolume(  encasing_solid, 
                            materials::FR4(), 
                            this->GetModel());

    this->SetLogicalVolume(encasing_logic);

    /// FR4 SUPPORT

    G4double aux =  this->GetThickness() 
                    -this->GetWindowThickness();
    G4Box* aux_solid =
      new G4Box(  "AUX", 
                  (this->GetTransverseDim())/2., 
                  aux/2., 
                  (this->GetTransverseDim())/2.);

    G4Box* sensarea_solid =
      new G4Box(  "MPPC_SENSAREA", 
                  (this->GetSensareaTransverseDim())/2., 
                  (this->GetSensareaThickness())/2., 
                  (this->GetSensareaTransverseDim())/2.);

    G4SubtractionSolid* support_solid = 
        new G4SubtractionSolid("MPPC_SUPPORT", aux_solid, sensarea_solid, 
                            new G4RotationMatrix{},
                            G4ThreeVector(  0., 
                                            (aux/2.) -((this->GetSensareaThickness())/2.), 
                                            0.));

    G4LogicalVolume* support_logic =
      new G4LogicalVolume(  support_solid, 
                            materials::FR4(), 
                            "MPPC_SUPPORT");

    if(reflective_support_){
        //VIKUITI coating for the support
        G4OpticalSurface* support_coating = 
        new G4OpticalSurface( "SUPPORT_COATING", 
                              unified, 
                              ground, 
                              dielectric_metal, 
                              1);

        support_coating->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());
        new G4LogicalSkinSurface( "SUPPORT_COATING", 
                                  support_logic, 
                                  support_coating); 
    }

    new G4PVPlacement(nullptr,
                      G4ThreeVector(  0., 
                                      (-1.*(this->GetThickness())/2.) +(aux/2.), 
                                      0.),
                      support_logic, 
                      this->GetModel(), 
                      encasing_logic, 
                      false, 
                      0, 
                      true);

    /// OPTICAL WINDOW
    G4Material* op_sil = materials::OpticalSilicone();
    op_sil->SetMaterialPropertiesTable(opticalprops::HamamatsuEpoxy());

    G4Box* window_solid =
      new G4Box(  "MPPC_WINDOW", 
                  (this->GetTransverseDim())/2., 
                  (this->GetWindowThickness())/2., 
                  (this->GetTransverseDim())/2.);

    G4LogicalVolume* window_logic =
      new G4LogicalVolume(  window_solid, 
                            op_sil, 
                            "MPPC_WINDOW");

    new G4PVPlacement(  nullptr, 
                        G4ThreeVector(  0., 
                                        (this->GetThickness())/2. -(this->GetWindowThickness())/2., 
                                         0.),
                        window_logic, 
                        "MPPC_WINDOW", 
                        encasing_logic, 
                        false, 
                        0, 
                        true);

    /// PHOTOSENSITIVE AREA
    G4LogicalVolume* sensarea_logic =
      new G4LogicalVolume(sensarea_solid,
                          G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"),
                          "MPPC_SENSAREA");

    new G4PVPlacement(  nullptr, 
                        G4ThreeVector(  0., 
                                        ((this->GetThickness())/2.) -(this->GetWindowThickness()) -((this->GetSensareaThickness())/2.), 
                                        0.),          
                        sensarea_logic, 
                        "MPPC_SENSAREA", 
                        encasing_logic, 
                        false, 
                        0, 
                        true);

    /// OPTICAL PROPERTIES
    std::pair<G4int, G4double*> sensarea_energy = GetSensareaEnergyArray();
    G4int energy_size = sensarea_energy.first;
    G4double* energy = sensarea_energy.second;

    std::pair<G4int, G4double*> sensarea_efficiency = GetSensareaEfficiencyArray();
    G4int efficiency_size = sensarea_efficiency.first;
    G4double* efficiency = sensarea_efficiency.second;
  
    std::pair<G4int, G4double*> sensarea_reflectivity = GetSensareaReflectivityArray();
    G4int reflectivity_size = sensarea_reflectivity.first;
    G4double* reflectivity = sensarea_reflectivity.second;

    if(energy_size!=efficiency_size || efficiency_size!=reflectivity_size){
      G4Exception("[SiPMMPPC]", "Construct()",
                  FatalException, "Sizes of energy, effciency and reflectivity arrays must match.");
    }

    G4double test_energy[]  = {opticalprops::optPhotMinE_, opticalprops::optPhotMaxE_};
    G4double test_efficiency[]    = {1., 1.};

    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();
    photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, energy_size);
    photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   energy_size);

    G4OpticalSurface* photosensor_opsurf =
      new G4OpticalSurface("MPPC_OPSURF", unified, polished, dielectric_metal);

    photosensor_opsurf->SetMaterialPropertiesTable(photosensor_mpt);
    new G4LogicalSkinSurface("MPPC_OPSURF", sensarea_logic, photosensor_opsurf); 
   
    /// SENSITIVE DETECTOR
    G4SDManager* SDMptr = G4SDManager::GetSDMpointer();
    SensorSD* sensdet = dynamic_cast<SensorSD*>(SDMptr->FindSensitiveDetector("/PHOTOSENSOR/"+(this->GetModel())));
    if(!sensdet){
        sensdet = new SensorSD("/PHOTOSENSOR/"+(this->GetModel()));
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