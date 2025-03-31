#include "LArScintillationGenerator.h"

#include "DetectorConstruction.h"
#include "GeometryBase.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include <G4OpticalPhoton.hh>
#include <G4RunManager.hh>
#include <G4PrimaryVertex.hh>
#include <G4Event.hh>
#include <G4RandomTools.hh>

#include "CLHEP/Units/SystemOfUnits.h"

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(LArScintillationGenerator, G4VPrimaryGenerator)

LArScintillationGenerator::LArScintillationGenerator():
G4VPrimaryGenerator(),
msg_(0),
geom_(0), 
lambertian_(true),
pn_x_(0.), pn_y_(1.), pn_z_(0.),
region_(""),
rd_{},
gen_(rd_()),
bin_edges_{},
emission_spectrum_{},
sampler_(0)
{
  
  LoadNCheckLArData();
  sampler_ = new std::piecewise_constant_distribution(
      bin_edges_.begin(), bin_edges_.end(), emission_spectrum_.begin()
        );

  msg_ = new G4GenericMessenger(
    this,
    "/Generator/LArScintillation/",
    "Control commands of LAr scintillation generator."
  );

  msg_->DeclareProperty(
    "lambertian",
    lambertian_,
    "Whether to use a lambertian emitter or a collimated photon emitter."
  );

  msg_->DeclareProperty(
    "pn_x",
    pn_x_,
    "X coordinate of the lambertian-emitter plane normal vector or the photon direction, depending on the value given to the 'lambertian_' attribute."
  );
  
  msg_->DeclareProperty(
    "pn_y",
    pn_y_,
    "Y coordinate of the lambertian-emitter plane normal vector or the photon direction, depending on the value given to the 'lambertian_' attribute."
  );

  msg_->DeclareProperty(
    "pn_z",
    pn_z_,
    "Z coordinate of the lambertian-emitter plane normal vector or the photon direction, depending on the value given to the 'lambertian_' attribute."
  );

  msg_->DeclareProperty(
    "region",
    region_,
    "Set the region of the geometry where the vertex will be generated."
  );

  DetectorConstruction* detconst = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();

}

LArScintillationGenerator::~LArScintillationGenerator()
{
  if(sampler_)  delete sampler_;
  if(msg_)      delete msg_;
}

void LArScintillationGenerator::GeneratePrimaryVertex(G4Event* event)
{

  G4PrimaryParticle* a_photon = new G4PrimaryParticle(G4OpticalPhoton::Definition());
  G4ThreeVector photon_momentum_dir = (lambertian_ ? G4LambertianRand(G4ThreeVector(pn_x_, pn_y_, pn_z_)) : G4ThreeVector(pn_x_, pn_y_, pn_z_));
  a_photon->SetMomentumDirection(photon_momentum_dir);
  a_photon->SetPolarization(G4PlaneVectorRand(photon_momentum_dir));
  a_photon->SetKineticEnergy(RandomEnergy());
  

  // Generate an initial position for the particle using the geometry
  G4ThreeVector position = geom_->GenerateVertex(region_);

  // Particle generated at start-of-event
  G4double time = 0.;

  // Create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, time);

    // Add particle to the vertex and this to the event
  vertex->SetPrimary(a_photon);
  event->AddPrimaryVertex(vertex);
  return;
}

G4double LArScintillationGenerator::RandomEnergy(){
    if(sampler_){
        return sampler_->operator()<std::mt19937>(gen_);
    }
    else{
        G4Exception("[LArScintillationGenerator]", "RandomEnergy()", JustWarning,
        "PTP histogram sampler was not set. Energy could not be sampled.");
        return -1.;
    }
}

void LArScintillationGenerator::LoadNCheckLArData(){

    // Data from blue curve (LAr excited by proton beam) of figure 4 of DOI: 10.1140/epjc/s10052-013-2618-0
    bin_edges_ = {
      h_Planck * c_light / (312.041 * nm), h_Planck * c_light / (306.538 * nm), h_Planck * c_light / (301.479 * nm), h_Planck * c_light / (296.598 * nm),
      h_Planck * c_light / (291.272 * nm), h_Planck * c_light / (286.657 * nm), h_Planck * c_light / (282.041 * nm), h_Planck * c_light / (275.828 * nm),
      h_Planck * c_light / (269.26 * nm), h_Planck * c_light / (262.87 * nm), h_Planck * c_light / (256.302 * nm), h_Planck * c_light / (249.467 * nm),
      h_Planck * c_light / (242.722 * nm), h_Planck * c_light / (236.509 * nm), h_Planck * c_light / (230.917 * nm), h_Planck * c_light / (225.858 * nm),
      h_Planck * c_light / (221.154 * nm), h_Planck * c_light / (216.272 * nm), h_Planck * c_light / (211.923 * nm), h_Planck * c_light / (208.018 * nm),
      h_Planck * c_light / (203.58 * nm), h_Planck * c_light / (198.254 * nm), h_Planck * c_light / (192.574 * nm), h_Planck * c_light / (187.781 * nm),
      h_Planck * c_light / (183.698 * nm), h_Planck * c_light / (180.68 * nm), h_Planck * c_light / (178.462 * nm), h_Planck * c_light / (175.71 * nm),
      h_Planck * c_light / (172.959 * nm), h_Planck * c_light / (170.917 * nm), h_Planck * c_light / (168.609 * nm), h_Planck * c_light / (165.858 * nm),
      h_Planck * c_light / (163.284 * nm), h_Planck * c_light / (161.509 * nm), h_Planck * c_light / (159.645 * nm), h_Planck * c_light / (157.249 * nm),
      h_Planck * c_light / (155.473 * nm), h_Planck * c_light / (154.32 * nm), h_Planck * c_light / (153.166 * nm), h_Planck * c_light / (152.189 * nm),
      h_Planck * c_light / (151.391 * nm), h_Planck * c_light / (150.68 * nm), h_Planck * c_light / (150.059 * nm), h_Planck * c_light / (149.083 * nm),
      h_Planck * c_light / (148.107 * nm), h_Planck * c_light / (147.574 * nm), h_Planck * c_light / (146.953 * nm), h_Planck * c_light / (145.71 * nm),
      h_Planck * c_light / (144.29 * nm), h_Planck * c_light / (142.959 * nm), h_Planck * c_light / (141.716 * nm), h_Planck * c_light / (140.118 * nm),
      h_Planck * c_light / (138.343 * nm), h_Planck * c_light / (136.834 * nm), h_Planck * c_light / (135.059 * nm), h_Planck * c_light / (133.284 * nm),
      h_Planck * c_light / (131.775 * nm), h_Planck * c_light / (129.911 * nm), h_Planck * c_light / (127.604 * nm), h_Planck * c_light / (125.385 * nm),
      h_Planck * c_light / (123.432 * nm), h_Planck * c_light / (122.101 * nm), h_Planck * c_light / (120.503 * nm), h_Planck * c_light / (118.195 * nm)
    };
    
    emission_spectrum_ = {
      3.7e-06, 4e-06, 5.8e-06, 6.3e-06, 8.9e-06, 1.24e-05, 1.46e-05, 1.84e-05, 1.53e-05, 1.17e-05, 1.03e-05, 9.9e-06, 8.1e-06, 7.8e-06, 7.6e-06, 
      7.9e-06, 1.02e-05, 1.28e-05, 1.67e-05, 1.64e-05, 1.54e-05, 1.39e-05, 1.55e-05, 1.75e-05, 2.29e-05, 2.54e-05, 2.15e-05, 1.9e-05, 1.94e-05, 1.31e-05,
      1.07e-05, 1.21e-05, 1.55e-05, 2.71e-05, 3.01e-05, 2.35e-05, 3.85e-05, 6.87e-05, 0.0001358, 0.0002741, 0.0005886, 0.0014306, 0.003406, 0.0059493, 0.0042758,
      0.0027714, 0.0014314, 0.0020338, 0.0031385, 0.0058327, 0.0106176, 0.0214331, 0.03822, 0.0710305, 0.1375806, 0.2504622, 0.4751756, 0.7802588, 1.0, 0.8837221,
      0.5970604, 0.3950928, 0.231022
    };

    if(DataIsIllFormed()){
        G4Exception("[LArScintillationGenerator]", "LoadNCheckLArData()", FatalException,
        "The provided data is ill-formed. Check LArScintillationGenerator::DataIsIllFormed for more info.");
    }

    return;
}

G4bool LArScintillationGenerator::DataIsIllFormed(){
    if(bin_edges_.size()!=1+emission_spectrum_.size()) return true;
    else return false;
}