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
pn_x_(0.), pn_y_(1.), pn_z_(0.),
region_(""),
rd_{},
gen_(rd_()),
bin_edges_{},
emission_spectrum_{},
sampler_(0)
{
  
  LoadNCheckPTPData();
  sampler_ = new std::piecewise_constant_distribution(
      bin_edges_.begin(), bin_edges_.end(), emission_spectrum_.begin()
        );

  msg_ = new G4GenericMessenger(this, "/Generator/LArScintillation/",
    "Control commands of LAr scintillation generator.");

  msg_->DeclareProperty("region", region_,
    "Set the region of the geometry where the vertex will be generated.");

  msg_->DeclareProperty("pn_x", pn_x_, "X coordinate of PTP plane normal vector.");
  msg_->DeclareProperty("pn_y", pn_y_, "Y coordinate of PTP plane normal vector.");
  msg_->DeclareProperty("pn_z", pn_z_, "Z coordinate of PTP plane normal vector.");

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

  G4ThreeVector photon_momentum_dir = G4LambertianRand(G4ThreeVector(pn_x_, pn_y_, pn_z_));
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

void LArScintillationGenerator::LoadNCheckPTPData(){

    bin_edges_ = {  h_Planck*c_light/(308.842*nm), h_Planck*c_light/(305.646*nm), h_Planck*c_light/(301.7045*nm), h_Planck*c_light/(297.6565*nm), h_Planck*c_light/(293.8215*nm), 
                    h_Planck*c_light/(289.9865*nm), h_Planck*c_light/(286.6845*nm), h_Planck*c_light/(282.6365*nm), h_Planck*c_light/(277.523*nm), h_Planck*c_light/(272.5165*nm), 
                    h_Planck*c_light/(267.4035*nm), h_Planck*c_light/(262.7165*nm), h_Planck*c_light/(257.603*nm), h_Planck*c_light/(251.744*nm), h_Planck*c_light/(246.844*nm), 
                    h_Planck*c_light/(242.796*nm), h_Planck*c_light/(238.0025*nm), h_Planck*c_light/(232.996*nm), h_Planck*c_light/(229.161*nm), h_Planck*c_light/(225.9655*nm), 
                    h_Planck*c_light/(222.024*nm), h_Planck*c_light/(217.8695*nm), h_Planck*c_light/(213.822*nm), h_Planck*c_light/(209.8805*nm), h_Planck*c_light/(207.43*nm), 
                    h_Planck*c_light/(204.3405*nm), h_Planck*c_light/(200.3995*nm), h_Planck*c_light/(197.0975*nm), h_Planck*c_light/(193.795*nm), h_Planck*c_light/(190.2795*nm), 
                    h_Planck*c_light/(186.551*nm), h_Planck*c_light/(182.823*nm), h_Planck*c_light/(179.521*nm), h_Planck*c_light/(175.686*nm), h_Planck*c_light/(171.9575*nm), 
                    h_Planck*c_light/(169.081*nm), h_Planck*c_light/(166.5245*nm), h_Planck*c_light/(163.8615*nm), h_Planck*c_light/(161.6245*nm), h_Planck*c_light/(159.707*nm), 
                    h_Planck*c_light/(157.051*nm), h_Planck*c_light/(154.068*nm), h_Planck*c_light/(151.9305*nm), h_Planck*c_light/(150.972*nm), h_Planck*c_light/(150.3325*nm), 
                    h_Planck*c_light/(149.8*nm), h_Planck*c_light/(149.2675*nm), h_Planck*c_light/(148.6285*nm), h_Planck*c_light/(148.096*nm), h_Planck*c_light/(147.67*nm), 
                    h_Planck*c_light/(147.244*nm), h_Planck*c_light/(146.711*nm), h_Planck*c_light/(145.965*nm), h_Planck*c_light/(145.0065*nm), h_Planck*c_light/(143.622*nm), 
                    h_Planck*c_light/(142.45*nm), h_Planck*c_light/(141.8405*nm), h_Planck*c_light/(141.2015*nm), h_Planck*c_light/(140.354*nm), h_Planck*c_light/(139.3955*nm), 
                    h_Planck*c_light/(138.4635*nm), h_Planck*c_light/(137.611*nm), h_Planck*c_light/(136.6785*nm), h_Planck*c_light/(135.8265*nm), h_Planck*c_light/(135.1*nm), 
                    h_Planck*c_light/(134.202*nm), h_Planck*c_light/(133.3495*nm), h_Planck*c_light/(132.4365*nm), h_Planck*c_light/(131.5155*nm), h_Planck*c_light/(130.8765*nm), 
                    h_Planck*c_light/(130.0935*nm), h_Planck*c_light/(128.9215*nm), h_Planck*c_light/(127.3235*nm), h_Planck*c_light/(125.5125*nm), h_Planck*c_light/(124.234*nm), 
                    h_Planck*c_light/(123.4885*nm), h_Planck*c_light/(122.8495*nm), h_Planck*c_light/(121.7345*nm), h_Planck*c_light/(120.1435*nm)};

    emission_spectrum_ = {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.001, 0.001, 0.001, 0.0, 0.0, 0.0, 0.001, 0.001, 0.002, 
                            0.002, 0.003, 0.005, 0.009, 0.014, 0.024, 0.036, 0.054, 0.082, 0.113, 0.176, 0.275, 0.36, 0.475, 0.624, 0.82, 1.0, 0.82, 0.624, 0.441, 0.336, 0.236}; 

    if(DataIsIllFormed()){
        G4Exception("[PTPPhotonGenerator]", "LoadNCheckPTPData()", FatalException,
        "The provided data is ill-formed. Check PTPPhotonGenerator::DataIsIllFormed for more info.");
    }

    return;
}

G4bool LArScintillationGenerator::DataIsIllFormed(){
    if(bin_edges_.size()!=1+emission_spectrum_.size()) return true;
    else return false;
}