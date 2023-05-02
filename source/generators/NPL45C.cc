#include "NPL45C.h"

#include "DetectorConstruction.h"
#include "GeometryBase.h"
#include "FactoryBase.h"
#include "RandomUtils.h"

#include <G4GenericMessenger.hh>
#include <G4OpticalPhoton.hh>
#include <G4RunManager.hh>
#include <G4PrimaryVertex.hh>
#include <G4Event.hh>
#include <G4RandomTools.hh>

#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(NPL45C, G4VPrimaryGenerator)

NPL45C::NPL45C():
G4VPrimaryGenerator(),
msg_(0),
geom_(0), 
ed_x_(0.), ed_y_(1.), ed_z_(0.),
region_(""),
rd_{},
gen_(rd_()),
wl_bin_edges_{},
emission_spectrum_{},
wl_sampler_(0)
{

  LoadNCheckSpectrumData();
  wl_sampler_ = new std::piecewise_constant_distribution(
      wl_bin_edges_.begin(), wl_bin_edges_.end(), emission_spectrum_.begin());

  msg_ = new G4GenericMessenger(this, "/Generator/NPL45C/",
    "Control commands of NPL45C generator.");

  msg_->DeclareProperty("region", region_,
    "Set the region of the geometry where the vertex will be generated.");

  msg_->DeclareProperty("ed_x", ed_x_, "X coordinate of laser emission direction.");
  msg_->DeclareProperty("ed_y", ed_y_, "Y coordinate of laser emission direction.");
  msg_->DeclareProperty("ed_z", ed_z_, "Z coordinate of laser emission direction.");

  DetectorConstruction* detconst = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();

}

NPL45C::~NPL45C()
{
  if(wl_sampler_)     delete wl_sampler_;
  if(msg_)            delete msg_;
}

void NPL45C::GeneratePrimaryVertex(G4Event* event)
{
  G4PrimaryParticle* a_photon = new G4PrimaryParticle(G4OpticalPhoton::Definition());

  // Generate photon momentum direction out of sampled angles
  G4double norm = sqrt((ed_x_*ed_x_)+(ed_y_*ed_y_)+(ed_z_*ed_z_));
  G4ThreeVector photon_momentum_dir = G4ThreeVector(ed_x_/norm, ed_y_/norm, ed_z_/norm);
  a_photon->SetMomentumDirection(photon_momentum_dir);
  a_photon->SetPolarization(G4PlaneVectorRand(photon_momentum_dir));

  // Sample photon energy
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

G4double NPL45C::RandomEnergy(){
    if(wl_sampler_){
        return wl_sampler_->operator()<std::mt19937>(gen_);
    }
    else{
        G4Exception("[NPL45C]", "RandomEnergy()", JustWarning,
        "Laser spectrum sampler was not set. Energy could not be sampled.");
        return -1.;
    }
}

void NPL45C::LoadNCheckSpectrumData(){

    wl_bin_edges_ = { h_Planck*c_light/(439.04*nm), h_Planck*c_light/(439.15*nm), h_Planck*c_light/(439.25*nm), h_Planck*c_light/(439.35*nm), 
                      h_Planck*c_light/(439.46*nm), h_Planck*c_light/(439.56*nm), h_Planck*c_light/(439.66*nm), h_Planck*c_light/(439.76*nm), 
                      h_Planck*c_light/(439.87*nm), h_Planck*c_light/(439.97*nm), h_Planck*c_light/(440.07*nm), h_Planck*c_light/(440.18*nm), 
                      h_Planck*c_light/(440.28*nm), h_Planck*c_light/(440.38*nm), h_Planck*c_light/(440.49*nm), h_Planck*c_light/(440.59*nm), 
                      h_Planck*c_light/(440.69*nm), h_Planck*c_light/(440.79*nm), h_Planck*c_light/(440.9*nm), h_Planck*c_light/(441.0*nm), 
                      h_Planck*c_light/(441.1*nm), h_Planck*c_light/(441.17*nm), h_Planck*c_light/(441.22*nm), h_Planck*c_light/(441.28*nm), 
                      h_Planck*c_light/(441.32*nm), h_Planck*c_light/(441.36*nm), h_Planck*c_light/(441.39*nm), h_Planck*c_light/(441.4*nm), 
                      h_Planck*c_light/(441.43*nm), h_Planck*c_light/(441.48*nm), h_Planck*c_light/(441.56*nm), h_Planck*c_light/(441.62*nm), 
                      h_Planck*c_light/(441.76*nm), h_Planck*c_light/(441.83*nm), h_Planck*c_light/(441.86*nm), h_Planck*c_light/(441.92*nm), 
                      h_Planck*c_light/(441.97*nm), h_Planck*c_light/(442.03*nm), h_Planck*c_light/(442.04*nm), h_Planck*c_light/(442.09*nm), 
                      h_Planck*c_light/(442.11*nm), h_Planck*c_light/(442.16*nm), h_Planck*c_light/(442.19*nm), h_Planck*c_light/(442.25*nm), 
                      h_Planck*c_light/(442.27*nm), h_Planck*c_light/(442.31*nm), h_Planck*c_light/(442.36*nm), h_Planck*c_light/(442.43*nm), 
                      h_Planck*c_light/(442.54*nm), h_Planck*c_light/(442.59*nm), h_Planck*c_light/(442.62*nm), h_Planck*c_light/(442.66*nm), 
                      h_Planck*c_light/(442.69*nm), h_Planck*c_light/(442.74*nm), h_Planck*c_light/(442.77*nm), h_Planck*c_light/(442.77*nm), 
                      h_Planck*c_light/(442.78*nm), h_Planck*c_light/(442.8*nm), h_Planck*c_light/(442.81*nm), h_Planck*c_light/(442.83*nm), 
                      h_Planck*c_light/(442.83*nm), h_Planck*c_light/(442.88*nm), h_Planck*c_light/(442.92*nm), h_Planck*c_light/(442.94*nm), 
                      h_Planck*c_light/(442.96*nm), h_Planck*c_light/(443.05*nm), h_Planck*c_light/(443.05*nm), h_Planck*c_light/(443.05*nm), 
                      h_Planck*c_light/(443.07*nm), h_Planck*c_light/(443.07*nm), h_Planck*c_light/(443.14*nm), h_Planck*c_light/(443.22*nm), 
                      h_Planck*c_light/(443.25*nm), h_Planck*c_light/(443.25*nm), h_Planck*c_light/(443.25*nm), h_Planck*c_light/(443.28*nm), 
                      h_Planck*c_light/(443.31*nm), h_Planck*c_light/(443.38*nm), h_Planck*c_light/(443.47*nm), h_Planck*c_light/(443.58*nm), 
                      h_Planck*c_light/(443.68*nm), h_Planck*c_light/(443.78*nm), h_Planck*c_light/(443.88*nm), h_Planck*c_light/(443.99*nm), 
                      h_Planck*c_light/(444.09*nm), h_Planck*c_light/(444.19*nm), h_Planck*c_light/(444.3*nm), h_Planck*c_light/(444.4*nm), 
                      h_Planck*c_light/(444.5*nm), h_Planck*c_light/(444.6*nm), h_Planck*c_light/(444.71*nm), h_Planck*c_light/(444.81*nm), 
                      h_Planck*c_light/(444.91*nm), h_Planck*c_light/(444.98*nm), h_Planck*c_light/(445.10*nm)};

    emission_spectrum_ = {0.01795, 0.01765, 0.0178, 0.01765, 0.0184, 0.01825, 0.01795, 0.01765, 0.01825, 
                          0.01795, 0.01765, 0.01765, 0.01735, 0.01765, 0.0172, 0.01765, 0.01825, 0.0181, 
                          0.0193, 0.02605, 0.04251, 0.06982, 0.0991, 0.13288, 0.15766, 0.18694, 0.22298, 
                          0.259, 0.30406, 0.27702, 0.25, 0.22973, 0.21171, 0.22973, 0.27027, 0.31306, 
                          0.36937, 0.3491, 0.31306, 0.27027, 0.30406, 0.36261, 0.42342, 0.38964, 0.36487, 
                          0.31982, 0.27477, 0.22748, 0.26802, 0.30631, 0.34009, 0.36261, 0.38513, 0.45721, 
                          0.54054, 0.5946, 0.68694, 0.79279, 0.86712, 1.0, 0.94594, 0.88739, 0.76802, 
                          0.68694, 0.5946, 0.44144, 0.39414, 0.3536, 0.30856, 0.27027, 0.26802, 0.23965, 
                          0.16439, 0.12245, 0.19959, 0.07772, 0.04119, 0.01943, 0.01705, 0.01675, 0.01645, 
                          0.016, 0.01585, 0.01645, 0.01615, 0.01645, 0.01555, 0.016, 0.01645, 0.0163, 
                          0.0163, 0.01555, 0.01555, 0.01622};

    if(SpectrumDataIsIllFormed()){
        G4Exception("[NPL45C]", "LoadNCheckSpectrumData()", FatalException,
        "The provided data is ill-formed. Check NPL45C::SpectrumDataIsIllFormed for more info.");
    }

    return;
}

G4bool NPL45C::SpectrumDataIsIllFormed(){
    if(wl_bin_edges_.size()!=1+emission_spectrum_.size()) return true;
    return false;
}