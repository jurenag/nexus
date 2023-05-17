#include "MonochromaticBeam.h"

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

REGISTER_CLASS(MonochromaticBeam, G4VPrimaryGenerator)

MonochromaticBeam::MonochromaticBeam():
G4VPrimaryGenerator(),
msg_              (nullptr  ),
geom_             (nullptr  ),
region_           (""       ), 
wavelength_mean_  (400. *nm ),
wavelength_std_   (4.   *nm ),
ed_x_             (0.       ), 
ed_y_             (1.       ), 
ed_z_             (0.       ),
generator_        {         },
distribution_ptr_ (nullptr  ) 
{

  msg_ = new G4GenericMessenger(this, "/Generator/MonochromaticBeam/",
    "Control commands of MonochromaticBeam generator.");

  msg_->DeclareProperty("region", region_,
    "Set the region of the geometry where the vertex will be generated.");

  G4GenericMessenger::Command& wm_cmd =
    msg_->DeclareProperty("wavelength_mean", wavelength_mean_,
		    "Mean of the normal distribution from which wavelengths will be sampled.");
  wm_cmd.SetUnitCategory("Length");
  wm_cmd.SetParameterName("wavelength_mean", false);
  wm_cmd.SetRange("wavelength_mean>0.");

  G4GenericMessenger::Command& ws_cmd =
    msg_->DeclareProperty("wavelength_std", wavelength_std_,
		    "Standard deviation of the normal from which wavelengths will be sampled.");
  ws_cmd.SetUnitCategory("Length");
  ws_cmd.SetParameterName("wavelength_std", false);
  ws_cmd.SetRange("wavelength_std>0.");

  msg_->DeclareProperty("ed_x", ed_x_, "X coordinate of the emission direction.");
  msg_->DeclareProperty("ed_y", ed_y_, "Y coordinate of the emission direction.");
  msg_->DeclareProperty("ed_z", ed_z_, "Z coordinate of the emission direction.");

  DetectorConstruction* detconst = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();
}

MonochromaticBeam::~MonochromaticBeam()
{
  if(msg_)              delete msg_;
  if(distribution_ptr_) delete distribution_ptr_;
}

void MonochromaticBeam::GeneratePrimaryVertex(G4Event* event)
{
  G4PrimaryParticle* a_photon = new G4PrimaryParticle(G4OpticalPhoton::Definition());

  // Generate photon momentum direction out of sampled angles
  G4double norm = sqrt((ed_x_*ed_x_)+(ed_y_*ed_y_)+(ed_z_*ed_z_));
  G4ThreeVector photon_momentum_dir = G4ThreeVector(ed_x_/norm, ed_y_/norm, ed_z_/norm);
  a_photon->SetMomentumDirection(photon_momentum_dir);
  a_photon->SetPolarization(G4PlaneVectorRand(photon_momentum_dir));

  // Sample photon energy
  if(!distribution_ptr_) SetUpSampler();      // Ideally, the sampler *distribution_ptr_ should
                                              // be initialized in the constructor of MonochromaticBeam.
                                              // However, in order to initialize the normal-distribution
                                              // sampler, one must give the mean and the std, which are 
                                              // taken from a macro file. You can check that the data 
                                              // from the macro file is read after the constructor 
                                              // execution. Hence, the distribution sampler cannot
                                              // be initialized within the constructor. An (inefficient)
                                              // solution is the one that is implemented right now.
                                              // It is inefficient because it has to evaluate a conditional
                                              // statement everytime a vertex is generated. An alternative
                                              // solution is to use the Box-Muller method, which gives a
                                              // close analytic formula to get normal-samples out of
                                              // uniformly-distributed random samples. These formulas could
                                              // be coded here, using wavelength_mean_ and wavelength_std_,
                                              // thus avoiding evaluation of conditional statements.
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

void MonochromaticBeam::SetUpSampler(){
  if(!distribution_ptr_){
    distribution_ptr_ = new std::normal_distribution<G4double>(wavelength_mean_, wavelength_std_);
  }
  return;
};

G4double MonochromaticBeam::RandomEnergy(G4int max_iter){
  G4int i=0;
  G4double random_wl;   // Depending on the values of wavelength_mean_
                        // and wavelength_std_, negative energies might
                        // be sampled. If a negative wavelength is sampled, 
                        // then it is dumped and another wavelength is
                        // sampled. If the sampled wavelengths happen to 
                        // be negative for max_iter times in a row, then
                        // the execution of the G4 application is interrupted.
                        // If that's the case, the values given to mean and std
                        // are probably unphysical.
  G4bool wl_is_ok = false;
  do{
    i += 1;
    random_wl = (*distribution_ptr_)(generator_);
    if(random_wl>0.) wl_is_ok = true;
  }while(!wl_is_ok && i<=max_iter);

  if(wl_is_ok) return h_Planck*c_light/(random_wl);
  else{
    G4Exception("[MonochromaticBeam]", "RandomEnergy()", JustWarning,
    "After the maximum number of sampling iterations, no physical energy could be sampled.");
  }
}