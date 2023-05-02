#include "LED375L.h"

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

REGISTER_CLASS(LED375L, G4VPrimaryGenerator)

LED375L::LED375L():
G4VPrimaryGenerator(),
msg_(0),
geom_(0), 
//pn_x_(0.), pn_y_(1.), pn_z_(0.),
along_positive_z_(true),
region_(""),
rd_{},
gen_(rd_()),
wl_bin_edges_{},
emission_spectrum_{},
wl_sampler_(0),
angle_bin_edges_{},
angular_distribution_{},
angle_sampler_(0)
{

  LoadNCheckSpectrumData();
  wl_sampler_ = new std::piecewise_constant_distribution(
      wl_bin_edges_.begin(), wl_bin_edges_.end(), emission_spectrum_.begin());

  LoadNCheckAngleData();
  angle_sampler_ = new std::piecewise_constant_distribution(
      angle_bin_edges_.begin(), angle_bin_edges_.end(), angular_distribution_.begin());

  msg_ = new G4GenericMessenger(this, "/Generator/LED375L/",
    "Control commands of LED375L generator.");

  msg_->DeclareProperty("region", region_,
    "Set the region of the geometry where the vertex will be generated.");

  /*
  msg_->DeclareProperty("pn_x", pn_x_, "X coordinate of LED emission direction.");
  msg_->DeclareProperty("pn_y", pn_y_, "Y coordinate of LED emission direction.");
  msg_->DeclareProperty("pn_z", pn_z_, "Z coordinate of LED emission direction.");
  */


  G4GenericMessenger::Command& apz_cmd =
    msg_->DeclareProperty("along_positive_z", along_positive_z_,
  		    "Whether the LED aims to positive values of the Z-axis");

  DetectorConstruction* detconst = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();

}

LED375L::~LED375L()
{
  if(wl_sampler_)     delete wl_sampler_;
  if(angle_sampler_)  delete angle_sampler_;
  if(msg_)            delete msg_;
}

void LED375L::GeneratePrimaryVertex(G4Event* event)
{
  G4PrimaryParticle* a_photon = new G4PrimaryParticle(G4OpticalPhoton::Definition());

  // Sample angles
  G4double theta = RandomAngle();
  G4double phi = nexus::UniformRandomInRange(2*pi, 0.);

  // Generate photon momentum direction out of smapled angles
  G4ThreeVector photon_momentum_dir = G4ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
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

G4double LED375L::RandomEnergy(){
    if(wl_sampler_){
        return wl_sampler_->operator()<std::mt19937>(gen_);
    }
    else{
        G4Exception("[LED375L]", "RandomEnergy()", JustWarning,
        "LED spectrum sampler was not set. Energy could not be sampled.");
        return -1.;
    }
}

G4double LED375L::RandomAngle(){
    if(angle_sampler_){
        return angle_sampler_->operator()<std::mt19937>(gen_);
    }
    else{
        G4Exception("[LED375L]", "RandomAngle()", JustWarning,
        "LED angular distribution sampler was not set. Angle could not be sampled.");
        return -1.;
    }
}

void LED375L::LoadNCheckSpectrumData(){

    wl_bin_edges_ = {   h_Planck*c_light/(351.19*nm), h_Planck*c_light/(353.08*nm), h_Planck*c_light/(354.97*nm), 
                        h_Planck*c_light/(356.86*nm), h_Planck*c_light/(358.75*nm), h_Planck*c_light/(360.64*nm), 
                        h_Planck*c_light/(362.05*nm), h_Planck*c_light/(362.99*nm), h_Planck*c_light/(363.77*nm), 
                        h_Planck*c_light/(364.45*nm), h_Planck*c_light/(364.96*nm), h_Planck*c_light/(365.33*nm), 
                        h_Planck*c_light/(365.64*nm), h_Planck*c_light/(366.1*nm), h_Planck*c_light/(366.26*nm), 
                        h_Planck*c_light/(366.41*nm), h_Planck*c_light/(366.56*nm), h_Planck*c_light/(366.87*nm), 
                        h_Planck*c_light/(367.44*nm), h_Planck*c_light/(367.65*nm), h_Planck*c_light/(367.85*nm), 
                        h_Planck*c_light/(367.95*nm), h_Planck*c_light/(368.1*nm), h_Planck*c_light/(368.35*nm), 
                        h_Planck*c_light/(369.03*nm), h_Planck*c_light/(369.5*nm), h_Planck*c_light/(369.51*nm), 
                        h_Planck*c_light/(369.69*nm), h_Planck*c_light/(369.71*nm), h_Planck*c_light/(369.83*nm), 
                        h_Planck*c_light/(370.23*nm), h_Planck*c_light/(370.64*nm), h_Planck*c_light/(370.78*nm), 
                        h_Planck*c_light/(371.06*nm), h_Planck*c_light/(371.73*nm), h_Planck*c_light/(371.86*nm), 
                        h_Planck*c_light/(372.27*nm), h_Planck*c_light/(372.68*nm), h_Planck*c_light/(372.82*nm), 
                        h_Planck*c_light/(373.24*nm), h_Planck*c_light/(373.38*nm), h_Planck*c_light/(373.52*nm), 
                        h_Planck*c_light/(374.47*nm), h_Planck*c_light/(374.74*nm), h_Planck*c_light/(375.16*nm), 
                        h_Planck*c_light/(375.49*nm), h_Planck*c_light/(375.5*nm), h_Planck*c_light/(375.66*nm), 
                        h_Planck*c_light/(375.98*nm), h_Planck*c_light/(376.55*nm), h_Planck*c_light/(377.23*nm), 
                        h_Planck*c_light/(377.64*nm), h_Planck*c_light/(378.04*nm), h_Planck*c_light/(378.15*nm), 
                        h_Planck*c_light/(378.52*nm), h_Planck*c_light/(379.0*nm), h_Planck*c_light/(379.57*nm), 
                        h_Planck*c_light/(380.39*nm), h_Planck*c_light/(380.94*nm), h_Planck*c_light/(381.38*nm), 
                        h_Planck*c_light/(382.17*nm), h_Planck*c_light/(383.0*nm), h_Planck*c_light/(383.81*nm), 
                        h_Planck*c_light/(384.7*nm), h_Planck*c_light/(385.46*nm), h_Planck*c_light/(387.1*nm), 
                        h_Planck*c_light/(387.65*nm), h_Planck*c_light/(388.65*nm), h_Planck*c_light/(390.07*nm), 
                        h_Planck*c_light/(391.65*nm), h_Planck*c_light/(393.39*nm), h_Planck*c_light/(395.28*nm), 
                        h_Planck*c_light/(397.17*nm), h_Planck*c_light/(399.06*nm), h_Planck*c_light/(400.95*nm), 
                        h_Planck*c_light/(402.84*nm), h_Planck*c_light/(404.72*nm), h_Planck*c_light/(406.61*nm), 
                        h_Planck*c_light/(408.5*nm), h_Planck*c_light/(410.39*nm), h_Planck*c_light/(412.28*nm), 
                        h_Planck*c_light/(414.17*nm), h_Planck*c_light/(416.06*nm), h_Planck*c_light/(417.95*nm), 
                        h_Planck*c_light/(419.84*nm), h_Planck*c_light/(421.73*nm), h_Planck*c_light/(423.62*nm), 
                        h_Planck*c_light/(425.51*nm), h_Planck*c_light/(427.4*nm), h_Planck*c_light/(429.29*nm), 
                        h_Planck*c_light/(431.18*nm), h_Planck*c_light/(433.07*nm), h_Planck*c_light/(434.96*nm), 
                        h_Planck*c_light/(436.84*nm), h_Planck*c_light/(438.73*nm), h_Planck*c_light/(439.99*nm),
                        h_Planck*c_light/(441.00*nm)};

    emission_spectrum_ = {  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.01112, 0.03918, 0.07816, 0.12238, 0.17324, 0.21889, 
                            0.265, 0.30352, 0.32919, 0.36627, 0.39480, 0.43331, 0.47182, 0.51556, 0.56025, 0.6059, 0.65724, 
                            0.69574, 0.74074, 0.72593, 0.81128, 0.77705, 0.88148, 0.84932, 0.91605, 0.95803, 0.92593, 0.89383, 
                            0.8716, 0.91852, 0.96049, 1.0, 0.97531, 0.91852, 0.85433, 0.81234, 0.78271, 0.82222, 0.84939, 0.79507, 
                            0.74852, 0.7143, 0.67864, 0.63014, 0.60247, 0.60987, 0.5926, 0.57166, 0.52317, 0.4761, 0.43473, 0.40247, 
                            0.40247, 0.38025, 0.34345, 0.29495, 0.27407, 0.26419, 0.22365, 0.19259, 0.17284, 0.14567, 0.11953, 
                            0.089, 0.05306, 0.02254, 0.0097, 0.00114, 0.00067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00067, 0.00019, 0.0, 
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if(SpectrumDataIsIllFormed()){
        G4Exception("[LED375L]", "LoadNCheckSpectrumData()", FatalException,
        "The provided data is ill-formed. Check LED375L::SpectrumDataIsIllFormed for more info.");
    }

    return;
}

void LED375L::LoadNCheckAngleData(){

    angle_bin_edges_ = {-24.615*deg, -23.568*deg, -22.521*deg, -21.474*deg, -20.429*deg, -19.646*deg, -19.213*deg, -18.838*deg, 
                        -18.522*deg, -18.206*deg, -18.006*deg, -17.834*deg, -17.634*deg, -17.318*deg, -17.002*deg, -16.628*deg, 
                        -16.282*deg, -15.936*deg, -15.59*deg, -15.243*deg, -14.897*deg, -14.551*deg, -14.205*deg, -13.771*deg, 
                        -12.989*deg, -12.115*deg, -11.501*deg, -10.887*deg, -10.012*deg, -8.9632*deg, -7.9165*deg, -6.8713*deg, 
                        -5.9141*deg, -5.0443*deg, -4.0871*deg, -3.0417*deg, -1.9944*deg, -0.9471*deg, 0.01093*deg, 0.70608*deg, 
                        1.314*deg, 2.1842*deg, 3.0577*deg, 3.671*deg, 4.1976*deg, 4.724*deg, 5.3378*deg, 6.0387*deg, 6.914*deg, 
                        7.9637*deg, 9.0131*deg, 10.062*deg, 11.109*deg, 12.157*deg, 13.205*deg, 14.079*deg, 14.693*deg, 
                        15.132*deg, 15.483*deg, 15.835*deg, 16.187*deg, 16.539*deg, 16.978*deg, 17.505*deg, 18.032*deg, 
                        18.558*deg, 18.998*deg, 19.349*deg, 19.701*deg, 20.14*deg, 20.754*deg, 21.629*deg, 22.677*deg, 
                        23.724*deg, 24.684*deg, 25.684*deg};
    
    angular_distribution_ = {   0.02058, 0.0234, 0.02766, 0.03383, 0.05008, 0.08313, 0.12293, 0.16559, 0.21499, 0.27013, 0.31425, 
                                0.34878, 0.39002, 0.44445, 0.49935, 0.55449, 0.59477, 0.62784, 0.66236, 0.69832, 0.73428, 0.77024, 
                                0.80619, 0.84646, 0.87855, 0.86075, 0.8233, 0.78513, 0.74598, 0.72578, 0.73003, 0.75396, 0.78729, 
                                0.82351, 0.85741, 0.87846, 0.87552, 0.87162, 0.8946, 0.93197, 0.96863, 1.0, 0.98651, 0.95387, 
                                0.91449, 0.87705, 0.83744, 0.80071, 0.76467, 0.7296, 0.69932, 0.6796, 0.66945, 0.66556, 0.65543, 
                                0.629, 0.59058, 0.55698, 0.52242, 0.48643, 0.44899, 0.41155, 0.37267, 0.33139, 0.29202, 0.25074, 
                                0.21138, 0.17395, 0.13939, 0.09907, 0.0573, 0.02847, 0.01881, 0.01492, 0.01246}; 

    if(AngularDataIsIllFormed()){
        G4Exception("[LED375L]", "LoadNCheckAngleData()", FatalException,
        "The provided data is ill-formed. Check LED375L::AngularDataIsIllFormed for more info.");
    }

    return;
}

G4bool LED375L::SpectrumDataIsIllFormed(){
    if(wl_bin_edges_.size()!=1+emission_spectrum_.size()) return true;
    return false;
}

G4bool LED375L::AngularDataIsIllFormed(){
    if(angle_bin_edges_.size()!=1+angular_distribution_.size()) return true;
    return false;
}