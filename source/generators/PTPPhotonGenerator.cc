#include "PTPPhotonGenerator.h"

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

REGISTER_CLASS(PTPPhotonGenerator, G4VPrimaryGenerator)

PTPPhotonGenerator::PTPPhotonGenerator():
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

  msg_ = new G4GenericMessenger(this, "/Generator/PTPPhoton/",
    "Control commands of PTP photon generator.");

  msg_->DeclareProperty("region", region_,
    "Set the region of the geometry where the vertex will be generated.");

  msg_->DeclareProperty("pn_x", pn_x_, "X coordinate of PTP plane normal vector.");
  msg_->DeclareProperty("pn_y", pn_y_, "Y coordinate of PTP plane normal vector.");
  msg_->DeclareProperty("pn_z", pn_z_, "Z coordinate of PTP plane normal vector.");

  DetectorConstruction* detconst = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();

}

PTPPhotonGenerator::~PTPPhotonGenerator()
{
  if(sampler_)  delete sampler_;
  if(msg_)      delete msg_;
}

void PTPPhotonGenerator::GeneratePrimaryVertex(G4Event* event)
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

G4double PTPPhotonGenerator::RandomEnergy(){
    if(sampler_){
        return sampler_->operator()<std::mt19937>(gen_);
    }
    else{
        G4Exception("[PTPPhotonGenerator]", "RandomEnergy()", JustWarning,
        "PTP histogram sampler was not set. Energy could not be sampled.");
        return -1.;
    }
}

void PTPPhotonGenerator::LoadNCheckPTPData(){

    bin_edges_ =    {3.104*eV, 3.1261*eV, 3.1372*eV, 3.1545*eV, 3.1656*eV, 
                    3.1781*eV, 3.1869*eV, 3.1946*eV, 3.2016*eV, 3.2102*eV, 
                    3.2136*eV, 3.222*eV, 3.2307*eV, 3.2389*eV, 3.2454*eV, 
                    3.2542*eV, 3.2628*eV, 3.2706*eV, 3.2796*eV, 3.2879*eV, 
                    3.2972*eV, 3.3052*eV, 3.3068*eV, 3.3129*eV, 3.3152*eV, 
                    3.3192*eV, 3.3243*eV, 3.3258*eV, 3.3284*eV, 3.3319*eV, 
                    3.3345*eV, 3.3376*eV, 3.3403*eV, 3.3422*eV, 3.3496*eV, 
                    3.3517*eV, 3.3555*eV, 3.3604*eV, 3.3659*eV, 3.3739*eV, 
                    3.3871*eV, 3.4076*eV, 3.4199*eV, 3.4287*eV, 3.4367*eV, 
                    3.4415*eV, 3.4471*eV, 3.4501*eV, 3.4519*eV, 3.4536*eV, 
                    3.4553*eV, 3.457*eV, 3.4588*eV, 3.4605*eV, 3.4622*eV, 
                    3.4663*eV, 3.4709*eV, 3.4723*eV, 3.4812*eV, 3.4967*eV, 
                    3.5099*eV, 3.5123*eV, 3.517*eV, 3.5183*eV, 3.5229*eV, 
                    3.5235*eV, 3.526*eV, 3.5298*eV, 3.5313*eV, 3.5362*eV, 
                    3.5373*eV, 3.5422*eV, 3.5503*eV, 3.5601*eV, 3.571*eV, 
                    3.5762*eV, 3.5779*eV, 3.587*eV, 3.5962*eV, 3.6038*eV, 
                    3.6092*eV, 3.6141*eV, 3.6155*eV, 3.6191*eV, 3.6226*eV, 
                    3.6262*eV, 3.6266*eV, 3.6273*eV, 3.6279*eV, 3.6288*eV, 
                    3.6323*eV, 3.6338*eV, 3.6357*eV, 3.6365*eV, 3.6403*eV, 
                    3.6409*eV, 3.6434*eV, 3.6443*eV, 3.6443*eV, 3.6461*eV, 
                    3.6469*eV, 3.6513*eV, 3.6513*eV, 3.652*eV, 3.6556*eV, 
                    3.6582*eV, 3.659*eV, 3.6616*eV, 3.6633*eV, 3.6651*eV, 
                    3.6668*eV, 3.6705*eV, 3.672*eV, 3.6737*eV, 3.6737*eV, 
                    3.6737*eV, 3.6764*eV, 3.6806*eV, 3.6829*eV, 3.6859*eV, 
                    3.6897*eV, 3.6932*eV, 3.6974*eV, 3.7021*eV, 3.7115*eV, 
                    3.7255*eV, 3.7368*eV};

    emission_spectrum_ =     {0.0509, 0.066, 0.0763, 0.092, 0.1075, 0.1226, 0.1382, 
                            0.153, 0.1681, 0.186, 0.1981, 0.2136, 0.231, 0.2493, 
                            0.2687, 0.2847, 0.2995, 0.3136, 0.3277, 0.3404, 0.3568, 
                            0.3803, 0.3977, 0.4071, 0.4253, 0.4473, 0.4708, 0.4875, 
                            0.5056, 0.5411, 0.528, 0.5526, 0.5738, 0.5896, 0.6021, 
                            0.6126, 0.6298, 0.6466, 0.6612, 0.6815, 0.6993, 0.7076, 
                            0.7215, 0.742, 0.7631, 0.7815, 0.7966, 0.8221, 0.8113, 
                            0.8519, 0.837, 0.8855, 0.8678, 0.9153, 0.8986, 0.9286, 
                            0.949, 0.9662, 0.986, 0.9979, 0.9762, 0.9588, 0.9407, 
                            0.9207, 0.8922, 0.9066, 0.8777, 0.8603, 0.8419, 0.8219, 
                            0.8014, 0.7826, 0.7579, 0.7402, 0.7552, 0.7746, 0.7881, 
                            0.8055, 0.8116, 0.7945, 0.7775, 0.757, 0.74, 0.7188, 
                            0.6967, 0.6588, 0.6773, 0.6433, 0.6089, 0.6252, 0.5701, 
                            0.5903, 0.5374, 0.5537, 0.5096, 0.5247, 0.465, 0.4933, 
                            0.4796, 0.4343, 0.4503, 0.4192, 0.3877, 0.4034, 0.3674, 
                            0.3484, 0.3276, 0.2898, 0.3061, 0.2575, 0.2739, 0.2216, 
                            0.2387, 0.2072, 0.1935, 0.1798, 0.1598, 0.1424, 0.1243, 
                            0.1048, 0.0821, 0.0607, 0.0404, 0.0222, 0.0068, 0.0054}; 

    if(DataIsIllFormed()){
        G4Exception("[PTPPhotonGenerator]", "LoadNCheckPTPData()", FatalException,
        "The provided data is ill-formed. Check PTPPhotonGenerator::DataIsIllFormed for more info.");
    }

    return;
}

G4bool PTPPhotonGenerator::DataIsIllFormed(){
    if(bin_edges_.size()!=1+emission_spectrum_.size()) return true;
    else return false;
}