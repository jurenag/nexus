
#ifndef LED375L_H
#define LED375L_H

#include <G4VPrimaryGenerator.hh>

#include <random>
#include <initializer_list>
#include <vector>

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;


namespace nexus {

  class GeometryBase;

  class LED375L: public G4VPrimaryGenerator
  {
  public:
    /// Constructor
    LED375L();
    /// Destructor
    ~LED375L();

    /// This method is invoked at the beginning of the event. It sets
    /// a primary vertex (that is, a photon in a given position and time)
    /// in the event.
    void GeneratePrimaryVertex(G4Event*);

  private:    
    G4GenericMessenger* msg_;

    const GeometryBase* geom_;      ///< Pointer to the detector geometry
    
    // For further implementation:
    //G4double pn_x_, pn_y_, pn_z_;   ///< Coordinates for the normal vector to the emitter plane 

    G4bool along_positive_z_;       ///< Whether the LED aims to positive values of the Z-axis
    G4String region_;               ///< Region label that is given to the geometry to generate a vertex

    std::random_device rd_;     ///< Weak random generator (used to seed mt19937 random gen.)
    std::mt19937 gen_;          ///< Random number generator

    std::vector<G4double> wl_bin_edges_;                          ///< Bin edges for the LED emission spectrum
    std::vector<G4double> emission_spectrum_;                     ///< LED emission spectrum
    std::piecewise_constant_distribution<G4double>* wl_sampler_;  ///< LED emission spectrum histogram sampler

    std::vector<G4double> angle_bin_edges_;                         ///< Bin edges for the LED angular distribution
    std::vector<G4double> angular_distribution_;                    ///< LED angular distribution
    std::piecewise_constant_distribution<G4double>* angle_sampler_; ///< LED emission spectrum histogram sampler

    void LoadNCheckSpectrumData();  ///< Loads LED emission spectrum data into wl_bin_edges_ and emission_spectrum_
    void LoadNCheckAngleData();     ///< Loads LED angular distribution data into angle_bin_edges_ and angular_distribution   
    G4bool SpectrumDataIsIllFormed();   ///< Checks whether the loaded data is ill-formed
    G4bool AngularDataIsIllFormed();
    G4double RandomEnergy();    ///< Generate a random kinetic energy according to LED emission spectrum.
    G4double RandomAngle();     ///< Generate a random emission angle according to the LED angular distribution
  };

} // end namespace nexus

#endif
