
#ifndef MONOCHROMATIC_BEAM_H
#define MONOCHROMATIC_BEAM_H

#include <G4VPrimaryGenerator.hh>

#include <random>

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

namespace nexus {

  class GeometryBase;
  // This class is aims to implement a monochromatic collimated light beam.
  // The wavelength of the generated optical photons is sampled from a normal
  // distribution whose mean and std can be tuned via attributes.
  class MonochromaticBeam: public G4VPrimaryGenerator
  {
  public:
    /// Constructor
    MonochromaticBeam();
    /// Destructor
    ~MonochromaticBeam();

    /// This method is invoked at the beginning of the event. It sets
    /// a primary vertex (that is, a photon in a given position and time)
    /// in the event.
    void GeneratePrimaryVertex(G4Event*);

  private:    
    G4GenericMessenger* msg_;
    const GeometryBase* geom_;    ///< Pointer to the detector geometry
    
    G4String region_;             ///< Region label that is given to the geometry to generate a vertex

    // For further implementation:
    G4double wavelength_mean_;    ///< Mean wavelength (in G4 units) of the normal 
                                  ///< distribution from which wavelengths will be sampled
    G4double wavelength_std_;     ///< Standard deviation (in G4 units) of the normal 
                                  ///< from which wavelengths will be sampled
    G4double ed_x_, ed_y_, ed_z_; ///< Coordinates of the emission direction.

    std::default_random_engine generator_;
    std::normal_distribution<G4double>* distribution_ptr_;

    void SetUpSampler();                        ///< Initialize the normal sampler
    G4double RandomEnergy(G4int max_iter = 3);  ///< Generate a random kinetic energy according to LED emission spectrum.
  };

} // end namespace nexus

#endif
