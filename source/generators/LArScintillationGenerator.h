
#ifndef LAR_PHOTON_GENERATOR_H
#define LAR_PHOTON_GENERATOR_H

#include <G4VPrimaryGenerator.hh>

#include <random>
#include <initializer_list>
#include <vector>

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;


namespace nexus {

  class GeometryBase;

  class LArScintillationGenerator: public G4VPrimaryGenerator
  {
  public:
    /// Constructor
    LArScintillationGenerator();
    /// Destructor
    ~LArScintillationGenerator();

    /// This method is invoked at the beginning of the event. It sets
    /// a primary vertex (that is, a photon in a given position and time)
    /// in the event.
    void GeneratePrimaryVertex(G4Event*);

  private:    
    G4GenericMessenger* msg_;

    const GeometryBase* geom_;      ///< Pointer to the detector geometry
    G4double pn_x_, pn_y_, pn_z_;   ///< Coordinates for the normal vector to the lambertian emitter plane 
    G4String region_;               ///< Region label that is given to the geometry to generate a vertex

    std::random_device rd_;     ///< Weak random generator (used to seed mt19937 random gen.)
    std::mt19937 gen_;          ///< Random number generator

    std::vector<G4double> bin_edges_;           ///< Bin edges for the PTP emission spectrum
    std::vector<G4double> emission_spectrum_;   ///< PTP emission spectrum
    std::piecewise_constant_distribution<G4double>* sampler_; ///< LAr scintillation emission spectrum histogram sampler

    void LoadNCheckPTPData();   ///< Loads LAr emission spectrum data into bin_edges_ and emission_spectrum_
    G4bool DataIsIllFormed();     ///< Checks whether the loaded data is ill-formed
    G4double RandomEnergy();    ///< Generate a random kinetic energy according to LAr emission spectrum.
  };

} // end namespace nexus

#endif
