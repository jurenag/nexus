#ifndef APEX_H
#define APEX_H

#include "GeometryBase.h"

class G4VPhysicalVolume;
class G4GenericMessenger;
class G4UserLimits;

namespace nexus {

  /// APEX

  class APEX: public GeometryBase
  {

  // This class models the APEX geometry (proposal for the PDS of the DUNE FD3) You can find its dimensions here:
  // You can find more information here: 
  // https://indico.fnal.gov/event/58097/contributions/276023/attachments/171487/231299/FD3-APEX-CollMtg-Sept27-23.pdf

  public:
    ///Constructor
    APEX();
    ///Destructor
    ~APEX();

    void Construct();                   ///< Constructs the geometry

  private:

    void ConstructWLSPlate(G4VPhysicalVolume*) const;               ///< Called by Construct(). Adds the WLS plate.
    void ConstructSiPMSAndBoard(G4VPhysicalVolume*) const;          ///< Alternative to ConstructBoard, which fixes the colliding-volumes problem.
    void ConstructReflectiveFoil(G4VPhysicalVolume*) const;         ///< Called by Construct(). Adds the reflective foil that encloses every APEX face but one.
                                                                    ///<                        This reflective foil has holes which match the SiPM positions.
    void ConstructDichroicFilter(G4VPhysicalVolume*) const;         ///< Called by Construct(). Adds the dichroic filter (DF).

    void ConstructBoard(G4VPhysicalVolume*) const;        ///< Deprecated. Constructs a SiPM board (a SiPMBoard object)
                                                          // The reason why this one is deprecated is the following one. 
                                                          // The board encasing volume collides into the reflective foil 
                                                          // volume. As a matter of fact, the order of construction of
                                                          // the reflective foil and the SiPMBoard matters. Although 
                                                          // constructing first the board, then the reflective foil, is 
                                                          // preferred, both construction orders lead to (unprobable but 
                                                          // possible) un-physical photon tracks (photons entering volumes
                                                          // which are prohibited due to the defective protrusions).
                                                          // As an alternative to this, I coded  APEX::ConstructSiPMSAndBoard,
                                                          // which constructs the board and the SiPMs one by one, without 
                                                          // the need for an encasing volume which may collide into the 
                                                          // reflective foil.
    
    G4ThreeVector GenerateVertex(const G4String&) const;
    G4bool GeometryIsIllFormed();                                   ///< Checks whether the specified geometry, up to the given parameters, is feasible

  private:

    ///---- General attributes ----///
    G4String surrounding_media_;                                    ///< Which media to place the APEX in
                                                                    ///< 'lar'  - The APEX is placed in Liquid ARgon
                                                                    ///< 'gar'  - The APEX is placed in Gaseous ARgon
                                                                    ///< 'air'  - The APEX is placed in air - ¡Note that the implemented air has no bulk-absorption length yet!
                                                                    ///< Default behaviour is that of surrounding_media_=='lar'.
    G4double MLS_thickn_;                                           ///< Thickness of the DF multilayer structure (MLS) which is deposited on top of the WLS plate
    G4double MLS_rindex_;                                           ///< Effective refractive index of the multi-layer structure
    G4double coating_thickn_;                                       ///< Thickness of the coating layer that is deposited over the MLS
    G4double coating_rindex_;                                       ///< Refractive index of the coating layer that is deposited over the MLS
    G4bool remove_coating_;                                         ///< Whether to remove the coating layer that is deposited over the MLS
    G4bool remove_MLS_;                                             ///< Whether to remove the DF (the MLS) together with the coating layer that is deposited on top of it
    G4double plate_length_, plate_thickn_, plate_width_;            ///< WLS plate dimensions
    G4double WLSp_rindex_;                                          ///< Refractive index of the wavelength shifting plate
    G4double secondary_wls_attlength_;                              ///< Attenuation length of the secondary WLShifter (the WLS plate), in case EJ286 is used.
    G4double cromophore_concentration_;                             ///< Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter (the WLS plate), in case G2P_FB118 is used.
    G4double reflective_foil_thickn_;                               ///< Reflective foil thickness
    G4int SiPM_code_;                                               ///< Integer signalling which SiPM to construct
                                                                    ///< 1                  -> Hamamatsu S13360-6050VE
                                                                    ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                                                    ///< 3                  -> FBK-NUV-HD-CRYO-TT
                                                                    ///< 4                  -> Broadcom AFBR-S4N44P044M (2x2 SiPM array)
                                                                    ///< Any other integer  -> PerfectSiPMMPPC (100% efficiency)
    G4int num_phsensors_;                                           ///< Number of photosensors per board
    G4int board_position_code_;                                     ///< Integer signalling where to place the SiPM board. In any case, the board is arranged parallely to the APEX sides which are plate_length_ long.
                                                                    ///< 1                  -> Facing the middle of one of the two largest WLS plate faces.
                                                                    ///< Any other integer  -> Facing one of the smallest WLS plate faces.
    G4double gap_;                                                  ///< Gap between the photosensors and the WLS plate. A negative gap can help modelate the immersion of the SiPMs into the dimples. Be careful not to collide the SiPMs into the plate.
    G4bool ref_phsensors_supports_;                                 ///< Whether photosensors supports are reflective (the FR4 box that supports the SiPM)
    /// Dimples may be used in the future, but not for now ----------------------------------------------------------------------------------------------------
    /// For the moment, APEX will use cryo glue so that SiPMs are in optical contact with the plate. ---------------------------------------------------------- 
    G4bool with_dimples_;                                           ///< Whether the plate has carved dimples on it.
    G4String dimple_type_;                                          ///< Dimple type. Might be 'flat', 'cylindrical' or 'spherical'.
    G4double flat_dimple_width_, flat_dimple_depth_;                ///< Used for flat dimples. The width of the dimple (along the board direction) and its depth, perpendicular to the plate surface.
    G4double curvy_dimple_radius_;                                  ///< Used for cylindrical or spherical dimples. Radius of the dimple.
    /// -------------------------------------------------------------------------------------------------------------------------------------------------------

    G4String generation_region_;                                    ///< Where to place the generation vertex (GV).
                                                                    ///< 'random'   - The GV is randomly sampled over the DF
                                                                    ///< 'custom'   - The GV is randomly sampled within a circle centered in (gen_x_, -, gen_z_)
                                                                    ///<              whose diameter is equal to gen_diameter_.
                                                                    ///< Default behaviour is that of generation_region_=='random'.
    G4double gen_x_, gen_z_;                                        ///< Average GV coordinates. It is only used if generation_region_=='custom' is True.
    G4double gen_diameter_;                                         ///< Diameter of the circle where the GV could be randomly sampled. It is only used if generation_region_=='custom' is True.
    G4String path_to_inwards_dichroic_data_;                        ///< Absolute path to the dichroic data file that is to be sampled for the light trying to enter the WLS plate. 
    G4String path_to_outwards_dichroic_data_;                       ///< Absolute path to the dichroic data file that is to be sampled for the light trying to escape the WLS plate.
                                                                    ///< WARNING: Unless you have re-compiled a custom version of Geant4, you lost the multi transmission curve functionality when resetting the laptop.
    G4double world_extra_thickn_;                                   ///< Extra thickness for the surrounding box world to wrap the APEX

    ///---- General internal attributes ----///                     ///< These attributes are internal. They must not be set by the user.
    G4double overall_length_, overall_width_, overall_thickn_;      ///< Overall APEX dimensions
    G4double board_length_;                                         ///< Length of the board which contains the SiPMs. 
                                                                    ///< It matches plate_length_ by default.
    G4GenericMessenger* msg_;                                       ///< Messenger for the definition of control commands
    G4UserLimits* ul_;                                              ///< Useful to set a maximum track length in the plate
  };
}

#endif