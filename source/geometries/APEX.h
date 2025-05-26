#ifndef APEX_H
#define APEX_H

#include "GeometryBase.h"
// Forward declaration of G4TwoVector (i.e. 'class G4TwoVector;')
// gives a compilation error which I don't know how to debug right now.
#include <G4TwoVector.hh>

class G4VPhysicalVolume;
class G4MaterialPropertiesTable;
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
    void ConstructAttachedDichroicFilter(G4VPhysicalVolume*) const; ///< Called by Construct() if !remove_MLS_ and !detach_DF_. Constructs a dichroic filter (DF)
                                                                    ///< which is attached to (in contact with) the WLS plate.
    void ConstructDetachedDichroicFilter(G4VPhysicalVolume*) const; ///< Called by Construct() if !remove_MLS_ and detach_DF_. Constructs a DF which is detached
                                                                    ///< from WLS plate.

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
    std::vector<G4TwoVector> GetTriangularPlatePrismBase() const;   ///< Returns the triangular prism base, as a function of the plate_length_ and 
                                                                    ///< plate_width_ attributes, which could be extruded to form the triangular
                                                                    ///< WLS plate, in case shape_code_==1.
    G4bool GeometryIsIllFormed();                                   ///< Checks whether the specified geometry, up to the given parameters, is feasible

  private:

    ///---- General attributes ----///
    G4String surrounding_media_;                                    ///< Which media to place the APEX in
                                                                    ///< 'lar'  - The APEX is placed in Liquid ARgon
                                                                    ///< 'gar'  - The APEX is placed in Gaseous ARgon
                                                                    ///< 'air'  - The APEX is placed in air - ¡Note that the implemented air has no bulk-absorption length yet!
                                                                    ///< Default behaviour is that of surrounding_media_=='lar'.
    G4int shape_code_;                                              ///< The shape of the built APEX (i.e. WLS plate and any upper layer, such as DF/no-DF substrates) depends
                                                                    ///< on this parameter. It can take the following values:
                                                                    ///< 0 -> Rectangular
                                                                    ///< 1 -> Triangular
    G4bool detach_DF_;                                              ///< Whether to detach the DF from the WLS plate. If so, as a consequence, a substrate that acts as a mechanical 
                                                                    ///< support for the MLS is added. Note that, if detach_DF_ is true, the DF is floating on top of the WLS plate,
                                                                    ///< which is an unphysical situation.
    G4double wlsp_DF_gap_;                                          ///< This parameter only makes a difference if detach_DF_ is set to true. Thickness of the gap which is left between the DF and the WLS plate.
    G4double DF_substrate_thickn_;                                  ///< This parameter only makes a difference if detach_DF_ is set to true. Thickness of the DF substrate.
    G4MaterialPropertiesTable* DF_substrate_mpt_;                   ///< This parameter only makes a difference if detach_DF_ is set to true. Material Properties Table of the DF substrate.
    G4double MLS_thickn_;                                           ///< Thickness of the DF multilayer structure (MLS) which is deposited on top of the WLS plate
    G4double MLS_rindex_;                                           ///< Effective refractive index of the multi-layer structure. Currently unused.
    G4double coating_thickn_;                                       ///< Thickness of the coating layer that is deposited over the MLS
    G4double coating_rindex_;                                       ///< Refractive index of the coating layer that is deposited over the MLS
    G4bool remove_coating_;                                         ///< Whether to remove the coating layer that is deposited over the MLS
    G4bool remove_MLS_;                                             ///< Whether to remove the DF (the MLS) together with the coating layer that is deposited on top of it
    G4double plate_length_, plate_thickn_, plate_width_;            ///< These are the WLS plate dimensions. In any case, plate_thickn_ is the span of the WLS plate along the y-axis. The rest of the
                                                                    ///< dimensions depend on the shape_code_ parameter. If shape_code_==0, then plate_length_ (resp. plate_width_) is the length of the
                                                                    ///< rectangular WLS plate along the x-axis (resp. z-axis). If shape_code_==1, then plate_length_ (resp. plate_width_) is the span
                                                                    ///< of the triangular WLS plate along the x-axis (resp. z-axis), which matches the base (resp. height) of the modelled isosceles
                                                                    ///< triangle. N.B.: Note that the X-Z parameterization of triangular case of the WLSPlate class is inverted with respect to the one
                                                                    ///< we use here. I.e. in the WLSPlate class, the height of the triangle is laid along the X axis, while here, it is laid along the
                                                                    ///< Z axis. The reason for this is that we want to keep the same orientation of the SiPM boards across the rectangular and
                                                                    ///< triangular cases. I.e. in both cases the SiPM board spans along the X-direction. Unifying this is convenient in the sense that
                                                                    ///< further code and offline analysis can be shared.
    G4double WLSp_rindex_;                                          ///< Refractive index of the wavelength shifting plate
    G4double secondary_wls_attlength_;                              ///< Constant (wavelength indepedent) attenuation length of the secondary WLShifter. For the particular case when the G2P_FB118() G4MaterialPropertiesTable is used, the G2P_FB118() function should take care of
                                                                    ///< of ignoring this input and setting the real (wavelength dependent) measured attenuation-length spectrum if a non positive (negative or null) is given to this parameter. For config_code_==1 (resp. 2), this
                                                                    ///< is the attenuation length for the WLS plate (WLS fibers).
    G4double cromophore_concentration_;                             ///< Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter (the WLS plate), in case G2P_FB118 is used.
    G4bool cryogenic_temperature_;                                  ///< Whether the secondary WLShifter is at cryogenic temperature or not. It only makes a difference if G2P_FB118 is used.
    G4double reflective_foil_thickn_;                               ///< Reflective foil thickness. For the case of shape_code_==0 (rectangular APEX), this foil thickness is exact. For the case of shape_code_==1
                                                                    ///< (triangular APEX), the foil thickness is only exact for the edge of the triangle which is aligned with the X-axis. For the other two (oblique)
                                                                    ///< edges, the foil thickness is approximately equal to this value. N.B.: The exact value of this parameter for such oblique sides is not important.
                                                                    ///< Exactly implementing it may not be worth it.
    G4bool remove_back_plane_foil_;                                 ///< If true, the WLS plate back plane (i.e. the plane which is parallel to the plane where the generated photons impinge) is not lined with reflective foil.
    G4bool remove_reflective_foil_;                                 ///< Whether to remove the reflective foil that covers the WLS plate.
    G4int SiPM_code_;                                               ///< Integer signalling which SiPM to construct
                                                                    ///< 1                  -> Hamamatsu S13360-6050VE
                                                                    ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                                                    ///< 3                  -> FBK-NUV-HD-CRYO-TT
                                                                    ///< 4                  -> Broadcom AFBR-S4N44P044M (2x2 SiPM array)
                                                                    ///< Any other integer  -> PerfectSiPMMPPC (100% efficiency)
    G4int num_phsensors_;                                           ///< Number of photosensors per board
    G4int board_position_code_;                                     ///< Integer signalling where to place the SiPM boards. This parameter only makes a difference if shape_code_ is set to 0.
                                                                    ///< In such case, it can take the following values:
                                                                    ///< 1                  -> One board facing the middle of one of the two largest WLS plate faces.
                                                                    ///< 2                  -> One board facing one of the smallest WLS plate faces.
                                                                    ///< Any other integer  -> One board facing each one of the smallest WLS plate faces (two boards in total).
                                                                    ///< In any case, the boards are arranged parallely to the APEX sides which are plate_length_ long. If shape_code_ is set
                                                                    ///< to 1, then only one board, facing the base of the triangular WLS plate, is placed.
    G4bool align_lower_edges_of_plate_and_SiPMs_;                   ///< If shape_code_ is set to 0, then this parameter only makes a difference if board_position_code_ is set to a value other
                                                                    ///< than 1. If shape_code_ is set to 1, then this parameter always applies. If set to true, the lower edges of the WLS plate
                                                                    ///< and the SiPMs are aligned. If set to false, then the center (along the plate_thickn_ dimension) of the SiPMs is aligned
                                                                    ///< with the center of the WLS plate
    G4double gap_;                                                  ///< Gap between the photosensors and the WLS plate. A negative gap can help modelate the immersion of the SiPMs into the dimples.
                                                                    ///< Be careful not to collide the SiPMs into the plate.
    G4bool ref_phsensors_supports_;                                 ///< Whether photosensors supports are reflective (the FR4 box that supports the SiPM)
    /// Dimples may be used in the future, but not for now ----------------------------------------------------------------------------------------------------
    /// For the moment, APEX will use cryo glue so that SiPMs are in optical contact with the plate. ---------------------------------------------------------- 
    G4bool with_dimples_;                                           ///< This parameter only makes a difference if shape_code_ is equal to 0. Whether the plate
                                                                    ///< has carved dimples on it.
    G4String dimple_type_;                                          ///< This parameter only makes a difference if shape_code_ is equal to 0. Dimple type. Might
                                                                    ///< be 'flat', 'cylindrical' or 'spherical'.
    G4double flat_dimple_width_, flat_dimple_depth_;                ///< This parameter only makes a difference if shape_code_ is equal to 0. Used for flat
                                                                    ///< dimples. The width of the dimple (along the board direction) and its depth, perpendicular
                                                                    ///< to the plate surface.
    G4double curvy_dimple_radius_;                                  ///< This parameter only makes a difference if shape_code_ is equal to 0. Used for cylindrical
                                                                    ///< or spherical dimples. Radius of the dimple.
    // N.B.: Note that, if at some point dimples should be simulated with a triangular APEX (i.e. with shape_code_==1), then such feature should be added first to
    // the WLSPlate class, which is called by APEX::ConstructWLSPlate. I.e. the WLSPlate class does not allow, at the moment of writing, to simulate dimples on a
    // triangular WLS plate.
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