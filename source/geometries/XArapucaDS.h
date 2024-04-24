#ifndef X_ARAPUCA_DS
#define X_ARAPUCA_DS

#include "GeometryBase.h"

class G4VPhysicalVolume;
class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4UserLimits;

namespace nexus {

  /// X-Arapuca Dual Scintillator (DS)

  class XArapucaDS: public GeometryBase
  {

  public:
    ///Constructor
    XArapucaDS();
    ///Destructor
    ~XArapucaDS();

    void Construct();                   ///< Constructs the geometry

  private:

    void ConstructWLSPlate(G4VPhysicalVolume*) const;               ///< Called by Construct(). Adds the WLS plate.
    void ConstructWLSRing(G4VPhysicalVolume*) const;                ///< Called by Construct(). 
    void ConstructPhotosensors(G4VPhysicalVolume*) const;           ///< DEPRECATED. Called by Construct(). Adds (floating) SiPMs.
    void ConstructBoards(G4VPhysicalVolume*) const;                 ///< Called by Construct(). Adds the board-mounted SiPMs.
    void ConstructReflectiveCase(G4VPhysicalVolume*) const;         ///< Called by Construct(). Adds the reflective case that encloses the X-ARAPUCA DS geometry.
    void ConstructSubstratesAssemblies(G4VPhysicalVolume*) const;   ///< Called by Construct(). Adds the coated substrates plane, i.e. coated substrates embedded in a frame
                                                                    ///< If double_sided_, the CSAs are placed on both sides of the X-ARAPUCA DS. If not, one reflective cover is added 
                                                                    ///< on the dead side.
    
    G4ThreeVector GenerateVertex(const G4String&) const;
    G4bool geometry_is_ill_formed();                                    ///<Check whether the WLS plate fits within the XArapuca cavity and the given number of photosensors fit along the XArapuca long side

  private:

  // Plate dimensions are used in ConstructCollectors. This is not ok. 
  // Plate dimensions should just refer to config_code_==1.

    ///---- General attributes ----///
    G4int config_code_;                                             ///< Value in (1,2) which labels which X-ARAPUCA DS configuration will be simulated
                                                                    ///< 1 -    X-ARAPUCA DS
    G4String surrounding_media_;                                    ///< Which media to place the XArapuca DS in
                                                                    ///< 'lar'  - The XArapuca DS is placed in Liquid ARgon
                                                                    ///< 'gar'  - The XArapuca DS is placed in Gaseous ARgon
                                                                    ///< 'air'  - The XArapuca DS is placed in air - ¡Note that the implemented air has no bulk-absorption length yet!
                                                                    ///< Default behaviour is that of surrounding_media_=='lar'.
    G4double internal_length_, internal_width_, internal_thickn_;   ///< Internal dimensions of the reflector cavity
    G4double CSA_thickn_;                                           ///< Stands for Coated Substrates Assembly. It is the frame thickness 
    G4double CS_thickn_;                                            ///< Thickness of the coated substrates.
    G4MaterialPropertiesTable* CS_mpt_;                             ///< Material Properties Table of the coated substrates.
    G4double CS_pos_wrt_CSA_pos_;                                   ///< Position (height) of the coated substrates with respect to the CSA position (height). This parameter can take values from 0 to 1. See * below
    G4bool substrates_are_coated_;                                  ///< Whether the substrates are coated with PTP or not
    G4double coating_thickn_;                                       ///< Thickness of the coating layer that is deposited over the substrates
    G4double outter_frame_width_along_wlsplength_;                  ///< Outter frame dimension along the length of the X-ARAPUCA DS
    G4double outter_frame_width_along_wlspwidth_;                   ///< Outter frame dimension along the width of the X-ARAPUCA DS
    G4double inner_frames_width_along_wlsplength_;                  ///< Inner frames dimension along the length of the X-ARAPUCA DS
    G4double inner_frames_width_along_wlspwidth_;                   ///< Inner frames dimension along the width of the X-ARAPUCA DS
    G4int cs_no_along_wlsplength_;                                  ///< Number of coated substrates along the length of the X-ARAPUCA DS
    G4int cs_no_along_wlspwidth_;                                   ///< Number of coated substrates along the width of the X-ARAPUCA DS
    G4bool CSA_frame_is_reflective_;                                ///< Whether the FR4 CSA frame is vikuiti-coated or not
    G4bool CSA_frame_is_specular_;                                  ///< Whether the vikuiti coating of the CSA frame is specular-spikely reflective or diffusively reflective. Only makes a difference if CSA_frame_is_reflective_==True.
    G4bool remove_CSs_;                                             ///< Whether to remove the coated substrates or not
    G4bool remove_CSA_frame_;                                       ///< Whether to remove the CSA frame
    G4double SS_cromophore_concentration_;                          ///< Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter, in case G2P_FB118 is used.
    G4double TS_cromophore_concentration_;                          ///< Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the tertiary WLShifter, in case FakeG2P_FB118 is used.
    G4double case_thickn_;                                          ///< Reflective foils thickness
    G4int SiPM_code_;                                               ///< Integer signalling which SiPM to construct
                                                                    ///< 1                  -> Hamamatsu S13360-6050VE
                                                                    ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                                                    ///< 3                  -> FBK-NUV-HD-CRYO-TT
                                                                    ///< Any other integer  -> PerfectSiPMMPPC (100% efficiency)
    G4int num_phsensors_;                                           ///< If config_code_==1, this is the number of SiPMs per long side
    G4double gap_;                                                  ///< Gap between photosensors and WLS ring. A negative gap can help modelate the immersion of the SiPMs into the flat dimples. Be careful not to collide the SiPMs into the plate.
    G4bool ref_phsensors_supports_;                                 ///< Whether photosensors supports are reflective (the FR4 box that supports the SiPM, regardless PS_config_code_)
    G4bool double_sided_;                                           ///< Whether there are coated substrates on both sides of the WLS plate
    G4bool collectors_are_reflective_;                              ///< Whether the collectors that replace the coated substrates are reflective or not
    G4String generation_region_;                                    ///< Where to place the generation vertex (GV).
                                                                    ///< 'random'   - The GV is randomly sampled over the CSA (including the frame)
                                                                    ///< 'substrates' - The GV is randomly sampled over the CSs (not including the frame)
                                                                    ///< 'custom'   - The GV is randomly sampled within a circle centered in (gen_x_, -, gen_z_)
                                                                    ///<              whose diameter is equal to gen_diameter_.
                                                                    ///< Default behaviour is that of generation_region_=='random'.
    G4double gen_x_, gen_z_;                                        ///< Average GV coordinates. It is only used if generation_region_=='custom' is True.
    G4double gen_diameter_;                                         ///< Diameter of the circle where the GV could be randomly sampled. It is only used if generation_region_=='custom' is True.
    G4double world_extra_thickn_;                                   ///< Extra thickness for the surrounding box world to wrap the XArapuca

    ///---- General internal attributes ----///                     ///< These attributes are internal. They must not be set by the user.
    G4double overall_length_, overall_width_, overall_thickn_;      ///< Overall X-ARAPUCA dimensions
    G4double CSA_length_, CSA_width_;                               ///< Coated substrates assembly transverse dimensions.
    G4double CS_length_, CS_width_;                                 ///< Coated substrate transverse dimensions. All of the filters of the CSA have the same dimensions.

    ///---- config_code_==1 parameters ----///
    G4double plate_length_, plate_width_, plate_thickn_;            ///< Secondary WLShifter plate dimensions
    G4double ring_thickness_;                                       ///< Tertiary WLShifter ring thickness
    G4double plate_ring_gap_;                                       ///< Gap in between the WLS plate and the WLS ring
    G4bool only_sipms_along_long_sides_;                            ///< Whether the photo sensors are installed only along the long sides of the X-ARAPUCA DS (i.e. the 'length' dimension), or along every side.
    G4bool with_boards_;                                            ///< Whether the photosensors are mounted on a board
    G4int PS_config_code_;                                          ///< Value in (1,2) which labels which photo sensors configuration will be simulated
                                                                    ///< 1 -    num_phsensors_ are evenly distributed along each chosen* side of the WLS plate.
                                                                    ///< 2 -    one big photosensor (whose transversal dimensions match those of the WLS plate side face) which matches the efficiency curve of HamamatsuS133606050VE is placed in front of each chosen* WLS 
                                                                    ///<        plate side face. In this case, which specific SiPM to implement cannot be chosen.  
                                                                    ///< WARNING: Currently, PS_config_code_==2 is only available if with_boards_==false, so setting with_boards_==true automatically entails that num_phsensors_ different SiPMs per side will be installed      
                                                                    ///< NOTE: chosen* up to only_sipms_along_long_sides_
    G4bool with_dimples_;                                           ///< Whether the ring has carved dimples on it.
    G4String dimple_type_;                                          ///< Type of the dimple, i.e. "cylindrical", "spherical" etc.
    G4double flat_dimple_width_, flat_dimple_depth_;                ///< Used for flat dimples. The width of the dimple (along the board direction) and its depth, perpendicular to the plate surface.
    G4double curvy_dimple_radius_;                                  ///< Used for cylindrical or spherical dimples. Radius of the dimple.

    G4GenericMessenger* msg_;                                       ///< Messenger for the definition of control commands
    G4UserLimits* ul_;                                              ///< Useful to set a maximum track length in the plate
  };
}

#endif

/*
________ ___________________ ________
        |                   |  
        |       CS          | <---CS_pos_wrt_CSA_pos_ = 1.  (Shallow)
        |___________________|
        |                   |
        |                   |
        |                   |       
        |                   |
  CSA   |                   |  CSA
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |___________________|
        |                   |
        |        CS         | <--- CS_pos_wrt_CSA_pos_ = 0. (Deepest position within CSA)
________|___________________|________


_____________________________________
                WLSPLATE
_____________________________________

________ ___________________ ________
        |                   |  
        |       CS          | <--- CS_pos_wrt_CSA_pos_ = 0. (Deepest position within CSA)
        |___________________|
        |                   |
        |                   |
        |                   |       
        |                   |
  CSA   |                   |  CSA
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |___________________|
        |                   |
        |        CS         | <--- CS_pos_wrt_CSA_pos_ = 1.  (Shallow)
________|___________________|________



*/
