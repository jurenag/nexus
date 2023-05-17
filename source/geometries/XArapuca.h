#ifndef X_ARAPUCA
#define X_ARAPUCA

#include "GeometryBase.h"

class G4VPhysicalVolume;
class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4UserLimits;

namespace nexus {

  /// X-Arapuca 

  class XArapuca: public GeometryBase
  {

  // This class is able to model the megacell geometry. You can find its dimensions here:
  // indico.fnal.gov/event/55302/contributions/245694/attachments/156811/204810/Preparation of VD-CBs WLS prototypes_220706.pdf

  public:
    ///Constructor
    XArapuca();
    ///Destructor
    ~XArapuca();

    void Construct();                   ///< Constructs the geometry

  private:

    void ConstructWLSPlate(G4VPhysicalVolume*) const;               ///< Called by Construct(). Adds the WLS plate.
    void ConstructFibers(G4VPhysicalVolume*) const;                 ///< Called by Construct(). 
    void ConstructPhotosensors(G4VPhysicalVolume*) const;           ///< DEPRECATED. Called by Construct(). Adds (floating) SiPMs.
    void ConstructBoards(G4VPhysicalVolume*) const;                 ///< Called by Construct(). Adds the board-mounted SiPMs.
    void ConstructReflectiveCase(G4VPhysicalVolume*) const;         ///< Called by Construct(). Adds the reflective case that encloses the X-ARAPUCA geometry.
    void ConstructCollectors(G4VPhysicalVolume*) const;             ///< Called by Construct(). Adds the artificial volumes that replace the dichroic filters.
    void ConstructDichroicAssemblies(G4VPhysicalVolume*) const;     ///< Called by Construct(). Adds the dichroic filters plane, i.e. dichroic filters embedded in a frame
                                                                    ///< If double_sided_, the DFAs are placed on both sides of the X-ARAPUCA. If not, one reflective cover is added 
                                                                    ///< on the dead side.
    
    
    G4ThreeVector GenerateVertex(const G4String&) const;
    G4bool geometry_is_ill_formed();                                    ///<Check whether the WLS plate fits within the XArapuca cavity and the given number of photosensors fit along the XArapuca long side

  private:

  // Plate dimensions are used in ConstructCollectors. This is not ok. 
  // Plate dimensions should just refer to config_code_==1.

    ///---- General attributes ----///
    G4int config_code_;                                             ///< Value in (1,2) which labels which X-ARAPUCA configuration will be simulated
                                                                    ///< 1 -    common X-ARAPUCA (with WLS plate)
                                                                    ///< 2 -    Replace WLS plate with optical fibers with SiPMs attached to its ends                                                                    
    G4String surrounding_media_;                                    ///< Which media to place the XArapuca in
                                                                    ///< 'lar'  - The XArapuca is placed in Liquid ARgon
                                                                    ///< 'gar'  - The XArapuca is placed in Gaseous ARgon
                                                                    ///< 'air'  - The XArapuca is placed in air - ¡Note that the implemented air has no bulk-absorption length yet!
                                                                    ///< Default behaviour is that of surrounding_media_=='lar'.
    G4double internal_length_, internal_width_, internal_thickn_;   ///< Internal dimensions of the reflector cavity
    G4double DFA_thickn_;                                           ///< Frame thickness 
    G4double DF_thickn_;                                            ///< Overall thickness of the dichroic filters (MLS+substrate). Must be smaller than the frame thickness.
    G4double DF_substrate_thickn_;                                  ///< Thickness of the dichroic filters substrate. Must be smaller than DF_thickn_.
    G4MaterialPropertiesTable* DF_substrate_mpt_;                   ///< Material Properties Table of the dichroic filters substrate.
    G4double DF_pos_wrt_DFA_pos_;                                   ///< Position (height) of the dichroic filters with respect to the DFA position (height). This parameter can take values from 0 to 1. See * below
    G4bool DF_are_coated_;                                          ///< Whether the filters are coated with PTP or not
    G4double coating_thickn_;                                       ///< Thickness of the coating layer that is deposited over the dichroic filters
    G4double outter_frame_width_along_wlsplength_;                  ///< Outter frame dimension along the length of the X-ARAPUCA
    G4double outter_frame_width_along_wlspwidth_;                   ///< Outter frame dimension along the width of the X-ARAPUCA
    G4double inner_frames_width_along_wlsplength_;                  ///< Inner frames dimension along the length of the X-ARAPUCA
    G4double inner_frames_width_along_wlspwidth_;                   ///< Inner frames dimension along the width of the X-ARAPUCA
    G4int df_no_along_wlsplength_;                                  ///< Number of dichroic filters along the length of the X-ARAPUCA
    G4int df_no_along_wlspwidth_;                                   ///< Number of dichroic filters along the width of the X-ARAPUCA
    G4bool DFA_frame_is_reflective_;                                ///< Whether the FR4 DFA frame is vikuiti-coated or not
    G4bool DFA_frame_is_specular_;                                  ///< Whether the vikuiti coating of the DFA frame is specular-spikely reflective or diffusively reflective. Only makes a difference if DFA_frame_is_reflective_==True.
    G4bool remove_DFs_;                                             ///< Whether to remove the dichroic filters or not
    G4bool remove_DFA_frame_;                                       ///< Whether to remove the dichroic filters assembly frame
    G4double secondary_wls_attlength_;                              ///< Attenuation length of the secondary WLShifter. For config_code_==1 (resp. 2), this is the attenuation length for the WLS plate (WLS fibers).
    G4double case_thickn_;                                          ///< Reflective foils thickness
    G4int SiPM_code_;                                               ///< Integer signalling which SiPM to construct
                                                                    ///< 1                  -> Hamamatsu S13360-6050VE
                                                                    ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                                                    ///< Any other integer  -> FBK-NUV-HD-CRYO-TT
    G4int num_phsensors_;                                           ///< If config_code_==1, this is the number of SiPMs per long side
                                                                    ///< If config_code_==2, this is the number of photosensors per side which is perpendicular to the fibers
    G4double gap_;                                                  ///< Gap between photosensors and WLS material. A negative gap can help modelate the immersion of the SiPMs into the flat dimples. Be careful not to collide the SiPMs into the plate.
    G4bool ref_phsensors_supports_;                                 ///< Whether photosensors supports are reflective (the FR4 box that supports the SiPM, regardless PS_config_code_)
    G4bool double_sided_;                                           ///< Whether there are dichroic filters on both sides of the WLS plate
    G4bool collectors_are_reflective_;                              ///< Whether the collectors that replace the dichroic filters are reflective or not                                      
    G4String generation_region_;                                    ///< Where to place the generation vertex (GV).
                                                                    ///< 'random'   - The GV is randomly sampled over the DFA (including the frame)
                                                                    ///< 'dichroic' - The GV is randomly sampled over the DFs (not including the frame)
                                                                    ///< 'custom'   - The GV is randomly sampled within a circle centered in (gen_x_, -, gen_z_)
                                                                    ///<              whose diameter is equal to gen_diameter_.
                                                                    ///< Default behaviour is that of generation_region_=='random'.
    G4double gen_x_, gen_z_;                                        ///< Average GV coordinates. It is only used if generation_region_=='custom' is True.
    G4double gen_diameter_;                                         ///< Diameter of the circle where the GV could be randomly sampled. It is only used if generation_region_=='custom' is True.
    G4String path_to_inwards_dichroic_data_;                        ///< Absolute path to dichroic data file that is to be sampled for the light trying to enter the XA cavity
    G4String path_to_outwards_dichroic_data_;                       ///< Absolute path to dichroic data file that is to be sampled for the light trying to escape the XA cavity
    G4double world_extra_thickn_;                                   ///< Extra thickness for the surrounding box world to wrap the XArapuca

    ///---- General internal attributes ----///                     ///< These attributes are internal. They must not be set by the user.
    G4double overall_length_, overall_width_, overall_thickn_;      ///< Overall X-ARAPUCA dimensions
    G4double DFA_length_, DFA_width_;                               ///< Dichroic filters assembly transverse dimensions.
    G4double DF_length_, DF_width_;                                 ///< Dichroic filter transverse dimensions. All of the filters of the DFA have the same dimensions.

    ///---- config_code_==1 parameters ----///
    G4double plate_length_, plate_width_, plate_thickn_;            ///< WLS plate thickness
    G4bool only_sipms_along_long_sides_;                            ///< Whether the photo sensors are installed only along the long sides of the X-ARAPUCA (i.e. the 'length' dimension), or along every side.
    G4bool with_boards_;                                            ///< Whether the photosensors are mounted on a board
    G4int PS_config_code_;                                          ///< Value in (1,2) which labels which photo sensors configuration will be simulated
                                                                    ///< 1 -    num_phsensors_ are evenly distributed along each chosen* side of the WLS plate.
                                                                    ///< 2 -    one big photosensor (whose transversal dimensions match those of the WLS plate side face) which matches the efficiency curve of HamamatsuS133606050VE is placed in front of each chosen* WLS 
                                                                    ///<        plate side face. In this case, which specific SiPM to implement cannot be chosen.  
                                                                    ///< WARNING: Currently, PS_config_code_==2 is only available if with_boards_==false, so setting with_boards_==true automatically entails that num_phsensors_ different SiPMs per side will be installed      
                                                                    ///< NOTE: chosen* up to only_sipms_along_long_sides_
    G4bool with_dimples_;                                           ///< Whether the plate has carved dimples on it.
    G4String dimple_type_;                                          ///< Type of the dimple, i.e. "cylindrical", "spherical" etc.
    G4double flat_dimple_width_, flat_dimple_depth_;                ///< Used for flat dimples. The width of the dimple (along the board direction) and its depth, perpendicular to the plate surface.
    G4double curvy_dimple_radius_;                                  ///< Used for cylindrical or spherical dimples. Radius of the dimple.
    ///---- config_code_==2 parameters ---///
    G4int fibers_no_;                                               ///< Total number of fibers are installed within the X-ARAPUCA. Must be divisible by fiber_planes_no_
    G4int fiber_planes_no_;                                         ///< How many fibers-horizontal-planes are installed within the X-ARAPUCA
    G4double fiber_radius_;                                         ///< Radius of each (cylindrical) fiber
    G4double fiber_length_;                                         ///< Length of each fiber
    G4bool along_long_side_;                                        ///< Whether the fibers are parallel to the long side of the X-ARAPUCA or not.
                                                                    ///< The two SiPMs board are always installed perpendicular to the fibers

    G4GenericMessenger* msg_;                                       ///< Messenger for the definition of control commands
    G4UserLimits* ul_;                                              ///< Useful to set a maximum track length in the plate
  };
}

#endif

/*
________ ___________________ ________
        |                   |  
        |       DF          | <--- DF_pos_wrt_DFA_pos_ = 1.  (Shallow)
        |___________________|
        |                   |
        |                   |
        |                   |       
        |                   |
  DFA   |                   |  DFA
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |___________________|
        |                   |
        |        DF         | <--- DF_pos_wrt_DFA_pos_ = 0. (Deepest position within DFA)
________|___________________|________


_____________________________________
                WLSPLATE
_____________________________________

________ ___________________ ________
        |                   |  
        |       DF          | <--- DF_pos_wrt_DFA_pos_ = 0. (Deepest position within DFA)
        |___________________|
        |                   |
        |                   |
        |                   |       
        |                   |
  DFA   |                   |  DFA
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |___________________|
        |                   |
        |        DF         | <--- DF_pos_wrt_DFA_pos_ = 1.  (Shallow)
________|___________________|________



*/
