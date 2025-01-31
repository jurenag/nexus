#ifndef X_ARAPUCA
#define X_ARAPUCA

#include "GeometryBase.h"

class G4VPhysicalVolume;
class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4UserLimits;

namespace nexus {

  class SiPMMPPC; // Declaring it outside of the namespace environment gives a compile-time error

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
    void ConstructPlateReflectiveWrap(G4VPhysicalVolume*) const;    ///< Called by Construct(). If config_code_ is 1, double_sided_ is False, with_boards_ is False and 
                                                                    ///< wrap_plate_in_ref_foil_ is True, then this function wraps the WLS plate with a reflective foil 
                                                                    ///< which covers all faces but the one facing the DFA.
    
    
    G4ThreeVector GenerateVertex(const G4String&) const;
    G4bool geometry_is_ill_formed();                                    ///< Check whether the WLS plate fits within the XArapuca cavity and the given number of photosensors fit along the XArapuca long side
    SiPMMPPC * ChooseSiPM(G4int) const;                                 ///< Returns a pointer to a new object of a class which inherits from the SiPMMPPC abstract class.
                                                                        ///< Such class depends on the integer value given to ChooseSiPM according to the following enumeration
                                                                        ///< 1: HamamatsuS133606050VE
                                                                        ///< 2: HamamatsuS133605075HQR
                                                                        ///< 3: FbkNuvHdCryoTT
                                                                        ///< 4: BroadcomAFBRS4N44P044M
                                                                        ///< any other integer: PerfectSiPMMPPC
    G4MaterialPropertiesTable * ChooseDFSubstrateMPT(G4String) const;   ///< Returns a pointer to a new G4MaterialPropertiesTable object, which implements the optical
                                                                        ///< properties of the dichroic filters substrate. The mapping encoded by this function is the
                                                                        ///< following:
                                                                        ///< "SCHOTT_B270" - opticalprops::SCHOTT_B270()
                                                                        ///< "SCHOTT_BOROFLOAT_33" - opticalprops::SCHOTT_BOROFLOAT_33()
                                                                        ///< "FUSED_SILICA" - opticalprops::FusedSilica()

  private:

  // Plate dimensions are used in ConstructCollectors. This is not ok. 
  // Plate dimensions should just refer to config_code_==1.

    ///---- General attributes ----///
    G4int config_code_;                                             ///< Value in (1,2) which labels which X-ARAPUCA configuration will be simulated
                                                                    ///< 1 -    common X-ARAPUCA (with WLS plate)
                                                                    ///< 2 -    Replace WLS plate with optical fibers with SiPMs attached to its ends                                                                    
    G4String surrounding_media_;                                    ///< Which media to place the XArapuca in
                                                                    ///< "lar"  - The XArapuca is placed in Liquid ARgon
                                                                    ///< "gar"  - The XArapuca is placed in Gaseous ARgon
                                                                    ///< "air"  - The XArapuca is placed in air - ¡Note that the implemented air has no bulk-absorption length yet!
                                                                    ///< Default behaviour is that of surrounding_media_=="lar".
    G4double internal_length_, internal_width_, internal_thickn_;   ///< Internal dimensions of the reflector cavity
    G4double DFA_thickn_;                                           ///< Frame thickness 
    G4double DF_thickn_;                                            ///< Overall thickness of the dichroic filters (MLS+substrate). Must be smaller than the frame thickness.
    G4double DF_substrate_thickn_;                                  ///< Thickness of the dichroic filters substrate. Must be smaller than DF_thickn_.
    G4String DF_substrate_mpt_;                                     ///< Gives the material properties table of the dichroic filters substrate. It can take the following values:
                                                                    ///< "SCHOTT_B270" - opticalprops::SCHOTT_B270()
                                                                    ///< "SCHOTT_BOROFLOAT_33" - opticalprops::SCHOTT_BOROFLOAT_33()
                                                                    ///< "FUSED_SILICA" - opticalprops::FusedSilica()
    G4double DF_pos_wrt_DFA_pos_;                                   ///< Position (height) of the dichroic filters with respect to the DFA position (height). This parameter can take values from 0 to 1. See * below
    G4bool DF_are_coated_;                                          ///< Whether the filters are coated with PTP or not
    G4double coating_thickn_;                                       ///< Thickness of the coating layer that is deposited over the dichroic filters
    G4double coating_rindex_;                                       ///< Refractive index of the coating layer
    G4double outter_frame_width_along_wlsplength_;                  ///< Outter frame dimension along the length of the X-ARAPUCA
    G4double outter_frame_width_along_wlspwidth_;                   ///< Outter frame dimension along the width of the X-ARAPUCA
    G4double inner_frames_width_along_wlsplength_;                  ///< Inner frames dimension along the length of the X-ARAPUCA
    G4double inner_frames_width_along_wlspwidth_;                   ///< Inner frames dimension along the width of the X-ARAPUCA
    G4int df_no_along_wlsplength_;                                  ///< Number of dichroic filters along the length of the X-ARAPUCA
    G4int df_no_along_wlspwidth_;                                   ///< Number of dichroic filters along the width of the X-ARAPUCA
    G4bool DFA_frame_is_reflective_;                                ///< Whether the FR4 DFA frame is vikuiti-coated or not
    G4double vikuiti_reflectivity_scale_factor_;                    ///< Scale factor for the vikuiti reflectivity curve. It must belong to the [0., 1.] range. Note that this affects every volume which implements the vikuiti optical
                                                                    ///< properties in the XArapuca geometry. P.e. it also affects the vikuiti that coats the SiPM boards.
    G4bool DFA_frame_is_specular_;                                  ///< Whether the vikuiti coating of the DFA frame is specular-spikely reflective or diffusively reflective. It only makes a difference if DFA_frame_is_reflective_==True.
    G4bool remove_DFs_;                                             ///< Whether to remove the dichroic filters or not
    G4bool remove_DFA_frame_;                                       ///< Whether to remove the dichroic filters assembly frame
    G4double secondary_wls_attlength_;                              ///< Attenuation length of the secondary WLShifter, in case EJ286 is used. For config_code_==1 (resp. 2), this is the attenuation length for the WLS plate (WLS fibers).
    G4double cromophore_concentration_;                             ///< Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter, in case G2P_FB118 is used. For config_code_==1 (resp. 2), this is the cromophore concentration for the WLS plate (WLS fibers).
    G4double case_thickn_;                                          ///< Reflective foils thickness
    G4int SiPM_code_;                                               ///< Integer signalling which SiPM to construct
                                                                    ///< 1                  -> Hamamatsu S13360-6050VE
                                                                    ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                                                    ///< 3                  -> FBK-NUV-HD-CRYO-TT
                                                                    ///< 4                  -> Broadcom AFBR-S4N44P044M (2x2 SiPM array)
                                                                    ///< Any other integer  -> PerfectSiPMMPPC (100% efficiency)
    G4int num_phsensors_;                                           ///< If config_code_==1, this is the number of SiPMs per long side
                                                                    ///< If config_code_==2, this is the number of photosensors per side which is perpendicular to the fibers
    G4double gap_;                                                  ///< Gap between photosensors and WLS material. A negative gap can help modelate the immersion of the SiPMs into the flat dimples. Be careful not to collide the SiPMs into the secondary-WLS material.
    G4bool ref_phsensors_supports_;                                 ///< Whether photosensors supports are reflective (the FR4 box that supports the SiPM, regardless PS_config_code_)
    G4bool double_sided_;                                           ///< Whether there are dichroic filters on both sides of the secondary-WLS material
    G4bool collectors_are_reflective_;                              ///< Whether the collectors that replace the dichroic filters are reflective or not                                      
    G4String generation_region_;                                    ///< Where to place the generation vertex (GV).
                                                                    ///< "random"   - The GV is randomly sampled over the DFA (including the frame)
                                                                    ///< "dichroic" - The GV is randomly sampled over the DFs (not including the frame)
                                                                    ///< "custom"   - The GV is randomly sampled within a circle centered in (gen_x_, -, gen_z_)
                                                                    ///<              whose diameter is equal to gen_diameter_.
                                                                    ///< Default behaviour is that of generation_region_=="random".
    G4double gen_x_, gen_z_;                                        ///< Average GV coordinates. It is only used if generation_region_=="custom" is True.
    G4double gen_diameter_;                                         ///< Diameter of the circle where the GV could be randomly sampled. It is only used if generation_region_=="custom" is True.
    G4String path_to_inwards_dichroic_data_;                        ///< Absolute path to dichroic data file that is to be sampled for the light trying to enter the XA cavity
    G4String path_to_outwards_dichroic_data_;                       ///< Absolute path to dichroic data file that is to be sampled for the light trying to escape the XA cavity
    G4double world_extra_thickn_;                                   ///< Extra thickness for the surrounding box world to wrap the XArapuca

    ///---- General internal attributes ----///                     ///< These attributes are internal. They must not be set by the user.
    G4double overall_length_, overall_width_, overall_thickn_;      ///< Overall X-ARAPUCA dimensions
    G4double DFA_length_, DFA_width_;                               ///< Dichroic filters assembly transverse dimensions.
    G4double DF_length_, DF_width_;                                 ///< Dichroic filter transverse dimensions. All of the filters of the DFA have the same dimensions.

    ///---- config_code_==1 parameters ----///
    G4double plate_length_, plate_width_, plate_thickn_;            ///< WLS plate dimensions (X, Z and Y, respectively)
    G4double tunneling_probability_;                                ///< Probability that a photon tunnels through the surface of the WLSPlate (inwards or outwards), meaning that it is 
                                                                    ///< straightfoward-ly transmitted (no deflection due to Frensel refraction). The goal of this parameter is to let 
                                                                    ///< us model surface imperfections of the WLS plate.
    G4bool sipms_at_x_plus_;                                        ///< Whether to place sipms at the XA side which is contained in x>0.0
    G4bool sipms_at_x_minus_;                                       ///< Whether to place sipms at the XA side which is contained in x<0.0
    G4bool sipms_at_z_plus_;                                        ///< Whether to place sipms at the XA side which is contained in z>0.0
    G4bool sipms_at_z_minus_;                                       ///< Whether to place sipms at the XA side which is contained in z<0.0
    G4bool wrap_plate_in_ref_foil_;                                 ///< This parameter only makes a difference if double_sided_ is False and with_boards_ is False. In such case, whether 
                                                                    ///< to wrap the WLS plate in a reflective foil a la APEX. In any other case, the reflective wrap is not constructed.
                                                                    ///< This reflective wrap is mean to cover:
                                                                    ///<    - the small faces which are not instrumented with SiPMs,
                                                                    ///<    - the gaps in between the SiPMs in the small faces which are instrumented 
                                                                    ///<    - and one of the big faces of the WLS plate.
                                                                    ///< Hence, it does not makes sense for a double sided XArapuca. On the other hand, the reason for requiring 
                                                                    ///< with_boards_==False, is a computational one. The SiPMBoard volume, which is placed if with_boards_ is True, collides 
                                                                    ///< into the reflective wrap volume, giving a non-desired behaviour for some photons. A workaround is to place individual 
                                                                    ///< SiPMs, which will fit in the wrap holes, and whose volume won't collide into the wrap one.
    G4bool with_boards_;                                            ///< Whether the photosensors are mounted on a board
    G4bool add_blocks_between_sipms_;                               ///< Whether to add reflective blocks filling the gaps in between the board-mounted SiPMs. This option is only available 
                                                                    ///< if with_boards_ is true. These blocks are nearly coplanar with the SiPM surface: Their surface is 0.1 millimeters
                                                                    ///< behind the SiPMs surface - i.e. the SiPMs protrude 0.1 mm from these blocks. There are two reasons for this offset:
                                                                    ///<    1)  These blocks are part of a SiPMBoard object, whose logical volume is wrapped with a G4LogicalSkinSurface.
                                                                    ///<        On the other hand, if tunneling_probability_ is different from 0.0, then the secondary-WLS material logical 
                                                                    ///<        volume will also be wrapped with a G4LogicalSkinSurface. If that's the case, and gap_ is 0.0, then two 
                                                                    ///<        G4LogicalSkinSurface's match in space, giving an undefined behaviour which is prone to bugs. P.e. last time 
                                                                    ///<        this surfaces-match was observed, photons arriving to that surface from within a WLSPlate were absorbed.
                                                                    ///<    2)  Setting the blocks to be perfectly coplanar with the SiPM surface is not realistic either way, so we are
                                                                    ///<        safe introducing this 0.1 mm gap in between the secondary-WLS material and the SiPMs surface.
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
    G4bool cut_plate_;                                              ///< Whether to split the WLS plate up into two pieces.
                                                                    ///< WARNING: Although the application works fine, the visualization may break
                                                                    ///<          down when with_dimples_ and cut_plate_ are simultaneously true.
    G4double cut_angle_;                                            ///< This parameter only makes a difference if cut_plate_ is true. Angle of the cut, with respect to the Z axis (the axis along which the WLS plate is dz_ long)
    G4double cut_thickness_;                                        ///< This parameter only makes a difference if cut_plate_ is true. Thickness of the cut/crack that is carved from the plate.
    G4bool place_foil_at_the_cut_;                                  ///< This parameter only makes a difference if cut_plate_ is true. Whether to place (two) pieces of foil covering the inner faces of the cut.
                                                                    ///< If that's the case, then the material properties table of such foils is opticalprops::Vikuiti() and their thickness match cut_thickness_/5.
                                                                    ///< These foils are in optical contact contact with the WLS plate volume (meaning that there's no gab in between a foil and the WLS plate). Hence, if 
                                                                    ///< tunneling_probability_ is different from 0.0, the foils G4LogicalSkinSurface's match (in space) that of the WLS plate. However, in this case (contrary
                                                                    ///< to the one explained in the add_blocks_between_sipms_ attribute documentation) the photons which arrive at the cut behave as expected (are reflected with
                                                                    ///< the vikuiti reflection probability).

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
