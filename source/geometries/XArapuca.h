#ifndef X_ARAPUCA
#define X_ARAPUCA

#include "GeometryBase.h"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4GenericMessenger;
class G4UserLimits;

namespace nexus {

  /// X-Arapuca 

  class XArapuca: public GeometryBase
  {
  public:
    ///Constructor
    XArapuca();
    ///Destructor
    ~XArapuca();

    void Construct();                   ///< Constructs the geometry

  private:

    void ConstructWLSPlate(G4VPhysicalVolume*) const;         ///< Called by Construct(). Adds the WLS plate.
    void ConstructPhotosensors(G4VPhysicalVolume*) const;     ///< DEPRECTAED. Called by Construct(). Adds (floating) SiPMs.
    void ConstructBoards(G4VPhysicalVolume*) const;           ///< Called by Construct(). Adds the board-mounted SiPMs.
    void ConstructReflectiveBox(G4VPhysicalVolume*) const;    ///< Called by Construct(). Adds the reflective box that encloses the X-ARAPUCA geometry.
    void ConstructCollectors(G4VPhysicalVolume*) const;       ///< Called by Construct(). Adds the artificial volumes that replace the dichoric filters.
    void ConstructDichroicFilters(G4VPhysicalVolume*) const;  ///< Called by Construct(). Adds the dichroic filters.
    
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:

    G4double internal_length_, internal_width_, internal_thickn_;   ///<Internal dimensions of the reflector cavity
    G4double plate_length_, plate_width_, plate_thickn_;            ///<WLS plate thickness
    G4double case_thickn_;                                          ///<Reflective foils thickness
    G4double df_thickn_;                                            ///<Dichroic filter thickness
    G4double gap_;                                                  ///<Gap between photosensors and WLS plate
    G4int num_phsensors_;                                           ///<Number of photosensors per long side
    G4bool ref_phsensors_supports_;                                 ///<Whether photosensors supports are reflective (the tiny FR4 box that supports the SiPM)
    G4bool with_boards_;                                            ///<Whether the photosensors are mounted on a board
    G4bool double_sided_;                                           ///<Whether there are dichroic filters on both sides of the WLS plate
    G4bool collectors_are_reflective_;                              ///<Whether the collectores that replace the dichroic filters are reflective or not                                      
    G4bool random_generation_vertex_;                               ///<Whether the generation vertex is randomly sample over the inner dichroic surface. If not, the generation vertex matches the geometric center of such surface.
    G4bool remove_filters_;                                         ///<Whether to remove the filters or not
    G4String path_to_dichroic_data_;                                ///<Absolute path to dichroic data file
    G4double world_extra_thickn_;                                   ///<Extra thickness for the surrounding box world to wrap the XArapuca

    G4GenericMessenger* msg_;                                           ///< Messenger for the definition of control commands
    G4UserLimits* ul_;

    G4bool geometry_is_ill_formed();                                    ///<Check whether the WLS plate fits within the XArapuca cavity and the given number of photosensors fit along the XArapuca long side
  };
}

#endif