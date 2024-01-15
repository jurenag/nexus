#include "SiPMBoard.h"
#include "SiPMMPPC.h"
#include "HamamatsuS133606050VE.h"
#include "HamamatsuS133605075HQR.h"
#include "FbkNuvHdCryoTT.h"
#include "BroadcomAFBRS4N44P044M.h"
#include "PerfectSiPMMPPC.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "Visibilities.h"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SystemOfUnits.hh>
#include <G4MultiUnion.hh>
#include <G4SubtractionSolid.hh>

namespace nexus {

  SiPMBoard::SiPMBoard():
    ///These are the DUNE SP TDR vol. IX dimensions
    board_length_(488.*mm), ///< X
    board_height_(8.*mm),   ///< Y
    board_thickn_(1.*mm),   ///< Z
    ref_phsensors_supports_(true),
    add_blocks_between_sipms_(false),
    base_id_(0),
    num_phsensors_(24),
    SiPM_code_(1)
  {
  }

  SiPMBoard::~SiPMBoard()
  {
  }

  G4int SiPMBoard::GetNumPhsensors()        const   { return num_phsensors_;}
  G4int SiPMBoard::GetSiPMCode()            const   { return SiPM_code_;    }
  G4double SiPMBoard::GetBoardLength()      const   { return board_length_; }
  G4double SiPMBoard::GetBoardHeight()      const   { return board_height_; }
  G4double SiPMBoard::GetOverallHeight()    const   
  { 
    SiPMMPPC* sipm_ptr = nullptr;
    if(SiPM_code_==1){
      sipm_ptr = new HamamatsuS133606050VE();
    }
    else if(SiPM_code_==2){
      sipm_ptr = new HamamatsuS133605075HQR();
    }
    else if(SiPM_code_==3){
      sipm_ptr = new FbkNuvHdCryoTT();
    }
    else if(SiPM_code_==4){
      sipm_ptr = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm_ptr = new PerfectSiPMMPPC();
    }
    G4double sipm_height = sipm_ptr->GetTransverseDim();
    delete sipm_ptr;
    return sipm_height>GetBoardHeight() ? sipm_height : GetBoardHeight();
  }

  G4double SiPMBoard::GetBoardThickness()   const   { return board_thickn_; }
  G4double SiPMBoard::GetOverallThickness() const
  {
    SiPMMPPC* sipm_ptr = nullptr;
    if(SiPM_code_==1){
      sipm_ptr = new HamamatsuS133606050VE();
    }
    else if(SiPM_code_==2){
      sipm_ptr = new HamamatsuS133605075HQR();
    }
    else if(SiPM_code_==3){
      sipm_ptr = new FbkNuvHdCryoTT();
    }
    else if(SiPM_code_==4){
      sipm_ptr = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm_ptr = new PerfectSiPMMPPC();
    }
    G4double sipm_thickn = sipm_ptr->GetThickness();
    delete sipm_ptr;
    return board_thickn_ + sipm_thickn;
  }
  G4double SiPMBoard::GetHasBlocks()        const   { return add_blocks_between_sipms_; }

  G4bool SiPMBoard::GeometryIsIllFormed() const
  {
      // Imagine the board in this configuration:
      //
      //    ----------------------------------------------
      //    |  | s |  | s |  | s |  | s |  | s |  | s |  |
      //    ----------------------------------------------  
      //
      // Where | s | depicts a SiPM. In this context, SiPMs height is allowed to be such that
      // SiPMs stand out over and under the board, but they cannot stand out the board via its 
      // lateral sides. I.e., this is allowed:
      //
      //           ___
      //        --|   |--
      //   ...    | s |   ...
      //        --|___|--
      //
      // But this is not:
      //
      //   ___------
      //  |_s_|      ...
      //      ------
      //
      // Therefore, there are not any constraints on the height nor the thickness of the board, but
      // only on its length depending the number of photosensors the board must allocate.


    SiPMMPPC* sipm_ptr = nullptr;
    if(SiPM_code_==1){
      sipm_ptr = new HamamatsuS133606050VE();
    }
    else if(SiPM_code_==2){
      sipm_ptr = new HamamatsuS133605075HQR();
    }
    else if(SiPM_code_==3){
      sipm_ptr = new FbkNuvHdCryoTT();
    }
    else if(SiPM_code_==4){
      sipm_ptr = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm_ptr = new PerfectSiPMMPPC();
    }
    G4double sipm_width = sipm_ptr->GetTransverseDim();
    delete sipm_ptr;

    if(sipm_width*num_phsensors_>GetBoardLength()){
        return true;
    }
    else{
        return false;
    }
  }


  void SiPMBoard::Construct()
  {

    if(GeometryIsIllFormed()){
        G4Exception("[SiPMBoard]", "Construct()", FatalException,
        "The provided parameters do not describe a feasible geometry.");
    }

    //ENCASING
    G4String encasing_name = "BOARD_ENCASING";
    const G4double encasing_height = GetOverallHeight();
    const G4double encasing_thickn = GetOverallThickness();
    G4Box* encasing_solid = 
        new G4Box(encasing_name, board_length_/2., encasing_height/2., encasing_thickn/2.);

    G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
    lAr->SetMaterialPropertiesTable(opticalprops::LAr());

    G4LogicalVolume* encasing_logic = 
        new G4LogicalVolume(encasing_solid, lAr, encasing_name);

    encasing_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

    this->SetLogicalVolume(encasing_logic);

    ConstructBoard(encasing_logic);
    ConstructSiPMs(encasing_logic);
    return;

  }

  void SiPMBoard::ConstructBoard(G4LogicalVolume* encasing_logic_vol)
  {

    //Check the availability of the mother volume
    if(!encasing_logic_vol){
        G4Exception("[SiPMBoard]", "ConstructBoard()", FatalException,
        "The encasing logical volume is not available.");
    }

    //FR4 piece (the board itself)
    G4String board_name = "BOARD";
    G4VSolid* board_solid = nullptr;

    if(!add_blocks_between_sipms_){
      board_solid = dynamic_cast<G4VSolid*>(new G4Box(board_name, board_length_/2., 
                                                                  board_height_/2., 
                                                                  board_thickn_/2.));
    }
    else{
      SiPMMPPC* sipm_ptr = nullptr;
      if(SiPM_code_==1){
        sipm_ptr = new HamamatsuS133606050VE();
      }
      else if(SiPM_code_==2){
        sipm_ptr = new HamamatsuS133605075HQR();
      }
      else if(SiPM_code_==3){
        sipm_ptr = new FbkNuvHdCryoTT();
      }
      else if(SiPM_code_==4){
        sipm_ptr = new BroadcomAFBRS4N44P044M();
      }
      else{
        sipm_ptr = new PerfectSiPMMPPC();
      }
      board_solid = dynamic_cast<G4VSolid*>(new G4Box(board_name, board_length_/2., 
                                                                  board_height_/2., 
                                                                  GetOverallThickness()/2.));

      sipm_ptr->Construct();  // I won't be able to get its logical volume unless the SiPMMPPC::Construct()
                              // method has been previosly called. Otherwise you get a segmentation fault.

      G4double tolerance = 1.*mm; // Some tolerance to prevent matching surfaces in boolean subtraction
      G4Box* sipm_wide_carvings = new G4Box(board_name, sipm_ptr->GetTransverseDim()/2., 
                                                        (board_height_+tolerance)/2.,             // Add tolerance at both, the height
                                                        (sipm_ptr->GetThickness()+tolerance)/2.); // dimension and the thickness dimension.

      G4MultiUnion* board_carvings = new G4MultiUnion("BOARD_CARVINGS");
      G4double x_pos, z_pos;
      G4Transform3D* transform_ptr = nullptr;
      G4RotationMatrix* rot = new G4RotationMatrix();

      for(G4int i=0; i<num_phsensors_; i++){
        x_pos = (-1.*board_length_/2.) + ((0.5 + i)*board_length_/num_phsensors_);
        transform_ptr = new G4Transform3D(*rot, G4ThreeVector(x_pos, 0., 0.));
        board_carvings->AddNode(*sipm_wide_carvings, *transform_ptr);
      }
      board_carvings->Voxelize();

      z_pos = (GetOverallThickness()/2.) -GetBoardThickness() -(sipm_ptr->GetThickness()/2.) -(tolerance/2.);
      board_solid = dynamic_cast<G4VSolid*>(new G4SubtractionSolid( board_name, board_solid,
                                                                    board_carvings, 
                                                                    nullptr, G4ThreeVector(0., 0., z_pos)));
      delete sipm_ptr;
    }

    G4LogicalVolume* board_logic = 
        new G4LogicalVolume(board_solid, materials::FR4(), board_name);

    G4VisAttributes board_col = nexus::White();
    //board_col.SetForceSolid(true);
    board_logic->SetVisAttributes(board_col);

    //VIKUITI coating for the board
    const G4String bc_name = "BOARD_COATING";
    G4OpticalSurface* board_coating = 
      new G4OpticalSurface(bc_name, unified, ground, dielectric_metal, 1);
    
    board_coating->SetMaterialPropertiesTable(opticalprops::specularspikeVIKUITI());
    new G4LogicalSkinSurface(bc_name, board_logic, board_coating); 

    G4double pos;
    if(!add_blocks_between_sipms_){
      pos = (GetOverallThickness()-board_thickn_)/2.;
    }
    else{         // In this case, the geometric
      pos = 0.0;  // center of board_solid matches
    }             // that of the encasing volume

    //Place it
    G4double aux = GetOverallThickness();
    new G4PVPlacement(nullptr, 
        G4ThreeVector(0., 0., pos),
        board_logic, "COATED_BOARD", 
        encasing_logic_vol,
        false, 0, true);

    return;
  }

  void SiPMBoard::ConstructSiPMs(G4LogicalVolume* encasing_logic_vol)
  {

    //Check the availability of the mother volume
    if(!encasing_logic_vol){
        G4Exception("[SiPMBoard]", "ConstructSiPMs()", FatalException,
        "The encasing logical volume is not available.");
    }

    SiPMMPPC * sipm = nullptr;
    if(SiPM_code_==1){
      sipm = dynamic_cast<HamamatsuS133606050VE*>(sipm);
      sipm = new HamamatsuS133606050VE();
    }
    else if(SiPM_code_==2){
      sipm = dynamic_cast<HamamatsuS133605075HQR*>(sipm);
      sipm = new HamamatsuS133605075HQR();
    }
    else if(SiPM_code_==3){
      sipm = dynamic_cast<FbkNuvHdCryoTT*>(sipm);
      sipm = new FbkNuvHdCryoTT();
    }
    else if(SiPM_code_==4){
      sipm = dynamic_cast<BroadcomAFBRS4N44P044M*>(sipm);
      sipm = new BroadcomAFBRS4N44P044M();
    }
    else{
      sipm = dynamic_cast<PerfectSiPMMPPC*>(sipm);
      sipm = new PerfectSiPMMPPC();
    }
    sipm->SetReflectiveSupports(ref_phsensors_supports_);  
    sipm->Construct();
    G4double sipm_thickn = sipm->GetThickness();
    G4LogicalVolume* sipm_logic_vol = sipm->GetLogicalVolume();

    if (!sipm_logic_vol) {
      G4Exception("[SiPMBoard]", "ConstructSiPMs()",
                  FatalException, "Null pointer to logical volume.");
    }

    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateX(90.*deg);

    G4int phsensor_id = base_id_;
    for (G4int i=0; i<num_phsensors_; ++i) {
      G4ThreeVector pos((-1.*board_length_/2.) + ((0.5 + i)*board_length_/num_phsensors_),
                        0.,
                        (GetOverallThickness()/2.) -GetBoardThickness() - (sipm_thickn/2.));

      new G4PVPlacement(rot, pos,
                        sipm_logic_vol, sipm->GetModel(),
                        encasing_logic_vol, false, phsensor_id, true);

      phsensor_id += 1;
    }
    
    return;

  }

}

