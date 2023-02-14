#include "SiPMBoard.h"
#include "SiPMMPPC.h"
#include "HamamatsuS133606050VE.h"
#include "HamamatsuS133605075HQR.h"
#include "FbkNuvHdCryoTT.h"

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


namespace nexus {

  SiPMBoard::SiPMBoard():
    ///These are the DUNE SP TDR vol. IX dimensions
    board_length_(488.*mm), ///< X
    board_height_(8.*mm),   ///< Y
    board_thickn_(1.*mm),   ///< Z
    ref_phsensors_supports_(true),
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
    else{
      sipm_ptr = new FbkNuvHdCryoTT();
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
    else{
      sipm_ptr = new FbkNuvHdCryoTT();
    }
    G4double sipm_thickn = sipm_ptr->GetThickness();
    delete sipm_ptr;
    return board_thickn_ + sipm_thickn;
  }

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
    else{
      sipm_ptr = new FbkNuvHdCryoTT();
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
    const G4double encasing_thickn = GetOverallThickness();
    G4Box* encasing_solid = 
        new G4Box(encasing_name, board_length_/2., board_height_/2., encasing_thickn/2.);

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
    G4Box* board_solid =
        new G4Box(board_name, board_length_/2., board_height_/2., board_thickn_/2.);
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

    //Place it
    G4double aux = GetOverallThickness();
    new G4PVPlacement(nullptr, 
        G4ThreeVector(0., 0., (aux-board_thickn_)/2.),
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
    else{
      sipm = dynamic_cast<FbkNuvHdCryoTT*>(sipm);
      sipm = new FbkNuvHdCryoTT();
    }
    sipm->SetReflectiveSupports(ref_phsensors_supports_);  
    sipm->Construct();
    G4double sipm_thickn = sipm->GetThickness();
    G4LogicalVolume* sipm_logic_vol = sipm->GetLogicalVolume();

    G4VisAttributes sipm_col = nexus::Red();
    sipm_col.SetForceSolid(true);
    sipm_logic_vol->SetVisAttributes(sipm_col);

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

