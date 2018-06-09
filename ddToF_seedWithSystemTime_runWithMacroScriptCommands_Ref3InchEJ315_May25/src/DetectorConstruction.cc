//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "CellParameterisation.hh"
#include "HodoscopeSD.hh"
#include "CollimatorSD.hh"
#include "G4SubtractionSolid.hh"
#include "waterSD.hh"
#include "waterTankSD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

///G4ThreadLocal G4FieldManager* DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fMessenger(0),
  fHodoscope1Logical(0), fHodoscope2Logical(0),
  fCollimatorLogical(0),fCollimatorContainerLogical(0),
  fWaterLogical(0), 
  fWaterTankLogical(0),
  fTwommStealDDLogical(0), fFourmmLeadDDLogical(0),
  fVisAttributes(),
  fArmAngle(0.*deg), fArmRotation(0), fSecondArmPhys(0)

{
    fArmRotation = new G4RotationMatrix();
    fArmRotation->rotateY(fArmAngle);
    
    // define commands for this class
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete fArmRotation;
    delete fMessenger;
    
    for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
    {
      delete fVisAttributes[i];
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Construct materials
    ConstructMaterials();
    G4Material* aluminium = G4Material::GetMaterial("G4_Al");
    G4Material* air = G4Material::GetMaterial("G4_AIR");
    G4Material* lead = G4Material::GetMaterial("G4_Pb");
    //G4Material* argonGas = G4Material::GetMaterial("_Ar");
    G4Material* liquidWater = G4Material::GetMaterial("G4_WATER");
    ///G4Material* scintillator 
    ///  = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    ///G4Material* polyethylene = G4Material::GetMaterial("G4_POLYETHYLENE");
    G4Material* steel_304 = G4Material::GetMaterial("stainless steel 304");
    ///G4Material* D2O = G4Material::GetMaterial("heavy water");
    G4Material* materialBC501A = G4Material::GetMaterial("materialBC501A");
    G4Material* EJ315 = G4Material::GetMaterial("EJ315");
    
    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = false; ///true;

    // geometries --------------------------------------------------------------
    // experimental hall (world volume)
    G4VSolid* worldSolid 
      = new G4Box("worldBox",10.*m,10.*m,10.*m);
    G4LogicalVolume* worldLogical
      = new G4LogicalVolume(worldSolid,air,"worldLogical");
    G4VPhysicalVolume* worldPhysical
      = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                          false,0,checkOverlaps);


    ///~~~~~~~~~~~~~~~~~~ first arm, reserved for neutron source ~~~~~~~~~~~~~~~~~
    G4VSolid* firstArmSolid 
      = new G4Box("firstArmBox",15.*cm,15.*cm,11.*cm);/// 11 is to avoid overlap of first & secondArmLogical 
    G4LogicalVolume* firstArmLogical
      = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,32.85*cm),firstArmLogical,/// 23.15 cm is the Z coordinate of the central point of the reflector in the G4 coordinate system.
                      "firstArmPhysical",worldLogical,
                      false,0,checkOverlaps);

   
    /// Reflector in first arm, named as collimator, material is EJ-315.
    G4VSolid* collimatorSolid 
      ///= new G4Tubs("collimatorSlice",0.0*cm,3.*cm,0.1*cm,0.,2.0*pi);/// phi = 6.0 cm
      = new G4Tubs("collimatorSlice",0.0*cm,3.8*cm,0.1*cm,0.,2.0*pi);/// 3", phi = 7.6 cm
    fCollimatorLogical
      = new G4LogicalVolume(collimatorSolid, /*D2O*/ EJ315,"collimatorLogical");
    ///for (G4int i=0;i<100;i++)/// Length: 100 * 0.2 cm = 20.0 cm
    for (G4int i=0;i<38;i++)/// Length: 38 * 0.2 cm = 7.6 cm
    {
        ///G4double z1 = ( (i-50)*0.2 + 0.1 )*cm;
        ///new G4PVPlacement(0,G4ThreeVector(0.5,0.,z1),fCollimatorLogical,
        G4double z1 = ( (i-19)*0.2 + 0.1 )*cm;
        new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),fCollimatorLogical,
                          "collimatorPhysical",firstArmLogical,
                          false,i,checkOverlaps);
    }
    
    /// Reflector container in first arm, material is polyethylene. 
    G4VSolid* collimatorContainerSolid 
      ///= new G4Tubs("collimatorSlice",3.0*cm,3.1*cm,0.1*cm,0.,2.0*pi);
      = new G4Tubs("collimatorSlice",3.8*cm,3.95*cm,0.1*cm,0.,2.0*pi);
    fCollimatorContainerLogical
      = new G4LogicalVolume(collimatorContainerSolid, aluminium /*polyethylene*/,"collimatorContainerLogical");
    ///for (G4int i=0;i<100;i++)/// Length: 100 * 0.2 cm = 20.0 cm
    for (G4int i=0;i<38;i++)/// Length: 38 * 0.2 cm = 7.6 cm
    {
        ///G4double z1 = ( (i-50)*0.2 + 0.1 )*cm;
        ///new G4PVPlacement(0,G4ThreeVector(0.5,0.,z1),fCollimatorContainerLogical,
        G4double z1 = ( (i-19)*0.2 + 0.1 )*cm;
        new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),fCollimatorContainerLogical,
                          "collimatorContainerPhysical",firstArmLogical,
                          false,i,checkOverlaps);
    }
    
    /// ~~~~~~~~~~~~~~~~~ second arm, water, water tank ~~~~~~~~~~~~~~~~~~~~~~
    
    G4RotationMatrix* rotD90X = new G4RotationMatrix();
    rotD90X->rotateX(90.*deg);
    G4ThreeVector zTrans(0, 0, 0);

    
    G4VSolid* secondArmSolid 
      = new G4Box("secondArmBox",1.5*m,1.5*m,146.*cm);
    G4LogicalVolume* secondArmLogical
      = new G4LogicalVolume(secondArmSolid,air,"secondArmLogical");
    fSecondArmPhys
      = new G4PVPlacement(fArmRotation,G4ThreeVector(0.,0.,190.*cm),secondArmLogical,
                          "fSecondArmPhys",worldLogical,
                          false,0,checkOverlaps);
    
    /// waterModerator : cylinder vertical placed with a hole going through its center radially
    G4VSolid* waterSolid 
      = new G4Tubs("waterModerator",0.0*cm,100.*cm,100.*cm,0.,2.0*pi);

    G4VSolid* couduitShapeSolid /// to be subtracted to get space for plactic conduit and scintillators
      = new G4Tubs("conduitShape",0.0*cm,10.*cm,101.*cm,0.,2.0*pi);/// 2cm shorter than 1 m to avoid overlap

    G4SubtractionSolid* subtraction =
    new G4SubtractionSolid("water Tank - conduitShape", waterSolid, couduitShapeSolid, rotD90X, zTrans);

    fWaterLogical
      = new G4LogicalVolume(subtraction,liquidWater,"waterLogical");
      ///= new G4LogicalVolume(subtraction,air,"waterLogical");
    new G4PVPlacement(rotD90X,G4ThreeVector(0.,0.,0.*cm),fWaterLogical,
                          "waterPhysical",secondArmLogical,
                          false,0,checkOverlaps);

    ///water Tank 
    G4VSolid* waterTankSolid 
      ///= new G4Tubs("water Tank",100.0*cm,101.8*cm,0.25*cm,0.,2.0*pi);/// water tank thick is 1.8 cm.
      = new G4Tubs("water Tank",100.0*cm,100.5*cm,0.25*cm,0.,2.0*pi);/// water tank thick is 0.5 cm.
    fWaterTankLogical
      ///= new G4LogicalVolume(waterTankSolid,  polyethylene ,"waterTankLogical");
      = new G4LogicalVolume(waterTankSolid,  steel_304 ,"waterTankLogical");
      ///= new G4LogicalVolume(waterTankSolid,  air ,"waterTankLogical");  
    for (G4int i=0;i<400;i++)/// 0.5 * 400 = 2 m height
    {
        G4double y1 = ( (i-200)*0.5 + 0.25 )*cm;
        new G4PVPlacement(rotD90X,G4ThreeVector(0.,y1,0.),fWaterTankLogical,
                          "waterTankPhysical",secondArmLogical,
                          false,i,checkOverlaps);
    }

    ///scintillators: BC501A
    G4VSolid* BC501ASolid_3Inch /// 3 inches BC501A  
      = new G4Tubs("BC501A",0.0*cm,4.0*cm,0.1*cm,0.,2.0*pi); 

    ///scintillator in the middle of water
    fHodoscope1Logical
      = new G4LogicalVolume(BC501ASolid_3Inch, materialBC501A ,"hodoscope1Logical");/// air is tmp
    for (G4int i=0;i<40;i++)/// 40 * 0.2 cm = 8 cm height
    {
        G4double y1 = ( (i-20)*0.2 + 0.1 )*cm;
        new G4PVPlacement(rotD90X,G4ThreeVector(0.,y1,0.),fHodoscope1Logical,
                          "hodoscope1Physical",secondArmLogical,
                          false,i,checkOverlaps);
    }


    G4VSolid* BC501ASolid_5Inch /// 5 inches BC501A
      ///= new G4Tubs("BC501A",0.0*cm,6.5*cm,13.*cm,0.,2.0*pi); 
      = new G4Tubs("BC501A",0.0*cm,6.5*cm,0.1*cm,0.,2.0*pi); 

    ///scintillator outside of water
    fHodoscope2Logical
      = new G4LogicalVolume(BC501ASolid_5Inch, materialBC501A ,"hodoscope2Logical");
    for (G4int i=0;i<64;i++)/// 0.2cm * 64 = 12.8 cm height
    {
        G4double y1 = ( (i-32)*0.2 + 0.1 )*cm;
        new G4PVPlacement(rotD90X,G4ThreeVector(0.,y1,120.*cm),fHodoscope2Logical,
                          "hodoscope2Physical",secondArmLogical,
                          false,i,checkOverlaps);
    }



    /// 2 mm steal of DD generator
    G4VSolid* twommStealDD 
      = new G4Tubs("2 mm steal DD",7.9*cm,8.1*cm,0.5*cm,0.,2.0*pi); 

    fTwommStealDDLogical
      = new G4LogicalVolume(twommStealDD,steel_304 ,"twommStealDDLogical");
    for (G4int i=0;i<30;i++)/// 1.0cm * 30 = 30.0 cm length
    {
        G4double y1 = ( (i-15)*1.0 + 0.5)*cm;
        ///new G4PVPlacement(rotD90X,G4ThreeVector(-12.*cm,y1,-136.*cm),fTwommStealDDLogical,
        /// To get 160 degree scattering angle, DD moves forward to the water tank 12 cm.
        new G4PVPlacement(rotD90X,G4ThreeVector(-12.*cm,y1,-124.*cm),fTwommStealDDLogical,
                          "twommStealDDPhysical",secondArmLogical,
                          false,i,checkOverlaps);
    }

    /// 2 mm lead of DD generator
    G4VSolid* fourmmleadDD 
      = new G4Tubs("2 mm lead DD",8.11*cm,8.5*cm,0.5*cm,0.,2.0*pi); 

    fFourmmLeadDDLogical
      = new G4LogicalVolume(fourmmleadDD, lead ,"fourmmLeadDDLogical");
    for (G4int i=0;i<30;i++)/// 1.0cm * 30 = 30.0 cm length
    {
        G4double y1 = ( (i-15)*1.0 + 0.5)*cm;
        ///new G4PVPlacement(rotD90X,G4ThreeVector(-12.*cm,y1,-136.*cm),fFourmmLeadDDLogical,
        /// To get 160 degree scattering angle, DD moves forward to the water tank 12 cm.
        new G4PVPlacement(rotD90X,G4ThreeVector(-12.*cm,y1,-124.*cm),fFourmmLeadDDLogical,
                          "fourmmLeadDDPhysical",secondArmLogical,
                          false,i,checkOverlaps);
    }


    //
    // visualization attributes ------------------------------------------------
    // Step 6: uncomment visualization attributes of the newly created volumes
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    visAttributes->SetVisibility(false);
    worldLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(1.0,0.5,0.5));
    visAttributes->SetVisibility(false);
    firstArmLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    visAttributes->SetVisibility(false);
    secondArmLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.8,0.30,0.20));
    ///visAttributes->SetVisibility(false);
    fHodoscope1Logical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.0,0.5,0.9));
    ///visAttributes->SetVisibility(false);
    fHodoscope2Logical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.0,0.8888,0.0));
    visAttributes->SetVisibility(false);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.5888,0.8888,0.5));
    ///visAttributes->SetVisibility(false);
    fCollimatorLogical->SetVisAttributes(visAttributes);
    fCollimatorContainerLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.0));
    fWaterLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
    ///visAttributes->SetVisibility(false);
    fWaterTankLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.3, 0.0, 0.9));
    ///visAttributes->SetVisibility(false);
    fTwommStealDDLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.0, 0.6, 0.9));
    ///visAttributes->SetVisibility(false);
    fFourmmLeadDDLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
     
    // return the world physical volume ----------------------------------------
    
    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
   
    // sensitive detectors -----------------------------------------------------
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    G4VSensitiveDetector* hodoscope1 
      = new HodoscopeSD(SDname="/hodoscope1");
    SDman->AddNewDetector(hodoscope1);
    fHodoscope1Logical->SetSensitiveDetector(hodoscope1);

    G4VSensitiveDetector* hodoscope2 
      = new HodoscopeSD(SDname="/hodoscope2");
    SDman->AddNewDetector(hodoscope2);
    fHodoscope2Logical->SetSensitiveDetector(hodoscope2);
     
    G4VSensitiveDetector* collimator 
      = new CollimatorSD(SDname="/collimator");///"CollimatorSD" typoed as "HodoscopeSD", wasted around one day.
    SDman->AddNewDetector(collimator);
    fCollimatorLogical->SetSensitiveDetector(collimator);
    
    G4VSensitiveDetector* water 
      = new waterSD(SDname="/water");
    SDman->AddNewDetector(water);
    fWaterLogical->SetSensitiveDetector(water);

    G4VSensitiveDetector* waterTank 
      = new waterTankSD(SDname="/waterTank");
    SDman->AddNewDetector(waterTank);
    fWaterTankLogical->SetSensitiveDetector(waterTank);

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();

    // Air 
    nistManager->FindOrBuildMaterial("G4_AIR");
  
    // Argon gas
    nistManager->FindOrBuildMaterial("G4_Ar");
    // With a density different from the one defined in NIST
    // G4double density = 1.782e-03*g/cm3; 
    // nistManager->BuildMaterialWithNewDensity("_Ar","G4_Ar",density);
    // !! cases segmentation fault

    // liquid water 
    nistManager->FindOrBuildMaterial("G4_WATER");
  
    // Scintillator
    // (PolyVinylToluene, C_9H_10)
    ///nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    // CsI
    nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
    
    // Lead
    nistManager->FindOrBuildMaterial("G4_Pb");
    
    // aluminium 
    nistManager->FindOrBuildMaterial("G4_Al");
    
    // polyethylene
    ///nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");

    //Define Stainless-steel-304 
    G4double z, Z, A, b, fractionmass, densitySteel,density;
    G4String name,  symbol;
    G4int ncomponents/*, natoms*/;
    A  =  54.94*g/mole;
    G4Element* elMn   =  new G4Element("Manganese","Mn",Z = 25.,A);
    A = 28.09*g/mole;
    G4Element* elSi  = new G4Element("Silicon","Si",Z = 14.,A);
    A = 52.00*g/mole;
    G4Element* elCr  = new G4Element("Chromium","Cr",Z = 24.,A);
    A = 58.70*g/mole;
    G4Element* elNi  = new G4Element("Nickel","Ni",Z = 28.,A);
    A = 55.85*g/mole;
    G4Element* elFe  = new G4Element("Iron","Fe",Z = 26.,A);
    A = 12.0*g/mole;
    G4Element* elC = new G4Element(name="Carbon",symbol="C" , z= 6., A);

    densitySteel = 7.999*g/cm3 ;
    ///steel_304 = new G4Material(name="steel_304",densitySteel,6);
    G4Material* steel_304 = new G4Material(name="stainless steel 304",densitySteel,ncomponents=6);
    steel_304->AddElement(elMn, 0.02);
    steel_304->AddElement(elSi, 0.01);
    steel_304->AddElement(elCr, 0.19);
    steel_304->AddElement(elNi, 0.10);
    steel_304->AddElement(elFe, 0.6792);
    steel_304->AddElement(elC, 0.0008);

    /// material BC501A
    b = 1.0*g/mole;
    G4Element* elH = new G4Element(name="Hydrogen" ,symbol="H" , z= 1., b);
    density = 0.874*g/cm3;
    G4Material* materialBC501A = new G4Material(name="materialBC501A",density,ncomponents=2);
    materialBC501A->AddElement(elC, fractionmass=45.2*perCent);
    materialBC501A->AddElement(elH, fractionmass=54.8*perCent);        
        
    /// material D2O
   /* G4Element* O  = new G4Element("Oxygen" ,"O" , 8., 16.00*g/mole);
    G4Isotope* H2 = new G4Isotope("H2",1,2,2.0*g/mole);
    G4Element* D  = new G4Element("Heavy Water", "D", 1);
    D->AddIsotope(H2, 100*perCent);  
    G4Material* D2O = new G4Material("heavy water", 1.107*g/cm3, ncomponents=2,
                          kStateLiquid, 293.15*kelvin, 1*atmosphere);
    D2O->AddElement(D, natoms=2);
    D2O->AddElement(O, natoms=1);
   */ 
    /// material EJ315 
    density = 0.953*g/cm3; ///cauculated by number densities of D,H and C.
    G4Isotope* H2 = new G4Isotope("H2",1,2,2.0*g/mole);
    G4Element* D  = new G4Element("Heavy Water", "D", 1);
    D->AddIsotope(H2, 100*perCent);  
    G4Material* EJ315 = new G4Material(name="EJ315",density,ncomponents=3);

    EJ315->AddElement(elC, fractionmass=85.8*perCent);
    EJ315->AddElement(elH, fractionmass=0.05*perCent);        
    EJ315->AddElement(D, fractionmass=14.15*perCent);        
    
   /* b = 2.0*g/mole;
    G4Element* elD = new G4Element(name="Deuterium" ,symbol="D" , z= 1., b);
    density = 1.107*mg/cm3;
    G4Material* D2O = new G4Material(name="heavy water",density,ncomponents=2);
    D2O->AddElement(elD, fractionmass=0.2*perCent);
    D2O->AddElement(elO, fractionmass=0.8*perCent);        
       */ 
    /// Vacuum "Galactic"
    ///nistManager->FindOrBuildMaterial("G4_Galactic");

    // Vacuum "Air with low density"
    // G4Material* air = G4Material::GetMaterial("G4_AIR");
    // G4double density = 1.0e-5*air->GetDensity();
    // nistManager
    //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

    ///G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetArmAngle(G4double val)
{
    if (!fSecondArmPhys)
    {
        G4cerr << "Detector has not yet been constructed." << G4endl;
        return;
    }
    
    fArmAngle = val;
    *fArmRotation = G4RotationMatrix();  // make it unit vector
    fArmRotation->rotateY(fArmAngle);
    ///G4double x = -5.*m * std::sin(fArmAngle);
    ///G4double z = 5.*m * std::cos(fArmAngle);
    G4double x = -1.*m * std::sin(fArmAngle);
    G4double z = 1.*m * std::cos(fArmAngle);
    fSecondArmPhys->SetTranslation(G4ThreeVector(x,0.,z));
    
    // tell G4RunManager that we change the geometry
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
    // Define //detector command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this, 
                                        "/tutorial/detector/",
                                        "Detector control");

    // armAngle command
    G4GenericMessenger::Command& armAngleCmd
      = fMessenger->DeclareMethodWithUnit("armAngle","deg",
                                  &DetectorConstruction::SetArmAngle, 
                                  "Set rotation angle of the second arm.");
    armAngleCmd.SetParameterName("angle", true);
    armAngleCmd.SetRange("angle>=0. && angle<180.");
    armAngleCmd.SetDefaultValue("30.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
