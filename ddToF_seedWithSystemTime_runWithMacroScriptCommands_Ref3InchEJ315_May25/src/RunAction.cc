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
// $Id: RunAction.cc 74204 2013-10-01 07:04:43Z ihrivnac $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"
#include "Run.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{
  
  //=================================
  // Exercise 3 Step 1:
  // Create an output file containing
  // the histograms and the ntuple
  // Histograms:
  // 1D -> Scintillator bar 1 containing number of hits (50 bins from 0 to 50)
  // 1D -> Scintillator bar 2 containing number of hits (50 bins from 0 to 50)
  // 2D -> Scintillator bar 1 X vs Y binning: 50[-1000,1000] x 50[-300,300]
  // 2D -> Scintillator bar 2 X vs Y binning: 50[-1000,1000] x 50[-300,300]
  // Ntuple:
  // Integer Column: Number Hits Scintillator bar 1
  // Integer Column: Number Hits Scintillator bar 2
  // Double Column : Energy in Scintillator bar 1
  // Double Column : Energy in Scintillator bar 2
  // Double Column : Time in Scintillator bar 1
  // Double Column : Time in Scintillator bar 2
    
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  ///G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetVerboseLevel(0);
  ///analysisManager->SetFileName("300keVNeutronCalibration");
  analysisManager->SetFileName("two45MeVIsoNeutronsCalibration");

  // Book histograms, ntuple
  //

/*
  
  // Creating 1D histograms
  analysisManager
    ->CreateH1("ScintillatorBar1","Scintillator bar 1 # Hits", 50, 0., 15); // h1 Id = 0
  analysisManager
    ->CreateH1("ScintillatorBar2","Scintillator bar 2 # Hits", 50, 0., 15); // h1 Id = 1
  
  analysisManager
    ->CreateH1("E of First hit in Bar1","Deposited energy of first hit in bar 1", 50, 0.,2000); // h1 Id = 2
  analysisManager
    ->CreateH1("E of First hit in Bar2","Deposited energy of first hit in bar 2 ", 50, 0.,2000); // h1 Id = 3
  
  analysisManager
    ->CreateH1("Total E for all hits in Bar1","Deposited energy of all hits in bar 1", 50, 0.,5000); // h1 Id = 4
  analysisManager
    ->CreateH1("Total E of all hits in Bar2","Deposited energy of all hits in bar 2", 50, 0.,5000); // h1 Id = 5
  
  analysisManager
    ->CreateH1("Time of First hit in Bar1","Timing of first hit in bar 1", 1000, 0.,5000); // h1 Id = 6
  analysisManager
    ->CreateH1("Time of First hit in Bar2","Timing of first hit in bar 2 ", 1000, 0.,5000); // h1 Id = 7
  
  analysisManager
    ->CreateH1("Reconstructed E in Bar1","Reconstructed energy in bar 1", 100, 0.,50); // h1 Id = 8
  analysisManager
    ->CreateH1("Reconstructed E in Bar2","Reconstructed energy in bar 2 ", 100, 0.,50); // h1 Id = 9
  
  // Creating 2D histograms
  analysisManager                                                
    ->CreateH2("ScintillatorBar1 XY","Scintillator bar 1 X vs Y",           // h2 Id = 0
               50, -1000., 1000, 50, -300., 300.); 
  analysisManager                                                
    ->CreateH2("ScintillatorBar2 XY","Scintillator bar 2 X vs Y",           // h2 Id = 1
               50, -1500., 1500, 50, -300., 300.);


*/

  // Creating ntuple
  //
  analysisManager->CreateNtuple("three00keVNeutronCalibration", "hitsInfo");
  analysisManager->CreateNtupleIColumn("Bar1_N_of_Hits");  // column Id = 0
  analysisManager->CreateNtupleIColumn("Bar2_N_of_Hits");  // column Id = 1
  analysisManager->CreateNtupleDColumn("Bar1_Energy"); // column Id = 2
  analysisManager->CreateNtupleDColumn("Bar2_Energy"); // column Id = 3
  analysisManager->CreateNtupleDColumn("Bar1_Energy_1stHit"); // column Id = 4
  analysisManager->CreateNtupleDColumn("Bar2_Energy_1stHit"); // column Id = 5
  analysisManager->CreateNtupleDColumn("Bar1_Time1_1stHit");    // column Id = 6
  analysisManager->CreateNtupleDColumn("Bar2_Time2_1stHit");    // column Id = 7
  analysisManager->CreateNtupleDColumn("Bar1_Reconstructed_E");    // column Id = 8
  analysisManager->CreateNtupleDColumn("Bar2_Reconstructed_E");    // column Id = 9
  analysisManager->CreateNtupleIColumn("Reflector_N_of_Hits");    // column Id = 10
  analysisManager->CreateNtupleDColumn("Reflector_Time_1stHit");    // column Id = 11
  analysisManager->CreateNtupleDColumn("Reflector_Energy_1stHit");    // column Id = 12
  analysisManager->CreateNtupleDColumn("Reflector_Energy");    // column Id = 13
  analysisManager->CreateNtupleIColumn("water_N_of_Hits");  // column Id = 14 
  analysisManager->CreateNtupleDColumn("water_Time_1stHit");    // column Id = 15
  analysisManager->CreateNtupleDColumn("water_Energy_1stHit"); // column Id = 16
  analysisManager->CreateNtupleDColumn("water_Energy"); // column Id = 17
  analysisManager->CreateNtupleIColumn("waterTank_N_of_Hits");  // column Id = 18
  analysisManager->CreateNtupleDColumn("waterTank_Time_1stHit");    // column Id = 19
  analysisManager->CreateNtupleDColumn("waterTank_Energy_1stHit"); // column Id = 20
  analysisManager->CreateNtupleDColumn("waterTank_Energy"); // column Id = 21
  analysisManager->CreateNtupleDColumn("seed0"); // column Id = 22
  analysisManager->CreateNtupleDColumn("seed1"); // column Id = 23
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{

/*    
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  ///CLHEP::HepRandom::setTheSeed(24534575684783);
  long seeds[2];
  ///seeds[0] = 534524575674523;
  ///seeds[1] = 526345623452457;

  CLHEP::HepRandom::setTheSeeds(seeds);
  ///G4cout << "seed= " << CLHEP::HepRandom::getTheSeed() << G4endl;
  ///G4cout << "seed[0]= " << CLHEP::HepRandom::getTheSeed(*seeds) << G4endl;
  ///G4cout << "seed[1]= " << CLHEP::HepRandom::getTheSeed(*(seeds+1)) << G4endl;
*/
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //=================================
  // Exercise 3 Step 2:
  // Open output file at each new run

    // Open an output file
  // The default file name is set in RunAction::RunAction(),
  // it can be overwritten in a macro
  analysisManager->OpenFile();
  ///G4cout << G4endl << "Begin_of_RunAction" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//=================================
// Exercise 2:
// Collect the energy accumulated in the local Run
// And dump on screen the needed physics quantities
// for this particular run.
void RunAction::EndOfRunAction(const G4Run* run)
{
    const Run* myrun = dynamic_cast<const Run*>(run);
    if ( myrun )
    {
        G4int nEvets = myrun->GetNumberOfEvent();
        if ( nEvets < 1 )
        {
            G4ExceptionDescription msg;
            msg << "Run consists of 0 events";
            G4Exception("RunAction::EndOfRunAction()",
                        "Code001", JustWarning, msg);
            nEvets=1;
        }
        ///G4double bar1_ene = myrun->GetBar1Energy();
        ///G4double bar2_ene = myrun->GetBar2Energy();
        /*
        G4cout<<"Run["<<myrun->GetRunID()<<"] With: "<<nEvets<<"Events\n"
        <<" <E_Bar1>="<<G4BestUnit(bar1_ene/nEvets,"Energy")<<"\n"
        <<" <E_Bar2>="<<G4BestUnit(bar2_ene/nEvets,"Energy")<<"\n" 
        << G4endl;
        */
    } else {
        G4ExceptionDescription msg;
        msg << "Run is not of correct type, skypping analysis via RunAction";
        G4Exception("RunAction::EndOfRunAction()",
                    "Code001", JustWarning, msg);
    }
    
  //=================================
  // Exercise 3 Step 3:
  // Write and close output file
    
  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  ///G4cout << G4endl << "End_of_RunAction" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//=================================
// Exercise 1 Step 4:
// Implement this method,
// instead of creating a generic G4Run
// create a user-defined Run
G4Run* RunAction::GenerateRun() {
    return new Run;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
