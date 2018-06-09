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
// $Id: EventAction.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "HodoscopeHit.hh"
#include "CollimatorHit.hh"
#include "waterHit.hh"
#include "waterTankHit.hh"
#include "Analysis.hh"
///#include "PrimaryGeneratorAction.hh"
#include "exGPSPrimaryGeneratorAction.hh"


#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleGun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction(), 
  ///fHHC1ID(-1),
  ///fHHC2ID(-1),
  fBAR1HCID(-1),///Bar1
  fBAR2HCID(-1),///Bar2
  fCLMHCID(-1),///Collimator Hit Collection ID
  fWATERHCID(-1),///water
  fWATERTANKHCID(-1),///water Tank
  genAction(0)
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(0);
  seeds=NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* )
{

  ///G4int evtId = evt->GetEventID();
  ///G4cout << G4endl << "event ID in EventAction.cc is: " << evtId << G4endl;


/*  if (!genAction)
    genAction = dynamic_cast<const exGPSPrimaryGeneratorAction*>
      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  /// grab seeds
  seeds = genAction->GetEventSeeds();

  long seed0 = *seeds;
  long seed1 = *(seeds+1);
  G4cout << "EventAction.cc, seed0 is: " << seed0 << ", seed2 is:"  << seed2 << G4endl;
*/


    if (fBAR1HCID==-1) {
      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
      fBAR1HCID = sdManager->GetCollectionID("hodoscope1/hodoscopeColl");
      fBAR2HCID = sdManager->GetCollectionID("hodoscope2/hodoscopeColl");
      fCLMHCID = sdManager->GetCollectionID("collimator/collimatorColl"); /// fCLMHCID = -1
      fWATERHCID = sdManager->GetCollectionID("water/waterColl");
      fWATERTANKHCID = sdManager->GetCollectionID("waterTank/waterTankColl");
    }
}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    // =============================================
    // Exercise 3
    // Print on screen the hits of the hodoscope
    // Step 1: Get the hits collection of this event

    ///G4cout << G4endl << "begin of EndOfEventAction()" << G4endl; 
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) 
    {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found.\n"; 
        G4Exception("EventAction::EndOfEventAction()",
                    "Code001", JustWarning, msg);
        return;
    }   

    // Step 2: Using the memorised IDs get the collections
    // corresponding to the two hodoscopes


    // Get hits collections 
    HodoscopeHitsCollection* hHC1 
      = static_cast<HodoscopeHitsCollection*>(hce->GetHC(fBAR1HCID));
      
    HodoscopeHitsCollection* hHC2 
      = static_cast<HodoscopeHitsCollection*>(hce->GetHC(fBAR2HCID));
      
    CollimatorHitsCollection* hHCLM 
      = static_cast<CollimatorHitsCollection*>(hce->GetHC(fCLMHCID));
      
    waterHitsCollection* wtrHC ///water, not "waterHC" to avoid same name in Run.cc
      = static_cast<waterHitsCollection*>(hce->GetHC(fWATERHCID));
      
    waterTankHitsCollection* wtrTkHC ///water Tank, not "waterTankHC" to avoid same name in Run.cc
      = static_cast<waterTankHitsCollection*>(hce->GetHC(fWATERTANKHCID));
      
    ///if ( (!hHC1) || (!hHC2)  )
    if ( (!hHC1) || (!hHC2) || (!hHCLM)  || (!wtrHC) || (!wtrTkHC) )
    {
        G4ExceptionDescription msg;
        msg << "Some of hits collections of this event not found.\n"; 
        G4Exception("EventAction::EndOfEventAction()",
                    "Code001", JustWarning, msg);
        return;
    }



    //=================================
    // Fill histograms & ntuple
    
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
    // Fill histograms
 
    G4int n_hit_1 = hHC1->entries();///n_hit_1 is the number of hits per incident particle in bar 1.  

    for (G4int i=0;i<n_hit_1;i++)
    {
       HodoscopeHit* hit = (*hHC1)[i];
       G4ThreeVector Pos = hit->GetPos();
    }
 
 
    G4int n_hit_2 = hHC2->entries();
    ///analysisManager->FillH1(1, n_hit_2);

    for (G4int i=0;i<n_hit_2;i++)
    {
       HodoscopeHit* hit = (*hHC2)[i];
       G4ThreeVector Pos = hit->GetPos();
       ///analysisManager->FillH2(1, Pos.x()/cm, Pos.y()/cm);
    }
    
 
        
    // Fill ntuple
    

    analysisManager->FillNtupleIColumn(0, n_hit_1);
    analysisManager->FillNtupleIColumn(1, n_hit_2);
    
    // Bar1 Energy
    G4int totalBar1Hit = 0;
    G4double totalBar1E = 0.;
    G4double E_1st_hit_Bar1 = 0.;
    G4double T_1st_hit_Bar1 = 0.;
    ///for (G4int i=0;i<25;i++)
    for (G4int i=0;i<n_hit_1;i++)
    {
        HodoscopeHit* hit = (*hHC1)[i];
        G4double eDep = hit->GetEdep();
        if (eDep>0.)
        {
        T_1st_hit_Bar1 = (*hHC1)[0]->GetTime();
        E_1st_hit_Bar1 = (*hHC1)[0]->GetEdep();
            totalBar1Hit++;
            totalBar1E += eDep;
          }
    }
    analysisManager->FillNtupleDColumn(2, totalBar1E/keV);


     //Bar2 Energy
    G4int totalBar2Hit = 0;
    G4double totalBar2E = 0.;
    G4double E_1st_hit_Bar2 = 0.;
    G4double T_1st_hit_Bar2 = 0.;
    for (G4int i=0;i<n_hit_2;i++)
    {
        HodoscopeHit* hit = (*hHC2)[i];
        G4double eDep = hit->GetEdep();
        if (eDep>0.)
        {
        T_1st_hit_Bar2 = (*hHC2)[0]->GetTime();
        E_1st_hit_Bar2 = (*hHC2)[0]->GetEdep();
            totalBar2Hit++;
            totalBar2E += eDep;
        }
    }
    analysisManager->FillNtupleDColumn(3, totalBar2E/keV);


    // Bar1, Energy deposited in 1st hit 
    analysisManager->FillNtupleDColumn(4,E_1st_hit_Bar1/keV);

    // Bar2, Energy deposited in 1st hit
    analysisManager->FillNtupleDColumn(5,E_1st_hit_Bar2/keV);

    // Time 1
     if (T_1st_hit_Bar1 < 1000./ns){ 
     analysisManager->FillNtupleDColumn(6,T_1st_hit_Bar1/ns);
     }
     else{
     analysisManager->FillNtupleDColumn(6,1001./ns);
     }    
    // Reconstructed energy for bar 1: 195 cm from the neutron source
     analysisManager->FillNtupleDColumn(8, 1.986*1.0E+7 / (T_1st_hit_Bar1*T_1st_hit_Bar1));

    // Time 2
     if (T_1st_hit_Bar2 < 1000./ns){
        analysisManager->FillNtupleDColumn(7,T_1st_hit_Bar2/ns);
     }    
     else{
        analysisManager->FillNtupleDColumn(7,1001./ns);
     }    
    // Reconstructed energy for bar 2: 315 cm from the neutron source
     analysisManager->FillNtupleDColumn(9, 5.182*1.0E+7 / (T_1st_hit_Bar2*T_1st_hit_Bar2));

/*
     /// The energy of incident neutron.
     G4RunManager* runManager = G4RunManager::GetRunManager();
    ///PrimaryGeneratorAction* generator = (PrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
     exGPSPrimaryGeneratorAction* generator = (exGPSPrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
     G4ParticleDefinition* particle = generator->GetParticleSource()->GetParticleDefinition();
     G4String Particle = particle->GetParticleName();
     G4double energy = generator->GetParticleSource()->GetParticleEnergy();

     analysisManager->FillNtupleDColumn(10,energy/keV);

*/

    // Collimator / reflector Hits
 
    G4int n_hit_r = hHCLM->entries();///n_hit_r is the number of hits collimator / reflector.  
    analysisManager->FillNtupleIColumn(10, n_hit_r);

    /// reflector Energy and Time
    G4int totalReflectorHit = 0;
    G4double totalReflectorE = 0.;
    G4double E_1st_hit_r = 0.;
    G4double T_1st_hit_r = 0.;

    for (G4int i=0;i<n_hit_r;i++)
    {
        CollimatorHit* hit = (*hHCLM)[i];
        G4double eDep = hit->GetEdep();
        if (eDep>0.)
        {
        T_1st_hit_r = (*hHCLM)[0]->GetTime();
        E_1st_hit_r = (*hHCLM)[0]->GetEdep();
        totalReflectorHit++;
        totalReflectorE += eDep;
        }
    }

    // Time recorded in reflector
     if (T_1st_hit_r < 1000./ns){
        analysisManager->FillNtupleDColumn(11,T_1st_hit_r/ns);
     }    
     else{
        analysisManager->FillNtupleDColumn(11,1001./ns);
     }    
    analysisManager->FillNtupleDColumn(12, E_1st_hit_r/keV);
    analysisManager->FillNtupleDColumn(13, totalReflectorE/keV);

    /*~~~~~~~~~~~~~~~~~~~ water ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /// water hit
    G4int n_hit_wtr = wtrHC->entries();///n_hit_water is the number of hits per incident particle in water.  
    analysisManager->FillNtupleIColumn(14, n_hit_wtr);

    /// water Energy and Time
    G4int totalWtrHit = 0;
    G4double totalWtrE = 0.;
    G4double E_1st_hit_wtr = 0.;
    G4double T_1st_hit_wtr = 0.;

    for (G4int i=0;i<n_hit_wtr;i++)
    {
        waterHit* hit = (*wtrHC)[i];
        G4double eDep = hit->GetEdep();
        if (eDep>0.)
        {
        T_1st_hit_wtr = (*wtrHC)[0]->GetTime();
        E_1st_hit_wtr = (*wtrHC)[0]->GetEdep();
        totalWtrHit++;
        totalWtrE += eDep;
        }
    }

    // Time recorded in water
     if (T_1st_hit_wtr < 1000./ns){
     analysisManager->FillNtupleDColumn(15,T_1st_hit_wtr/ns);
     }    
     else{
     analysisManager->FillNtupleDColumn(15,1001./ns);
     }    

    analysisManager->FillNtupleDColumn(16, E_1st_hit_wtr/keV);
    analysisManager->FillNtupleDColumn(17, totalWtrE/keV);

    /*~~~~~~~~~~~~~~~~~~~ water Tank ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /// water Tank hit
    G4int n_hit_wtrTk = wtrTkHC->entries();///n_hit_wtrTk is the # of hits per incident particle.  
    analysisManager->FillNtupleIColumn(18, n_hit_wtrTk);

    /// water Tank Energy and Time
    G4int totalWtrTkHit = 0;
    G4double totalWtrTkE = 0.;
    G4double E_1st_hit_wtrTk = 0.;
    G4double T_1st_hit_wtrTk = 0.;

    for (G4int i=0;i<n_hit_wtrTk;i++)
    {
        waterTankHit* hit = (*wtrTkHC)[i];
        G4double eDep = hit->GetEdep();
        if (eDep>0.)
        {
        T_1st_hit_wtrTk = (*wtrTkHC)[0]->GetTime();///wtrTkHC was wrongly as wtrHC !!!
        E_1st_hit_wtrTk = (*wtrTkHC)[0]->GetEdep();///wtrTkHC was wrongly as wtrHC !
        totalWtrTkHit++;
        totalWtrTkE += eDep;
        }
    }

     if (T_1st_hit_wtrTk < 1000./ns){
     analysisManager->FillNtupleDColumn(19,T_1st_hit_wtrTk/ns);
     }    
     else{
     analysisManager->FillNtupleDColumn(19,1001./ns);
     }    

    analysisManager->FillNtupleDColumn(20, E_1st_hit_wtrTk/keV);
    analysisManager->FillNtupleDColumn(21, totalWtrTkE/keV);


    /// grab seeds
    if (!genAction)
    genAction = dynamic_cast<const exGPSPrimaryGeneratorAction*>
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    seeds = genAction->GetEventSeeds();
    ///long index = genAction->GetEventSeedIndex();
    long seed0 = *seeds;
    long seed1 = *(seeds+1);
    ///G4cout << "EventAction.cc, EndofEvent, seed0 is: " << seed0 << ", seed1 is:"  << seed1 << G4endl;
    analysisManager->FillNtupleDColumn(22, seed0);
    analysisManager->FillNtupleDColumn(23, seed1);

    analysisManager->AddNtupleRow();

    
    G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( printModulo==0 || event->GetEventID() % printModulo != 0) return;
    
 /*   G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary(0);
    
    G4cout << G4endl
           << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
           << primary->GetG4code()->GetParticleName()
           << " " << primary->GetMomentum() << G4endl;
*/      
    // Step 3: Loop on the two collections and dump on screen hits
     /// Print info on the monitor
    // Bar1 Energy
     ///  G4cout << "Bar1 has " << totalBar1Hit << " hits. " << G4endl;
    ///<< totalBar1E/keV << " (keV)" << ", 1st hit in bar1 : " << E_1st_hit_Bar1/keV << "(keV)" << G4endl;

    // Bar2 Energy
       ///G4cout << "waterTank has " << totalWtrTkHit << " hits." << G4endl; ///Total Edep is "
     ///  G4cout << "water has " << totalWtrHit << " hits." << G4endl; ///Total Edep is "
    ///<< totalBar2E/keV << " (keV)" << ", 1st hit in bar2 : " << E_1st_hit_Bar2/keV << "(keV)" << G4endl;
 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
