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

#include "Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "HodoscopeHit.hh"
#include "HodoscopeSD.hh"
#include "CollimatorHit.hh"
#include "CollimatorSD.hh"
#include "waterHit.hh"
#include "waterSD.hh"
#include "waterTankHit.hh"
#include "waterTankSD.hh"
///#include "HadCalorimeterHit.hh"
///#include "EmCalorimeterHit.hh"


Run::Run() :
    G4Run(),
    bar1_ene(0),
    bar2_ene(0),
    clm_ene(0),
    water_ene(0),
    waterTank_ene(0),
    BAR1HCID(-1), ///Bar 1 Hit counter ID
    BAR2HCID(-1),
    CLMHCID(-1),
    WATERHCID(-1),
    WATERTANKHCID(-1)
{ }


//=================================
// Exercise 1 Step 2:
// Get the hits collection from the event
// For the EM and HAD calorimeter sensitive detectors
// Accumulate on them
// Do not forget to forward the call to the base class RecordEvent
// method

void Run::RecordEvent(const G4Event* evt)
{
    //Forward call to base class
    G4Run::RecordEvent(evt);
  ///G4int evtId = evt->GetEventID();
  ///G4cout << G4endl << "event ID in Run.cc is: " << evtId << G4endl;

    if ( BAR1HCID == -1 || BAR2HCID == -1 || CLMHCID == -1 || WATERHCID == -1 || WATERTANKHCID == -1 ) {
      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
      BAR1HCID = sdManager->GetCollectionID("hodoscope1/hodoscopeColl");
      BAR2HCID = sdManager->GetCollectionID("hodoscope2/hodoscopeColl");
      CLMHCID = sdManager->GetCollectionID("collimator/collimatorColl");
      WATERHCID = sdManager->GetCollectionID("water/waterColl");
      WATERTANKHCID = sdManager->GetCollectionID("waterTank/waterTankColl");
    }
    G4HCofThisEvent* hce = evt->GetHCofThisEvent();
    if (!hce) {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found.\n";
        G4Exception("Run::RecordEvent()",
                    "Code001", JustWarning, msg);
        return;
  
    }


    ///const EmCalorimeterHitsCollection* emHC =
    ///    static_cast<const EmCalorimeterHitsCollection*>(hce->GetHC(BAR1HCID));
    const HodoscopeHitsCollection* Bar1HC =
        static_cast<const HodoscopeHitsCollection*>(hce->GetHC(BAR1HCID));
    const HodoscopeHitsCollection* Bar2HC =
        static_cast<const HodoscopeHitsCollection*>(hce->GetHC(BAR2HCID));
    const CollimatorHitsCollection* CHC =
        static_cast<const CollimatorHitsCollection*>(hce->GetHC(CLMHCID));
    const waterHitsCollection* waterHC =
        static_cast<const waterHitsCollection*>(hce->GetHC(WATERHCID));
    const waterTankHitsCollection* waterTankHC =
        static_cast<const waterTankHitsCollection*>(hce->GetHC(WATERTANKHCID));

    if ( !Bar1HC || !Bar2HC || !CHC || !waterHC || !waterTankHC)
    {
        G4ExceptionDescription msg;
        msg << "Some of hits collections of this event not found.\n";
        G4Exception("Run::RecordEvent()",
                    "Code001", JustWarning, msg);
        return;

    }
    ///G4double em = 0;
    ///G4double had = 0;
    G4double e_bar1 = 0;
    G4double e_bar2 = 0;
    G4double e_clm = 0;
    G4double e_water = 0;
    G4double e_waterTank = 0;

    for (size_t i=0;i<Bar1HC->GetSize();i++)
    {
        ///EmCalorimeterHit* hit = (*emHC)[i];
        HodoscopeHit* hit = (*Bar1HC)[i];
        e_bar1 += hit->GetEdep();
    }
    ///G4cout << "Bar1 has " << Bar1HC->GetSize() << " events. " << G4endl;

    for (size_t i=0;i<Bar2HC->GetSize();i++)
    {
        ///HadCalorimeterHit* hit = (*hadHC)[i];
        HodoscopeHit* hit = (*Bar2HC)[i];
        e_bar2 += hit->GetEdep();
    }

    ///G4cout << "Bar2 has " << Bar2HC->GetSize() << " events. " << G4endl;

    for (size_t i=0;i<CHC->GetSize();i++)
    {
        ///HadCalorimeterHit* hit = (*hadHC)[i];
        CollimatorHit* hit = (*CHC)[i];
        e_clm += hit->GetEdep();
    }

    ///G4cout << "Collimator has " << CHC->GetSize() << " events. " << G4endl;

    for (size_t i=0;i<waterHC->GetSize();i++)
    {
        waterHit* hit = (*waterHC)[i];
        e_water += hit->GetEdep();
    }

    ///G4cout << "water has " << waterHC->GetSize() << " events. " << G4endl;

    for (size_t i=0;i<waterTankHC->GetSize();i++)
    {
        waterTankHit* hit = (*waterTankHC)[i];
        e_waterTank += hit->GetEdep();
    }

    ///G4cout << "water Tank has " << waterTankHC->GetSize() << " events. " << G4endl;

    clm_ene += e_clm;
    bar2_ene += e_bar2;
    bar1_ene += e_bar1; 
    water_ene += e_water;
    waterTank_ene += e_waterTank; 

    ///if ( had+em > 0 )
        ///shower_shape += ( em/(had+em) );
}

//=================================
// Exercise 1 Step 3:
// Implement a Merge method.
// The meaning of this will become clear
// When Multi-threading is discussed.
// This method of a MASTER instance is called
// By all worker threads passing WORKER run as
// parameter. Get the energy from WORKER run and
// Accumulate in the MASTER
// Do not forget to forward the call to the base class
// method to accumulate the number of processed events!

/*
void Run::Merge(const G4Run* aRun)
{
    const Run* localRun = static_cast<const Run*>(aRun);
    bar1_ene += localRun->GetBar1Energy();
    bar2_ene += localRun->GetBar2Energy();
    ///shower_shape += localRun->GetShowerShape();
    //Forward call to base-class
    G4Run::Merge(localRun);
}
*/


