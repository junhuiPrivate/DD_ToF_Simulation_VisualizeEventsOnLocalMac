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
// Original script, $Id: HodoscopeSD.cc 76474 2013-11-11 10:36:34Z gcosmo $
/// Modified by Junhui Liao, Dec 2017.
/// \file waterSD.cc
/// \brief Implementation of the waterSD class

#include "waterSD.hh"
#include "waterHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

waterSD::waterSD(G4String name)
: G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    G4String HCname = "waterColl";
    collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

waterSD::~waterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void waterSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new waterHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    if (fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    hce->AddHitsCollection(fHCID,fHitsCollection);

        ///fHitsCollection->insert(hit);

   // fill scintillator bar hits with zero energy deposition
   // The followine lines result weird # of events, very small energy and accumulated timing.
   /*  G4int NCellsOneBar = 25;
    for (G4int i=0;i<NCellsOneBar;i++)
    {
        HodoscopeHit* hit = new HodoscopeHit(i);
        fHitsCollection->insert(hit);
    }
   */ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TotalEdepwaterSD = 0.;

G4bool waterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep==0.) return true;

    G4StepPoint* preStepPoint = step->GetPreStepPoint();

    G4TouchableHistory* touchable
      = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4int copyNo = touchable->GetVolume()->GetCopyNo();
    G4double hitTime = preStepPoint->GetGlobalTime();
    ///G4double depositeE = preStepPoint->GetTotalEnergy();
    
    // check if this finger already has a hit
    G4int ix = -1;
    for (G4int i=0;i<fHitsCollection->entries();i++)
    {
        if ((*fHitsCollection)[i]->GetID()==copyNo)
        {
            ix = i;
            break;
        }
    }

    if (ix>=0)
        // if it has, then take the earlier time
    {
        (*fHitsCollection)[ix]->SetEdep(edep); ///set energy to the collection.Very important,othewise energy spectrum is weird.
        if ((*fHitsCollection)[ix]->GetTime()>hitTime)
        { (*fHitsCollection)[ix]->SetTime(hitTime); }
    }
    else
        // if not, create a new hit and set it to the collection
    {
        waterHit* hit = new waterHit(copyNo,hitTime);
        ///HodoscopeHit* hit = new HodoscopeHit(copyNo,hitTime,depositeE);///this method shows the E 
        G4VPhysicalVolume* physical = touchable->GetVolume();             ///deposited by all incident n. why?
        hit->SetLogV(physical->GetLogicalVolume());
        G4AffineTransform transform 
          = touchable->GetHistory()->GetTopTransform();
        transform.Invert();
        hit->SetRot(transform.NetRotation());
        hit->SetPos(transform.NetTranslation());
        hit->SetEdep(edep);///set energy to the collection. Very important,othewise energy spectrum is weird.
        ///G4double TimeOfThisHit = hit->GetTime(); ///It's not necessary to do this
        ///hit->SetTime(TimeOfThisHit);
        fHitsCollection->insert(hit);
    // add energy deposition
        ///hit->AddEdep(edep);
    
    }    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
