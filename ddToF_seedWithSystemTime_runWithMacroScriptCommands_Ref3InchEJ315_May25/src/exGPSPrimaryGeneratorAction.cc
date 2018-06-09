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
// $Id: exGPSPrimaryGeneratorAction.cc 68020 2013-03-13 13:34:15Z gcosmo $
//
/// \file eventgenerator/exgps/src/exGPSPrimaryGeneratorAction.cc
/// \brief Implementation of the exGPSPrimaryGeneratorAction class
//

#include "exGPSPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exGPSPrimaryGeneratorAction::exGPSPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
   fParticleGun = new G4GeneralParticleSource();
  //  seeds=NULL;
  seeds[0] =-1;
  seeds[1] =-1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exGPSPrimaryGeneratorAction::~exGPSPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exGPSPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
 
  ///time_t systime = time(NULL);
  ///seeds[0] = (long) systime;
  ///seeds[1] = (long) (systime*G4UniformRand());
  seeds[0] = 1526459892; 
  seeds[1] = 468818116;
  CLHEP::HepRandom::setTheSeeds(seeds);

/*
  G4int evtId = event->GetEventID();
  seeds[0] = evtId + 1000000000;///1 G
  ///seeds[1] = (long) (seeds[0]*G4UniformRand());
  seeds[1] = evtId + 2000000000;///2 G
  CLHEP::HepRandom::setTheSeeds(seeds);
*/

  ///print the seed number on screen.
  ///seeds[0] = *CLHEP::HepRandom::getTheSeeds();
  ///seeds[1] = *(CLHEP::HepRandom::getTheSeeds()+1);
  ///G4cout << G4endl << "exGPSPrimaryGeneratorAction.cc, seed[0]= " << seeds[0] << ", seed[1]= " << seeds[1] <<  G4endl;
 
  fParticleGun->GeneratePrimaryVertex(event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
