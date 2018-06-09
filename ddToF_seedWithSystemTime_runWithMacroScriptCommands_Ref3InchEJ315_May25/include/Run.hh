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
// \file  Run.hh
// \brief Definition of Run class

#ifndef RUN_HH
#define RUN_HH

#include "G4Run.hh"

class G4Event;

//=================================
// Exercise 1 Step 1:
// Create a user-defined Run class
// This class accumulates the energy
// in Hadronic and Em calorimeters for the
// entire run. Energy is taken from hits.
// Add two datamembers that will hold the accumulated energy
// As an optimization, instead of using the hits collection
// names to search them, let's use their ID. Add two integers
// that holds IDs.

class Run : public G4Run {
public:
    Run();
  virtual ~Run() {};

  virtual void RecordEvent(const G4Event*);
  ///virtual void Merge(const G4Run*);
    ///G4double GetEmEnergy() const { return em_ene; }
    ///G4double GetHadEnergy() const { return had_ene; }
    ///G4double GetShowerShape() const { return shower_shape; }
    G4double GetBar1Energy() const { return bar1_ene; }
    G4double GetBar2Energy() const { return bar2_ene; }
    G4double GetCollimatorEnergy() const { return clm_ene; }
    G4double GetwaterEnergy() const { return water_ene; }
    G4double GetwaterTankEnergy() const { return waterTank_ene; }
private:
    ///G4double em_ene;
    G4double bar1_ene; ///accumulated energy in bar1.
    G4double bar2_ene;///accumulated energy in bar2.
    G4double clm_ene;///accumulated energy in collimator.
    G4double water_ene; ///accumulated energy in water.
    G4double waterTank_ene;///accumulated energy in waterTank.
    ///G4double shower_shape;
    ///G4int ECHCID;
    ///G4int HCHCID;
    G4int BAR1HCID; ///ID for bar1 hits collection.
    G4int BAR2HCID; ///ID for bar2 hits collection.
    G4int CLMHCID;
    G4int WATERHCID; ///ID for water hits collection.
    G4int WATERTANKHCID; ///ID for waterTank hits collection.
};
#endif
