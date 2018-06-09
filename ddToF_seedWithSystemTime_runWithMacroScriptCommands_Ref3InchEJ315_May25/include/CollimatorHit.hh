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
// Original : $Id: HodoscopeHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
// MOdified to : Collimator.hh 2014-08-11 17:00  Junhui Liao $
//
/// \file CollimatorHit.hh
/// \brief Definition of the CollimatorHit class

#ifndef CollimatorHit_h
#define CollimatorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Collimator hit
///
/// It records:
/// - the strip ID
/// - the particle time
/// - the strip logical volume, its position and energy

class CollimatorHit : public G4VHit
{
public:
    CollimatorHit(G4int i,G4double t);
    ///CollimatorHit(G4int i,G4double t,G4double e);
    CollimatorHit(G4int i);
    CollimatorHit(const CollimatorHit &right);
    virtual ~CollimatorHit();

    const CollimatorHit& operator=(const CollimatorHit &right);
    int operator==(const CollimatorHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void* aHit);
    
    void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;///Refer to "CollimatorHit.cc"
    virtual std::vector<G4AttValue>* CreateAttValues() const;///Refer to "CollimatorHit.cc"
    void Print();
    
    G4int GetID() const { return fId; } /// fId is strip ID.

    void SetCellID(G4int z) { fCellID = z; }
    G4int GetCellID() const { return fCellID; }

    void SetEdep(G4double de) { fEdep = de; }
    G4double GetEdep() const { return fEdep; }

    ///void SetEkin(G4double ek) { fEkin = ek; }
    ///G4double GetEkin() const { return fEkin; }

    void SetTime(G4double val) { fTime = val; }
    G4double GetTime() const { return fTime; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    ///void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    ///G4ThreeVector GetLocalPos() const { return fLocalPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }
    
private:
    G4int fId;
    G4int fCellID;
    G4double fEdep;
    G4double fTime;
    ///G4ThreeVector fLocalPos;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    const G4LogicalVolume* fPLogV;
};

typedef G4THitsCollection<CollimatorHit> CollimatorHitsCollection;

extern G4ThreadLocal G4Allocator<CollimatorHit>* CollimatorHitAllocator;

inline void* CollimatorHit::operator new(size_t)
{
    if (!CollimatorHitAllocator)
        CollimatorHitAllocator = new G4Allocator<CollimatorHit>;
    return (void*)CollimatorHitAllocator->MallocSingle();
}

inline void CollimatorHit::operator delete(void* aHit)
{
    CollimatorHitAllocator->FreeSingle((CollimatorHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
