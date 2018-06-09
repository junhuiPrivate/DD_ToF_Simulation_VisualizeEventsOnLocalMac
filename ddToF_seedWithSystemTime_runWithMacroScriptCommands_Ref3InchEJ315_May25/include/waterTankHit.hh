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

#ifndef waterTankHit_h
#define waterTankHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// waterTank hit
///
/// It records:
/// - the strip ID
/// - the particle time
/// - the strip logical volume, its position and rotation

class waterTankHit : public G4VHit
{
public:
    waterTankHit(G4int i,G4double t);
    ///waterTankHit(G4int i,G4double t,G4double e);
    waterTankHit(G4int i);
    waterTankHit(const waterTankHit &right);
    virtual ~waterTankHit();

    const waterTankHit& operator=(const waterTankHit &right);
    int operator==(const waterTankHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void* aHit);
    
    void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;///Refer to "HodoscopeHit.cc"
    virtual std::vector<G4AttValue>* CreateAttValues() const;///Refer to "HodoscopeHit.cc"
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

typedef G4THitsCollection<waterTankHit> waterTankHitsCollection;

extern G4ThreadLocal G4Allocator<waterTankHit>* waterTankHitAllocator;

inline void* waterTankHit::operator new(size_t)
{
    if (!waterTankHitAllocator)
        waterTankHitAllocator = new G4Allocator<waterTankHit>;
    return (void*)waterTankHitAllocator->MallocSingle();
}

inline void waterTankHit::operator delete(void* aHit)
{
    waterTankHitAllocator->FreeSingle((waterTankHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
