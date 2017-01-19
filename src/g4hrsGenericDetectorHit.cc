#include "g4hrsGenericDetectorHit.hh"

G4Allocator<g4hrsGenericDetectorHit> g4hrsGenericDetectorHitAllocator;

g4hrsGenericDetectorHit::g4hrsGenericDetectorHit(G4int det, G4int copy){
    fDetID  = det;
    fCopyID = copy;

    f3X = G4ThreeVector(-1e9, -1e9, -1e9);
    f3P = G4ThreeVector(-1e9, -1e9, -1e9);
    f3V = G4ThreeVector(-1e9, -1e9, -1e9);

    fP  = -1.0;
    fE  = -1.0;
    fM  = -1.0;

    fTrID  = -1;
    fPID   = (G4int) 1e9;
    fmTrID = -1;

    fGen   = 1;

}

g4hrsGenericDetectorHit::~g4hrsGenericDetectorHit(){
}

g4hrsGenericDetectorHit::g4hrsGenericDetectorHit(const g4hrsGenericDetectorHit &right) : G4VHit(){
    // copy constructor

    fDetID  = right.fDetID;
    fCopyID = right.fCopyID;
    f3X     = right.f3X;
    f3P     = right.f3P;
    f3V     = right.f3V;

    fP      = right.fP;
    fE      = right.fE;
    fM      = right.fM;

    fTrID   = right.fTrID;
    fPID    = right.fPID;
    fmTrID  = right.fmTrID;
    fGen    = right.fGen;

}

const g4hrsGenericDetectorHit& g4hrsGenericDetectorHit::operator =(const g4hrsGenericDetectorHit &right){
    (*this) = right;
    return *this;
}

G4int g4hrsGenericDetectorHit::operator==(const g4hrsGenericDetectorHit &right ) const {
    return (this==&right) ? 1 : 0;
}
