#include "g4hrsGenericDetectorSum.hh"

G4Allocator<g4hrsGenericDetectorSum> g4hrsGenericDetectorSumAllocator;

g4hrsGenericDetectorSum::g4hrsGenericDetectorSum(int detid, int copyid){
    fDetID  = detid;
    fCopyID = copyid;
    fEdep   = 0.0;
}

g4hrsGenericDetectorSum::~g4hrsGenericDetectorSum(){
}

g4hrsGenericDetectorSum::g4hrsGenericDetectorSum(const g4hrsGenericDetectorSum &right) : G4VHit(){
    // copy constructor
    fDetID  = right.fDetID;
    fCopyID = right.fCopyID;
    fEdep   = right.fEdep;
}

const g4hrsGenericDetectorSum& g4hrsGenericDetectorSum::operator =(const g4hrsGenericDetectorSum &right){
    (*this) = right;
    return *this;
}

G4int g4hrsGenericDetectorSum::operator==(const g4hrsGenericDetectorSum &right ) const {
    return (this==&right) ? 1 : 0;
}
