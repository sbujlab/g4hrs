#include "g4hrsVertex.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

g4hrsVertex::g4hrsVertex(){
    // Some default material
    fMaterial = NULL;
    fBeamE = 0.0*GeV;
    fRadLen = 0.0;
}

g4hrsVertex::~g4hrsVertex(){
}
