#ifndef __REMOLLVERTEX_HH
#define __REMOLLVERTEX_HH

#include "G4ThreeVector.hh"

/*!
  Vertex information that only the
  user defined generators will see
*/

class G4Material;

class g4hrsVertex {
    public:
	 g4hrsVertex();
	~g4hrsVertex();

	G4double    GetBeamE(){ return fBeamE; }
	G4double    GetRadLen(){ return fRadLen; }
	G4Material *GetMaterial(){ return fMaterial; }

    public:
	G4double fBeamE;
	G4double fRadLen;
	G4Material *fMaterial;
};

#endif//__REMOLLVERTEX_HH
