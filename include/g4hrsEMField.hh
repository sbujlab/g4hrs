// ********************************************************************
//
// $Id: g4hrsEMField.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//

#ifndef g4hrsEMField_H
#define g4hrsEMField_H

#include "G4ThreeVector.hh"
#include "G4ElectroMagneticField.hh"
#include "fields/BField_Septum.hh" //Septum Field class

class g4hrsEMField : public G4ElectroMagneticField
{
public:

	g4hrsEMField() ;                
	~g4hrsEMField() ;  

	inline void GetFieldValue(const G4double Point[4], G4double *Bfield ) const;
	//  Point[4] x,y,z,time
	//  Return as Bfield[0], [1], [2] the magnetic field x, y & z components
	//   and   as Bfield[3], [4], [5] the electric field x, y & z components

	G4bool DoesFieldChangeEnergy() const { return false; }

	inline void SetBField3V(G4ThreeVector v) { BField3V = v; bUseUniformBField=true;}
	inline G4ThreeVector GetBField3V() const { return BField3V; }

private:
	g4hrsEMFieldMessenger* messenger;

	bool bUseUniformBField;

	BField_Septum*  mBField_Septum;

	G4ThreeVector BField3V;

};

#endif
