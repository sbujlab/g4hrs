#ifndef __REMOLLGLOBALFIELD_HH
#define __REMOLLGLOBALFIELD_HH

/*!
   \class g4hrsGlobalField
   \brief Global field interface
*/

#include "G4MagneticField.hh"
#include "G4UImanager.hh"

class g4hrsMagneticField;

class g4hrsGlobalField : public G4MagneticField {
    public: 
	 g4hrsGlobalField();
	~g4hrsGlobalField();

	void AddNewField( G4String file );
	void SetFieldScale( G4String file, G4double scale  );
	void SetMagnetCurrent( G4String file, G4double scale  );

	void GetFieldValue( const G4double[], G4double *) const;

    private:
	std::vector<g4hrsMagneticField *> fFields;

	g4hrsMagneticField *GetFieldByName( G4String file );

};


#endif//__REMOLLGLOBALFIELD_HH
