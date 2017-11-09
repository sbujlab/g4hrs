// ********************************************************************
//
// $Id: g4hrsEMField.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//   User Field class Setup implementation.
//
//  
#include "g4hrsEMField.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:
g4hrsEMField::g4hrsEMField(): mBField_Septum(0)
{
	fSeptumMomentum = 1.063;
	fSeptumMapFile = "PREX_septumFieldMap.dat";
//	fSeptumMapFile = "PREX_juliette_May2017ERR_1300Am2.dat";	

        // FIXME - this should not be hardcoded
	if(fSeptumMomentum != 0 && fSeptumMapFile != "") {

//		mBField_Septum = new BField_Septum(pHRSMomentum,"PREX_septumFieldMap.dat");
		mBField_Septum = new BField_Septum(fSeptumMomentum,fSeptumMapFile);
	}

	bUseUniformBField=false;
	BField3V.set(0,0,0);
}

//////////////////////////////////////////////////////////////////////////
//
//  Deconstructors:
g4hrsEMField::~g4hrsEMField()
{
	if(mBField_Septum){ delete mBField_Septum; mBField_Septum = NULL;}
}


////////////////////////////////////////////////////////////////////////////
//input Point[4] (x,y,z,t) 
//
inline void g4hrsEMField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{  
	//////////////////////////////////////////////////////////
	//get BField
	if(this->bUseUniformBField) 
	{
		Bfield[0]=BField3V.x();
		Bfield[1]=BField3V.y();
		Bfield[2]=BField3V.z();
	}
	else
	{
		double pB[3],pPos[3]={Point[0]/cm,Point[1]/cm,Point[2]/cm};  //turn into cm
		for(int i=0;i<3;i++) Bfield[i]=0.0;  //reset


		//septum field read from map
		if(mBField_Septum)
		{
			for(int i=0;i<3;i++) pB[i]=0.0;  //reset
			if (! mBField_Septum->IsUniformField() )  mBField_Septum->GetBField(pPos,pB); 
			else  mBField_Septum->GetUniformField(pB); 
			for(int i=0;i<3;i++) Bfield[i] =pB[i]*tesla;
		}

	}

}
