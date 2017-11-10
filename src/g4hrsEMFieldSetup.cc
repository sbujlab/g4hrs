// ********************************************************************
//
// $Id: g4hrsEMField.hh,v 1.0, 2010/12/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//
//   User Field class Setup implementation.
//
#include "g4hrsEMFieldSetup.hh"
#include "g4hrsEMField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "BField_Dipole.hh"
#include "BField_Quad.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ios.hh"
//#include "QuadFringe.hh"

#include <iostream>
using namespace std;


//////////////////////////////////////////////////////////////////////////
g4hrsEMFieldSetup* g4hrsEMFieldSetup::fg4hrsEMFieldSetup=0;
g4hrsEMFieldSetup* g4hrsEMFieldSetup::Getg4hrsEMFieldSetup()
{ 
	if(!fg4hrsEMFieldSetup)  
	{
		G4cout<<"g4hrsEMFieldSetup is not initialized yet...exit...\n";
		exit(-99);
	}
	return fg4hrsEMFieldSetup; 
}

//////////////////////////////////////////////////////////////////////////
//
g4hrsEMFieldSetup::g4hrsEMFieldSetup()
: fChordFinder(0), fStepper(0), fIntgrDriver(0)
{

	G4ThreeVector pivotOffset(0.0, 0.0,  1053.79*mm );

	fHRSMomentum = 1.063*GeV;
	fSnakeModel = 49;
	fHRSAngle = 12.5*deg;
	fSeptumAngle = 5.0*deg;

	G4cout << "HRS angles: " << fHRSAngle <<  G4endl;

	fg4hrsEMFieldSetup=this;
	
	fEMfield = new g4hrsEMField();
	
	fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	UpdateField();

	G4double snakemagnumber = (-4.77577734*0.83762) / fHRSMomentum;

	G4int    quads_on = 1;
	G4int sos = 1;

	KAPPA1 =  0.;
	KAPPA2 =  0.;
	KAPPA3 =  0.;
	dipoleField = 0.;

	double a_sup = 0.1492;                    //Bore radius, in meters
	double l_sup = 0.9413;                    //Length of quad, in meters
	double a_sos = 0.12827;
	double l_sos = 0.70;

	//BELOW ARE SNAKE VALUES, NOT NIM VALUES                                                                                  
	double pTarget =             0.0  * cm;
	double pQ1en   = pTarget + 160.0  * cm;
	double pQ1ex   = pQ1en   +  94.13 * cm;
	double pQ2en   = pQ1ex   + 115.58 * cm;
	double pQ2ex   = pQ2en   + 182.66 * cm;
	//ABOVE ARE SNAKE VALUES, NOT NIM VALUES  

	//Local field  FZB2, Q1
	//double pHallCenter2Q1Face=1.69*m;//NIM
	double pHallCenter2Q1Face=pQ1en;//SNAKE
	//double pQ1Length= 80*cm;//NIM
	double pQ1Length= sos ? 70.    * cm : 94.13*cm;//SNAKE
	double pQ1Radius= sos ? 12.827 * cm : 0.1492 * m;//SNAKE
	double q1shift = sos ? 0.0 * m : 0.0 * m;
	double pQ1Pos_Z=(pHallCenter2Q1Face+94.13*cm/2.0+q1shift);//NIM
	//double pHallCenter2Q2Face=3.74*m;//NIM
	double pHallCenter2Q2Face=pQ2en;//SNAKE
	//double pQ2Length=180*cm;//NIM
	double pQ2Length=182.66*cm;//SNAKE
	double pQ2Radius=  0.3 * m;//SNAKE
	double pQ2Pos_Z=(pHallCenter2Q2Face+pQ2Length/2.0);
	//double pQ3Length=180*cm;//NIM
	double pQ3Length=182.68*cm;//SNAKE
	double pQ3Radius=  pQ2Radius;//SNAKE
	

	/////////////////////////
	// QUADRUPOLE SETTINGS //
	/////////////////////////
	
	// Default settings
	// PREX, 5 deg
	KAPPA1 = sos ? 0.260387 * tesla / 1.063 * fHRSMomentum/GeV: -0.8476  * tesla / snakemagnumber; //test3
	KAPPA2 =  0.93528 * tesla / snakemagnumber;//so don't worry!
	KAPPA3 =  1.15762 * tesla / snakemagnumber;

	// Alternate settings
/*
	// CREX, 4 deg
	KAPPA1 = -0.7146 * tesla / snakemagnumber * l_sup / l_sos * a_sos / a_sup;
	KAPPA1 *= 0.85;//this is acceptable as a first edit of Seamus
	//KAPPA1 *= 0.8;//Dustin wants me to tighten the transverse part
	KAPPA2 =  0.8680 * tesla / snakemagnumber;
	KAPPA3 =  1.1748 * tesla / snakemagnumber;
	// 6 deg
	KAPPA1 = 1.16070 * tesla * pQ1Radius / 1000.;
	KAPPA2 =-1.37117 * tesla * pQ2Radius / 1000.;
	KAPPA3 =-1.72440 * tesla * pQ3Radius / 1000.;
	// 4.5 deg
	KAPPA1 = 1.02555 * tesla * pQ1Radius / 1000. / 2.2 * fHRSMomentum/GeV;
	KAPPA2 =-1.33393 * tesla * pQ2Radius / 1000. / 2.2 * fHRSMomentum/GeV;;
	KAPPA3 =-1.70387 * tesla * pQ3Radius / 1000. / 2.2 * fHRSMomentum/GeV;;
	// 5.5 deg
	KAPPA1 = 1.10973 * tesla * pQ1Radius / 1000. / 2.2 * fHRSMomentum/GeV;
	KAPPA2 =-1.35688 * tesla * pQ2Radius / 1000. / 2.2 * fHRSMomentum/GeV;
	KAPPA3 =-1.71687 * tesla * pQ3Radius / 1000. / 2.2 * fHRSMomentum/GeV;
	// 6.5 deg
	KAPPA1 = 1.20291 * tesla * pQ1Radius / 1000. / 2.2 * fHRSMomentum/GeV;
	KAPPA2 =-1.38139 * tesla * pQ2Radius / 1000. / 2.2 * fHRSMomentum/GeV;
	KAPPA3 =-1.73056 * tesla * pQ3Radius / 1000. / 2.2 * fHRSMomentum/GeV;
	// PREX tune
	KAPPA1 = -0.8476  * tesla / snakemagnumber / .1492 / m; 
	KAPPA2 =  0.93528 * tesla / snakemagnumber / .300  / m; 
	KAPPA3 =  1.15762 * tesla / snakemagnumber / .300  / m;
	// STD tune
	KAPPA1 = -0.2445  * tesla / snakemagnumber / .1492 / m;
	KAPPA2 =  0.1939  * tesla / snakemagnumber / .300  / m;
	KAPPA3 =  0.1794  * tesla / snakemagnumber / .300  / m;
*/

	/////////////////////
	// DIPOLE SETTINGS //
	/////////////////////
	
	//G4double dipoleField = -0.4192 * tesla; //using dipolert2 from SNAKE, thin vb, this one is reconciled with SNAKE
	//G4double dipoleField = -0.4205 * tesla; //matches nicely with DATA (not snake comparison) 
	dipoleField = -0.4205 * tesla; //matches nicely with DATA (not snake comparison) 

	G4cout << "Magnets: " << KAPPA1 << " " << KAPPA2 << " " << dipoleField << " " << KAPPA3 << G4endl;

	G4RotationMatrix* LROTATED = new G4RotationMatrix;
	G4RotationMatrix* RROTATED = new G4RotationMatrix;
	LROTATED->rotateY( fHRSAngle );
	RROTATED->rotateY( -fHRSAngle );
	double pDipoleRCenterY=8.4  * m;
	double pDipoleRCenterZ=9.961 * m;//SNAKE

	G4ThreeVector     LORIGIND(pDipoleRCenterZ * sin( fHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( fHRSAngle ));
	G4ThreeVector     RORIGIND(pDipoleRCenterZ * sin( -fHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( -fHRSAngle ));

	G4ThreeVector     LORIGINQ1(pQ1Pos_Z * sin(fHRSAngle), 0., pQ1Pos_Z * cos(fHRSAngle));
	G4ThreeVector     RORIGINQ1(pQ1Pos_Z * sin(-fHRSAngle), 0., pQ1Pos_Z * cos(-fHRSAngle));
	G4RotationMatrix* LROTATEQ1 = new G4RotationMatrix;
	G4RotationMatrix* RROTATEQ1 = new G4RotationMatrix;
	LROTATEQ1->rotateY( fHRSAngle);
	RROTATEQ1->rotateY( -fHRSAngle);

	G4ThreeVector     LORIGINQ2(pQ2Pos_Z * sin(fHRSAngle), 0., pQ2Pos_Z * cos(fHRSAngle));
	G4ThreeVector     RORIGINQ2(pQ2Pos_Z * sin(-fHRSAngle), 0., pQ2Pos_Z * cos(-fHRSAngle));
	G4RotationMatrix* LROTATEQ2 = new G4RotationMatrix;
	G4RotationMatrix* RROTATEQ2 = new G4RotationMatrix;

	LROTATEQ2->rotateY( fHRSAngle);
	RROTATEQ2->rotateY( -fHRSAngle);

	//double pFPR       = 8.4  * m;//radius of curvature of dipole
	//double pFPA       = 9.96 * m;//distance from pivot to entrance of dipole//NIM
	//double pFPA       = pDipoleRCenterZ ;//distance from pivot to entrance of dipole//SNAKE
	//double pFPH       = pFPR * tan ( 22.5 * deg );//height
	//double pFPCenterX = ( pFPA + pFPH + ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2) ) *-sin( mRHRSAngle ); //not including half of Q3//NIM
	//double pFPCenterZ = ( pFPA + pFPH + ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2) ) * cos( mRHRSAngle ); //not including half of Q3//NIM
	//double pFPCenterY = ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2); //including half of Q3//NIM

	double pLQ3CenterX = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( fHRSAngle );//SNAKE
	double pLQ3CenterZ = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( fHRSAngle );//SNAKE
	double pLQ3CenterY = ( 3.5853101 * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	
	double pRQ3CenterX = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( -fHRSAngle );//SNAKE
	double pRQ3CenterZ = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( -fHRSAngle );//SNAKE
	double pRQ3CenterY = ( 3.5853101 * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	

	G4ThreeVector     LORIGINQ3(pLQ3CenterX, pLQ3CenterY, pLQ3CenterZ);//try this one
	G4ThreeVector     RORIGINQ3(pRQ3CenterX, pRQ3CenterY, pRQ3CenterZ);//try this one
	G4RotationMatrix* LROTATEQ3 = new G4RotationMatrix;
	G4RotationMatrix* RROTATEQ3 = new G4RotationMatrix;

	LROTATEQ3->rotateX(-45.0 * deg);
	LROTATEQ3->rotateY( fHRSAngle);
	RROTATEQ3->rotateX(-45.0 * deg);
	RROTATEQ3->rotateY( -fHRSAngle);

	////////////////////////////////////
	// CREATE THE HRS MAGNETIC FIELDS //
	////////////////////////////////////

	// LHRS Q1
	fMagFieldFZBL1 = new BField_Quad(KAPPA1, pivotOffset, LORIGINQ1, LROTATEQ1, pQ1Length, pQ1Radius, 1);
	fLocalFieldManagerFZBL1 = new G4FieldManager();
	UpdateFieldFZBL1();
	// RHRS Q1
	fMagFieldFZBR1 = new BField_Quad(KAPPA1, pivotOffset,RORIGINQ1, RROTATEQ1, pQ1Length, pQ1Radius, 1);
	fLocalFieldManagerFZBR1 = new G4FieldManager();
	UpdateFieldFZBR1();
	// LHRS Q2
	fMagFieldFZBL2 = new BField_Quad(KAPPA2,pivotOffset, LORIGINQ2, LROTATEQ2, pQ2Length, pQ2Radius, 2);
	fLocalFieldManagerFZBL2 = new G4FieldManager();
	UpdateFieldFZBL2();
	// RHRS Q2
	fMagFieldFZBR2 = new BField_Quad(KAPPA2,pivotOffset, RORIGINQ2, RROTATEQ2, pQ2Length, pQ2Radius, 2);
	fLocalFieldManagerFZBR2 = new G4FieldManager();
	UpdateFieldFZBR2();
	// LHRS D
	fMagFieldFZBL3 = new BField_Dipole( dipoleField,pivotOffset, LORIGIND, LROTATED );
	fLocalFieldManagerFZBL3 = new G4FieldManager();
	UpdateFieldFZBL3();
	// RHRS D
	fMagFieldFZBR3 = new BField_Dipole( dipoleField,pivotOffset, RORIGIND, RROTATED );
	fLocalFieldManagerFZBR3 = new G4FieldManager();
	UpdateFieldFZBR3();
	// LHRS Q3
	fMagFieldFZBL4 = new BField_Quad(KAPPA3, LORIGINQ3,pivotOffset, LROTATEQ3, pQ3Length, pQ3Radius, 3);
	fLocalFieldManagerFZBL4 = new G4FieldManager();
	UpdateFieldFZBL4();
	// RHRS Q3
	fMagFieldFZBR4 = new BField_Quad(KAPPA3, RORIGINQ3,pivotOffset, RROTATEQ3, pQ3Length, pQ3Radius, 3);
	fLocalFieldManagerFZBR4 = new G4FieldManager();
	UpdateFieldFZBR4();

}

/////////////////////////////////////////////////////////////////////////////////
//

g4hrsEMFieldSetup::~g4hrsEMFieldSetup()
{
	if(fChordFinder) delete fChordFinder;
	if(fStepper)     delete fStepper;
	if(fEquation)    delete fEquation;
	if(fEMfield)     delete fEMfield;
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void g4hrsEMFieldSetup::UpdateField()
{
  fFieldManager->SetDetectorField(fEMfield);
  fFieldManager->CreateChordFinder( fEMfield );
}


/////////////////////////////////////////////////////////////////////////////
void g4hrsEMFieldSetup::UpdateFieldFZBL1()
{
	fLocalFieldManagerFZBL1->SetDetectorField(fMagFieldFZBL1);
	fLocalFieldManagerFZBL1->CreateChordFinder(fMagFieldFZBL1);
}
void g4hrsEMFieldSetup::UpdateFieldFZBR1()
{
	fLocalFieldManagerFZBR1->SetDetectorField(fMagFieldFZBR1);
	fLocalFieldManagerFZBR1->CreateChordFinder(fMagFieldFZBR1);
}

/////////////////////////////////////////////////////////////////////////////
void g4hrsEMFieldSetup::UpdateFieldFZBL2()
{
	fLocalFieldManagerFZBL2->SetDetectorField(fMagFieldFZBL2);
	fLocalFieldManagerFZBL2->CreateChordFinder(fMagFieldFZBL2);
}
void g4hrsEMFieldSetup::UpdateFieldFZBR2()
{
	fLocalFieldManagerFZBR2->SetDetectorField(fMagFieldFZBR2);
	fLocalFieldManagerFZBR2->CreateChordFinder(fMagFieldFZBR2);
}

void g4hrsEMFieldSetup::UpdateFieldFZBL3()
{
  fLocalFieldManagerFZBL3->SetDetectorField(fMagFieldFZBL3);
  fLocalFieldManagerFZBL3->CreateChordFinder(fMagFieldFZBL3);
}
void g4hrsEMFieldSetup::UpdateFieldFZBR3()
{
  fLocalFieldManagerFZBR3->SetDetectorField(fMagFieldFZBR3);
  fLocalFieldManagerFZBR3->CreateChordFinder(fMagFieldFZBR3);
}


/////////////////////////////////////////////////////////////////////////////
void g4hrsEMFieldSetup::UpdateFieldFZBL4()
{
	fLocalFieldManagerFZBL4->SetDetectorField(fMagFieldFZBL4);
	fLocalFieldManagerFZBL4->CreateChordFinder(fMagFieldFZBL4);
}
void g4hrsEMFieldSetup::UpdateFieldFZBR4()
{
	fLocalFieldManagerFZBR4->SetDetectorField(fMagFieldFZBR4);
	fLocalFieldManagerFZBR4->CreateChordFinder(fMagFieldFZBR4);
}


///////////////////////////////////////////////////////////////////////////////
