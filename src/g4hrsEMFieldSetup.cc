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
#include "g4hrsTune.hh"

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

	fTune = g4hrsTune::GetTune(); 

	fMomentumScale = fTune->GetMomentumScale(); 		// assume PREX default
  	fHRSMomentum = fTune->HRSMomentum;
  	fSnakeModel = 49;
  	fHRSAngle = fTune->HRSAngle;
  	fSeptumAngle = fTune->septumAngle;

	KAPPA1 = (fTune->kappaQuad1)*(fHRSMomentum/fMomentumScale);
	KAPPA2 = (fTune->kappaQuad2)*(fHRSMomentum/fMomentumScale);
	DPBEND = (fTune->bendDipole)*(fHRSMomentum/fMomentumScale);
	KAPPA3 = (fTune->kappaQuad3)*(fHRSMomentum/fMomentumScale);
	quads_on = fTune->quadsOn;
	sos = fTune->sosQuad;
		

//  G4cout << "HRS angles: " << fHRSAngle <<  G4endl;

  //G4cout << "Quad fringe?" << G4endl;
  //QuadFringe* fringe = new QuadFringe();
  //G4cout << "Quad fringe!" << G4endl;
  fg4hrsEMFieldSetup=this;
  //global EM field
  fEMfield = new g4hrsEMField();
//  fEquation = new G4EqMagElectricField(fEMfield);
  fMinStep  = 5.*mm; 
//  fStepperType = 4 ;     // ClassicalRK4 -- the default stepper
  
// 	fFieldManager = new G4FieldManager();
//	UpdateField();

  fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
//  fChordFinder = 0;   //will be set in UpdateField()
  UpdateField();
  
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

  G4cout << "Magnets: " << KAPPA1 << " " << KAPPA2 << " " << DPBEND << " " << KAPPA3 << G4endl;

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
  
  //cout << pQ3CenterX << " " << pFPCenterX << endl;
  //cout << pQ3CenterY << " " << pFPCenterY << endl;
  //cout << pQ3CenterZ << " " << pFPCenterZ << endl;
  G4ThreeVector     LORIGINQ3(pLQ3CenterX, pLQ3CenterY, pLQ3CenterZ);//try this one
  G4ThreeVector     RORIGINQ3(pRQ3CenterX, pRQ3CenterY, pRQ3CenterZ);//try this one
  G4RotationMatrix* LROTATEQ3 = new G4RotationMatrix;
  G4RotationMatrix* RROTATEQ3 = new G4RotationMatrix;
  
  LROTATEQ3->rotateX(-45.0 * deg);
  LROTATEQ3->rotateY( fHRSAngle);
  RROTATEQ3->rotateX(-45.0 * deg);
  RROTATEQ3->rotateY( -fHRSAngle);
 
  fMagFieldFZBL1 = new BField_Quad(KAPPA1, pivotOffset, LORIGINQ1, LROTATEQ1, pQ1Length, pQ1Radius, 1);
  fLocalFieldManagerFZBL1 = new G4FieldManager();
  UpdateFieldFZBL1();
  fMagFieldFZBR1 = new BField_Quad(KAPPA1, pivotOffset,RORIGINQ1, RROTATEQ1, pQ1Length, pQ1Radius, 1);
  fLocalFieldManagerFZBR1 = new G4FieldManager();
  UpdateFieldFZBR1();
  
  //Local field  FZB2, Q2
  fMagFieldFZBL2 = new BField_Quad(KAPPA2,pivotOffset, LORIGINQ2, LROTATEQ2, pQ2Length, pQ2Radius, 2);
  fLocalFieldManagerFZBL2 = new G4FieldManager();
  UpdateFieldFZBL2();
  fMagFieldFZBR2 = new BField_Quad(KAPPA2,pivotOffset, RORIGINQ2, RROTATEQ2, pQ2Length, pQ2Radius, 2);
  fLocalFieldManagerFZBR2 = new G4FieldManager();
  UpdateFieldFZBR2();
  
  fMagFieldFZBL3 = new BField_Dipole(DPBEND,pivotOffset, LORIGIND, LROTATED );
  fLocalFieldManagerFZBL3 = new G4FieldManager();
  UpdateFieldFZBL3();
  fMagFieldFZBR3 = new BField_Dipole(DPBEND,pivotOffset, RORIGIND, RROTATED );
  fLocalFieldManagerFZBR3 = new G4FieldManager();
  UpdateFieldFZBR3();
  
  //Local field  FZB4, Q3
  fMagFieldFZBL4 = new BField_Quad(KAPPA3, LORIGINQ3,pivotOffset, LROTATEQ3, pQ3Length, pQ3Radius, 3);
  fLocalFieldManagerFZBL4 = new G4FieldManager();
  UpdateFieldFZBL4();
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
  //SetStepper();
//  fStepper = new G4ClassicalRK4( fEquation, 8 );
  //G4cout<<"g4hrsEMFieldSetup:: The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fFieldManager->SetDetectorField(fEMfield);
  fFieldManager->CreateChordFinder( fEMfield );
  /*
  if(fChordFinder) delete fChordFinder;
  fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  fFieldManager->SetChordFinder( fChordFinder );
  */
}


/////////////////////////////////////////////////////////////////////////////
void g4hrsEMFieldSetup::UpdateFieldFZBL1()
{
  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBL1 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL1->SetDetectorField(fMagFieldFZBL1);
	fLocalFieldManagerFZBL1->CreateChordFinder(fMagFieldFZBL1);
/*
	if(fChordFinderFZBL1) delete fChordFinderFZBL1;
	fIntgrDriverFZBL1 = new G4MagInt_Driver(fMinStep,fStepperFZBL1,fStepperFZBL1->GetNumberOfVariables());
	fChordFinderFZBL1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL1, fMinStep, fStepperFZBL1);
	fLocalFieldManagerFZBL1->SetChordFinder( fChordFinderFZBL1 );
*/	
}
void g4hrsEMFieldSetup::UpdateFieldFZBR1()
{
  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBR1 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR1->SetDetectorField(fMagFieldFZBR1);
	fLocalFieldManagerFZBR1->CreateChordFinder(fMagFieldFZBR1);
/*
	if(fChordFinderFZBR1) delete fChordFinderFZBR1;
	fIntgrDriverFZBR1 = new G4MagInt_Driver(fMinStep,fStepperFZBR1,fStepperFZBR1->GetNumberOfVariables());
	fChordFinderFZBR1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR1, fMinStep, fStepperFZBR1);
	fLocalFieldManagerFZBR1->SetChordFinder( fChordFinderFZBR1 );
*/	
}

/////////////////////////////////////////////////////////////////////////////
void g4hrsEMFieldSetup::UpdateFieldFZBL2()
{
  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBL2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL2->SetDetectorField(fMagFieldFZBL2);
	fLocalFieldManagerFZBL2->CreateChordFinder(fMagFieldFZBL2);
/*
	if(fChordFinderFZBL2) delete fChordFinderFZBL2;
	fIntgrDriverFZBL2 = new G4MagInt_Driver(fMinStep,fStepperFZBL2,fStepperFZBL2->GetNumberOfVariables());
	fChordFinderFZBL2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL2, fMinStep, fStepperFZBL2);
	fLocalFieldManagerFZBL2->SetChordFinder( fChordFinderFZBL2 );
*/
}
void g4hrsEMFieldSetup::UpdateFieldFZBR2()
{
  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBR2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR2->SetDetectorField(fMagFieldFZBR2);
	fLocalFieldManagerFZBR2->CreateChordFinder(fMagFieldFZBR2);
/*
	if(fChordFinderFZBR2) delete fChordFinderFZBR2;
	fIntgrDriverFZBR2 = new G4MagInt_Driver(fMinStep,fStepperFZBR2,fStepperFZBR2->GetNumberOfVariables());
	fChordFinderFZBR2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR2, fMinStep, fStepperFZBR2);
	fLocalFieldManagerFZBR2->SetChordFinder( fChordFinderFZBR2 );
*/
}

void g4hrsEMFieldSetup::UpdateFieldFZBL3()
{
//  fStepperFZBL3 = new G4ClassicalRK4( fEquationFZBL3, 8 );
  //G4cout<<"g4hrsEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;

  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBL3 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fLocalFieldManagerFZBL3->SetDetectorField(fMagFieldFZBL3);
  fLocalFieldManagerFZBL3->CreateChordFinder(fMagFieldFZBL3);

	fLocalFieldManagerFZBL3->GetChordFinder()->GetIntegrationDriver()->SetHmin(fMinStep);

/*  
  if(fChordFinderFZBL3) delete fChordFinderFZBL3;
  fChordFinderFZBL3 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL3, fMinStep, fStepperFZBL3);
  fLocalFieldManagerFZBL3->SetChordFinder( fChordFinderFZBL3 );

  fLocalFieldManagerFZBL3->SetMinimumEpsilonStep( 1.0e-5 );
  fLocalFieldManagerFZBL3->SetMaximumEpsilonStep( 1.0e-4 );
  fLocalFieldManagerFZBL3->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer
*/
}
void g4hrsEMFieldSetup::UpdateFieldFZBR3()
{
//  fStepperFZBR3 = new G4ClassicalRK4( fEquationFZBR3, 8 );
  //G4cout<<"g4hrsEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;

  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBR3 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fLocalFieldManagerFZBR3->SetDetectorField(fMagFieldFZBR3);
  fLocalFieldManagerFZBR3->CreateChordFinder(fMagFieldFZBR3);
	
	fLocalFieldManagerFZBR3->GetChordFinder()->GetIntegrationDriver()->SetHmin(fMinStep);
 
/* 
  if(fChordFinderFZBR3) delete fChordFinderFZBR3;
  fChordFinderFZBR3 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR3, fMinStep, fStepperFZBR3);
  fLocalFieldManagerFZBR3->SetChordFinder( fChordFinderFZBR3 );

  fLocalFieldManagerFZBR3->SetMinimumEpsilonStep( 1.0e-5 );
  fLocalFieldManagerFZBR3->SetMaximumEpsilonStep( 1.0e-4 );
  fLocalFieldManagerFZBR3->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer
*/
}


/////////////////////////////////////////////////////////////////////////////
void g4hrsEMFieldSetup::UpdateFieldFZBL4()
{
  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBL4 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL4->SetDetectorField(fMagFieldFZBL4);
	fLocalFieldManagerFZBL4->CreateChordFinder(fMagFieldFZBL4);
/*
	if(fChordFinderFZBL4) delete fChordFinderFZBL4;
	fIntgrDriverFZBL4 = new G4MagInt_Driver(fMinStep,fStepperFZBL4,fStepperFZBL4->GetNumberOfVariables());
	fChordFinderFZBL4 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL4, fMinStep, fStepperFZBL4);
	fLocalFieldManagerFZBL4->SetChordFinder( fChordFinderFZBL4 );
*/
}
void g4hrsEMFieldSetup::UpdateFieldFZBR4()
{
  //G4cout<<"g4hrsEMFieldSetup:: The minimal step for FZBR4 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR4->SetDetectorField(fMagFieldFZBR4);
	fLocalFieldManagerFZBR4->CreateChordFinder(fMagFieldFZBR4);
/*
	if(fChordFinderFZBR4) delete fChordFinderFZBR4;
	fIntgrDriverFZBR4 = new G4MagInt_Driver(fMinStep,fStepperFZBR4,fStepperFZBR4->GetNumberOfVariables());
	fChordFinderFZBR4 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR4, fMinStep, fStepperFZBR4);
	fLocalFieldManagerFZBR4->SetChordFinder( fChordFinderFZBR4 );
*/
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//


void g4hrsEMFieldSetup::SetStepper()
{
	G4int nvar = 8;

	if(fStepper) delete fStepper;

	switch ( fStepperType )
	{
	case 0:
		fStepper = new G4ExplicitEuler( fEquation, nvar );
		//G4cout<<"g4hrsEMFieldSetup:: G4ExplicitEuler is calledS"<<G4endl;
		break;
	case 1:
		fStepper = new G4ImplicitEuler( fEquation, nvar );
		//G4cout<<"g4hrsEMFieldSetup:: G4ImplicitEuler is called"<<G4endl;
		break;
	case 2:
		fStepper = new G4SimpleRunge( fEquation, nvar );
		//G4cout<<"g4hrsEMFieldSetup:: G4SimpleRunge is called"<<G4endl;
		break;
	case 3:
		fStepper = new G4SimpleHeum( fEquation, nvar );
		//G4cout<<"g4hrsEMFieldSetup:: G4SimpleHeum is called"<<G4endl;
		break;
	case 4:
		fStepper = new G4ClassicalRK4( fEquation, nvar );
		//G4cout<<"g4hrsEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;
		break;
	case 5:
		fStepper = new G4CashKarpRKF45( fEquation, nvar );
		//G4cout<<"g4hrsEMFieldSetup:: G4CashKarpRKF45 is called"<<G4endl;
		break;

	//The following not working for electric field
	case 6:
		fStepper = 0 ; // new G4RKG3_Stepper( fMagEquation );
		//G4cout<<"g4hrsEMFieldSetup:: G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;
		break;
	case 7:
		fStepper = 0 ; // new G4HelixExplicitEuler( fMagEquation );
		//G4cout<<"g4hrsEMFieldSetup:: G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 8:
		fStepper = 0 ; //  new G4HelixImplicitEuler( fMagEquation );
		//G4cout<<"g4hrsEMFieldSetup:: G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 9:
		fStepper = 0 ; //  new G4HelixSimpleRunge( fMagEquation );
		//G4cout<<"g4hrsEMFieldSetup:: G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
		break;
	default: fStepper = 0;
	}
}


///////////////////////////////////////////////////////////////////////////////
