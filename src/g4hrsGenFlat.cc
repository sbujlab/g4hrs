#include "g4hrsGenFlat.hh"

#include "CLHEP/Random/RandFlat.h"

#include "G4Material.hh"
#include "G4PhysicalConstants.hh"

#include "g4hrsEvent.hh"
#include "g4hrsVertex.hh"
#include "g4hrstypes.hh"

g4hrsGenFlat::g4hrsGenFlat(){
	fTh_min = 0.0*deg;
    	fTh_max = 5.0*deg;
	fPh_min = -90.0*deg;
	fPh_max = 90.0*deg;
	fE_min = 1.063*GeV;
	fE_max = 1.063*GeV;

    fApplyMultScatt = false;
}

g4hrsGenFlat::~g4hrsGenFlat(){
}

void g4hrsGenFlat::SamplePhysics(g4hrsVertex *vert, g4hrsEvent *evt){
    // Generate event flat in phase space

    double beamE = vert->GetBeamE();

    double mp = 0.938*GeV;

    double th = acos(CLHEP::RandFlat::shoot(cos(fTh_max), cos(fTh_min)));
    double ph = CLHEP::RandFlat::shoot(fPh_min, fPh_max);
    double ef = CLHEP::RandFlat::shoot(fE_min, fE_max);
	
    evt->SetEffCrossSection(1);

    if( vert->GetMaterial()->GetNumberOfElements() != 1 ){
	G4cerr << __FILE__ << " line " << __LINE__ << 
	    ": Error!  Some lazy programmer didn't account for complex materials in the moller process!" << G4endl;
	exit(1);
    }


    double Q2 = 2.0*beamE*ef*(1.0-cos(th));
    evt->SetQ2( Q2 );

    G4double APV = Q2*0.8e-4/GeV/GeV; // Empirical APV value, 
                                      // stolen from mollerClass.C in mollersim

    evt->SetAsymmetry(APV);


    evt->SetW2( mp*mp + 2.0*mp*(beamE-ef) - Q2 );

    evt->ProduceNewParticle( G4ThreeVector(0.0, 0.0, 0.0), 
	                     G4ThreeVector(ef*sin(th)*cos(ph), ef*sin(th)*sin(ph), ef*cos(th) ), 
			     "e-" );

    return;
}
