#include "g4hrsVEventGen.hh"

#include "G4RotationMatrix.hh"

#include "g4hrsBeamTarget.hh"
#include "g4hrsVertex.hh"
#include "g4hrsEvent.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"

g4hrsVEventGen::g4hrsVEventGen() {
    fBeamTarg = g4hrsBeamTarget::GetBeamTarget();
    fRunData  = g4hrsRun::GetRun()->GetData();

    fSampType       = kMainTarget;
    fApplyMultScatt = false;
}

g4hrsVEventGen::~g4hrsVEventGen() {
}

g4hrsEvent *g4hrsVEventGen::GenerateEvent() {
    // Set up beam/target vertex
    g4hrsVertex vert   = fBeamTarg->SampleVertex(fSampType);

    /////////////////////////////////////////////////////////////////////
    // Create and initialize values for event
    g4hrsEvent *thisev = new g4hrsEvent();
    thisev->fVertexPos    = fBeamTarg->fVer;
    if( fApplyMultScatt ) {
        thisev->fBeamMomentum = fBeamTarg->fSampE*(fBeamTarg->fDir.unit());
    } else {
        thisev->fBeamMomentum = fBeamTarg->fSampE*G4ThreeVector(0.0, 0.0, 1.0);
    }
    /////////////////////////////////////////////////////////////////////

    SamplePhysics(&vert, thisev);

    PolishEvent(thisev);

    return thisev;
}
	
void g4hrsVEventGen::PolishEvent(g4hrsEvent *ev) {
    /*!
       Here it's our job to:
          Make sure the event is sane
          Apply multiple scattering effects to the final
        products if applicable
      Calculate rates from our given luminosity
      Calculate measured asymmetry from polarization
      Calculate vertex offsets
     */

    if( !ev->EventIsSane() ) {
        G4cerr << __FILE__ << " line " << __LINE__ << ":  Event check failed for generator " << fName << ".  Aborting" << G4endl;
        ev->Print();
        exit(1);
    }

    G4ThreeVector rotax      = (fBeamTarg->fDir.cross(G4ThreeVector(0.0, 0.0, 1.0))).unit();
    G4RotationMatrix msrot;
    msrot.rotate(fBeamTarg->fDir.theta(), rotax);

    std::vector<G4ThreeVector>::iterator iter;

    if( fApplyMultScatt ) {
        for( iter = ev->fPartRealMom.begin(); iter != ev->fPartRealMom.end(); iter++ ) {
            //  rotate direction vectors based on multiple scattering
            (*iter) *= msrot;
        }

        // Rotate position offsets due to multiple scattering
        for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
            //  rotate direction vectors based on multiple scattering
            (*iter) *= msrot;
        }
    }

    // Add base vertex
    for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
        (*iter) += ev->fVertexPos;
    }
    

    if ( ev->fRate == 0 ){// If the rate is set to 0 then calculate it using the cross section
    	ev->fRate  = ev->fEffXs*fBeamTarg->GetEffLumin()/((G4double) fRunData->GetNthrown());
    }
    else { // For LUND - calculate rate and cross section	
    	ev->fEffXs = ev->fRate*((G4double) fRunData->GetNthrown())/(fBeamTarg->GetEffLumin());
    	ev->fRate = ev->fRate/((G4double) fRunData->GetNthrown());
    }

    ev->fmAsym = ev->fAsym*fBeamTarg->fBeamPol;

    return;
}

