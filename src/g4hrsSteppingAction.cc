#include "g4hrsSteppingAction.hh"
//#include "g4hrsSteppingActionMessenger.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"

g4hrsSteppingAction::g4hrsSteppingAction()
:drawFlag(false)
{
///  new g4hrsSteppingActionMessenger(this);

    fEnableKryptonite = true;
}

void g4hrsSteppingAction::UserSteppingAction(const G4Step *aStep) {
    G4Track* fTrack = aStep->GetTrack();
    G4Material* material = fTrack->GetMaterial();


    // Don't continue in these materials
    if( (   material->GetName()=="Tungsten" 
//        ||  material->GetName()=="Pb"
	||  material->GetName()=="Copper" )
	    && fEnableKryptonite
	){
      fTrack->SetTrackStatus(fStopAndKill); // kill the current track
      // fTrack->SetTrackStatus(fKillTrackAndSecondaries); // kill the current track and also associated secondaries
    }


}


