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
	rad = 1.;
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

	
	G4String volName = fTrack->GetVolume()->GetName();
	
	G4double x = fTrack->GetPosition().x();
	G4double y = fTrack->GetPosition().y();
	G4double z = fTrack->GetPosition().z();
	G4ThreeVector p = fTrack->GetMomentum();
	
	if(volName.find("virtualBoundaryPhys_sen") != G4String::npos) {
		fX_sen = x;
		fY_sen = y;
		fZ_sen = z;
		fP_sen = p.mag();	
		fTheta_sen = p.theta()/rad;
		fPhi_sen = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_sm") != G4String::npos) {
		fX_sm = x;
		fY_sm = y;
		fZ_sm = z;
		fP_sm = p.mag();
		fTheta_sm = p.theta()/rad;
		fPhi_sm = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_sex") != G4String::npos) {
		fX_sex = x;
		fY_sex = y;
		fZ_sex = z;
		fP_sex = p.mag();
		fTheta_sex = p.theta()/rad;
		fPhi_sex = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_coil") != G4String::npos) {
		fX_coil = x;
		fY_coil = y;
		fZ_coil = z;
		fP_coil = p.mag();
		fTheta_coil = p.theta()/rad;
		fPhi_coil = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_mid") != G4String::npos) {
		fX_mid = x;
		fY_mid = y;
		fZ_mid = z;
		fP_mid = p.mag();
		fTheta_mid = p.theta()/rad;
		fPhi_mid = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_col") != G4String::npos) {
		fX_col = x;
		fY_col = y;
		fZ_col = z;
		fP_col = p.mag();
		fTheta_col = p.theta()/rad;
		fPhi_col = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q1en") != G4String::npos) {
		fX_q1en = x;
		fY_q1en = y;
		fZ_q1en = z;
		fP_q1en = p.mag();
		fTheta_q1en = p.theta()/rad;
		fPhi_q1en = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q1ex") != G4String::npos) {
		fX_q1ex = x;
		fY_q1ex = y;
		fZ_q1ex = z;
		fP_q1ex = p.mag();
		fTheta_q1ex = p.theta()/rad;
		fPhi_q1ex = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q2en") != G4String::npos) {
		fX_q2en = x;
		fY_q2en = y;
		fZ_q2en = z;
		fP_q2en = p.mag();
		fTheta_q2en = p.theta()/rad;
		fPhi_q2en = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q2ex") != G4String::npos) {
		fX_q2ex = x;
		fY_q2ex = y;
		fZ_q2ex = z;
		fP_q2ex = p.mag();
		fTheta_q2ex = p.theta()/rad;
		fPhi_q2ex = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_den") != G4String::npos) {
		fX_den = x;
		fY_den = y;
		fZ_den = z;
		fP_den = p.mag();
		fTheta_den = p.theta()/rad;
		fPhi_den = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_dex") != G4String::npos) {
		fX_dex = x;
		fY_dex = y;
		fZ_dex = z;
		fP_dex = p.mag();
		fTheta_dex = p.theta()/rad;
		fPhi_dex = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q3en") != G4String::npos) {
		fX_q3en = x;
		fY_q3en = y;
		fZ_q3en = z;
		fP_q3en = p.mag();
		fTheta_q3en = p.theta()/rad;
		fPhi_q3en = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q3ex") != G4String::npos) {
		fX_q3ex = x;
		fY_q3ex = y;
		fZ_q3ex = z;
		fP_q3ex = p.mag();
		fTheta_q3ex = p.theta()/rad;
		fPhi_q3ex = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_vdc") != G4String::npos) {
		fX_vdc = x;
		fY_vdc = y;
		fZ_vdc = z;
		fP_vdc = p.mag();
		fTheta_vdc = p.theta()/rad;
		fPhi_vdc = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_qz1") != G4String::npos) {
		fX_qz1 = x;
		fY_qz1 = y;
		fZ_qz1 = z;
		fP_qz1 = p.mag();
		fTheta_qz1 = p.theta()/rad;
		fPhi_qz1 = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_qz2") != G4String::npos) {
		fX_qz2 = x;
		fY_qz2 = y;
		fZ_qz2 = z;
		fP_qz2 = p.mag();
		fTheta_qz2 = p.theta()/rad;
		fPhi_qz2 = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_fp") != G4String::npos) {
		fX_fp = x;
		fY_fp = y;
		fZ_fp = z;
		fP_fp = p.mag();
		fTheta_fp = p.theta()/rad;
		fPhi_fp = p.phi()/rad;
	} 

}


