#include "g4hrsSteppingAction.hh"
//#include "g4hrsSteppingActionMessenger.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4TouchableHandle.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"

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
	
	G4double x = fTrack->GetPosition().x();
	G4double y = fTrack->GetPosition().y();
	G4double z = fTrack->GetPosition().z();
	G4ThreeVector p = fTrack->GetMomentum();
	
	G4Navigator* fParallelNavigator = G4TransportationManager::GetTransportationManager()->GetNavigator("g4hrsparallel");
	G4String volName = fParallelNavigator->LocateGlobalPointAndSetup(G4ThreeVector(x,y,z))->GetName();


	G4AffineTransform trans = fParallelNavigator->GetGlobalToLocalTransform();
	G4ThreeVector localPos = trans.TransformPoint(G4ThreeVector(x,y,z));	

//	G4cout << volName << "\n";
//	G4cout << "GLOBAL POSITION " << x << " " << y << " " << z << "\n";
//	G4cout << "LOCAL  POSITION " << localPos.x() << " " << localPos.y() << " " << localPos.z() << "\n\n";
/*	
	G4TouchableHandle touch = new G4TouchableHistory();
	fParallelNavigator->LocateGlobalPointAndUpdateTouchableHandle(G4ThreeVector(x,y,z),G4ThreeVector(0.,0.,0.),touch);	
	G4cout << touch->GetVolume()->GetName() << "\n";	
	G4ThreeVector localPos = touch->GetHistory()->GetTopTransform().TransformPoint(G4ThreeVector(x,y,z));	
	G4cout << "GLOBAL POSITION " << x << " " << y << " " << z << "\n";
	G4cout << "LOCAL  POSITION " << localPos.x() << " " << localPos.y() << " " << localPos.z() << "\n\n";
*/

	G4cout << "Begin transform test\n\n";

	G4RotationMatrix transRot = G4RotationMatrix();
	transRot.rotateY(3.14159/2.);
	G4ThreeVector transTrans = G4ThreeVector(0.,0.,5.);
	G4AffineTransform theTrans = G4AffineTransform(transRot,transTrans);

	G4ThreeVector myPoint = G4ThreeVector(5.,0.,5.);
	G4cout << "The original point is " << myPoint << "\n";

	G4ThreeVector newPoint = theTrans.TransformPoint(myPoint);
	G4cout << "The new point is " << newPoint << "\n";



	if(volName.find("virtualBoundaryPhys_sen") != G4String::npos) {
		fX_sen = x/1000.;
		fY_sen = y/1000.;
		fZ_sen = z/1000.;
		fP_sen = p.mag();	
		fTheta_sen = p.theta()/rad;
		fPhi_sen = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_sm") != G4String::npos) {
		fX_sm = x/1000.;
		fY_sm = y/1000.;
		fZ_sm = z/1000.;
		fP_sm = p.mag();
		fTheta_sm = p.theta()/rad;
		fPhi_sm = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_sex") != G4String::npos) {
		fX_sex = x/1000.;
		fY_sex = y/1000.;
		fZ_sex = z/1000.;
		fP_sex = p.mag();
		fTheta_sex = p.theta()/rad;
		fPhi_sex = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_coil") != G4String::npos) {
		fX_coil = x/1000.;
		fY_coil = y/1000.;
		fZ_coil = z/1000.;
		fP_coil = p.mag();
		fTheta_coil = p.theta()/rad;
		fPhi_coil = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_mid") != G4String::npos) {
		fX_mid = x/1000.;
		fY_mid = y/1000.;
		fZ_mid = z/1000.;
		fP_mid = p.mag();
		fTheta_mid = p.theta()/rad;
		fPhi_mid = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_col") != G4String::npos) {
		fX_col = x/1000.;
		fY_col = y/1000.;
		fZ_col = z/1000.;
		fP_col = p.mag();
		fTheta_col = p.theta()/rad;
		fPhi_col = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q1en_LHRS") != G4String::npos) {
		fX_q1en_L = x/1000.;
		fY_q1en_L = y/1000.;
		fZ_q1en_L = z/1000.;
		fP_q1en_L = p.mag();
		fTheta_q1en_L = p.theta()/rad;
		fPhi_q1en_L = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q1ex_LHRS") != G4String::npos) {
		fX_q1ex_L = x/1000.;
		fY_q1ex_L = y/1000.;
		fZ_q1ex_L = z/1000.;
		fP_q1ex_L = p.mag();
		fTheta_q1ex_L = p.theta()/rad;
		fPhi_q1ex_L = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q1en_RHRS") != G4String::npos) {
		fX_q1en_R = x/1000.;
		fY_q1en_R = y/1000.;
		fZ_q1en_R = z/1000.;
		fP_q1en_R = p.mag();
		fTheta_q1en_R = p.theta()/rad;
		fPhi_q1en_L = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q1ex_RHRS") != G4String::npos) {
		fX_q1ex_R = x/1000.;
		fY_q1ex_R = y/1000.;
		fZ_q1ex_R = z/1000.;
		fP_q1ex_R = p.mag();
		fTheta_q1ex_R = p.theta()/rad;
		fPhi_q1ex_R = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q2en_LHRS") != G4String::npos) {
		fX_q2en_L = x/1000.;
		fY_q2en_L = y/1000.;
		fZ_q2en_L = z/1000.;
		fP_q2en_L = p.mag();
		fTheta_q2en_L = p.theta()/rad;
		fPhi_q2en_L = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q2ex_LHRS") != G4String::npos) {
		fX_q2ex_L = x/1000.;
		fY_q2ex_L = y/1000.;
		fZ_q2ex_L = z/1000.;
		fP_q2ex_L = p.mag();
		fTheta_q2ex_L = p.theta()/rad;
		fPhi_q2ex_L = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q2en_RHRS") != G4String::npos) {
		fX_q2en_R = x/1000.;
		fY_q2en_R = y/1000.;
		fZ_q2en_R = z/1000.;
		fP_q2en_R = p.mag();
		fTheta_q2en_R = p.theta()/rad;
		fPhi_q2en_L = p.phi()/rad;
	} else if(volName.find("virtualBoundaryPhys_q2ex_RHRS") != G4String::npos) {
		fX_q2ex_R = x/1000.;
		fY_q2ex_R = y/1000.;
		fZ_q2ex_R = z/1000.;
		fP_q2ex_R = p.mag();
		fTheta_q2ex_R = p.theta()/rad;
		fPhi_q2ex_R = p.phi()/rad;
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
	} else if(volName.find("parallelDetPhys") != G4String::npos) {
//		G4cout << "in DetPhys vb...\n";
		fX_par = x;
		fY_par = y;
		fZ_par = z;
		fP_par = p.mag();
		fTheta_par = p.theta()/rad;
		fPhi_par = p.phi()/rad;
	} 
 

 

}


