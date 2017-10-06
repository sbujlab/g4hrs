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
//#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//#include "hamcPREXTrans.hh"
#include "prex_forward.hh"

g4hrsSteppingAction::g4hrsSteppingAction()
:drawFlag(false)
{

///  new g4hrsSteppingActionMessenger(this);
	rad = 1.;
    fEnableKryptonite = true;
//	fTransportFunction = new hamcPREXTrans();
	//FIXME set in macro
	fSeptumAngle = 5.*deg;
	fHRSAngle = 12.5*deg;
	fHRSMom = 1.063*GeV;

	int n = 5;
	np = &n;
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



	////////////////////////////////////////////////
	// Virtual boundaries and transport functions //
	////////////////////////////////////////////////
	
	if(fTrack->GetParentID()==0) {


		G4double x = fTrack->GetPosition().x();
		G4double y = fTrack->GetPosition().y();
		G4double z = fTrack->GetPosition().z();
		G4ThreeVector p = fTrack->GetMomentum();
		
		// Navigator of parallel world to get volume name (to identify virtual boundaries) and hall to transport coordinate transformation
		G4Navigator* fParallelNavigator = G4TransportationManager::GetTransportationManager()->GetNavigator("g4hrsparallel");

	
		G4String volName = fParallelNavigator->LocateGlobalPointAndSetup(G4ThreeVector(x,y,z))->GetName();
		G4AffineTransform trans = fParallelNavigator->GetGlobalToLocalTransform();
		G4ThreeVector localPos = trans.TransformPoint(G4ThreeVector(x,y,z));	

		G4StepPoint* prePoint = aStep->GetPreStepPoint(); 

		G4RotationMatrix rotate_sen, rotate_sm, rotate_sex;
		rotate_sen.rotateY(-fSeptumAngle);
		rotate_sen.rotateZ(90.*deg);
		rotate_sm.rotateY((fSeptumAngle+fHRSAngle)/2.);
		rotate_sm.rotateZ(-90.*deg);
		rotate_sex.rotateY(fHRSAngle);
		rotate_sex.rotateZ(-90.*deg);
			
		// Axis transformation from HCS to TCS in each septum region
		G4AffineTransform transportAxis_sen = G4AffineTransform(rotate_sen);
		G4AffineTransform transportAxis_sm = G4AffineTransform(rotate_sm);
		G4AffineTransform transportAxis_sex = G4AffineTransform(rotate_sex);
		
		// Since we will be transforming POINTS (not AXES), we want the INVERSE of the axis transformation
		G4AffineTransform transport_sen = transportAxis_sen.Inverse();
		G4AffineTransform transport_sm = transportAxis_sm.Inverse();
		G4AffineTransform transport_sex = transportAxis_sex.Inverse();
		

		if(fTrack->GetCurrentStepNumber()<=1) {
				
			fX0 = prePoint->GetPosition().x();
			fY0 = prePoint->GetPosition().y();
			fZ0 = prePoint->GetPosition().z();
			fP0 = prePoint->GetMomentum().mag();
			fTh0 = prePoint->GetMomentum().getTheta();
			fPh0 = prePoint->GetMomentum().getPhi();	
			
	/*
			fX0 = 2.*mm;
			fY0 = 2.*mm;
			fZ0 = 0.;
			fP0 = 1.063*GeV;
			fTh0 = 5.*deg;
			fPh0 = 0.*deg;
	*/

			// Define position, angle vectors to be transformed
			G4ThreeVector position = G4ThreeVector(fX0, fY0, fZ0);
			G4ThreeVector momentum = prePoint->GetMomentum();
			//Transform
			
			G4ThreeVector position_tr = transport_sen.TransformPoint(position);
			G4ThreeVector momentum_tr = transport_sen.TransformPoint(momentum);
		
			// Get transport coordinates
			fX0_tr = position_tr.x();
			fY0_tr = position_tr.y();
			fZ0_tr = position_tr.z();
		
			// Convert from polar & azimuthal angles to dx/dz and dy/dz
			double dx = sin(momentum_tr.theta())*cos(momentum_tr.phi());
			double dy = sin(momentum_tr.theta())*sin(momentum_tr.phi());
			double dz = cos(momentum_tr.theta());
			fTh0_tr = dx/dz;
			fPh0_tr = dy/dz;

//			G4cout << "HCS\n" << fX0 << "\t" << fY0 << "\t" << fZ0 << "\t" << fTh0 << "\t" << fPh0 << "\n";
//			G4cout << "TCS\n" << fX0_tr << "\t" << fY0_tr << "\t" << fZ0_tr << "\t" << fTh0_tr << "\t" << fPh0_tr << "\n";
			
			r0[0] = (float)fX0;
			r0[1] = (float)fTh0;
			r0[2] = (float)fY0;
			r0[3] = (float)fPh0;
			r0[4] = (float)(fP0/fHRSMom);	

//			fTransportFunction->Acceptance_C(r0, x_tf, th_tf, y_tf, ph_tf, 0);
//			float xtest = x_sp_sen_(r0,np);	
//			G4cout << xtest << "\n\n";
			
			for (int i = 0; i<5; i++) {
				G4cout << r0[i] << "\n";
			}


		}

		if(volName.find("virtualBoundaryPhys_sen") != G4String::npos) {

			// HCS
			fX_sen = x/1000.;
			fY_sen = y/1000.;
			fZ_sen = z/1000.;
			fP_sen = p.mag();	
			fTheta_sen = p.theta()/rad;
			fPhi_sen = p.phi()/rad;
			// TCS
			G4ThreeVector pos_tr = transport_sen.TransformPoint(G4ThreeVector(fX_sen, fY_sen, fZ_sen));
			G4ThreeVector mom_tr = transport_sen.TransformPoint(p);
			fX_sen_tr = pos_tr.x(); 			
			fY_sen_tr = pos_tr.y(); 			
			fZ_sen_tr = pos_tr.z(); 			
			fP_sen_tr = mom_tr.mag();
			fTheta_sen_tr = (sin(mom_tr.theta())*cos(mom_tr.phi()))/cos(mom_tr.theta());	
			fPhi_sen_tr = (sin(mom_tr.theta())*sin(mom_tr.phi()))/cos(mom_tr.theta());	
			// Transport function
			fX_sen_tf = x_sp_sen_(r0, np);
			fY_sen_tf = y_sp_sen_(r0, np);
			fTheta_sen_tf = t_sp_sen_(r0, np);
			fPhi_sen_tf = p_sp_sen_(r0, np);

		} else if(volName.find("virtualBoundaryPhys_sm") != G4String::npos) {
			
			fX_sm = x/1000.;
			fY_sm = y/1000.;
			fZ_sm = z/1000.;
			fP_sm = p.mag();
			fTheta_sm = p.theta()/rad;
			fPhi_sm = p.phi()/rad;
			
			G4ThreeVector pos_tr = transport_sm.TransformPoint(G4ThreeVector(fX_sm, fY_sm, fZ_sm));
			G4ThreeVector mom_tr = transport_sm.TransformPoint(p);
			fX_sm_tr = pos_tr.x(); 			
			fY_sm_tr = pos_tr.y(); 			
			fZ_sm_tr = pos_tr.z(); 			
			fP_sm_tr = mom_tr.mag();
			fTheta_sm_tr = (sin(mom_tr.theta())*cos(mom_tr.phi()))/cos(mom_tr.theta());	
			fPhi_sm_tr = (sin(mom_tr.theta())*sin(mom_tr.phi()))/cos(mom_tr.theta());	

			fX_sm_tf = x_sp_sm_(r0, np);
			fY_sm_tf = y_sp_sm_(r0, np);
			fTheta_sm_tf = t_sp_sm_(r0, np);
			fPhi_sm_tf = p_sp_sm_(r0, np);

		} else if(volName.find("virtualBoundaryPhys_sex") != G4String::npos) {

			fX_sex = x/1000.;
			fY_sex = y/1000.;
			fZ_sex = z/1000.;
			fP_sex = p.mag();
			fTheta_sex = p.theta()/rad;
			fPhi_sex = p.phi()/rad;
			
			G4ThreeVector pos_tr = transport_sex.TransformPoint(G4ThreeVector(fX_sex, fY_sex, fZ_sex));
			G4ThreeVector mom_tr = transport_sex.TransformPoint(p);
			fX_sex_tr = pos_tr.x(); 			
			fY_sex_tr = pos_tr.y(); 			
			fZ_sex_tr = pos_tr.z(); 			
			fP_sex_tr = mom_tr.mag();
			fTheta_sex_tr = (sin(mom_tr.theta())*cos(mom_tr.phi()))/cos(mom_tr.theta());	
			fPhi_sex_tr = (sin(mom_tr.theta())*sin(mom_tr.phi()))/cos(mom_tr.theta());	

			fX_sex_tf = x_sp_sex_(r0, np);
			fY_sex_tf = y_sp_sex_(r0, np);
			fTheta_sex_tf = t_sp_sex_(r0, np);
			fPhi_sex_tf = p_sp_sex_(r0, np);

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
 
	}  // end if ParentID == 0
 

}


