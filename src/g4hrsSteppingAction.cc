#include "g4hrsSteppingAction.hh"
#include "g4hrsTransportFunction.hh"
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
#include "G4SystemOfUnits.hh"

g4hrsSteppingAction::g4hrsSteppingAction()
:drawFlag(false)
{

///  new g4hrsSteppingActionMessenger(this);
	rad = 1.;
    fEnableKryptonite = true;
	fTransportFunction = new g4hrsTransportFunction();
	//FIXME set in macro
	fSeptumAngle = 5.*deg;
	fHRSAngle = 12.5*deg;
	fHRSMom = 1.063*GeV;

	nelements = 5;

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
		
	G4StepPoint* prePoint = aStep->GetPreStepPoint(); 
	
	if(fTrack->GetParentID()==0) {
	
		G4double x = fTrack->GetPosition().x();
		G4double y = fTrack->GetPosition().y();
		G4double z = fTrack->GetPosition().z();
		G4ThreeVector position = fTrack->GetPosition();
		G4ThreeVector momentum = fTrack->GetMomentum();

		if(fTrack->GetCurrentStepNumber()<=1) {
			fX0 = prePoint->GetPosition().x()/1000.;
			fY0 = prePoint->GetPosition().y()/1000.;
			fZ0 = prePoint->GetPosition().z()/1000.;
			fP0 = prePoint->GetMomentum().mag();
			fTh0 = prePoint->GetMomentum().getTheta();
			fPh0 = prePoint->GetMomentum().getPhi();	

			G4ThreeVector position0 = prePoint->GetPosition();	
			G4ThreeVector momentum0 = prePoint->GetMomentum();
			
			// Test transform with special position/angle		
	/*
			fX0 = 0.*mm;
			fY0 = 0.*mm;
			fZ0 = 0.;
			fP0 = 1.063*GeV;
			fTh0 = 5.*deg;
			fPh0 = 0.*deg;
			G4ThreeVector position0 = G4ThreeVector(fX0, fY0, fZ0);
			G4ThreeVector momentum0 = G4ThreeVector(fP0*sin(fTh0)*cos(fPh0),fP0*sin(fTh0)*sin(fPh0),fP0*cos(fTh0));
*/


			if(momentum0.x() > 0.){
				fLHRS = 1;
				septum_angle = -fSeptumAngle; 
				hrs_angle = -fHRSAngle;
				sign = -1.;
				
			}
			if(momentum0.x() < 0.){
				fRHRS = 1;
				septum_angle = fSeptumAngle;
				hrs_angle = fHRSAngle;
				sign = +1.;
			}

			G4RotationMatrix rotate_targ;
			rotate_targ.rotateY(septum_angle);	
			rotate_targ.rotateZ(90.*deg);
							
			G4AffineTransform transportAxis_targ = G4AffineTransform(rotate_targ);
			G4AffineTransform transport_targ = transportAxis_targ.Inverse();

			//Transform
			G4ThreeVector position0_tr = transport_targ.TransformPoint(position0);
			G4ThreeVector momentum0_tr = transport_targ.TransformPoint(momentum0);
		
			fX0_tr = position0_tr.x();
			fY0_tr = position0_tr.y();
			fZ0_tr = position0_tr.z();
			double dx = sin(momentum0_tr.theta())*cos(momentum0_tr.phi());
			double dy = sin(momentum0_tr.theta())*sin(momentum0_tr.phi());
			double dz = cos(momentum0_tr.theta());
			fTh0_tr = dx/dz;
			fPh0_tr = dy/dz;

			r0[0] = (float)fX0_tr;
			r0[1] = (float)(tan(fTh0_tr));
			r0[2] = (float)(fY0_tr*sign);			
			r0[3] = (float)(tan(fPh0_tr)*sign);
			r0[4] = (float)(fP0/fHRSMom);	

			// Only need to call transport functions once for event
			goodParticle = fTransportFunction->CallTransportFunction(r0, x_tf, th_tf, y_tf, ph_tf); 	
		
		}

		

	
		// Navigator of parallel world to get volume name (to identify virtual boundaries) and hall to transport coordinate transformation
		G4Navigator* fParallelNavigator = G4TransportationManager::GetTransportationManager()->GetNavigator("g4hrsparallel");

		G4String volName = fParallelNavigator->LocateGlobalPointAndSetup(G4ThreeVector(x,y,z))->GetName();
		G4AffineTransform trans = fParallelNavigator->GetGlobalToLocalTransform();
		G4ThreeVector position_tr = trans.TransformPoint(position);	
		G4ThreeVector momentum_tr = trans.TransformPoint(momentum);

		// Must explicitly define transform for septum, as these virtual boundaries are in hall coordinate system
		G4RotationMatrix rotate_sen, rotate_sm, rotate_sex;
		rotate_sen.rotateY(fSeptumAngle);
		rotate_sen.rotateZ(90.*deg);
		rotate_sm.rotateY((fSeptumAngle+fHRSAngle)/2.);
		rotate_sm.rotateZ(90.*deg);
		rotate_sex.rotateY(fHRSAngle);
		rotate_sex.rotateZ(90.*deg);
			
		// Axis transformation from HCS to TCS in each septum region
		G4AffineTransform transportAxis_sen = G4AffineTransform(rotate_sen);
		G4AffineTransform transportAxis_sm = G4AffineTransform(rotate_sm);
		G4AffineTransform transportAxis_sex = G4AffineTransform(rotate_sex);
		
		// Since we will be transforming POINTS (not AXES), we want the INVERSE of the axis transformation
		G4AffineTransform transport_sen = transportAxis_sen.Inverse();
		G4AffineTransform transport_sm = transportAxis_sm.Inverse();
		G4AffineTransform transport_sex = transportAxis_sex.Inverse();

		if(volName.find("virtualBoundaryPhys_sen") != G4String::npos) {
			// HCS
			fX_sen = x/1000.;
			fY_sen = y/1000.;
			fZ_sen = z/1000.;
			fTh_sen = momentum.theta()/rad;
			fPh_sen = momentum.phi()/rad;
			// TCS
			G4ThreeVector pos_tr = transport_sen.TransformPoint(G4ThreeVector(fX_sen, fY_sen, fZ_sen));
			G4ThreeVector mom_tr = transport_sen.TransformPoint(momentum);
			fX_sen_tr = pos_tr.x(); 			
			fY_sen_tr = pos_tr.y(); 			
			fZ_sen_tr = pos_tr.z(); 			
			fTh_sen_tr = (sin(mom_tr.theta())*cos(mom_tr.phi()))/cos(mom_tr.theta());	
			fPh_sen_tr = (sin(mom_tr.theta())*sin(mom_tr.phi()))/cos(mom_tr.theta());	
			// Transport function
			if(goodParticle) {
				fX_sen_tf = x_tf[sen];				
				fTh_sen_tf = th_tf[sen];				
				fY_sen_tf = y_tf[sen]*sign;	//				
				fPh_sen_tf = ph_tf[sen]*sign;	// Swap y/phi back if in left arm			
			}
		} else if(volName.find("virtualBoundaryPhys_sm") != G4String::npos) {
			// HCS
			fX_sm = x/1000.;
			fY_sm = y/1000.;
			fZ_sm = z/1000.;
			fTh_sm = momentum.theta()/rad;
			fPh_sm = momentum.phi()/rad;
			// TCS
			G4ThreeVector pos_tr = transport_sm.TransformPoint(G4ThreeVector(fX_sm, fY_sm, fZ_sm));
			G4ThreeVector mom_tr = transport_sm.TransformPoint(momentum);
			fX_sm_tr = pos_tr.x(); 			
			fY_sm_tr = pos_tr.y(); 			
			fZ_sm_tr = pos_tr.z(); 			
			fTh_sm_tr = (sin(mom_tr.theta())*cos(mom_tr.phi()))/cos(mom_tr.theta());	
			fPh_sm_tr = (sin(mom_tr.theta())*sin(mom_tr.phi()))/cos(mom_tr.theta());	
			// Transport function
			if(goodParticle) {
				fX_sm_tf = x_tf[sm];				
				fTh_sm_tf = th_tf[sm];				
				fY_sm_tf = y_tf[sm]*sign;		 				
				fPh_sm_tf = ph_tf[sm]*sign;				
			}
		} else if(volName.find("virtualBoundaryPhys_sex") != G4String::npos) {
			// HCS
			fX_sex = x/1000.;
			fY_sex = y/1000.;
			fZ_sex = z/1000.;
			fTh_sex = momentum.theta()/rad;
			fPh_sex = momentum.phi()/rad;
			// TCS
			G4ThreeVector pos_tr = transport_sex.TransformPoint(G4ThreeVector(fX_sex, fY_sex, fZ_sex));
			G4ThreeVector mom_tr = transport_sex.TransformPoint(momentum);
			fX_sex_tr = pos_tr.x(); 			
			fY_sex_tr = pos_tr.y(); 			
			fZ_sex_tr = pos_tr.z(); 			
			fTh_sex_tr = (sin(mom_tr.theta())*cos(mom_tr.phi()))/cos(mom_tr.theta());	
			fPh_sex_tr = (sin(mom_tr.theta())*sin(mom_tr.phi()))/cos(mom_tr.theta());	
			// Transport function
			if(goodParticle) {
				fX_sex_tf = x_tf[sex];				
				fTh_sex_tf = th_tf[sex];				
				fY_sex_tf = y_tf[sex]*sign;				
				fPh_sex_tf = ph_tf[sex]*sign;				
			}
		} else if(volName.find("virtualBoundaryPhys_col") != G4String::npos) {
			fX_col = x/1000.;
			fY_col = y/1000.;
			fZ_col = z/1000.;
			fTh_col = momentum.theta()/rad;
			fPh_col = momentum.phi()/rad;
		} else if(volName.find("virtualBoundaryPhys_q1en_LHRS") != G4String::npos) {
			// HCS
			fX_q1en = x/1000.;
			fY_q1en = y/1000.;
			fZ_q1en = z/1000.;
			fTh_q1en = momentum.theta()/rad;
			fPh_q1en = momentum.phi()/rad;
			// TCS
			fX_q1en_tr = position_tr.x();
			fY_q1en_tr = position_tr.y();
			fTh_q1en_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q1en_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
		} else if(volName.find("virtualBoundaryPhys_q1ex_LHRS") != G4String::npos) {
			// HCS
			fX_q1ex = x/1000.;
			fY_q1ex = y/1000.;
			fZ_q1ex = z/1000.;
			fTh_q1ex = momentum.theta()/rad;
			fPh_q1ex = momentum.phi()/rad;
			// TCS
			fX_q1ex_tr = position_tr.x();
			fY_q1ex_tr = position_tr.y();
			fTh_q1ex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q1ex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_q1ex_tf = x_tf[q1ex];
				fTh_q1ex_tf = th_tf[q1ex];
				fY_q1ex_tf = -y_tf[q1ex];
				fPh_q1ex_tf = -ph_tf[q1ex];
			}
		} else if(volName.find("virtualBoundaryPhys_q1en_RHRS") != G4String::npos) {
			// HCS
			fX_q1en = x/1000.;
			fY_q1en = y/1000.;
			fZ_q1en = z/1000.;
			fTh_q1en = momentum.theta()/rad;
			fPh_q1en = momentum.phi()/rad;
			// TCS
			fX_q1en_tr = position_tr.x();
			fY_q1en_tr = position_tr.y();
			fTh_q1en_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q1en_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
		} else if(volName.find("virtualBoundaryPhys_q1ex_RHRS") != G4String::npos) {
			// HCS
			fX_q1ex = x/1000.;
			fY_q1ex = y/1000.;
			fZ_q1ex = z/1000.;
			fTh_q1ex = momentum.theta()/rad;
			fPh_q1ex = momentum.phi()/rad;
			// TCS
			fX_q1ex_tr = position_tr.x();
			fY_q1ex_tr = position_tr.y();
			fTh_q1ex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q1ex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_q1ex_tf = x_tf[q1ex];
				fTh_q1ex_tf = th_tf[q1ex];
				fY_q1ex_tf = y_tf[q1ex];
				fPh_q1ex_tf = ph_tf[q1ex];
			}
		} else if(volName.find("virtualBoundaryPhys_q2en_LHRS") != G4String::npos) {
			// HCS
			fX_q2en = x/1000.;
			fY_q2en = y/1000.;
			fZ_q2en = z/1000.;
			fTh_q2en = momentum.theta()/rad;
			fPh_q2en = momentum.phi()/rad;
			// TCS
			fX_q2en_tr = position_tr.x();
			fY_q2en_tr = position_tr.y();
			fTh_q2en_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q2en_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
		} else if(volName.find("virtualBoundaryPhys_q2ex_LHRS") != G4String::npos) {
			// HCS
			fX_q2ex = x/1000.;
			fY_q2ex = y/1000.;
			fZ_q2ex = z/1000.;
			fTh_q2ex = momentum.theta()/rad;
			fPh_q2ex = momentum.phi()/rad;
			// TCS
			fX_q2ex_tr = position_tr.x();
			fY_q2ex_tr = position_tr.y();
			fTh_q2ex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q2ex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_q2ex_tf = x_tf[q2ex];
				fTh_q2ex_tf = th_tf[q2ex];
				fY_q2ex_tf = -y_tf[q2ex];
				fPh_q2ex_tf = -ph_tf[q2ex];
			}
		} else if(volName.find("virtualBoundaryPhys_q2en_RHRS") != G4String::npos) {
			// HCS
			fX_q2en = x/1000.;
			fY_q2en = y/1000.;
			fZ_q2en = z/1000.;
			fTh_q2en = momentum.theta()/rad;
			fPh_q2en = momentum.phi()/rad;
			// TCS
			fX_q2en_tr = position_tr.x();
			fY_q2en_tr = position_tr.y();
			fTh_q2en_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q2en_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
		} else if(volName.find("virtualBoundaryPhys_q2ex_RHRS") != G4String::npos) {
			// HCS
			fX_q2ex = x/1000.;
			fY_q2ex = y/1000.;
			fZ_q2ex = z/1000.;
			fTh_q2ex = momentum.theta()/rad;
			fPh_q2ex = momentum.phi()/rad;
			// TCS
			fX_q2ex_tr = position_tr.x();
			fY_q2ex_tr = position_tr.y();
			fTh_q2ex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q2ex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_q2ex_tf = x_tf[q2ex];
				fTh_q2ex_tf = th_tf[q2ex];
				fY_q2ex_tf = y_tf[q2ex];
				fPh_q2ex_tf = ph_tf[q2ex];
			}
		} else if(volName.find("virtualBoundaryPhys_den_LHRS") != G4String::npos) {
			// HCS
			fX_den = x/1000.;
			fY_den = y/1000.;
			fZ_den = z/1000.;
			fTh_den = momentum.theta()/rad;
			fPh_den = momentum.phi()/rad;
			// TCS
			fX_den_tr = position_tr.x();
			fY_den_tr = position_tr.y();
			fTh_den_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_den_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			if(goodParticle) {
				fX_den_tf = x_tf[den];
				fTh_den_tf = th_tf[den];
				fY_den_tf = -y_tf[den];
				fPh_den_tf = -ph_tf[den];
			}
		} else if(volName.find("virtualBoundaryPhys_dex_LHRS") != G4String::npos) {
			// HCS
			fX_dex = x/1000.;
			fY_dex = y/1000.;
			fZ_dex = z/1000.;
			fTh_dex = momentum.theta()/rad;
			fPh_dex = momentum.phi()/rad;
			// TCS
			fX_dex_tr = position_tr.x();
			fY_dex_tr = position_tr.y();
			fTh_dex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_dex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_dex_tf = x_tf[dex];
				fTh_dex_tf = th_tf[dex];
				fY_dex_tf = -y_tf[dex];
				fPh_dex_tf = -ph_tf[dex];
			}
		} else if(volName.find("virtualBoundaryPhys_den_RHRS") != G4String::npos) {
			// HCS
			fX_den = x/1000.;
			fY_den = y/1000.;
			fZ_den = z/1000.;
			fTh_den = momentum.theta()/rad;
			fPh_den = momentum.phi()/rad;
			// TCS
			fX_den_tr = position_tr.x();
			fY_den_tr = position_tr.y();
			fTh_den_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_den_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			if(goodParticle) {
				fX_den_tf = x_tf[den];
				fTh_den_tf = th_tf[den];
				fY_den_tf = y_tf[den];
				fPh_den_tf = ph_tf[den];
			}
		} else if(volName.find("virtualBoundaryPhys_dex_RHRS") != G4String::npos) {
			// HCS
			fX_dex = x/1000.;
			fY_dex = y/1000.;
			fZ_dex = z/1000.;
			fTh_dex = momentum.theta()/rad;
			fPh_dex = momentum.phi()/rad;
			// TCS
			fX_dex_tr = position_tr.x();
			fY_dex_tr = position_tr.y();
			fTh_dex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_dex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_dex_tf = x_tf[dex];
				fTh_dex_tf = th_tf[dex];
				fY_dex_tf = y_tf[dex];
				fPh_dex_tf = ph_tf[dex];
			}
		} else if(volName.find("virtualBoundaryPhys_q3en_LHRS") != G4String::npos) {
			// HCS
			fX_q3en = x/1000.;
			fY_q3en = y/1000.;
			fZ_q3en = z/1000.;
			fTh_q3en = momentum.theta()/rad;
			fPh_q3en = momentum.phi()/rad;
			// TCS
			fX_q3en_tr = position_tr.x();
			fY_q3en_tr = position_tr.y();
			fTh_q3en_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q3en_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			if(goodParticle) {
				fX_q3en_tf = x_tf[q3en];
				fTh_q3en_tf = th_tf[q3en];
				fY_q3en_tf = -y_tf[q3en];
				fPh_q3en_tf = -ph_tf[q3en];
			}
		} else if(volName.find("virtualBoundaryPhys_q3ex_LHRS") != G4String::npos) {
			// HCS
			fX_q3ex = x/1000.;
			fY_q3ex = y/1000.;
			fZ_q3ex = z/1000.;
			fTh_q3ex = momentum.theta()/rad;
			fPh_q3ex = momentum.phi()/rad;
			// TCS
			fX_q3ex_tr = position_tr.x();
			fY_q3ex_tr = position_tr.y();
			fTh_q3ex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q3ex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_q3ex_tf = x_tf[q3ex];
				fTh_q3ex_tf = th_tf[q3ex];
				fY_q3ex_tf = -y_tf[q3ex];
				fPh_q3ex_tf = -ph_tf[q3ex];
			}
		} else if(volName.find("virtualBoundaryPhys_q3en_RHRS") != G4String::npos) {
			// HCS
			fX_q3en = x/1000.;
			fY_q3en = y/1000.;
			fZ_q3en = z/1000.;
			fTh_q3en = momentum.theta()/rad;
			fPh_q3en = momentum.phi()/rad;
			// TCS
			fX_q3en_tr = position_tr.x();
			fY_q3en_tr = position_tr.y();
			fTh_q3en_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q3en_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			if(goodParticle) {
				fX_q3en_tf = x_tf[q3en];
				fTh_q3en_tf = th_tf[q3en];
				fY_q3en_tf = y_tf[q3en];
				fPh_q3en_tf = ph_tf[q3en];
			}
		} else if(volName.find("virtualBoundaryPhys_q3ex_RHRS") != G4String::npos) {
			// HCS
			fX_q3ex = x/1000.;
			fY_q3ex = y/1000.;
			fZ_q3ex = z/1000.;
			fTh_q3ex = momentum.theta()/rad;
			fPh_q3ex = momentum.phi()/rad;
			// TCS
			fX_q3ex_tr = position_tr.x();
			fY_q3ex_tr = position_tr.y();
			fTh_q3ex_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_q3ex_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_q3ex_tf = x_tf[q3ex];
				fTh_q3ex_tf = th_tf[q3ex];
				fY_q3ex_tf = y_tf[q3ex];
				fPh_q3ex_tf = ph_tf[q3ex];
			}
		} else if(volName.find("virtualBoundaryPhys_vdc_LHRS") != G4String::npos) {
			// HCS
			fX_vdc = x/1000.;
			fY_vdc = y/1000.;
			fZ_vdc = z/1000.;
			fTh_vdc = momentum.theta()/rad;
			fPh_vdc = momentum.phi()/rad;
			// TCS
			fX_vdc_tr = position_tr.x();
			fY_vdc_tr = position_tr.y();
			fTh_vdc_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_vdc_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
 		} else if(volName.find("virtualBoundaryPhys_vdc_RHRS") != G4String::npos) {
			// HCS
			fX_vdc = x/1000.;
			fY_vdc = y/1000.;
			fZ_vdc = z/1000.;
			fTh_vdc = momentum.theta()/rad;
			fPh_vdc = momentum.phi()/rad;
			// TCS
			fX_vdc_tr = position_tr.x();
			fY_vdc_tr = position_tr.y();
			fTh_vdc_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_vdc_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
 		} else if(volName.find("virtualBoundaryPhys_fp_LHRS") != G4String::npos) {
			// HCS
			fX_fp = x/1000.;
			fY_fp = y/1000.;
			fZ_fp = z/1000.;
			fTh_fp = momentum.theta()/rad;
			fPh_fp = momentum.phi()/rad;
			// TCS
			fX_fp_tr = position_tr.x();
			fY_fp_tr = position_tr.y();
			fTh_fp_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_fp_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_fp_tf = x_tf[fp];
				fTh_fp_tf = th_tf[fp];
				fY_fp_tf = -y_tf[fp];
				fPh_fp_tf = -ph_tf[fp];
			}
		} else if(volName.find("virtualBoundaryPhys_fp_LHRS") != G4String::npos) {
			// HCS
			fX_fp = x/1000.;
			fY_fp = y/1000.;
			fZ_fp = z/1000.;
			fTh_fp = momentum.theta()/rad;
			fPh_fp = momentum.phi()/rad;
			// TCS
			fX_fp_tr = position_tr.x();
			fY_fp_tr = position_tr.y();
			fTh_fp_tr = (sin(momentum_tr.theta())*cos(momentum_tr.phi()))/cos(momentum_tr.theta());	
			fPh_fp_tr = (sin(momentum_tr.theta())*sin(momentum_tr.phi()))/cos(momentum_tr.theta());
			// Transport function
			if(goodParticle) {
				fX_fp_tf = x_tf[fp];
				fTh_fp_tf = th_tf[fp];
				fY_fp_tf = y_tf[fp];
				fPh_fp_tf = ph_tf[fp];
			}
		}
 
	}  // end if ParentID == 0
 

}


