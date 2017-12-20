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
	
	fSeptumAngle = 5.*deg;
	fHRSAngle = 12.5*deg;
	fHRSMomentum = 1.063*GeV;
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

		// Get hall coordinates in meters to output to ROOT	
		G4double x = fTrack->GetPosition().x()/1000.;
		G4double y = fTrack->GetPosition().y()/1000.;
		G4double z = fTrack->GetPosition().z()/1000.;
		// Get position 3-vector in millimeters for transforms
		G4ThreeVector position = fTrack->GetPosition();
		// Get momentum 3-vector
		G4ThreeVector momentum = fTrack->GetMomentum();

		if(fTrack->GetCurrentStepNumber()<=1) {
		
			G4ThreeVector position0 = prePoint->GetPosition()/1000.;	
			G4ThreeVector momentum0 = prePoint->GetMomentum();
		
			fX0 = position0.x();
			fY0 = position0.y();
			fZ0 = position0.z();
			fP0 = momentum0.mag();
			fTh0 = momentum0.getTheta();
			fPh0 = momentum0.getPhi();	

			
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
			fTh0_tr = momentum0_tr.x()/momentum0_tr.z();
			fPh0_tr = momentum0_tr.y()/momentum0_tr.z();

			r0[0] = (float)fX0_tr;
			r0[1] = (float)(fTh0_tr);
			r0[2] = (float)(fY0_tr*sign);			
			r0[3] = (float)(fPh0_tr*sign);
			r0[4] = (float)((fP0-fHRSMomentum)/fHRSMomentum);	

			// Only need to call transport functions once for event
			goodParticle = fTransportFunction->CallTransportFunction(r0, x_tf, th_tf, y_tf, ph_tf); 	
		
		}
	
/* TROUBLESHOOTING 
		// Center of LHRS focal plane in HCS for coordinate transform test
		G4ThreeVector fphc = G4ThreeVector(21770.815*sin(fHRSAngle),8330.421,1053.79+21770.815*cos(fHRSAngle));
		// Unit momentum or central trajectory exiting dipole in HCS for coordinate transform test
		G4ThreeVector cthc = G4ThreeVector(cos(45.*deg)*sin(fHRSAngle),sin(45.*deg),cos(45.*deg)*cos(fHRSAngle));
*/
		
		
		// Navigator of parallel world to get volume name (to identify virtual boundaries) and hall to transport coordinate transformation
		G4Navigator* fParallelNavigator = G4TransportationManager::GetTransportationManager()->GetNavigator("g4hrsparallel");
		G4String volName = fParallelNavigator->LocateGlobalPointAndSetup(position)->GetName();

		if(volName.find("virtualBoundaryPhys") != G4String::npos) { 
		

			//This gets the full transform (rotation + translation) that is required for transforming POSITION		
			G4AffineTransform position_transform = fParallelNavigator->GetGlobalToLocalTransform();
			//However, the MOMENTUM transform requires rotation ONLY!! 
			G4RotationMatrix momentum_rotation = position_transform.NetRotation();
			G4AffineTransform momentum_transform = G4AffineTransform(momentum_rotation);

			G4ThreeVector position_tr = position_transform.TransformPoint(position)/1000.;	
			G4ThreeVector momentum_tr = momentum_transform.TransformPoint(momentum);
	
			// Must explicitly define transform for septum, as these virtual boundaries are in hall coordinate system
			//Rotations
			G4RotationMatrix rotate_sen, rotate_sm, rotate_sex;
			rotate_sen.rotateY(septum_angle);
			rotate_sen.rotateZ(90.*deg);
			rotate_sm.rotateY((septum_angle + hrs_angle)/2.);
			rotate_sm.rotateZ(90.*deg);
			rotate_sex.rotateY(hrs_angle);
			rotate_sex.rotateZ(90.*deg);
			//Translations
			double pivotZOffset = 105.379*cm;
			double septumZPosition = 69.99937*cm; //from snake, matches values in world and parallel world construction
			septumZPosition+=pivotZOffset; //put origin at target center, not Hall A pivot
			double septumLength = 74.*cm;
			double z_sen = septumZPosition - septumLength/2.;
			double x_sen = z_sen*tan(fSeptumAngle);
			//Length of chord connecting central trajectory at septum entrance and septum exit
			double chord = septumLength/cos((fHRSAngle+fSeptumAngle)/2.);
			//Radius of circle for central trajectory in ideal septum
			double R = sqrt((chord*chord)/(2.*(1-cos(fHRSAngle-fSeptumAngle))));
			//Center of aforementioned circle
			double xc = R*cos(fSeptumAngle) + x_sen;
			double zc = z_sen - R*sin(fSeptumAngle);
			//x,z position of central trajectory at septum middle
			double z_sm = septumZPosition;
			double x_sm = xc - sqrt(R*R - (z_sm - zc)*(z_sm - zc));
			//x,z position of central trajectory at septum exit
			double z_sex = septumZPosition + septumLength/2.;
			double x_sex = xc - sqrt(R*R - (z_sex - zc)*(z_sex - zc));
		
			if(fRHRS) {
				x_sen*=-1.;
				x_sm*=-1.;
				x_sex*=-1.;
			}
	
			G4ThreeVector translate_sen = G4ThreeVector(x_sen,0.,z_sen);
			G4ThreeVector translate_sm = G4ThreeVector(x_sm,0.,z_sm);
			G4ThreeVector translate_sex = G4ThreeVector(x_sex,0.,z_sex);

			G4AffineTransform position_transform_sen = G4AffineTransform(rotate_sen,translate_sen);
			position_transform_sen.Invert();
			G4AffineTransform position_transform_sm  = G4AffineTransform(rotate_sm, translate_sm );
			position_transform_sm.Invert();
			G4AffineTransform position_transform_sex = G4AffineTransform(rotate_sex,translate_sex);
			position_transform_sex.Invert();

			G4AffineTransform momentum_transform_sen = G4AffineTransform(rotate_sen);
			momentum_transform_sen.Invert();
			G4AffineTransform momentum_transform_sm  = G4AffineTransform(rotate_sm);
			momentum_transform_sm.Invert();
			G4AffineTransform momentum_transform_sex = G4AffineTransform(rotate_sex);
			momentum_transform_sex.Invert();
			
			//TROUBLESHOOTING
			G4ThreeVector senhc = G4ThreeVector(x_sen,-2.,z_sen);
			G4ThreeVector semhc = G4ThreeVector(x_sm,0.,z_sm);
			G4ThreeVector sexhc = G4ThreeVector(x_sex,0.,z_sex);
			G4ThreeVector senct = G4ThreeVector(sin(fSeptumAngle),0.,cos(fSeptumAngle));
			G4ThreeVector semct = G4ThreeVector(sin((fSeptumAngle+fHRSAngle)/2.),0.,cos((fSeptumAngle+fHRSAngle)/2.));
			G4ThreeVector sexct = G4ThreeVector(sin(fHRSAngle),0.,cos(fHRSAngle));
			G4ThreeVector sentr = position_transform_sen.TransformPoint(senhc);
			G4ThreeVector sencttr = momentum_transform_sen.TransformPoint(senct);
			G4ThreeVector semtr = position_transform_sm.TransformPoint(semhc);
			G4ThreeVector semcttr = momentum_transform_sm.TransformPoint(semct);
			G4ThreeVector sextr = position_transform_sex.TransformPoint(sexhc);
			G4ThreeVector sexcttr = momentum_transform_sex.TransformPoint(sexct);
	
			if(volName == "virtualBoundaryPhys_sen") {
				// HCS
				fX_sen = x;
				fY_sen = y;
				fZ_sen = z;
				fTh_sen = momentum.theta()/rad;
				fPh_sen = momentum.phi()/rad;
				// TCS
				G4ThreeVector pos_tr = position_transform_sen.TransformPoint(position)/1000.;
				G4ThreeVector mom_tr = momentum_transform_sen.TransformPoint(momentum);
				fX_sen_tr = pos_tr.x(); 			
				fY_sen_tr = pos_tr.y(); 			
				fZ_sen_tr = pos_tr.z(); 			
				fTh_sen_tr = mom_tr.x()/mom_tr.z();	
				fPh_sen_tr = mom_tr.y()/mom_tr.z();	
				// Transport function
				if(goodParticle) {
					fX_sen_tf = x_tf[sen];				
					fTh_sen_tf = th_tf[sen];				
					fY_sen_tf = y_tf[sen]*sign;	//				
					fPh_sen_tf = ph_tf[sen]*sign;	// Swap y/phi back if in left arm			
				}
			} else if(volName == "virtualBoundaryPhys_sm") {
				// HCS
				fX_sm = x;
				fY_sm = y;
				fZ_sm = z;
				fTh_sm = momentum.theta()/rad;
				fPh_sm = momentum.phi()/rad;
				// TCS
				G4ThreeVector pos_tr = position_transform_sm.TransformPoint(position)/1000.;
				G4ThreeVector mom_tr = momentum_transform_sm.TransformPoint(momentum);
				fX_sm_tr = pos_tr.x(); 			
				fY_sm_tr = pos_tr.y(); 			
				fZ_sm_tr = pos_tr.z(); 			
				fTh_sm_tr = mom_tr.x()/mom_tr.z();	
				fPh_sm_tr = mom_tr.y()/mom_tr.z();	
				// Transport function
				if(goodParticle) {
					fX_sm_tf = x_tf[sm];				
					fTh_sm_tf = th_tf[sm];				
					fY_sm_tf = y_tf[sm]*sign;		 				
					fPh_sm_tf = ph_tf[sm]*sign;				
				}
			} else if(volName == "virtualBoundaryPhys_sex") {
				// HCS
				fX_sex = x;
				fY_sex = y;
				fZ_sex = z;
				fTh_sex = momentum.theta()/rad;
				fPh_sex = momentum.phi()/rad;
				// TCS
				G4ThreeVector pos_tr = position_transform_sex.TransformPoint(position)/1000.;
				G4ThreeVector mom_tr = momentum_transform_sex.TransformPoint(momentum);
				fX_sex_tr = pos_tr.x(); 			
				fY_sex_tr = pos_tr.y(); 			
				fZ_sex_tr = pos_tr.z(); 			
				fTh_sex_tr = mom_tr.x()/mom_tr.z();	
				fPh_sex_tr = mom_tr.y()/mom_tr.z();	
				// Transport function
				if(goodParticle) {
					fX_sex_tf = x_tf[sex];				
					fTh_sex_tf = th_tf[sex];				
					fY_sex_tf = y_tf[sex]*sign;				
					fPh_sex_tf = ph_tf[sex]*sign;				
				}
			} else if(volName == "virtualBoundaryPhys_col") {
				fX_col = x;
				fY_col = y;
				fZ_col = z;
				fTh_col = momentum.theta()/rad;
				fPh_col = momentum.phi()/rad;
			} else if(volName == "virtualBoundaryPhys_q1en_LHRS") {
				// HCS
				fX_q1en = x;
				fY_q1en = y;
				fZ_q1en = z;
				fTh_q1en = momentum.theta()/rad;
				fPh_q1en = momentum.phi()/rad;
				// TCS
				fX_q1en_tr = position_tr.x();
				fY_q1en_tr = position_tr.y();
				fTh_q1en_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q1en_tr = momentum_tr.y()/momentum_tr.z();
			} else if(volName == "virtualBoundaryPhys_q1ex_LHRS") {
				// HCS
				fX_q1ex = x;
				fY_q1ex = y;
				fZ_q1ex = z;
				fTh_q1ex = momentum.theta()/rad;
				fPh_q1ex = momentum.phi()/rad;
				// TCS
				fX_q1ex_tr = position_tr.x();
				fY_q1ex_tr = position_tr.y();
				fTh_q1ex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q1ex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q1ex_tf = x_tf[q1ex];
					fTh_q1ex_tf = th_tf[q1ex];
					fY_q1ex_tf = -y_tf[q1ex];
					fPh_q1ex_tf = -ph_tf[q1ex];
				}
			} else if(volName == "virtualBoundaryPhys_q1en_RHRS") {
				// HCS
				fX_q1en = x;
				fY_q1en = y;
				fZ_q1en = z;
				fTh_q1en = momentum.theta()/rad;
				fPh_q1en = momentum.phi()/rad;
				// TCS
				fX_q1en_tr = position_tr.x();
				fY_q1en_tr = position_tr.y();
				fTh_q1en_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q1en_tr = momentum_tr.y()/momentum_tr.z();
			} else if(volName == "virtualBoundaryPhys_q1ex_RHRS") {
				// HCS
				fX_q1ex = x;
				fY_q1ex = y;
				fZ_q1ex = z;
				fTh_q1ex = momentum.theta()/rad;
				fPh_q1ex = momentum.phi()/rad;
				// TCS
				fX_q1ex_tr = position_tr.x();
				fY_q1ex_tr = position_tr.y();
				fTh_q1ex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q1ex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q1ex_tf = x_tf[q1ex];
					fTh_q1ex_tf = th_tf[q1ex];
					fY_q1ex_tf = y_tf[q1ex];
					fPh_q1ex_tf = ph_tf[q1ex];
				}
			} else if(volName == "virtualBoundaryPhys_q2en_LHRS") {
				// HCS
				fX_q2en = x;
				fY_q2en = y;
				fZ_q2en = z;
				fTh_q2en = momentum.theta()/rad;
				fPh_q2en = momentum.phi()/rad;
				// TCS
				fX_q2en_tr = position_tr.x();
				fY_q2en_tr = position_tr.y();
				fTh_q2en_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q2en_tr = momentum_tr.y()/momentum_tr.z();
			} else if(volName == "virtualBoundaryPhys_q2ex_LHRS") {
				// HCS
				fX_q2ex = x;
				fY_q2ex = y;
				fZ_q2ex = z;
				fTh_q2ex = momentum.theta()/rad;
				fPh_q2ex = momentum.phi()/rad;
				// TCS
				fX_q2ex_tr = position_tr.x();
				fY_q2ex_tr = position_tr.y();
				fTh_q2ex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q2ex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q2ex_tf = x_tf[q2ex];
					fTh_q2ex_tf = th_tf[q2ex];
					fY_q2ex_tf = -y_tf[q2ex];
					fPh_q2ex_tf = -ph_tf[q2ex];
				}
			} else if(volName == "virtualBoundaryPhys_q2en_RHRS") {
				// HCS
				fX_q2en = x;
				fY_q2en = y;
				fZ_q2en = z;
				fTh_q2en = momentum.theta()/rad;
				fPh_q2en = momentum.phi()/rad;
				// TCS
				fX_q2en_tr = position_tr.x();
				fY_q2en_tr = position_tr.y();
				fTh_q2en_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q2en_tr = momentum_tr.y()/momentum_tr.z();
			} else if(volName == "virtualBoundaryPhys_q2ex_RHRS") {
				// HCS
				fX_q2ex = x;
				fY_q2ex = y;
				fZ_q2ex = z;
				fTh_q2ex = momentum.theta()/rad;
				fPh_q2ex = momentum.phi()/rad;
				// TCS
				fX_q2ex_tr = position_tr.x();
				fY_q2ex_tr = position_tr.y();
				fTh_q2ex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q2ex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q2ex_tf = x_tf[q2ex];
					fTh_q2ex_tf = th_tf[q2ex];
					fY_q2ex_tf = y_tf[q2ex];
					fPh_q2ex_tf = ph_tf[q2ex];
				}
			} else if(volName == "virtualBoundaryPhys_den_LHRS") {
				// HCS
				fX_den = x;
				fY_den = y;
				fZ_den = z;
				fTh_den = momentum.theta()/rad;
				fPh_den = momentum.phi()/rad;
				// TCS
				fX_den_tr = position_tr.x();
				fY_den_tr = position_tr.y();
				fTh_den_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_den_tr = momentum_tr.y()/momentum_tr.z();
				if(goodParticle) {
					fX_den_tf = x_tf[den];
					fTh_den_tf = th_tf[den];
					fY_den_tf = -y_tf[den];
					fPh_den_tf = -ph_tf[den];
				}
			} else if(volName == "virtualBoundaryPhys_dex_LHRS") {
				// HCS
				fX_dex = x;
				fY_dex = y;
				fZ_dex = z;
				fTh_dex = momentum.theta()/rad;
				fPh_dex = momentum.phi()/rad;
				// TCS
				fX_dex_tr = position_tr.x();
				fY_dex_tr = position_tr.y();
				fTh_dex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_dex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_dex_tf = x_tf[dex];
					fTh_dex_tf = th_tf[dex];
					fY_dex_tf = -y_tf[dex];
					fPh_dex_tf = -ph_tf[dex];
				}
			} else if(volName == "virtualBoundaryPhys_den_RHRS") {
				// HCS
				fX_den = x;
				fY_den = y;
				fZ_den = z;
				fTh_den = momentum.theta()/rad;
				fPh_den = momentum.phi()/rad;
				// TCS
				fX_den_tr = position_tr.x();
				fY_den_tr = position_tr.y();
				fTh_den_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_den_tr = momentum_tr.y()/momentum_tr.z();
				if(goodParticle) {
					fX_den_tf = x_tf[den];
					fTh_den_tf = th_tf[den];
					fY_den_tf = y_tf[den];
					fPh_den_tf = ph_tf[den];
				}
			} else if(volName == "virtualBoundaryPhys_dex_RHRS") {
				// HCS
				fX_dex = x;
				fY_dex = y;
				fZ_dex = z;
				fTh_dex = momentum.theta()/rad;
				fPh_dex = momentum.phi()/rad;
				// TCS
				fX_dex_tr = position_tr.x();
				fY_dex_tr = position_tr.y();
				fTh_dex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_dex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_dex_tf = x_tf[dex];
					fTh_dex_tf = th_tf[dex];
					fY_dex_tf = y_tf[dex];
					fPh_dex_tf = ph_tf[dex];
				}
			} else if(volName == "virtualBoundaryPhys_q3en_LHRS") {
				// HCS
				fX_q3en = x;
				fY_q3en = y;
				fZ_q3en = z;
				fTh_q3en = momentum.theta()/rad;
				fPh_q3en = momentum.phi()/rad;
				// TCS
				fX_q3en_tr = position_tr.x();
				fY_q3en_tr = position_tr.y();
				fTh_q3en_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q3en_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q3en_tf = x_tf[q3en];
					fTh_q3en_tf = th_tf[q3en];
					fY_q3en_tf = -y_tf[q3en];
					fPh_q3en_tf = -ph_tf[q3en];
				}
			} else if(volName == "virtualBoundaryPhys_q3ex_LHRS") {
				// HCS
				fX_q3ex = x;
				fY_q3ex = y;
				fZ_q3ex = z;
				fTh_q3ex = momentum.theta()/rad;
				fPh_q3ex = momentum.phi()/rad;
				// TCS
				fX_q3ex_tr = position_tr.x();
				fY_q3ex_tr = position_tr.y();
				fTh_q3ex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q3ex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q3ex_tf = x_tf[q3ex];
					fTh_q3ex_tf = th_tf[q3ex];
					fY_q3ex_tf = -y_tf[q3ex];
					fPh_q3ex_tf = -ph_tf[q3ex];
				}
			} else if(volName == "virtualBoundaryPhys_q3en_RHRS") {
				// HCS
				fX_q3en = x;
				fY_q3en = y;
				fZ_q3en = z;
				fTh_q3en = momentum.theta()/rad;
				fPh_q3en = momentum.phi()/rad;
				// TCS
				fX_q3en_tr = position_tr.x();
				fY_q3en_tr = position_tr.y();
				fTh_q3en_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q3en_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q3en_tf = x_tf[q3en];
					fTh_q3en_tf = th_tf[q3en];
					fY_q3en_tf = y_tf[q3en];
					fPh_q3en_tf = ph_tf[q3en];
				}
			} else if(volName == "virtualBoundaryPhys_q3ex_RHRS") {
				// HCS
				fX_q3ex = x;
				fY_q3ex = y;
				fZ_q3ex = z;
				fTh_q3ex = momentum.theta()/rad;
				fPh_q3ex = momentum.phi()/rad;
				// TCS
				fX_q3ex_tr = position_tr.x();
				fY_q3ex_tr = position_tr.y();
				fTh_q3ex_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_q3ex_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_q3ex_tf = x_tf[q3ex];
					fTh_q3ex_tf = th_tf[q3ex];
					fY_q3ex_tf = y_tf[q3ex];
					fPh_q3ex_tf = ph_tf[q3ex];
				}
			} else if(volName == "virtualBoundaryPhys_vdc_LHRS") {
				// HCS
				fX_vdc = x;
				fY_vdc = y;
				fZ_vdc = z;
				fTh_vdc = momentum.theta()/rad;
				fPh_vdc = momentum.phi()/rad;
				// TCS
				fX_vdc_tr = position_tr.x();
				fY_vdc_tr = position_tr.y();
				fTh_vdc_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_vdc_tr = momentum_tr.y()/momentum_tr.z();
			} else if(volName == "virtualBoundaryPhys_vdc_RHRS") {
				// HCS
				fX_vdc = x;
				fY_vdc = y;
				fZ_vdc = z;
				fTh_vdc = momentum.theta()/rad;
				fPh_vdc = momentum.phi()/rad;
				// TCS
				fX_vdc_tr = position_tr.x();
				fY_vdc_tr = position_tr.y();
				fTh_vdc_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_vdc_tr = momentum_tr.y()/momentum_tr.z();
			} else if(volName == "virtualBoundaryPhys_fp_LHRS") {
				// HCS
				fX_fp = x;
				fY_fp = y;
				fZ_fp = z;
				fTh_fp = momentum.theta()/rad;
				fPh_fp = momentum.phi()/rad;
				// TCS
				fX_fp_tr = position_tr.x();
				fY_fp_tr = position_tr.y();
				fTh_fp_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_fp_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_fp_tf = x_tf[fp];
					fTh_fp_tf = th_tf[fp];
					fY_fp_tf = -y_tf[fp];
					fPh_fp_tf = -ph_tf[fp];
				}
			} else if(volName == "virtualBoundaryPhys_fp_RHRS") {
				// HCS
				fX_fp = x;
				fY_fp = y;
				fZ_fp = z;
				fTh_fp = momentum.theta()/rad;
				fPh_fp = momentum.phi()/rad;
				// TCS
				fX_fp_tr = position_tr.x();
				fY_fp_tr = position_tr.y();
				fTh_fp_tr = momentum_tr.x()/momentum_tr.z();	
				fPh_fp_tr = momentum_tr.y()/momentum_tr.z();
				// Transport function
				if(goodParticle) {
					fX_fp_tf = x_tf[fp];
					fTh_fp_tf = th_tf[fp];
					fY_fp_tf = y_tf[fp];
					fPh_fp_tf = ph_tf[fp];
				}
			}
	
		} //end if volName contains virtualBoundaryPhys
 
	}  // end if ParentID == 0
 

}
