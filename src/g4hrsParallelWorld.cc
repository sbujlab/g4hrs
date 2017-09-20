#include "g4hrsParallelWorld.hh"
#include "g4hrsGenericDetector.hh"
#include "g4hrsBeamTarget.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"
#include "g4hrsIO.hh"
#include "g4hrsMaterial.hh"
#include "g4hrsEMFieldSetup.hh"


#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "globals.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "G4ios.hh"

#include "G4UnitsTable.hh"

#include "G4Polycone.hh"
#include "G4Tubs.hh"


//visual
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Box.hh"

#define __DET_STRLEN 200
#define __MAX_DETS 5000

g4hrsParallelWorld::g4hrsParallelWorld(G4String parallelWorldName) :G4VUserParallelWorld(parallelWorldName) {

    // Default geometry file


    fHRSAngle=12.5*deg;
    fSeptumAngle=5.0*deg;

    fTargetW = 2.0*2.54*cm;
    fTargetH = 2.0*2.54*cm;
    fTargetL = 5*mm;

    fTargetX = 0.0;
    fTargetY = 0.0;
    fTargetZ = 0.0;

    fSeptum_UseUniformB=0;

    fPivot2SieveFace = 800*mm;
    fPivotXOffset =  0.0*deg;
    fPivotYOffset =  0.0*deg;
    fPivotZOffset =  1053.79*mm;

    fSetupSieveSlit = false;

    fUseSeptumPlusStdHRS=0;
    fSnakeModel = 49;

    fSetupStdScatChamber = 0;

    //#SetupHRS have the following candidates: 
    //# 0: do nothing; 1: will build septum, sieve and VB; 
    //# 2: add Q1; 3: add Q2 ; 4: add Dipole and Q3  
    fSetupHRS = 4;


    fScatChamberRin   =  484.025*mm;
    fScatChamberRout  =  509.425*mm;
    fScatChamberL     =  692.15*mm;

    fScatChamberXOffset = 0.0;
    fScatChamberYOffset = 0.0;
    fScatChamberZOffset = -1053.79*mm;
    fSetupCREXGeometry = true;

    //fEMFieldSetup = NULL;

    fFieldX = 45000*mm;
    fFieldY = 45000*mm;
    fFieldZ = 45000*mm;

}

g4hrsParallelWorld::~g4hrsParallelWorld() {
    }

void g4hrsParallelWorld::Construct() {
    
	G4VPhysicalVolume *worldVolume = GetWorld();
//	CreateSeptum(worldVolume->GetLogicalVolume());
//	CreateHRS(worldVolume->GetLogicalVolume());
	ConstructSD(worldVolume->GetLogicalVolume());

	return;
}


void g4hrsParallelWorld::CreateHRS(G4LogicalVolume* motherLogical)
{
	int IS_NIM = 0;//if it's 0, then we go by SNAKE, if it is 1, we go by NIM.
	//As it turns out, I think we want to go by SNAKE, so that is why it is zero right now.
	G4VPhysicalVolume* theHRSPhys=0;

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	G4RotationMatrix *pRotX45deg=new G4RotationMatrix();
	pRotX45deg->rotateX(-45*deg); 
	G4RotationMatrix *pRotX30deg=new G4RotationMatrix();
	pRotX30deg->rotateX(60*deg); 
	G4RotationMatrix *pRotX105deg=new G4RotationMatrix();
	pRotX105deg->rotateX(165*deg); 
	G4RotationMatrix *pRotXLHRSdeg=new G4RotationMatrix();
	pRotXLHRSdeg->rotateY(fHRSAngle); 
	G4RotationMatrix *pRotXRHRSdeg=new G4RotationMatrix();
	pRotXRHRSdeg->rotateY(-fHRSAngle); 

	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String SDname;
	G4VSensitiveDetector* Q1WindowSD=new g4hrsGenericDetector(SDname="Q1Window", 2);
	//Build the HRS, positions were taken fron the NIM paper
	//The Apertures are defined as:
	//Q1 exit is a circle of radius 0.15m
	//Q3 entrance and exit are circles of radius 0.30m
	//The dipole entrance and exit are trapezoids:
	//-0.40 m<x<0.40 m (+x is down) 
	//y=+-(0.125*(1-(1.25*x/8.40)) (smallest gap at positive x, x in m)

	/////////////////////////
	// HRS QQDQ Containner
	/////////////////////////
	//Build a container using polycone, covering 270+/-12 degrees, 1.46m to 20.76m
	//3.05m below the beam line is the ground, need to subtract everything beneath that
	//looks like this:
	/*
	HRS container:covering 24 degrees in X-Z plane, 10.05 m height.
	Stuff inside this containner will also use the pivot as the origin, but do not 
	need to worry about the rotation of the HRS.
	//                         7.0m above beam line 
	//                       -----------------------------| 20.76m from pivot,
	//                      /                             |
	//                     /                        Q3    |
	//                    /                               |
	//                   /                       E        /
	//                  /                     L           |
	//          --------                   O              |
	//         /                        P                 |
	//    -----                     I                     |
	//----|----- Q1 --- Q2 ---- D     --------------------|------ beam line -----
	//    |                                               |
	//    ------------------------------------------------|
	//    1.46m from povot, 3.05 m below beam line

	*/

	//double pHRSContainerRin=1.46*m,pHRSContainerRout=20.76*m;//JixieMode
	double pHRSContainerRin=1.46*m,pHRSContainerRout=25*m;//NickieMode it just has to be big enough to fit my VB at FP
	//double pBeamLine2Ground;
	//pBeamLine2Ground=-3.05*m;
	//build the container with polycone

	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {pHRSContainerRin,pHRSContainerRin,2.5*m,
		3.7*m,9.0*m,pHRSContainerRout-3.0*m,pHRSContainerRout};
	double rOuter_HRSContainer[] = {pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,
		pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,pHRSContainerRout};
	//double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
	//2.0*m,2.0*m,7.0*m,7.0*m};
	double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
					2.0*m,2.0*m,11.0*m,11.0*m};
	G4Polycone* HRSContainerSolid = new G4Polycone("HRSContainer",258.0*deg,24.0*deg,
		kNPlane_HRSContainer,zPlane_HRSContainer,rInner_HRSContainer,rOuter_HRSContainer);

	////build the container using tube
	//double pHRSContainerHeight=14.0*m;
	//G4VSolid* HRSContainerTub = new G4Tubs("HRSContainerTub",pHRSContainerRin,
	//	pHRSContainerRout,pHRSContainerHeight/2.0,258.0*deg,24.0*deg);
	////ground
	//double pGoundHeight=10*m;
	//G4VSolid* GroundTub = new G4Tubs("GroundTub",0,30*m,pGoundHeight/2,0*deg,360.0*deg);
	////tube subtract the ground
	//G4SubtractionSolid* HRSContainerSolid = new G4SubtractionSolid("HRSContainer",
	//		HRSContainerTub,GroundTub,
	//		0,G4ThreeVector(0,0,pBeamLine2Ground-pGoundHeight/2));

	G4LogicalVolume* LHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		0,"LHRSContainerLogical",0,0,0);
	G4LogicalVolume* RHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		0,"RHRSContainerLogical",0,0,0);

	G4VisAttributes* MagFieldVisAtt = new G4VisAttributes(G4Colour(1., 0., 1.));	

        G4VisAttributes *HallVisAtt = new G4VisAttributes(G4Colour(0.,1.0,0.));
//        HallVisAtt->SetVisibility(false);


	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt);
	//RHRSContainerLogical->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotLHRSContainer=new G4RotationMatrix();
	pRotLHRSContainer->rotateX(-270*deg);
	pRotLHRSContainer->rotateZ(-fHRSAngle);  

	G4RotationMatrix *pRotRHRSContainer=new G4RotationMatrix();
	pRotRHRSContainer->rotateX(-270*deg);
	pRotRHRSContainer->rotateZ(fHRSAngle);  

	if(fSetupHRS>=2)
	{
		new G4PVPlacement(pRotLHRSContainer,G4ThreeVector(0,0,0.+fPivotZOffset),
			LHRSContainerLogical,"LHRSContainerPhys",motherLogical,0,0,0);
	}
	if(fSetupHRS>=2)
	{
		new G4PVPlacement(pRotRHRSContainer,G4ThreeVector(0,0,0.+fPivotZOffset),
			RHRSContainerLogical,"RHRSContainerPhys",motherLogical,0,0,0);
	}


	//////////////////////////
	//Q1 Entrance Window //
	//////////////////////////
	bool pSetupQ1EntranceWindow=true;
	//G2P and CREX have their own VB defined at a different location
	if( fSetupCREXGeometry)  pSetupQ1EntranceWindow=false;
	if(pSetupQ1EntranceWindow) 
	{

		//this part is trying to place a virtual boundary at the Q1 entrance

		//Place both left and right VB for HRS, which is pHRSContainerRin+4*mm away from the 
		//hall center(1.462m). This aperture is a round disk of 29.8 cm diameter and 2 mm thick
		//The real Q1 vacumn entrance to hall center is 1.312m, 

		double pHRSQ1WinThick = 2*mm;
		G4VSolid* HRSQ1WinSolid = new G4Tubs("HRSQ1WinTub",0.0,14.9*cm,
			pHRSQ1WinThick/2,0.0,360.0*deg);
		G4LogicalVolume* HRSQ1WinLogical = new G4LogicalVolume(HRSQ1WinSolid,
			0,"HRSQ1WinLogical",0,0,0);
		SDman->AddNewDetector(Q1WindowSD);
		HRSQ1WinLogical->SetSensitiveDetector(Q1WindowSD);

                G4VisAttributes *LightYellowVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.2)); 
		HRSQ1WinLogical->SetVisAttributes(LightYellowVisAtt); 

		//since the container has been rotated by 90 deg about x axis,
		//y'=z  z'=-y ==> I have to put this offset as -y
		double pHallCenter2Q1Win=pHRSContainerRin+4*mm;  //place it at the first 1.464 m

		if(fSetupHRS)
		{
			new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pHallCenter2Q1Win,0),
				HRSQ1WinLogical,"virtualBoundaryPhys_LHRSQ1Win",LHRSContainerLogical,0,0,0);
		}
		if(fSetupHRS)
		{
			new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pHallCenter2Q1Win,0),
				HRSQ1WinLogical,"virtualBoundaryPhys_RHRSQ1Win",RHRSContainerLogical,0,0,0);
		}
	}

	int    sos   = 1;
	if( fSnakeModel >= 52){
	  sos = 1;
	}

        double fringe_extension = 1.00;



	//BELOW ARE SNAKE VALUES, NOT NIM VALUES
	double pTarget =             0.0  * cm;
	double pQ1en   = pTarget + 160.0  * cm;
	double pQ1ex   = pQ1en   +  94.13 * cm;
	double pQ2en   = pQ1ex   + 115.58 * cm;
	double pQ2ex   = pQ2en   + 182.66 * cm;
	//ABOVE ARE SNAKE VALUES, NOT NIM VALUES

        double pQ1Length = sos ?  70. * cm : pQ1ex - pQ1en;


        double pHallCenter2LQ1Face = pQ1en;
        double pHallCenter2RQ1Face = pQ1en;


        double pQ2Rin=30.0*cm;
        double pQ2Rout=75.0*cm;
        double pQ2Length = pQ2ex - pQ2en;
        double pQ2LengthMag=pQ2Length * fringe_extension;

        double pHallCenter2LQ2Face = pQ2en;
        double  pHallCenter2RQ2Face = pQ2en;


        double pQ3Length = 182.68 * cm;
        double pQ3CenterY =  3.5853101 * m  + pQ3Length / sqrt(2.) / 2.;
        double pQ3CenterZ = 17.0257042 * m  + pQ3Length / sqrt(2.) / 2;
        double pQ3Rin=30.0*cm;
        double pQ3Rout=75.0*cm;

        double pLQ1Pos_Z_en=(pHallCenter2LQ1Face);//NIM
        double pLQ1Pos_Z_ex=(pHallCenter2LQ1Face + pQ1Length);//NIM
        double pLQ2Pos_Z_en=(pHallCenter2LQ2Face);//NIM
        double pLQ2Pos_Z_ex=(pHallCenter2LQ2Face + pQ2Length);//NIM
        double pRQ1Pos_Z_en=(pHallCenter2RQ1Face);//NIM
        double pRQ1Pos_Z_ex=(pHallCenter2RQ1Face + pQ1Length);//NIM
        double pRQ2Pos_Z_en=(pHallCenter2RQ2Face);//NIM
        double pRQ2Pos_Z_ex=(pHallCenter2RQ2Face + pQ2Length);//NIM
        double pLDPos_Z_en = pLQ2Pos_Z_ex + 4.42 * m;
        double pRDPos_Z_en = pRQ2Pos_Z_ex + 4.42 * m;
        double pLDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
        double pRDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
        double pLDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;
        double pRDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;


	/////////////////////////
	// HRS Q1              // Nickie is also adding the Q1 collimator here
	/////////////////////////
	//flag to trigger sos quad instead of old quad for the proper runs
	//K. Allada and B. Schmookler:
	//Magnetic Length = 70 cm
	//Radius to Pole Tip = 12.827  cm
	
	double q1shift = sos ? 0.0 * m : 0.0 * m;
	double pQ1Rin  = sos ? 12.827 * cm : 15.0*cm;
	double pQ1Rout = sos ? 35.0   * cm : 35.0*cm;//for now, keep the outer same for either
	double pQ1PosAct = pQ1ex - pQ1en;
	double pQ1LengthMag=pQ1Length * fringe_extension;
	//double pQ1LengthMag=94.1*cm;
	//double pQ1Length=(1.698M-1.36*m)*2+0.8*m;
        
        double pFPCenterZ , pFPCenterY ;
	double pVDCCenterZ, pVDCCenterY;
	double pQZ1CenterZ, pQZ1CenterY;
	double pQZ2CenterZ, pQZ2CenterY;
	//if( IS_NIM == 1 ){
	//double pFPR      = 8.4  * m;//radius of curvature of dipole
	//double pFPA      = ( IS_NIM == 1 ) ? 9.96 * m : pDipoleRCenterZ;//distance from pivot to entrance of dipole
	//double pFPH      = pFPR * tan ( 45. / 2. * deg );
	//pFPCenterZ=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) );
	//pFPCenterY=(                 pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2);
	//}else{
	pVDCCenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pVDCCenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pQZ1CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ1CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ2CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pQZ2CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pFPCenterZ  = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	pFPCenterY  = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	//}
	
	G4double vb_thickness = 0.5 * mm;
	G4VSolid* FPSolid = new G4Tubs("FPTub",0,pQ3Rout * 2,vb_thickness,0.0,360.0*deg);
//	G4VSolid* PlaneSolid1 = new G4Tubs("PlaneTub",0,pQ1Rout,vb_thickness,0.0,360.0*deg); //circles
//	G4VSolid* PlaneSolid2 = new G4Tubs("PlaneTub",0,pQ2Rout,vb_thickness,0.0,360.0*deg); //circles
	G4VSolid* PlaneSolid1 = new G4Tubs("PlaneTub",0,200.*cm,vb_thickness,0.0,360.0*deg); //circles
	G4VSolid* PlaneSolid2 = new G4Tubs("PlaneTub",0,200.*cm,vb_thickness,0.0,360.0*deg); //circles

	double X_S0 = 29.5 * cm;
	double Y_S0 = 35.5 * cm;
	double Z_S0 = vb_thickness;

	double X_S1 = 29.5 * cm;
	double Y_S1 = 35.5 * cm;
	double Z_S1 = vb_thickness;

	double X_Q1 = 29.5 * cm;
	double Y_Q1 = 35.5 * cm;
	double Z_Q1 = vb_thickness;

	double X_Q2 = 29.5 * cm;
	double Y_Q2 = 35.5 * cm;
	double Z_Q2 = vb_thickness;

	G4VSolid* Scint0  = new G4Box("Scint0" , X_S0 / 2.0, Y_S0 / 2.0, Z_S0 / 2.0);
	G4VSolid* Scint1  = new G4Box("Scint1" , X_S1 / 2.0, Y_S1 / 2.0, Z_S1 / 2.0);
	G4VSolid* Quartz1 = new G4Box("Quartz1", X_Q1 / 2.0, Y_Q1 / 2.0, Z_Q1 / 2.0);
	G4VSolid* Quartz2 = new G4Box("Quartz2", X_Q2 / 2.0, Y_Q2 / 2.0, Z_Q2 / 2.0);

	G4VSolid* magneticSolid = new G4Box("magneticBox",fFieldX/2.0,fFieldY/2.0,fFieldZ/2.0);


	G4LogicalVolume* LPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		0,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		0,"RPlaneLogical1",0,0,0);
	G4LogicalVolume* LPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		0,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		0,"RPlaneLogical1",0,0,0);

        G4LogicalVolume* LFPLogical = new G4LogicalVolume(FPSolid,
                0,"LFPLogical",0,0,0);
        G4LogicalVolume* RFPLogical = new G4LogicalVolume(FPSolid,
                0,"RFPLogical",0,0,0);


	LPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	LPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotVDCInContainer=new G4RotationMatrix();
	pRotVDCInContainer->rotateX(0.*deg); 
	G4RotationMatrix *pRotFPInContainer=new G4RotationMatrix();
	pRotFPInContainer->rotateX(-45*deg); 

	
        double pHallCenter2Col = 1.38 * m;
	double pPaulColT        = 0.01 * m;
	double pPaulX = ( - pHallCenter2Col - pPaulColT * 2. ) * cos(fHRSAngle) ;
	double pPaulY = ( - pHallCenter2Col - pPaulColT * 2. ) * sin(fHRSAngle);
	if(fSnakeModel == 49 || fSnakeModel == 48 || fSnakeModel > 50 ){
	  //double pLQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/1.0);//SNAKE
	  new G4PVPlacement(pRotXLHRSdeg,G4ThreeVector(pPaulY,0,-pPaulX),
			    LPlaneLogical1,"virtualBoundaryPhys_col_LHRS",motherLogical,0,0,0);//col
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_en,0),
			    LPlaneLogical1,"virtualBoundaryPhys_q1en_LHRS",LHRSContainerLogical,0,0,0);//q1en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_ex,0),
			    LPlaneLogical1,"virtualBoundaryPhys_q1ex_LHRS",LHRSContainerLogical,0,0,0);//q1ex
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_en,0),
			    LPlaneLogical2,"virtualBoundaryPhys_q2en_LHRS",LHRSContainerLogical,0,0,0);//q2en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_ex,0),
			    LPlaneLogical2,"virtualBoundaryPhys_q2ex_LHRS",LHRSContainerLogical,0,0,0);//q2ex
	  new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pLDPos_Z_en,0),
	  LPlaneLogical2,"virtualBoundaryPhys_den_LHRS",LHRSContainerLogical,0,0,0);//den
	  new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pLDPos_Z_ex,pLDPos_X_ex),
	  LPlaneLogical2,"virtualBoundaryPhys_dex_LHRS",LHRSContainerLogical,0,0,0);//dex
	  
	  
	  new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ + pQ3Length / sqrt(2) / 2,pQ3CenterY - pQ3Length / sqrt(2) / 2),
	  LPlaneLogical2,"virtualBoundaryPhys_q3en_LHRS",LHRSContainerLogical,0,0,0);//q3en
          new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ - pQ3Length / sqrt(2) / 2,pQ3CenterY + pQ3Length / sqrt(2) / 2),
                  LPlaneLogical2,"virtualBoundaryPhys_q3ex_LHRS",LHRSContainerLogical,0,0,0);//q3ex
          new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pVDCCenterZ,pVDCCenterY),
                  LFPLogical,"virtualBoundaryPhys_vdc_LHRS",LHRSContainerLogical,0,0,0);
          new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pFPCenterZ,pFPCenterY),
                  LFPLogical,"virtualBoundaryPhys_fp_LHRS",LHRSContainerLogical,0,0,0);
          new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pQZ1CenterZ,pQZ1CenterY),
                  LFPLogical,"virtualBoundaryPhys_qz1_LHRS",LHRSContainerLogical,0,0,0);
          new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pQZ2CenterZ,pQZ2CenterY),
                  LFPLogical,"virtualBoundaryPhys_qz2_LHRS",LHRSContainerLogical,0,0,0);


	}
	if(fSnakeModel == 49 || fSnakeModel == 48 || fSnakeModel > 50 ){
	  new G4PVPlacement(pRotXRHRSdeg,G4ThreeVector(-pPaulY, 0, -pPaulX),
			    LPlaneLogical1,"virtualBoundaryPhys_col_RHRS",motherLogical,0,0,0);//q1en

	  //double pRQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/1.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z_en,0),
			    RPlaneLogical1,"virtualBoundaryPhys_q1en_RHRS",RHRSContainerLogical,0,0,0);//q1en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z_ex,0),
			    RPlaneLogical1,"virtualBoundaryPhys_q1ex_RHRS",RHRSContainerLogical,0,0,0);//q1ex
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_en,0),
			    RPlaneLogical2,"virtualBoundaryPhys_q2en_RHRS",RHRSContainerLogical,0,0,0);//q2en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_ex,0),
			    RPlaneLogical2,"virtualBoundaryPhys_q2ex_RHRS",RHRSContainerLogical,0,0,0);//q2ex
	  new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pRDPos_Z_en,0),
	  RPlaneLogical2,"virtualBoundaryPhys_den_RHRS",RHRSContainerLogical,0,0,0);//den
	  new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pRDPos_Z_ex,pRDPos_X_ex),
	  RPlaneLogical2,"virtualBoundaryPhys_dex_RHRS",RHRSContainerLogical,0,0,0);//dex
	  
	  
	  new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ + pQ3Length / sqrt(2) / 2,pQ3CenterY - pQ3Length / sqrt(2) / 2),
	  RPlaneLogical2,"virtualBoundaryPhys_q3en_RHRS",RHRSContainerLogical,0,0,0);//q3en
	  new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ - pQ3Length / sqrt(2) / 2,pQ3CenterY + pQ3Length / sqrt(2) / 2),
	  RPlaneLogical2,"virtualBoundaryPhys_q3ex_RHRS",RHRSContainerLogical,0,0,0);//q3ex
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pVDCCenterZ,pVDCCenterY),
			    RFPLogical,"virtualBoundaryPhys_vdc_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pFPCenterZ,pFPCenterY),
			    RFPLogical,"virtualBoundaryPhys_fp_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pQZ1CenterZ,pQZ1CenterY),
			    RFPLogical,"virtualBoundaryPhys_qz1_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pQZ2CenterZ,pQZ2CenterY),
			    RFPLogical,"virtualBoundaryPhys_qz2_RHRS",RHRSContainerLogical,0,0,0);

	}
	
	//#endif
	//////////////////////////////////////////////////////////

	return;

}

void g4hrsParallelWorld::ConstructSD(G4LogicalVolume* motherLogical) {

    // Add detectors here
/*	
	G4SDManager *SDMan = G4SDManager::GetSDMpointer();
	G4VSensitiveDetector* virtual45 = new g4hrsGenericDetector("virtual45",45);
	SDMan->AddNewDetector(virtual45);

	G4VSolid* parallelDetSolid = new G4Box("parallelDetSolid", 500., 500., 0.05);
        G4LogicalVolume* parallelDetLogical = new G4LogicalVolume(parallelDetSolid,NULL,"parallelDetLogical",0,0,0);		
	parallelDetLogical->SetSensitiveDetector(virtual45);
	G4VPhysicalVolume* parallelDetPhys = new G4PVPlacement(0,G4ThreeVector(fTargetX,fTargetY,fTargetZ + 69.99937 - 74./2. + fPivotZOffset - 10.*cm),
		parallelDetLogical,"parallelDetPhys",motherLogical,0,0);
*/	
	///////////////////////////////////
	// Solids for virtual boundaries //
	///////////////////////////////////
	G4double vb_thickness = 0.5 * mm;
	G4VSolid* FPSolid = new G4Tubs("FPTub",0,75.*cm * 2,vb_thickness,0.0,360.0*deg);
	G4VSolid* PlaneSolid1 = new G4Tubs("PlaneTub",0,35.*cm,vb_thickness,0.0,360.0*deg); //circles
	G4VSolid* PlaneSolid2 = new G4Tubs("PlaneTub",0,75.*cm,vb_thickness,0.0,360.0*deg); //circles

	G4LogicalVolume* LPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,0,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,0,"RPlaneLogical1",0,0,0);
	G4LogicalVolume* LPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,0,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,0,"RPlaneLogical1",0,0,0);

        G4LogicalVolume* LFPLogical = new G4LogicalVolume(FPSolid,0,"LFPLogical",0,0,0);
        G4LogicalVolume* RFPLogical = new G4LogicalVolume(FPSolid,0,"RFPLogical",0,0,0);

	G4VisAttributes* virtualBoundaryVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));	

	LPlaneLogical1->SetVisAttributes(virtualBoundaryVisAtt);	
	LPlaneLogical2->SetVisAttributes(virtualBoundaryVisAtt);
	RPlaneLogical1->SetVisAttributes(virtualBoundaryVisAtt);	
	RPlaneLogical2->SetVisAttributes(virtualBoundaryVisAtt);



	///////////////////////////////
	// Septum virtual boundaries //
	///////////////////////////////

	// Septum virtual boundaries will be placed directly in the (parallel) world volume
	
	double pSeptumX      = 140.0  * cm;
	double pSeptumY      = 84.4   * cm;
	double pSeptumZ      = 74.0   * cm;
	//double pSeptumPlaceZ = 70.414 * cm;
	double pSeptumPlaceZ = 69.99937 * cm;
	new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ - 0.5 * pSeptumZ + 2 * vb_thickness + fPivotZOffset),
		LPlaneLogical2,"virtualBoundaryPhys_sen",motherLogical,0,0);//sen
	new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ + fPivotZOffset),
		LPlaneLogical2,"virtualBoundaryPhys_sm",motherLogical,0,0);//sm
	new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ + 0.5 * pSeptumZ + fPivotZOffset),
	    	LPlaneLogical2,"virtualBoundaryPhys_sex",motherLogical,0,0);//sex
	new G4PVPlacement(0,G4ThreeVector(0,0,36. * cm + fPivotZOffset),
	    	LPlaneLogical2,"virtualBoundaryPhys_coil",motherLogical,0,0);//coil
	new G4PVPlacement(0,G4ThreeVector(0,0,-50. * cm + fPivotZOffset),
	    	LPlaneLogical2,"virtualBoundaryPhys_mid",motherLogical,0,0);//mid

	
	////////////////////////////
	// HRS virtual boundaries //
	////////////////////////////
	
	// LHRS (RHRS) virutal boundaries will be placed in a logical LHRS container (RHRS container) which will be placed in the world volume
	// Note that the HRS containers are created along -y axis, rotated to along the z axis, then rotated +/- HRS angle
	// Therefore the placement of virtual boundaries in the HRS containers will not be with respect to the hall axes 


	////////////////////
	// HRS containers //
	////////////////////

	double pHRSContainerRin=1.46*m,pHRSContainerRout=25*m;//NickieMode it just has to be big enough to fit my VB at FP
	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {	pHRSContainerRin,
						pHRSContainerRin,
						2.5*m,
						3.7*m,
						9.0*m,
						pHRSContainerRout-3.0*m,
						pHRSContainerRout};
	
	double rOuter_HRSContainer[] = {	pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout};
	
	double zPlane_HRSContainer[] = {	-2.0*m,
						1.0*m,
						1.0*m,
						2.0*m,
						2.0*m,
						11.0*m,
						11.0*m};
	
	
	G4Polycone* HRSContainerSolid = new G4Polycone("HRSContainer",258.0*deg,24.0*deg,
		kNPlane_HRSContainer,zPlane_HRSContainer,rInner_HRSContainer,rOuter_HRSContainer);

	G4LogicalVolume* LHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		0,"LHRSContainerLogical",0,0);
	G4LogicalVolume* RHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		0,"RHRSContainerLogical",0,0);

	
	G4VisAttributes *HallVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
	HallVisAtt->SetVisibility(false);

	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt);
	
	G4RotationMatrix *pRotLHRSContainer=new G4RotationMatrix();
	pRotLHRSContainer->rotateX(-270*deg);
	pRotLHRSContainer->rotateZ(-fHRSAngle);  

	G4RotationMatrix *pRotRHRSContainer=new G4RotationMatrix();
	pRotRHRSContainer->rotateX(-270*deg);
	pRotRHRSContainer->rotateZ(fHRSAngle);  

	new G4PVPlacement(pRotLHRSContainer,G4ThreeVector(0.,0.,0.+fPivotZOffset),
		LHRSContainerLogical,"LHRSContainerPhys",motherLogical,0,0,0);
	new G4PVPlacement(pRotRHRSContainer,G4ThreeVector(0,0,0.+fPivotZOffset),
		RHRSContainerLogical,"RHRSContainerPhys",motherLogical,0,0,0);

	
	G4RotationMatrix* LHRSRot = new G4RotationMatrix();
	LHRSRot->rotateY(fHRSAngle); 
	G4RotationMatrix* RHRSRot = new G4RotationMatrix();
	RHRSRot->rotateY(-fHRSAngle); 
	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	G4RotationMatrix *pRotX45deg=new G4RotationMatrix();
	pRotX45deg->rotateX(-45*deg); 
	G4RotationMatrix *pRotX30deg=new G4RotationMatrix();
	pRotX30deg->rotateX(60*deg); 
	G4RotationMatrix *pRotX105deg=new G4RotationMatrix();
	pRotX105deg->rotateX(165*deg); 
	G4RotationMatrix *pRotXLHRSdeg=new G4RotationMatrix();
	pRotXLHRSdeg->rotateY(fHRSAngle); 
	G4RotationMatrix *pRotXRHRSdeg=new G4RotationMatrix();
	pRotXRHRSdeg->rotateY(-fHRSAngle); 


	//BELOW ARE SNAKE VALUES, NOT NIM VALUES
	double pTarget = fTargetZ;
	double pQ1en   = pTarget + 160.0  * cm;
	double pQ1ex   = pQ1en   +  94.13 * cm;
	double pQ2en   = pQ1ex   + 115.58 * cm;
	double pQ2ex   = pQ2en   + 182.66 * cm;
	//ABOVE ARE SNAKE VALUES, NOT NIM VALUES


	///////////////////////////
	// Q1 virtual boundaries //
	///////////////////////////
	double pQ1Length = pQ1ex - pQ1en;	
        double pHallCenter2LQ1Face = pQ1en;
        double pLQ1Pos_Z_en=(pHallCenter2LQ1Face);//NIM
        double pLQ1Pos_Z_ex=(pHallCenter2LQ1Face + pQ1Length);//NIM
	//LHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_en,0),
		LPlaneLogical1,"virtualBoundaryPhys_q1en_LHRS",LHRSContainerLogical,0,0,0);//q1en	
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_ex,0),
		LPlaneLogical1,"virtualBoundaryPhys_q1ex_LHRS",LHRSContainerLogical,0,0,0);//q1ex
	//RHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_en,0),
		LPlaneLogical1,"virtualBoundaryPhys_q1en_RHRS",RHRSContainerLogical,0,0,0);//q1en	
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_ex,0),
		LPlaneLogical1,"virtualBoundaryPhys_q1ex_RHRS",RHRSContainerLogical,0,0,0);//q1ex
	 


	///////////////////////////
	// Q2 virtual boundaries //
	///////////////////////////
        double pHallCenter2LQ2Face = pQ2en;
        double  pHallCenter2RQ2Face = pQ2en;
        double pQ2Length = pQ2ex - pQ2en;
        double pLQ2Pos_Z_en=(pHallCenter2LQ2Face);//NIM
        double pLQ2Pos_Z_ex=(pHallCenter2LQ2Face + pQ2Length);//NIM
        double pRQ2Pos_Z_en=(pHallCenter2RQ2Face);//NIM
        double pRQ2Pos_Z_ex=(pHallCenter2RQ2Face + pQ2Length);//NIM
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_en,0),		    
		LPlaneLogical2,"virtualBoundaryPhys_q2en_LHRS",LHRSContainerLogical,0,0,0);//q2en
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_ex,0),
		LPlaneLogical2,"virtualBoundaryPhys_q2ex_LHRS",LHRSContainerLogical,0,0,0);//q2ex
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_en,0),
		RPlaneLogical2,"virtualBoundaryPhys_q2en_RHRS",RHRSContainerLogical,0,0,0);//q2en
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_ex,0),
		RPlaneLogical2,"virtualBoundaryPhys_q2ex_RHRS",RHRSContainerLogical,0,0,0);//q2ex
	 


	//////////////////////////
	// D virtual boundaries //
	//////////////////////////
	double pLDPos_Z_en = pLQ2Pos_Z_ex + 4.42 * m;
	double pRDPos_Z_en = pRQ2Pos_Z_ex + 4.42 * m;
	double pLDPos_X_ex = 2.4603032 * m;
	double pRDPos_X_ex = 2.4603032 * m;
	double pLDPos_Z_ex = -15.9006973 * m;        
	double pRDPos_Z_ex = -15.9006973 * m;
	new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pLDPos_Z_en,0),
		LPlaneLogical2,"virtualBoundaryPhys_den_LHRS",LHRSContainerLogical,0,0,0);//den
	new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pLDPos_Z_ex,pLDPos_X_ex),
		LPlaneLogical2,"virtualBoundaryPhys_dex_LHRS",LHRSContainerLogical,0,0,0);//dex
	new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pRDPos_Z_en,0),
		RPlaneLogical2,"virtualBoundaryPhys_den_RHRS",RHRSContainerLogical,0,0,0);//den
	new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pRDPos_Z_ex,pRDPos_X_ex),
		RPlaneLogical2,"virtualBoundaryPhys_dex_RHRS",RHRSContainerLogical,0,0,0);//dex
	  
  




	return;
}




















