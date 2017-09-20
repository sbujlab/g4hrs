#ifndef g4hrsIO_HH
#define g4hrsIO_HH

#include "TROOT.h"
#include "TObject.h"
#include "G4Run.hh"
#include "g4hrstypes.hh"

#include "G4String.hh"

class TFile;
class TTree;

class g4hrsGenericDetectorHit;
class g4hrsGenericDetectorSum;
class g4hrsEvent;
class g4hrsSteppingAction;

#include <xercesc/dom/DOMElement.hpp>


#define __IO_MAXHIT 10000
#define __FILENAMELEN 255

// Units for output
#define __E_UNIT GeV
#define __L_UNIT m
#define __T_UNIT ns
#define __ANG_UNIT rad
#define __ASYMM_SCALE 1e-9 // ppb

class g4hrsIO {
    public:
	 g4hrsIO();
	~g4hrsIO();

	void SetFilename(G4String  fn);
	G4String GetFilename(){return fFilename;}

	void FillTree();
	void Flush();
	void WriteTree();

	void WriteRun();

	void InitializeTree();

	void GrabGDMLFiles( G4String fn );

    private:
	TFile *fFile;
	TTree *fTree;

	char fFilename[__FILENAMELEN];

	std::vector<G4String>       fGDMLFileNames;

	void SearchGDMLforFiles(G4String );
	void TraverseChildren(  xercesc::DOMElement * );

	//  Interfaces and buffers to the tree
	//  This is potentially going to get very long...

	// Event data
    public:
	void SetEventData(g4hrsEvent *);
    private:
	Int_t fNEvPart;

	Double_t fEvRate;
	Double_t fEvEffXS;
	Double_t fEvAsym;
	Double_t fEvmAsym;
	Double_t fEvBeamP;
	Double_t fEvQ2;
	Double_t fEvW2;
	Double_t fEvThCoM;

	Double_t fBmX;
	Double_t fBmY;
	Double_t fBmZ;
	Double_t fBmdX;
	Double_t fBmdY;
	Double_t fBmdZ;
	Double_t fBmth;
	Double_t fBmph;

	Int_t fEvPID[__IO_MAXHIT];

	Double_t fEvPart_X[__IO_MAXHIT];
	Double_t fEvPart_Y[__IO_MAXHIT];
	Double_t fEvPart_Z[__IO_MAXHIT];
	Double_t fEvPart_Px[__IO_MAXHIT];
	Double_t fEvPart_Py[__IO_MAXHIT];
	Double_t fEvPart_Pz[__IO_MAXHIT];
	Double_t fEvPart_Th[__IO_MAXHIT];
	Double_t fEvPart_Ph[__IO_MAXHIT];
	Double_t fEvPart_P[__IO_MAXHIT];
	Double_t fEvPart_tPx[__IO_MAXHIT];
	Double_t fEvPart_tPy[__IO_MAXHIT];
	Double_t fEvPart_tPz[__IO_MAXHIT];


	//  GenericDetectorHit
    public:
	void AddGenericDetectorHit(g4hrsGenericDetectorHit *);
    private:
	Int_t fNGenDetHit;
	Int_t fGenDetHit_det[__IO_MAXHIT];
	Int_t fGenDetHit_id[__IO_MAXHIT];

	Int_t fGenDetHit_trid[__IO_MAXHIT];
	Int_t fGenDetHit_pid[__IO_MAXHIT];
	Int_t fGenDetHit_gen[__IO_MAXHIT];
	Int_t fGenDetHit_mtrid[__IO_MAXHIT];

	Double_t fGenDetHit_X[__IO_MAXHIT];
	Double_t fGenDetHit_Y[__IO_MAXHIT];
	Double_t fGenDetHit_Z[__IO_MAXHIT];
	Double_t fGenDetHit_R[__IO_MAXHIT];
	Double_t fGenDetHit_Ph[__IO_MAXHIT];

	Double_t fGenDetHit_Px[__IO_MAXHIT];
	Double_t fGenDetHit_Py[__IO_MAXHIT];
	Double_t fGenDetHit_Pz[__IO_MAXHIT];
	Double_t fGenDetHit_P[__IO_MAXHIT];
	Double_t fGenDetHit_E[__IO_MAXHIT];
	Double_t fGenDetHit_M[__IO_MAXHIT];

	Double_t fGenDetHit_Vx[__IO_MAXHIT];
	Double_t fGenDetHit_Vy[__IO_MAXHIT];
	Double_t fGenDetHit_Vz[__IO_MAXHIT];

	Int_t fCollCut;

	//  GenericDetectorSum
    public:
	void AddGenericDetectorSum(g4hrsGenericDetectorSum *);
    private:
	Int_t fNGenDetSum;
	Int_t fGenDetSum_det[__IO_MAXHIT];
	Int_t fGenDetSum_id[__IO_MAXHIT];
	Double_t fGenDetSum_edep[__IO_MAXHIT];

	// Virtual boundary data 
    public:
	void SetSteppingAction(g4hrsSteppingAction *stepper) { fSteppingAction = stepper;}
	void ClearVirtualBoundaryData();
	void SetVirtualBoundaryData();
    private:
	g4hrsSteppingAction* fSteppingAction;

	G4double fX_sen;
	G4double fY_sen;
	G4double fZ_sen;
	G4double fP_sen;
	G4double fTheta_sen;
	G4double fPhi_sen;

	G4double fX_sm;
	G4double fY_sm;
	G4double fZ_sm;
	G4double fP_sm;
	G4double fTheta_sm;
	G4double fPhi_sm;

	G4double fX_sex;
	G4double fY_sex;
	G4double fZ_sex;
	G4double fP_sex;
	G4double fTheta_sex;
	G4double fPhi_sex;
	
	G4double fX_coil;
	G4double fY_coil;
	G4double fZ_coil;
	G4double fP_coil;
	G4double fTheta_coil;
	G4double fPhi_coil;

	G4double fX_mid;
	G4double fY_mid;
	G4double fZ_mid;
	G4double fP_mid;
	G4double fTheta_mid;
	G4double fPhi_mid;

	G4double fX_col;
	G4double fY_col;
	G4double fZ_col;
	G4double fP_col;
	G4double fTheta_col;
	G4double fPhi_col;
	
	G4double fX_q1en_L;
	G4double fY_q1en_L;
	G4double fZ_q1en_L;
	G4double fP_q1en_L;
	G4double fTheta_q1en_L;
	G4double fPhi_q1en_L;

	G4double fX_q1ex_L;
	G4double fY_q1ex_L;
	G4double fZ_q1ex_L;
	G4double fP_q1ex_L;
	G4double fTheta_q1ex_L;
	G4double fPhi_q1ex_L;
	
	G4double fX_q1en_R;
	G4double fY_q1en_R;
	G4double fZ_q1en_R;
	G4double fP_q1en_R;
	G4double fTheta_q1en_R;
	G4double fPhi_q1en_R;

	G4double fX_q1ex_R;
	G4double fY_q1ex_R;
	G4double fZ_q1ex_R;
	G4double fP_q1ex_R;
	G4double fTheta_q1ex_R;
	G4double fPhi_q1ex_R;

	G4double fX_q2en_L;
	G4double fY_q2en_L;
	G4double fZ_q2en_L;
	G4double fP_q2en_L;
	G4double fTheta_q2en_L;
	G4double fPhi_q2en_L;

	G4double fX_q2ex_L;
	G4double fY_q2ex_L;
	G4double fZ_q2ex_L;
	G4double fP_q2ex_L;
	G4double fTheta_q2ex_L;
	G4double fPhi_q2ex_L;
	
	G4double fX_q2en_R;
	G4double fY_q2en_R;
	G4double fZ_q2en_R;
	G4double fP_q2en_R;
	G4double fTheta_q2en_R;
	G4double fPhi_q2en_R;

	G4double fX_q2ex_R;
	G4double fY_q2ex_R;
	G4double fZ_q2ex_R;
	G4double fP_q2ex_R;
	G4double fTheta_q2ex_R;
	G4double fPhi_q2ex_R;

	G4double fX_den;
	G4double fY_den;
	G4double fZ_den;
	G4double fP_den;
	G4double fTheta_den;
	G4double fPhi_den;

	G4double fX_dex;
	G4double fY_dex;
	G4double fZ_dex;
	G4double fP_dex;
	G4double fTheta_dex;
	G4double fPhi_dex;

	G4double fX_q3en;
	G4double fY_q3en;
	G4double fZ_q3en;
	G4double fP_q3en;
	G4double fTheta_q3en;
	G4double fPhi_q3en;

	G4double fX_q3ex;
	G4double fY_q3ex;
	G4double fZ_q3ex;
	G4double fP_q3ex;
	G4double fTheta_q3ex;
	G4double fPhi_q3ex;

	G4double fX_vdc;
	G4double fY_vdc;
	G4double fZ_vdc;
	G4double fP_vdc;
	G4double fTheta_vdc;
	G4double fPhi_vdc;

	G4double fX_qz1;
	G4double fY_qz1;
	G4double fZ_qz1;
	G4double fP_qz1;
	G4double fTheta_qz1;
	G4double fPhi_qz1;

	G4double fX_qz2;
	G4double fY_qz2;
	G4double fZ_qz2;
	G4double fP_qz2;
	G4double fTheta_qz2;
	G4double fPhi_qz2;

	G4double fX_fp;
	G4double fY_fp;
	G4double fZ_fp;
	G4double fP_fp;
	G4double fTheta_fp;
	G4double fPhi_fp;
	
	G4double fX_par;
	G4double fY_par;
	G4double fZ_par;
	G4double fP_par;
	G4double fTheta_par;
	G4double fPhi_par;




};

#endif//g4hrsIO_HH
