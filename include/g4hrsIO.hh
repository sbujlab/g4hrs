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

	int numVB;
	int numVar;

	G4double VBdata[14][12];

	char VBName[14][6];
	char VarName[12][5];

	g4hrsSteppingAction* fSteppingAction;

	G4int fLHRS;
	G4int fRHRS;

	G4double fX0;
	G4double fY0;
	G4double fZ0;
	G4double fP0;
	G4double fTh0;
	G4double fPh0;

	G4double fX0_tr;
	G4double fY0_tr;
	G4double fZ0_tr;
	G4double fP0_tr;
	G4double fTh0_tr;
	G4double fPh0_tr;
	
	G4double fX_sen;
	G4double fY_sen;
	G4double fZ_sen;
	G4double fTh_sen;
	G4double fPh_sen;

	G4double fX_sen_tr;
	G4double fY_sen_tr;
	G4double fZ_sen_tr;
	G4double fTh_sen_tr;
	G4double fPh_sen_tr;

	G4double fX_sen_tf;
	G4double fY_sen_tf;
	G4double fZ_sen_tf;
	G4double fTh_sen_tf;
	G4double fPh_sen_tf;
	
	G4double fX_sm;
	G4double fY_sm;
	G4double fZ_sm;
	G4double fTh_sm;
	G4double fPh_sm;

	G4double fX_sm_tr;
	G4double fY_sm_tr;
	G4double fZ_sm_tr;
	G4double fTh_sm_tr;
	G4double fPh_sm_tr;

	G4double fX_sm_tf;
	G4double fY_sm_tf;
	G4double fZ_sm_tf;
	G4double fTh_sm_tf;
	G4double fPh_sm_tf;
	
	G4double fX_sex;
	G4double fY_sex;
	G4double fZ_sex;
	G4double fTh_sex;
	G4double fPh_sex;

	G4double fX_sex_tr;
	G4double fY_sex_tr;
	G4double fZ_sex_tr;
	G4double fTh_sex_tr;
	G4double fPh_sex_tr;

	G4double fX_sex_tf;
	G4double fY_sex_tf;
	G4double fZ_sex_tf;
	G4double fTh_sex_tf;
	G4double fPh_sex_tf;
	
	G4double fX_col;
	G4double fY_col;
	G4double fZ_col;
	G4double fTh_col;
	G4double fPh_col;

	G4double fX_col_tr;
	G4double fY_col_tr;
	G4double fZ_col_tr;
	G4double fTh_col_tr;
	G4double fPh_col_tr;

	G4double fX_col_tf;
	G4double fY_col_tf;
	G4double fZ_col_tf;
	G4double fTh_col_tf;
	G4double fPh_col_tf;
	
	G4double fX_q1en;
	G4double fY_q1en;
	G4double fZ_q1en;
	G4double fTh_q1en;
	G4double fPh_q1en;

	G4double fX_q1en_tr;
	G4double fY_q1en_tr;
	G4double fZ_q1en_tr;
	G4double fTh_q1en_tr;
	G4double fPh_q1en_tr;

	G4double fX_q1ex;
	G4double fY_q1ex;
	G4double fZ_q1ex;
	G4double fTh_q1ex;
	G4double fPh_q1ex;

	G4double fX_q1ex_tr;
	G4double fY_q1ex_tr;
	G4double fZ_q1ex_tr;
	G4double fTh_q1ex_tr;
	G4double fPh_q1ex_tr;

	G4double fX_q1ex_tf;
	G4double fY_q1ex_tf;
	G4double fZ_q1ex_tf;
	G4double fTh_q1ex_tf;
	G4double fPh_q1ex_tf;
	
	G4double fX_q2en;
	G4double fY_q2en;
	G4double fZ_q2en;
	G4double fTh_q2en;
	G4double fPh_q2en;

	G4double fX_q2en_tr;
	G4double fY_q2en_tr;
	G4double fZ_q2en_tr;
	G4double fTh_q2en_tr;
	G4double fPh_q2en_tr;

	G4double fX_q2ex;
	G4double fY_q2ex;
	G4double fZ_q2ex;
	G4double fTh_q2ex;
	G4double fPh_q2ex;

	G4double fX_q2ex_tr;
	G4double fY_q2ex_tr;
	G4double fZ_q2ex_tr;
	G4double fTh_q2ex_tr;
	G4double fPh_q2ex_tr;

	G4double fX_q2ex_tf;
	G4double fY_q2ex_tf;
	G4double fZ_q2ex_tf;
	G4double fTh_q2ex_tf;
	G4double fPh_q2ex_tf;
	
	G4double fX_den;
	G4double fY_den;
	G4double fZ_den;
	G4double fTh_den;
	G4double fPh_den;

	G4double fX_den_tr;
	G4double fY_den_tr;
	G4double fZ_den_tr;
	G4double fTh_den_tr;
	G4double fPh_den_tr;

	G4double fX_den_tf;
	G4double fY_den_tf;
	G4double fZ_den_tf;
	G4double fTh_den_tf;
	G4double fPh_den_tf;
	
	G4double fX_dex;
	G4double fY_dex;
	G4double fZ_dex;
	G4double fTh_dex;
	G4double fPh_dex;

	G4double fX_dex_tr;
	G4double fY_dex_tr;
	G4double fZ_dex_tr;
	G4double fTh_dex_tr;
	G4double fPh_dex_tr;

	G4double fX_dex_tf;
	G4double fY_dex_tf;
	G4double fZ_dex_tf;
	G4double fTh_dex_tf;
	G4double fPh_dex_tf;
	
	G4double fX_q3en;
	G4double fY_q3en;
	G4double fZ_q3en;
	G4double fTh_q3en;
	G4double fPh_q3en;

	G4double fX_q3en_tr;
	G4double fY_q3en_tr;
	G4double fZ_q3en_tr;
	G4double fTh_q3en_tr;
	G4double fPh_q3en_tr;

	G4double fX_q3en_tf;
	G4double fY_q3en_tf;
	G4double fZ_q3en_tf;
	G4double fTh_q3en_tf;
	G4double fPh_q3en_tf;
	
	G4double fX_q3ex;
	G4double fY_q3ex;
	G4double fZ_q3ex;
	G4double fTh_q3ex;
	G4double fPh_q3ex;

	G4double fX_q3ex_tr;
	G4double fY_q3ex_tr;
	G4double fZ_q3ex_tr;
	G4double fTh_q3ex_tr;
	G4double fPh_q3ex_tr;

	G4double fX_q3ex_tf;
	G4double fY_q3ex_tf;
	G4double fZ_q3ex_tf;
	G4double fTh_q3ex_tf;
	G4double fPh_q3ex_tf;
		
	G4double fX_vdc;
	G4double fY_vdc;
	G4double fZ_vdc;
	G4double fTh_vdc;
	G4double fPh_vdc;

	G4double fX_vdc_tr;
	G4double fY_vdc_tr;
	G4double fZ_vdc_tr;
	G4double fTh_vdc_tr;
	G4double fPh_vdc_tr;

	G4double fX_fp;
	G4double fY_fp;
	G4double fZ_fp;
	G4double fTh_fp;
	G4double fPh_fp;

	G4double fX_fp_tr;
	G4double fY_fp_tr;
	G4double fZ_fp_tr;
	G4double fTh_fp_tr;
	G4double fPh_fp_tr;

	G4double fX_fp_tf;
	G4double fY_fp_tf;
	G4double fZ_fp_tf;
	G4double fTh_fp_tf;
	G4double fPh_fp_tf;


};

#endif//g4hrsIO_HH
