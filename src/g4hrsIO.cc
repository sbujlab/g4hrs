#include "g4hrsIO.hh"

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "G4ParticleDefinition.hh"

#include "g4hrsGenericDetectorHit.hh"
#include "g4hrsGenericDetectorSum.hh"
#include "g4hrsEvent.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"
#include "g4hrsBeamTarget.hh"
#include "g4hrsSteppingAction.hh"

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMNode.hpp>


g4hrsIO::g4hrsIO(){
    
    fTree = NULL;
    InitializeTree();
    // Default filename
    strcpy(fFilename, "g4hrsout.root");

    fFile = NULL;
}

g4hrsIO::~g4hrsIO(){
    if( fTree ){ delete fTree; }
    fTree = NULL;
    if( fFile ){ delete fFile; }
    fFile = NULL;
}

void g4hrsIO::SetFilename(G4String fn){
    G4cout << "Setting output file to " << fn << G4endl;
    strcpy(fFilename, fn.data());
}

void g4hrsIO::InitializeTree(){
    if( fFile ){
	fFile->Close();
	delete fFile;
    }

    fFile = new TFile(fFilename, "RECREATE");

    if( fTree ){ delete fTree; }

    fTree = new TTree("T", "Geant4 Moller Simulation");

    fTree->SetMaxTreeSize(1900000000); // 1.9GB

    // Event information
    fTree->Branch("rate",     &fEvRate,   "rate/D");
    fTree->Branch("ev.A",     &fEvAsym,   "ev.A/D");
    fTree->Branch("ev.Am",    &fEvmAsym,  "ev.Am/D");
    fTree->Branch("ev.xs",    &fEvEffXS,  "ev.xs/D");
    fTree->Branch("ev.Q2",    &fEvQ2,     "ev.Q2/D");
    fTree->Branch("ev.W2",    &fEvW2,     "ev.W2/D");
    fTree->Branch("ev.thcom", &fEvThCoM,  "ev.thcom/D");
    fTree->Branch("ev.beamp",  &fEvBeamP,   "ev.beamp/D");

    fTree->Branch("ev.npart", &fNEvPart   ,     "ev.npart/I");
    fTree->Branch("ev.pid",   &fEvPID,      "ev.pid[ev.npart]/I");
    fTree->Branch("ev.vx",    &fEvPart_X,   "ev.vx[ev.npart]/D");
    fTree->Branch("ev.vy",    &fEvPart_Y,   "ev.vy[ev.npart]/D");
    fTree->Branch("ev.vz",    &fEvPart_Z,   "ev.vz[ev.npart]/D");
    fTree->Branch("ev.p",     &fEvPart_P,   "ev.p[ev.npart]/D");
    fTree->Branch("ev.px",    &fEvPart_Px,  "ev.px[ev.npart]/D");
    fTree->Branch("ev.py",    &fEvPart_Py,  "ev.py[ev.npart]/D");
    fTree->Branch("ev.pz",    &fEvPart_Pz,  "ev.pz[ev.npart]/D");
    fTree->Branch("ev.th",    &fEvPart_Th,     "ev.th[ev.npart]/D");
    fTree->Branch("ev.ph",    &fEvPart_Ph,     "ev.ph[ev.npart]/D");
    fTree->Branch("ev.tpx",    &fEvPart_tPx,  "ev.tpx[ev.npart]/D");
    fTree->Branch("ev.tpy",    &fEvPart_tPy,  "ev.tpy[ev.npart]/D");
    fTree->Branch("ev.tpz",    &fEvPart_tPz,  "ev.tpz[ev.npart]/D");

    fTree->Branch("bm.x",    &fBmX,  "bm.x/D");
    fTree->Branch("bm.y",    &fBmY,  "bm.y/D");
    fTree->Branch("bm.z",    &fBmZ,  "bm.z/D");
    fTree->Branch("bm.dx",    &fBmdX,  "bm.dx/D");
    fTree->Branch("bm.dy",    &fBmdY,  "bm.dy/D");
    fTree->Branch("bm.dz",    &fBmdZ,  "bm.dz/D");
    fTree->Branch("bm.th",    &fBmth,  "bm.th/D");
    fTree->Branch("bm.ph",    &fBmph,  "bm.ph/D");

    // GenericDetectorHit
    fTree->Branch("hit.n",    &fNGenDetHit,     "hit.n/I");
    fTree->Branch("hit.det",  &fGenDetHit_det,  "hit.det[hit.n]/I");
    fTree->Branch("hit.vid",  &fGenDetHit_id,   "hit.vid[hit.n]/I");

    fTree->Branch("hit.pid",  &fGenDetHit_pid,   "hit.pid[hit.n]/I");
    fTree->Branch("hit.trid", &fGenDetHit_trid,  "hit.trid[hit.n]/I");
    fTree->Branch("hit.mtrid",&fGenDetHit_mtrid, "hit.mtrid[hit.n]/I");
    fTree->Branch("hit.gen",  &fGenDetHit_gen,   "hit.gen[hit.n]/I");

    fTree->Branch("hit.x",    &fGenDetHit_X,   "hit.x[hit.n]/D");
    fTree->Branch("hit.y",    &fGenDetHit_Y,   "hit.y[hit.n]/D");
    fTree->Branch("hit.z",    &fGenDetHit_Z,   "hit.z[hit.n]/D");
    fTree->Branch("hit.r",    &fGenDetHit_R,   "hit.r[hit.n]/D");
    fTree->Branch("hit.ph",   &fGenDetHit_Ph,  "hit.ph[hit.n]/D");

    fTree->Branch("hit.px",   &fGenDetHit_Px,   "hit.px[hit.n]/D");
    fTree->Branch("hit.py",   &fGenDetHit_Py,   "hit.py[hit.n]/D");
    fTree->Branch("hit.pz",   &fGenDetHit_Pz,   "hit.pz[hit.n]/D");

    fTree->Branch("hit.vx",   &fGenDetHit_Vx,   "hit.vx[hit.n]/D");
    fTree->Branch("hit.vy",   &fGenDetHit_Vy,   "hit.vy[hit.n]/D");
    fTree->Branch("hit.vz",   &fGenDetHit_Vz,   "hit.vz[hit.n]/D");

    fTree->Branch("hit.p",    &fGenDetHit_P,   "hit.p[hit.n]/D");
    fTree->Branch("hit.e",    &fGenDetHit_E,   "hit.e[hit.n]/D");
    fTree->Branch("hit.m",    &fGenDetHit_M,   "hit.m[hit.n]/D");

    fTree->Branch("hit.colCut",    &fCollCut,     "hit.colCut/I");

    // GenericDetectorSum
    fTree->Branch("sum.n",    &fNGenDetSum,     "sum.n/I");
    fTree->Branch("sum.det",  &fGenDetSum_det,  "sum.det[sum.n]/I");
    fTree->Branch("sum.vid",  &fGenDetSum_id,   "sum.vid[sum.n]/I");
    fTree->Branch("sum.edep", &fGenDetSum_edep, "sum.edep[sum.n]/D");

	//Virtual boundary data
	
	fTree->Branch("lhrs",	&fLHRS,		"lhrs/I");
	fTree->Branch("rhrs",	&fRHRS,		"rhrs/I");
	
	fTree->Branch("x0",	&fX0,		"x0/D");
	fTree->Branch("th0",	&fTh0,		"th0/D");
	fTree->Branch("p0",	&fP0,		"p0/D");
	fTree->Branch("y0",	&fY0,		"y0/D");
	fTree->Branch("ph0",	&fPh0,		"ph0/D");

	fTree->Branch("x0_tr",	&fX0_tr,	"x0_tr/D");
	fTree->Branch("th0_tr",	&fTh0_tr,		"th0_tr/D");
	fTree->Branch("p0_tr",	&fP0_tr,		"p0_tr/D");
	fTree->Branch("y0_tr",	&fY0_tr,		"y0_tr/D");
	fTree->Branch("ph0_tr",	&fPh0_tr,		"ph0_tr/D");

	fTree->Branch("x_sen", 	&fX_sen,	"x_sen/D");
	fTree->Branch("y_sen", 	&fY_sen,	"y_sen/D");
	fTree->Branch("z_sen", 	&fZ_sen,	"z_sen/D");
	fTree->Branch("th_sen", &fTh_sen,	"th_sen/D");
	fTree->Branch("ph_sen", &fPh_sen,	"ph_sen/D");

	fTree->Branch("x_sen_tr", 	&fX_sen_tr,	"x_sen_tr/D");
	fTree->Branch("y_sen_tr", 	&fY_sen_tr,	"y_sen_tr/D");
	fTree->Branch("z_sen_tr", 	&fZ_sen_tr,	"z_sen_tr/D");
	fTree->Branch("th_sen_tr", &fTh_sen_tr,	"th_sen_tr/D");
	fTree->Branch("ph_sen_tr", &fPh_sen_tr,	"ph_sen_tr/D");

	fTree->Branch("x_sen_tf", 	&fX_sen_tf,	"x_sen_tf/D");
	fTree->Branch("y_sen_tf", 	&fY_sen_tf,	"y_sen_tf/D");
	fTree->Branch("z_sen_tf", 	&fZ_sen_tf,	"z_sen_tf/D");
	fTree->Branch("th_sen_tf", &fTh_sen_tf,	"th_sen_tf/D");
	fTree->Branch("ph_sen_tf", &fPh_sen_tf,	"ph_sen_tf/D");

	fTree->Branch("x_sm", 	&fX_sm,	"x_sm/D");
	fTree->Branch("y_sm", 	&fY_sm,	"y_sm/D");
	fTree->Branch("z_sm", 	&fZ_sm,	"z_sm/D");
	fTree->Branch("th_sm", &fTh_sm,	"th_sm/D");
	fTree->Branch("ph_sm", &fPh_sm,	"ph_sm/D");

	fTree->Branch("x_sm_tr", 	&fX_sm_tr,	"x_sm_tr/D");
	fTree->Branch("y_sm_tr", 	&fY_sm_tr,	"y_sm_tr/D");
	fTree->Branch("z_sm_tr", 	&fZ_sm_tr,	"z_sm_tr/D");
	fTree->Branch("th_sm_tr", &fTh_sm_tr,	"th_sm_tr/D");
	fTree->Branch("ph_sm_tr", &fPh_sm_tr,	"ph_sm_tr/D");

	fTree->Branch("x_sm_tf", 	&fX_sm_tf,	"x_sm_tf/D");
	fTree->Branch("y_sm_tf", 	&fY_sm_tf,	"y_sm_tf/D");
	fTree->Branch("z_sm_tf", 	&fZ_sm_tf,	"z_sm_tf/D");
	fTree->Branch("th_sm_tf", &fTh_sm_tf,	"th_sm_tf/D");
	fTree->Branch("ph_sm_tf", &fPh_sm_tf,	"ph_sm_tf/D");

	fTree->Branch("x_sex", 	&fX_sex,	"x_sex/D");
	fTree->Branch("y_sex", 	&fY_sex,	"y_sex/D");
	fTree->Branch("z_sex", 	&fZ_sex,	"z_sex/D");
	fTree->Branch("th_sex", &fTh_sex,	"th_sex/D");
	fTree->Branch("ph_sex", &fPh_sex,	"ph_sex/D");

	fTree->Branch("x_sex_tr", 	&fX_sex_tr,	"x_sex_tr/D");
	fTree->Branch("y_sex_tr", 	&fY_sex_tr,	"y_sex_tr/D");
	fTree->Branch("z_sex_tr", 	&fZ_sex_tr,	"z_sex_tr/D");
	fTree->Branch("th_sex_tr", &fTh_sex_tr,	"th_sex_tr/D");
	fTree->Branch("ph_sex_tr", &fPh_sex_tr,	"ph_sex_tr/D");

	fTree->Branch("x_sex_tf", 	&fX_sex_tf,	"x_sex_tf/D");
	fTree->Branch("y_sex_tf", 	&fY_sex_tf,	"y_sex_tf/D");
	fTree->Branch("z_sex_tf", 	&fZ_sex_tf,	"z_sex_tf/D");
	fTree->Branch("th_sex_tf", &fTh_sex_tf,	"th_sex_tf/D");
	fTree->Branch("ph_sex_tf", &fPh_sex_tf,	"ph_sex_tf/D");

	fTree->Branch("x_col", 	&fX_col,	"x_col/D");
	fTree->Branch("y_col", 	&fY_col,	"y_col/D");
	fTree->Branch("z_col", 	&fZ_col,	"z_col/D");
	fTree->Branch("th_col", &fTh_col,	"th_col/D");
	fTree->Branch("ph_col", &fPh_col,	"ph_col/D");

	fTree->Branch("x_col_tr", 	&fX_col_tr,	"x_col_tr/D");
	fTree->Branch("y_col_tr", 	&fY_col_tr,	"y_col_tr/D");
	fTree->Branch("z_col_tr", 	&fZ_col_tr,	"z_col_tr/D");
	fTree->Branch("th_col_tr", &fTh_col_tr,	"th_col_tr/D");
	fTree->Branch("ph_col_tr", &fPh_col_tr,	"ph_col_tr/D");

	fTree->Branch("x_col_tf", 	&fX_col_tf,	"x_col_tf/D");
	fTree->Branch("y_col_tf", 	&fY_col_tf,	"y_col_tf/D");
	fTree->Branch("z_col_tf", 	&fZ_col_tf,	"z_col_tf/D");
	fTree->Branch("th_col_tf", &fTh_col_tf,	"th_col_tf/D");
	fTree->Branch("ph_col_tf", &fPh_col_tf,	"ph_col_tf/D");
	
	fTree->Branch("x_col", 	&fX_col,	"x_col/D");
	fTree->Branch("y_col", 	&fY_col,	"y_col/D");
	fTree->Branch("z_col", 	&fZ_col,	"z_col/D");
	fTree->Branch("th_col", &fTh_col,	"th_col/D");
	fTree->Branch("ph_col", &fPh_col,	"ph_col/D");

	fTree->Branch("x_col_tr", 	&fX_col_tr,	"x_col_tr/D");
	fTree->Branch("y_col_tr", 	&fY_col_tr,	"y_col_tr/D");
	fTree->Branch("z_col_tr", 	&fZ_col_tr,	"z_col_tr/D");
	fTree->Branch("th_col_tr", &fTh_col_tr,	"th_col_tr/D");
	fTree->Branch("ph_col_tr", &fPh_col_tr,	"ph_col_tr/D");

	fTree->Branch("x_col_tf", 	&fX_col_tf,	"x_col_tf/D");
	fTree->Branch("y_col_tf", 	&fY_col_tf,	"y_col_tf/D");
	fTree->Branch("z_col_tf", 	&fZ_col_tf,	"z_col_tf/D");
	fTree->Branch("th_col_tf", &fTh_col_tf,	"th_col_tf/D");
	fTree->Branch("ph_col_tf", &fPh_col_tf,	"ph_col_tf/D");
	
	fTree->Branch("x_q1en", 	&fX_q1en,	"x_q1en/D");
	fTree->Branch("y_q1en", 	&fY_q1en,	"y_q1en/D");
	fTree->Branch("z_q1en", 	&fZ_q1en,	"z_q1en/D");
	fTree->Branch("th_q1en", &fTh_q1en,	"th_q1en/D");
	fTree->Branch("ph_q1en", &fPh_q1en,	"ph_q1en/D");

	fTree->Branch("x_q1en_tr", 	&fX_q1en_tr,	"x_q1en_tr/D");
	fTree->Branch("y_q1en_tr", 	&fY_q1en_tr,	"y_q1en_tr/D");
	fTree->Branch("z_q1en_tr", 	&fZ_q1en_tr,	"z_q1en_tr/D");
	fTree->Branch("th_q1en_tr", &fTh_q1en_tr,	"th_q1en_tr/D");
	fTree->Branch("ph_q1en_tr", &fPh_q1en_tr,	"ph_q1en_tr/D");
	
	fTree->Branch("x_q1ex", 	&fX_q1ex,	"x_q1ex/D");
	fTree->Branch("y_q1ex", 	&fY_q1ex,	"y_q1ex/D");
	fTree->Branch("z_q1ex", 	&fZ_q1ex,	"z_q1ex/D");
	fTree->Branch("th_q1ex", &fTh_q1ex,	"th_q1ex/D");
	fTree->Branch("ph_q1ex", &fPh_q1ex,	"ph_q1ex/D");

	fTree->Branch("x_q1ex_tr", 	&fX_q1ex_tr,	"x_q1ex_tr/D");
	fTree->Branch("y_q1ex_tr", 	&fY_q1ex_tr,	"y_q1ex_tr/D");
	fTree->Branch("z_q1ex_tr", 	&fZ_q1ex_tr,	"z_q1ex_tr/D");
	fTree->Branch("th_q1ex_tr", &fTh_q1ex_tr,	"th_q1ex_tr/D");
	fTree->Branch("ph_q1ex_tr", &fPh_q1ex_tr,	"ph_q1ex_tr/D");

	fTree->Branch("x_q1ex_tf", 	&fX_q1ex_tf,	"x_q1ex_tf/D");
	fTree->Branch("y_q1ex_tf", 	&fY_q1ex_tf,	"y_q1ex_tf/D");
	fTree->Branch("z_q1ex_tf", 	&fZ_q1ex_tf,	"z_q1ex_tf/D");
	fTree->Branch("th_q1ex_tf", &fTh_q1ex_tf,	"th_q1ex_tf/D");
	fTree->Branch("ph_q1ex_tf", &fPh_q1ex_tf,	"ph_q1ex_tf/D");
	
	fTree->Branch("x_q2en", 	&fX_q2en,	"x_q2en/D");
	fTree->Branch("y_q2en", 	&fY_q2en,	"y_q2en/D");
	fTree->Branch("z_q2en", 	&fZ_q2en,	"z_q2en/D");
	fTree->Branch("th_q2en", &fTh_q2en,	"th_q2en/D");
	fTree->Branch("ph_q2en", &fPh_q2en,	"ph_q2en/D");

	fTree->Branch("x_q2en_tr", 	&fX_q2en_tr,	"x_q2en_tr/D");
	fTree->Branch("y_q2en_tr", 	&fY_q2en_tr,	"y_q2en_tr/D");
	fTree->Branch("z_q2en_tr", 	&fZ_q2en_tr,	"z_q2en_tr/D");
	fTree->Branch("th_q2en_tr", &fTh_q2en_tr,	"th_q2en_tr/D");
	fTree->Branch("ph_q2en_tr", &fPh_q2en_tr,	"ph_q2en_tr/D");

	fTree->Branch("x_q2ex", 	&fX_q2ex,	"x_q2ex/D");
	fTree->Branch("y_q2ex", 	&fY_q2ex,	"y_q2ex/D");
	fTree->Branch("z_q2ex", 	&fZ_q2ex,	"z_q2ex/D");
	fTree->Branch("th_q2ex", &fTh_q2ex,	"th_q2ex/D");
	fTree->Branch("ph_q2ex", &fPh_q2ex,	"ph_q2ex/D");

	fTree->Branch("x_q2ex_tr", 	&fX_q2ex_tr,	"x_q2ex_tr/D");
	fTree->Branch("y_q2ex_tr", 	&fY_q2ex_tr,	"y_q2ex_tr/D");
	fTree->Branch("z_q2ex_tr", 	&fZ_q2ex_tr,	"z_q2ex_tr/D");
	fTree->Branch("th_q2ex_tr", &fTh_q2ex_tr,	"th_q2ex_tr/D");
	fTree->Branch("ph_q2ex_tr", &fPh_q2ex_tr,	"ph_q2ex_tr/D");

	fTree->Branch("x_q2ex_tf", 	&fX_q2ex_tf,	"x_q2ex_tf/D");
	fTree->Branch("y_q2ex_tf", 	&fY_q2ex_tf,	"y_q2ex_tf/D");
	fTree->Branch("z_q2ex_tf", 	&fZ_q2ex_tf,	"z_q2ex_tf/D");
	fTree->Branch("th_q2ex_tf", &fTh_q2ex_tf,	"th_q2ex_tf/D");
	fTree->Branch("ph_q2ex_tf", &fPh_q2ex_tf,	"ph_q2ex_tf/D");
	
	fTree->Branch("x_den", 	&fX_den,	"x_den/D");
	fTree->Branch("y_den", 	&fY_den,	"y_den/D");
	fTree->Branch("z_den", 	&fZ_den,	"z_den/D");
	fTree->Branch("th_den", &fTh_den,	"th_den/D");
	fTree->Branch("ph_den", &fPh_den,	"ph_den/D");

	fTree->Branch("x_den_tr", 	&fX_den_tr,	"x_den_tr/D");
	fTree->Branch("y_den_tr", 	&fY_den_tr,	"y_den_tr/D");
	fTree->Branch("z_den_tr", 	&fZ_den_tr,	"z_den_tr/D");
	fTree->Branch("th_den_tr", &fTh_den_tr,	"th_den_tr/D");
	fTree->Branch("ph_den_tr", &fPh_den_tr,	"ph_den_tr/D");

	fTree->Branch("x_den_tf", 	&fX_den_tf,	"x_den_tf/D");
	fTree->Branch("y_den_tf", 	&fY_den_tf,	"y_den_tf/D");
	fTree->Branch("z_den_tf", 	&fZ_den_tf,	"z_den_tf/D");
	fTree->Branch("th_den_tf", &fTh_den_tf,	"th_den_tf/D");
	fTree->Branch("ph_den_tf", &fPh_den_tf,	"ph_den_tf/D");
	
	fTree->Branch("x_dex", 	&fX_dex,	"x_dex/D");
	fTree->Branch("y_dex", 	&fY_dex,	"y_dex/D");
	fTree->Branch("z_dex", 	&fZ_dex,	"z_dex/D");
	fTree->Branch("th_dex", &fTh_dex,	"th_dex/D");
	fTree->Branch("ph_dex", &fPh_dex,	"ph_dex/D");

	fTree->Branch("x_dex_tr", 	&fX_dex_tr,	"x_dex_tr/D");
	fTree->Branch("y_dex_tr", 	&fY_dex_tr,	"y_dex_tr/D");
	fTree->Branch("z_dex_tr", 	&fZ_dex_tr,	"z_dex_tr/D");
	fTree->Branch("th_dex_tr", &fTh_dex_tr,	"th_dex_tr/D");
	fTree->Branch("ph_dex_tr", &fPh_dex_tr,	"ph_dex_tr/D");

	fTree->Branch("x_dex_tf", 	&fX_dex_tf,	"x_dex_tf/D");
	fTree->Branch("y_dex_tf", 	&fY_dex_tf,	"y_dex_tf/D");
	fTree->Branch("z_dex_tf", 	&fZ_dex_tf,	"z_dex_tf/D");
	fTree->Branch("th_dex_tf", &fTh_dex_tf,	"th_dex_tf/D");
	fTree->Branch("ph_dex_tf", &fPh_dex_tf,	"ph_dex_tf/D");
	
	fTree->Branch("x_q3en", 	&fX_q3en,	"x_q3en/D");
	fTree->Branch("y_q3en", 	&fY_q3en,	"y_q3en/D");
	fTree->Branch("z_q3en", 	&fZ_q3en,	"z_q3en/D");
	fTree->Branch("th_q3en", &fTh_q3en,	"th_q3en/D");
	fTree->Branch("ph_q3en", &fPh_q3en,	"ph_q3en/D");

	fTree->Branch("x_q3en_tr", 	&fX_q3en_tr,	"x_q3en_tr/D");
	fTree->Branch("y_q3en_tr", 	&fY_q3en_tr,	"y_q3en_tr/D");
	fTree->Branch("z_q3en_tr", 	&fZ_q3en_tr,	"z_q3en_tr/D");
	fTree->Branch("th_q3en_tr", &fTh_q3en_tr,	"th_q3en_tr/D");
	fTree->Branch("ph_q3en_tr", &fPh_q3en_tr,	"ph_q3en_tr/D");

	fTree->Branch("x_q3en_tf", 	&fX_q3en_tf,	"x_q3en_tf/D");
	fTree->Branch("y_q3en_tf", 	&fY_q3en_tf,	"y_q3en_tf/D");
	fTree->Branch("z_q3en_tf", 	&fZ_q3en_tf,	"z_q3en_tf/D");
	fTree->Branch("th_q3en_tf", &fTh_q3en_tf,	"th_q3en_tf/D");
	fTree->Branch("ph_q3en_tf", &fPh_q3en_tf,	"ph_q3en_tf/D");
	
	fTree->Branch("x_q3ex", 	&fX_q3ex,	"x_q3ex/D");
	fTree->Branch("y_q3ex", 	&fY_q3ex,	"y_q3ex/D");
	fTree->Branch("z_q3ex", 	&fZ_q3ex,	"z_q3ex/D");
	fTree->Branch("th_q3ex", &fTh_q3ex,	"th_q3ex/D");
	fTree->Branch("ph_q3ex", &fPh_q3ex,	"ph_q3ex/D");

	fTree->Branch("x_q3ex_tr", 	&fX_q3ex_tr,	"x_q3ex_tr/D");
	fTree->Branch("y_q3ex_tr", 	&fY_q3ex_tr,	"y_q3ex_tr/D");
	fTree->Branch("z_q3ex_tr", 	&fZ_q3ex_tr,	"z_q3ex_tr/D");
	fTree->Branch("th_q3ex_tr", &fTh_q3ex_tr,	"th_q3ex_tr/D");
	fTree->Branch("ph_q3ex_tr", &fPh_q3ex_tr,	"ph_q3ex_tr/D");

	fTree->Branch("x_q3ex_tf", 	&fX_q3ex_tf,	"x_q3ex_tf/D");
	fTree->Branch("y_q3ex_tf", 	&fY_q3ex_tf,	"y_q3ex_tf/D");
	fTree->Branch("z_q3ex_tf", 	&fZ_q3ex_tf,	"z_q3ex_tf/D");
	fTree->Branch("th_q3ex_tf", &fTh_q3ex_tf,	"th_q3ex_tf/D");
	fTree->Branch("ph_q3ex_tf", &fPh_q3ex_tf,	"ph_q3ex_tf/D");
		
	fTree->Branch("x_vdc", 	&fX_vdc,	"x_vdc/D");
	fTree->Branch("y_vdc", 	&fY_vdc,	"y_vdc/D");
	fTree->Branch("z_vdc", 	&fZ_vdc,	"z_vdc/D");
	fTree->Branch("th_vdc", &fTh_vdc,	"th_vdc/D");
	fTree->Branch("ph_vdc", &fPh_vdc,	"ph_vdc/D");

	fTree->Branch("x_vdc_tr", 	&fX_vdc_tr,	"x_vdc_tr/D");
	fTree->Branch("y_vdc_tr", 	&fY_vdc_tr,	"y_vdc_tr/D");
	fTree->Branch("z_vdc_tr", 	&fZ_vdc_tr,	"z_vdc_tr/D");
	fTree->Branch("th_vdc_tr", &fTh_vdc_tr,	"th_vdc_tr/D");
	fTree->Branch("ph_vdc_tr", &fPh_vdc_tr,	"ph_vdc_tr/D");

	fTree->Branch("x_fp", 	&fX_fp,	"x_fp/D");
	fTree->Branch("y_fp", 	&fY_fp,	"y_fp/D");
	fTree->Branch("z_fp", 	&fZ_fp,	"z_fp/D");
	fTree->Branch("th_fp", &fTh_fp,	"th_fp/D");
	fTree->Branch("ph_fp", &fPh_fp,	"ph_fp/D");

	fTree->Branch("x_fp_tr", 	&fX_fp_tr,	"x_fp_tr/D");
	fTree->Branch("y_fp_tr", 	&fY_fp_tr,	"y_fp_tr/D");
	fTree->Branch("z_fp_tr", 	&fZ_fp_tr,	"z_fp_tr/D");
	fTree->Branch("th_fp_tr", &fTh_fp_tr,	"th_fp_tr/D");
	fTree->Branch("ph_fp_tr", &fPh_fp_tr,	"ph_fp_tr/D");

	fTree->Branch("x_fp_tf", 	&fX_fp_tf,	"x_fp_tf/D");
	fTree->Branch("y_fp_tf", 	&fY_fp_tf,	"y_fp_tf/D");
	fTree->Branch("z_fp_tf", 	&fZ_fp_tf,	"z_fp_tf/D");
	fTree->Branch("th_fp_tf", &fTh_fp_tf,	"th_fp_tf/D");
	fTree->Branch("ph_fp_tf", &fPh_fp_tf,	"ph_fp_tf/D");

    return;
}

void g4hrsIO::FillTree(){
    if( !fTree ){ 
	fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
	return; 
    }

    fTree->Fill();
    fTree->GetCurrentFile();
}

void g4hrsIO::Flush(){
    //  Set arrays to 0
    fNGenDetHit = 0;
    fNGenDetSum = 0;
    fCollCut = 1; // default
}

void g4hrsIO::WriteTree(){
    assert( fFile );
    assert( fTree );

    if( !fFile->IsOpen() ){
	G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << ": TFile not open" << G4endl;
	exit(1);
    }

    G4cout << "Writing output to " << fFile->GetName() << "... ";

    fFile->cd();

    fTree->Write("T", TObject::kOverwrite);
    g4hrsRun::GetRun()->GetData()->Write("run_data", TObject::kOverwrite); 

    fTree->ResetBranchAddresses();
    delete fTree;
    fTree = NULL;

    fFile->Close();

    delete fFile;
    fFile = NULL;

    G4cout << "written" << G4endl;

    return;
}

///////////////////////////////////////////////////////////////////////////////
// Interfaces to output section ///////////////////////////////////////////////

// Event Data

void g4hrsIO::SetEventData(g4hrsEvent *ev){
    int n = ev->fPartType.size();
    if( n > __IO_MAXHIT ){
	G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ":  Buffer size exceeded!" << G4endl;
	return;
    }

    fNEvPart = n;

    fEvRate   = ev->fRate*s;
    fEvEffXS  = ev->fEffXs/microbarn;
    fEvAsym   = ev->fAsym/__ASYMM_SCALE;
    fEvmAsym  = ev->fmAsym/__ASYMM_SCALE;
    fEvBeamP  = ev->fBeamMomentum.mag()/__E_UNIT;

    fEvQ2     = ev->fQ2/__E_UNIT/__E_UNIT;
    fEvW2     = ev->fW2/__E_UNIT/__E_UNIT;
    fEvThCoM  = ev->fThCoM/deg; // specify this in degrees over anything else

    int idx;
    for( idx = 0; idx < n; idx++ ){
	fEvPID[idx] = ev->fPartType[idx]->GetPDGEncoding();

	fEvPart_X[idx] = ev->fPartPos[idx].x()/__L_UNIT;
	fEvPart_Y[idx] = ev->fPartPos[idx].y()/__L_UNIT;
	fEvPart_Z[idx] = ev->fPartPos[idx].z()/__L_UNIT;

	fEvPart_Px[idx] = ev->fPartRealMom[idx].x()/__E_UNIT;
	fEvPart_Py[idx] = ev->fPartRealMom[idx].y()/__E_UNIT;
	fEvPart_Pz[idx] = ev->fPartRealMom[idx].z()/__E_UNIT;
	fEvPart_Th[idx] = ev->fPartRealMom[idx].theta();
	fEvPart_Ph[idx] = ev->fPartRealMom[idx].phi();

	fEvPart_P[idx] = ev->fPartRealMom[idx].mag()/__E_UNIT;

	fEvPart_tPx[idx] = ev->fPartMom[idx].x()/__E_UNIT;
	fEvPart_tPy[idx] = ev->fPartMom[idx].y()/__E_UNIT;
	fEvPart_tPz[idx] = ev->fPartMom[idx].z()/__E_UNIT;
    }

    /////////////////////////////////////////////////
    //  Set beam data as well

    g4hrsBeamTarget *bt = g4hrsBeamTarget::GetBeamTarget();

    fBmX = bt->fVer.x()/__L_UNIT;
    fBmY = bt->fVer.y()/__L_UNIT;
    fBmZ = bt->fVer.z()/__L_UNIT;
    
    fBmdX = bt->fDir.x();
    fBmdY = bt->fDir.y();
    fBmdZ = bt->fDir.z();
    fBmth = bt->fDir.theta();
    fBmph = bt->fDir.phi()/deg;

    //    G4cout << "** fDir:: " << bt->fDir.x()/deg << "  " << bt->fDir.y()/deg << "  " << bt->fVer.z()/mm << G4endl;

    return;
}

// GenericDetectorHit

void g4hrsIO::AddGenericDetectorHit(g4hrsGenericDetectorHit *hit){
    int n = fNGenDetHit;
    if( n > __IO_MAXHIT ){
	G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ":  Buffer size exceeded!" << G4endl;
	return;
    }

    fGenDetHit_det[n]  = hit->fDetID;
    fGenDetHit_id[n]   = hit->fCopyID;

    fGenDetHit_trid[n] = hit->fTrID;
    fGenDetHit_mtrid[n]= hit->fmTrID;
    fGenDetHit_pid[n]  = hit->fPID;
    fGenDetHit_gen[n]  = hit->fGen;

    fGenDetHit_X[n]  = hit->f3X.x()/__L_UNIT;
    fGenDetHit_Y[n]  = hit->f3X.y()/__L_UNIT;
    fGenDetHit_Z[n]  = hit->f3X.z()/__L_UNIT;
    fGenDetHit_R[n]  = sqrt(hit->f3X.x()*hit->f3X.x()+hit->f3X.y()*hit->f3X.y())/__L_UNIT;
    fGenDetHit_Ph[n] = hit->f3X.phi()/deg;

    fGenDetHit_Px[n]  = hit->f3P.x()/__E_UNIT;
    fGenDetHit_Py[n]  = hit->f3P.y()/__E_UNIT;
    fGenDetHit_Pz[n]  = hit->f3P.z()/__E_UNIT;

    fGenDetHit_Vx[n]  = hit->f3V.x()/__L_UNIT;
    fGenDetHit_Vy[n]  = hit->f3V.y()/__L_UNIT;
    fGenDetHit_Vz[n]  = hit->f3V.z()/__L_UNIT;

    fGenDetHit_P[n]  = hit->fP/__E_UNIT;
    fGenDetHit_E[n]  = hit->fE/__E_UNIT;
    fGenDetHit_M[n]  = hit->fM/__E_UNIT;

    fNGenDetHit++;

    // for collimator cut
    if( (hit->fDetID==200 && hit->f3X.perp()/__L_UNIT < 0.03) || 
      	(hit->fDetID==201 && hit->f3X.perp()/__L_UNIT < 0.05) )
      fCollCut=0;
}


// GenericDetectorSum

void g4hrsIO::AddGenericDetectorSum(g4hrsGenericDetectorSum *hit){
    int n = fNGenDetSum;
    if( n > __IO_MAXHIT ){
	G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ":  Buffer size exceeded!" << G4endl;
	return;
    }

    fGenDetSum_edep[n] = hit->fEdep/__E_UNIT;
    fGenDetSum_det[n]  = hit->fDetID;
    fGenDetSum_id[n]   = hit->fCopyID;

    fNGenDetSum++;
}

// Virtual boundary data

void g4hrsIO::ClearVirtualBoundaryData() {

	fSteppingAction->fLHRS = 0;
	fSteppingAction->fRHRS = 0;

	fSteppingAction->fX0 = -333.;
	fSteppingAction->fY0 = -333.;
	fSteppingAction->fZ0 = -333.;
	fSteppingAction->fP0 = -333.;
	fSteppingAction->fTh0 = -333.;
	fSteppingAction->fPh0 = -333.;
	
	fSteppingAction->fX0_tr = -333.;
	fSteppingAction->fY0_tr = -333.;
	fSteppingAction->fZ0_tr = -333.;
	fSteppingAction->fP0_tr = -333.;
	fSteppingAction->fTh0_tr = -333.;
	fSteppingAction->fPh0_tr = -333.;

	fSteppingAction->fX_sen = -333.;
	fSteppingAction->fY_sen = -333.;
	fSteppingAction->fZ_sen = -333.;
	fSteppingAction->fTh_sen = -333.;
	fSteppingAction->fPh_sen = -333.;
	
	fSteppingAction->fX_sen_tr = -333.;
	fSteppingAction->fY_sen_tr = -333.;
	fSteppingAction->fZ_sen_tr = -333.;
	fSteppingAction->fTh_sen_tr = -333.;
	fSteppingAction->fPh_sen_tr = -333.;
		
	fSteppingAction->fX_sen_tf = -333.;
	fSteppingAction->fY_sen_tf = -333.;
	fSteppingAction->fZ_sen_tf = -333.;
	fSteppingAction->fTh_sen_tf = -333.;
	fSteppingAction->fPh_sen_tf = -333.;
	
	fSteppingAction->fX_sm = -333.;
	fSteppingAction->fY_sm = -333.;
	fSteppingAction->fZ_sm = -333.;
	fSteppingAction->fTh_sm = -333.;
	fSteppingAction->fPh_sm = -333.;
	
	fSteppingAction->fX_sm_tr = -333.;
	fSteppingAction->fY_sm_tr = -333.;
	fSteppingAction->fZ_sm_tr = -333.;
	fSteppingAction->fTh_sm_tr = -333.;
	fSteppingAction->fPh_sm_tr = -333.;
		
	fSteppingAction->fX_sm_tf = -333.;
	fSteppingAction->fY_sm_tf = -333.;
	fSteppingAction->fZ_sm_tf = -333.;
	fSteppingAction->fTh_sm_tf = -333.;
	fSteppingAction->fPh_sm_tf = -333.;

	fSteppingAction->fX_sex = -333.;
	fSteppingAction->fY_sex = -333.;
	fSteppingAction->fZ_sex = -333.;
	fSteppingAction->fTh_sex = -333.;
	fSteppingAction->fPh_sex = -333.;
	
	fSteppingAction->fX_sex_tr = -333.;
	fSteppingAction->fY_sex_tr = -333.;
	fSteppingAction->fZ_sex_tr = -333.;
	fSteppingAction->fTh_sex_tr = -333.;
	fSteppingAction->fPh_sex_tr = -333.;
		
	fSteppingAction->fX_sex_tf = -333.;
	fSteppingAction->fY_sex_tf = -333.;
	fSteppingAction->fZ_sex_tf = -333.;
	fSteppingAction->fTh_sex_tf = -333.;
	fSteppingAction->fPh_sex_tf = -333.;

	fSteppingAction->fX_col = -333.;
	fSteppingAction->fY_col = -333.;
	fSteppingAction->fZ_col = -333.;
	fSteppingAction->fTh_col = -333.;
	fSteppingAction->fPh_col = -333.;
	
	fSteppingAction->fX_col_tr = -333.;
	fSteppingAction->fY_col_tr = -333.;
	fSteppingAction->fZ_col_tr = -333.;
	fSteppingAction->fTh_col_tr = -333.;
	fSteppingAction->fPh_col_tr = -333.;
		
	fSteppingAction->fX_col_tf = -333.;
	fSteppingAction->fY_col_tf = -333.;
	fSteppingAction->fZ_col_tf = -333.;
	fSteppingAction->fTh_col_tf = -333.;
	fSteppingAction->fPh_col_tf = -333.;

	fSteppingAction->fX_q1en = -333.;
	fSteppingAction->fY_q1en = -333.;
	fSteppingAction->fZ_q1en = -333.;
	fSteppingAction->fTh_q1en = -333.;
	fSteppingAction->fPh_q1en = -333.;
	
	fSteppingAction->fX_q1en_tr = -333.;
	fSteppingAction->fY_q1en_tr = -333.;
	fSteppingAction->fZ_q1en_tr = -333.;
	fSteppingAction->fTh_q1en_tr = -333.;
	fSteppingAction->fPh_q1en_tr = -333.;
		
	fSteppingAction->fX_q1ex = -333.;
	fSteppingAction->fY_q1ex = -333.;
	fSteppingAction->fZ_q1ex = -333.;
	fSteppingAction->fTh_q1ex = -333.;
	fSteppingAction->fPh_q1ex = -333.;
	
	fSteppingAction->fX_q1ex_tr = -333.;
	fSteppingAction->fY_q1ex_tr = -333.;
	fSteppingAction->fZ_q1ex_tr = -333.;
	fSteppingAction->fTh_q1ex_tr = -333.;
	fSteppingAction->fPh_q1ex_tr = -333.;
		
	fSteppingAction->fX_q1ex_tf = -333.;
	fSteppingAction->fY_q1ex_tf = -333.;
	fSteppingAction->fZ_q1ex_tf = -333.;
	fSteppingAction->fTh_q1ex_tf = -333.;
	fSteppingAction->fPh_q1ex_tf = -333.;

	fSteppingAction->fX_q2en = -333.;
	fSteppingAction->fY_q2en = -333.;
	fSteppingAction->fZ_q2en = -333.;
	fSteppingAction->fTh_q2en = -333.;
	fSteppingAction->fPh_q2en = -333.;
	
	fSteppingAction->fX_q2en_tr = -333.;
	fSteppingAction->fY_q2en_tr = -333.;
	fSteppingAction->fZ_q2en_tr = -333.;
	fSteppingAction->fTh_q2en_tr = -333.;
	fSteppingAction->fPh_q2en_tr = -333.;
		
	fSteppingAction->fX_q2ex = -333.;
	fSteppingAction->fY_q2ex = -333.;
	fSteppingAction->fZ_q2ex = -333.;
	fSteppingAction->fTh_q2ex = -333.;
	fSteppingAction->fPh_q2ex = -333.;
	
	fSteppingAction->fX_q2ex_tr = -333.;
	fSteppingAction->fY_q2ex_tr = -333.;
	fSteppingAction->fZ_q2ex_tr = -333.;
	fSteppingAction->fTh_q2ex_tr = -333.;
	fSteppingAction->fPh_q2ex_tr = -333.;
		
	fSteppingAction->fX_q2ex_tf = -333.;
	fSteppingAction->fY_q2ex_tf = -333.;
	fSteppingAction->fZ_q2ex_tf = -333.;
	fSteppingAction->fTh_q2ex_tf = -333.;
	fSteppingAction->fPh_q2ex_tf = -333.;

	fSteppingAction->fX_den = -333.;
	fSteppingAction->fY_den = -333.;
	fSteppingAction->fZ_den = -333.;
	fSteppingAction->fTh_den = -333.;
	fSteppingAction->fPh_den = -333.;
	
	fSteppingAction->fX_den_tr = -333.;
	fSteppingAction->fY_den_tr = -333.;
	fSteppingAction->fZ_den_tr = -333.;
	fSteppingAction->fTh_den_tr = -333.;
	fSteppingAction->fPh_den_tr = -333.;
		
	fSteppingAction->fX_den_tf = -333.;
	fSteppingAction->fY_den_tf = -333.;
	fSteppingAction->fZ_den_tf = -333.;
	fSteppingAction->fTh_den_tf = -333.;
	fSteppingAction->fPh_den_tf = -333.;

	fSteppingAction->fX_dex = -333.;
	fSteppingAction->fY_dex = -333.;
	fSteppingAction->fZ_dex = -333.;
	fSteppingAction->fTh_dex = -333.;
	fSteppingAction->fPh_dex = -333.;
	
	fSteppingAction->fX_dex_tr = -333.;
	fSteppingAction->fY_dex_tr = -333.;
	fSteppingAction->fZ_dex_tr = -333.;
	fSteppingAction->fTh_dex_tr = -333.;
	fSteppingAction->fPh_dex_tr = -333.;
		
	fSteppingAction->fX_dex_tf = -333.;
	fSteppingAction->fY_dex_tf = -333.;
	fSteppingAction->fZ_dex_tf = -333.;
	fSteppingAction->fTh_dex_tf = -333.;
	fSteppingAction->fPh_dex_tf = -333.;

	fSteppingAction->fX_q3en = -333.;
	fSteppingAction->fY_q3en = -333.;
	fSteppingAction->fZ_q3en = -333.;
	fSteppingAction->fTh_q3en = -333.;
	fSteppingAction->fPh_q3en = -333.;
	
	fSteppingAction->fX_q3en_tr = -333.;
	fSteppingAction->fY_q3en_tr = -333.;
	fSteppingAction->fZ_q3en_tr = -333.;
	fSteppingAction->fTh_q3en_tr = -333.;
	fSteppingAction->fPh_q3en_tr = -333.;
		
	fSteppingAction->fX_q3en_tf = -333.;
	fSteppingAction->fY_q3en_tf = -333.;
	fSteppingAction->fZ_q3en_tf = -333.;
	fSteppingAction->fTh_q3en_tf = -333.;
	fSteppingAction->fPh_q3en_tf = -333.;

	fSteppingAction->fX_q3ex = -333.;
	fSteppingAction->fY_q3ex = -333.;
	fSteppingAction->fZ_q3ex = -333.;
	fSteppingAction->fTh_q3ex = -333.;
	fSteppingAction->fPh_q3ex = -333.;
	
	fSteppingAction->fX_q3ex_tr = -333.;
	fSteppingAction->fY_q3ex_tr = -333.;
	fSteppingAction->fZ_q3ex_tr = -333.;
	fSteppingAction->fTh_q3ex_tr = -333.;
	fSteppingAction->fPh_q3ex_tr = -333.;
		
	fSteppingAction->fX_q3ex_tf = -333.;
	fSteppingAction->fY_q3ex_tf = -333.;
	fSteppingAction->fZ_q3ex_tf = -333.;
	fSteppingAction->fTh_q3ex_tf = -333.;
	fSteppingAction->fPh_q3ex_tf = -333.;

	fSteppingAction->fX_vdc = -333.;
	fSteppingAction->fY_vdc = -333.;
	fSteppingAction->fZ_vdc = -333.;
	fSteppingAction->fTh_vdc = -333.;
	fSteppingAction->fPh_vdc = -333.;
	
	fSteppingAction->fX_vdc_tr = -333.;
	fSteppingAction->fY_vdc_tr = -333.;
	fSteppingAction->fZ_vdc_tr = -333.;
	fSteppingAction->fTh_vdc_tr = -333.;
	fSteppingAction->fPh_vdc_tr = -333.;

	fSteppingAction->fX_fp = -333.;
	fSteppingAction->fY_fp = -333.;
	fSteppingAction->fZ_fp = -333.;
	fSteppingAction->fTh_fp = -333.;
	fSteppingAction->fPh_fp = -333.;
	
	fSteppingAction->fX_fp_tr = -333.;
	fSteppingAction->fY_fp_tr = -333.;
	fSteppingAction->fZ_fp_tr = -333.;
	fSteppingAction->fTh_fp_tr = -333.;
	fSteppingAction->fPh_fp_tr = -333.;
		
	fSteppingAction->fX_fp_tf = -333.;
	fSteppingAction->fY_fp_tf = -333.;
	fSteppingAction->fZ_fp_tf = -333.;
	fSteppingAction->fTh_fp_tf = -333.;
	fSteppingAction->fPh_fp_tf = -333.;

}

void g4hrsIO::SetVirtualBoundaryData() {

	fLHRS = fSteppingAction->fLHRS;
	fRHRS = fSteppingAction->fRHRS;

	fX0 = fSteppingAction->fX0;
	fY0 = fSteppingAction->fY0;
	fZ0 = fSteppingAction->fZ0;
	fP0 = fSteppingAction->fP0;
	fTh0 = fSteppingAction->fTh0;
	fPh0 = fSteppingAction->fPh0;

	fX0_tr = fSteppingAction->fX0_tr;
	fY0_tr = fSteppingAction->fY0_tr;
	fZ0_tr = fSteppingAction->fZ0_tr;
	fP0_tr = fSteppingAction->fP0_tr;
	fTh0_tr = fSteppingAction->fTh0_tr;
	fPh0_tr = fSteppingAction->fPh0_tr;

	fX_sen = fSteppingAction->fX_sen;
	fY_sen = fSteppingAction->fY_sen;
	fZ_sen = fSteppingAction->fZ_sen;
	fTh_sen = fSteppingAction->fTh_sen;
	fPh_sen = fSteppingAction->fPh_sen;

	fX_sen_tr = fSteppingAction->fX_sen_tr;
	fY_sen_tr = fSteppingAction->fY_sen_tr;
	fZ_sen_tr = fSteppingAction->fZ_sen_tr;
	fTh_sen_tr = fSteppingAction->fTh_sen_tr;
	fPh_sen_tr = fSteppingAction->fPh_sen_tr;

	fX_sen_tf = fSteppingAction->fX_sen_tf;
	fY_sen_tf = fSteppingAction->fY_sen_tf;
	fZ_sen_tf = fSteppingAction->fZ_sen_tf;
	fTh_sen_tf = fSteppingAction->fTh_sen_tf;
	fPh_sen_tf = fSteppingAction->fPh_sen_tf;

	fX_sm = fSteppingAction->fX_sm;
	fY_sm = fSteppingAction->fY_sm;
	fZ_sm = fSteppingAction->fZ_sm;
	fTh_sm = fSteppingAction->fTh_sm;
	fPh_sm = fSteppingAction->fPh_sm;

	fX_sm_tr = fSteppingAction->fX_sm_tr;
	fY_sm_tr = fSteppingAction->fY_sm_tr;
	fZ_sm_tr = fSteppingAction->fZ_sm_tr;
	fTh_sm_tr = fSteppingAction->fTh_sm_tr;
	fPh_sm_tr = fSteppingAction->fPh_sm_tr;

	fX_sm_tf = fSteppingAction->fX_sm_tf;
	fY_sm_tf = fSteppingAction->fY_sm_tf;
	fZ_sm_tf = fSteppingAction->fZ_sm_tf;
	fTh_sm_tf = fSteppingAction->fTh_sm_tf;
	fPh_sm_tf = fSteppingAction->fPh_sm_tf;

	fX_sex = fSteppingAction->fX_sex;
	fY_sex = fSteppingAction->fY_sex;
	fZ_sex = fSteppingAction->fZ_sex;
	fTh_sex = fSteppingAction->fTh_sex;
	fPh_sex = fSteppingAction->fPh_sex;

	fX_sex_tr = fSteppingAction->fX_sex_tr;
	fY_sex_tr = fSteppingAction->fY_sex_tr;
	fZ_sex_tr = fSteppingAction->fZ_sex_tr;
	fTh_sex_tr = fSteppingAction->fTh_sex_tr;
	fPh_sex_tr = fSteppingAction->fPh_sex_tr;

	fX_sex_tf = fSteppingAction->fX_sex_tf;
	fY_sex_tf = fSteppingAction->fY_sex_tf;
	fZ_sex_tf = fSteppingAction->fZ_sex_tf;
	fTh_sex_tf = fSteppingAction->fTh_sex_tf;
	fPh_sex_tf = fSteppingAction->fPh_sex_tf;

	fX_col = fSteppingAction->fX_col;
	fY_col = fSteppingAction->fY_col;
	fZ_col = fSteppingAction->fZ_col;
	fTh_col = fSteppingAction->fTh_col;
	fPh_col = fSteppingAction->fPh_col;

	fX_col_tr = fSteppingAction->fX_col_tr;
	fY_col_tr = fSteppingAction->fY_col_tr;
	fZ_col_tr = fSteppingAction->fZ_col_tr;
	fTh_col_tr = fSteppingAction->fTh_col_tr;
	fPh_col_tr = fSteppingAction->fPh_col_tr;

	fX_col_tf = fSteppingAction->fX_col_tf;
	fY_col_tf = fSteppingAction->fY_col_tf;
	fZ_col_tf = fSteppingAction->fZ_col_tf;
	fTh_col_tf = fSteppingAction->fTh_col_tf;
	fPh_col_tf = fSteppingAction->fPh_col_tf;

	fX_q1en = fSteppingAction->fX_q1en;
	fY_q1en = fSteppingAction->fY_q1en;
	fZ_q1en = fSteppingAction->fZ_q1en;
	fTh_q1en = fSteppingAction->fTh_q1en;
	fPh_q1en = fSteppingAction->fPh_q1en;

	fX_q1en_tr = fSteppingAction->fX_q1en_tr;
	fY_q1en_tr = fSteppingAction->fY_q1en_tr;
	fZ_q1en_tr = fSteppingAction->fZ_q1en_tr;
	fTh_q1en_tr = fSteppingAction->fTh_q1en_tr;
	fPh_q1en_tr = fSteppingAction->fPh_q1en_tr;

	fX_q1ex = fSteppingAction->fX_q1ex;
	fY_q1ex = fSteppingAction->fY_q1ex;
	fZ_q1ex = fSteppingAction->fZ_q1ex;
	fTh_q1ex = fSteppingAction->fTh_q1ex;
	fPh_q1ex = fSteppingAction->fPh_q1ex;

	fX_q1ex_tr = fSteppingAction->fX_q1ex_tr;
	fY_q1ex_tr = fSteppingAction->fY_q1ex_tr;
	fZ_q1ex_tr = fSteppingAction->fZ_q1ex_tr;
	fTh_q1ex_tr = fSteppingAction->fTh_q1ex_tr;
	fPh_q1ex_tr = fSteppingAction->fPh_q1ex_tr;

	fX_q1ex_tf = fSteppingAction->fX_q1ex_tf;
	fY_q1ex_tf = fSteppingAction->fY_q1ex_tf;
	fZ_q1ex_tf = fSteppingAction->fZ_q1ex_tf;
	fTh_q1ex_tf = fSteppingAction->fTh_q1ex_tf;
	fPh_q1ex_tf = fSteppingAction->fPh_q1ex_tf;

	fX_q2en = fSteppingAction->fX_q2en;
	fY_q2en = fSteppingAction->fY_q2en;
	fZ_q2en = fSteppingAction->fZ_q2en;
	fTh_q2en = fSteppingAction->fTh_q2en;
	fPh_q2en = fSteppingAction->fPh_q2en;

	fX_q2en_tr = fSteppingAction->fX_q2en_tr;
	fY_q2en_tr = fSteppingAction->fY_q2en_tr;
	fZ_q2en_tr = fSteppingAction->fZ_q2en_tr;
	fTh_q2en_tr = fSteppingAction->fTh_q2en_tr;
	fPh_q2en_tr = fSteppingAction->fPh_q2en_tr;

	fX_q2ex = fSteppingAction->fX_q2ex;
	fY_q2ex = fSteppingAction->fY_q2ex;
	fZ_q2ex = fSteppingAction->fZ_q2ex;
	fTh_q2ex = fSteppingAction->fTh_q2ex;
	fPh_q2ex = fSteppingAction->fPh_q2ex;

	fX_q2ex_tr = fSteppingAction->fX_q2ex_tr;
	fY_q2ex_tr = fSteppingAction->fY_q2ex_tr;
	fZ_q2ex_tr = fSteppingAction->fZ_q2ex_tr;
	fTh_q2ex_tr = fSteppingAction->fTh_q2ex_tr;
	fPh_q2ex_tr = fSteppingAction->fPh_q2ex_tr;

	fX_q2ex_tf = fSteppingAction->fX_q2ex_tf;
	fY_q2ex_tf = fSteppingAction->fY_q2ex_tf;
	fZ_q2ex_tf = fSteppingAction->fZ_q2ex_tf;
	fTh_q2ex_tf = fSteppingAction->fTh_q2ex_tf;
	fPh_q2ex_tf = fSteppingAction->fPh_q2ex_tf;

	fX_den = fSteppingAction->fX_den;
	fY_den = fSteppingAction->fY_den;
	fZ_den = fSteppingAction->fZ_den;
	fTh_den = fSteppingAction->fTh_den;
	fPh_den = fSteppingAction->fPh_den;

	fX_den_tr = fSteppingAction->fX_den_tr;
	fY_den_tr = fSteppingAction->fY_den_tr;
	fZ_den_tr = fSteppingAction->fZ_den_tr;
	fTh_den_tr = fSteppingAction->fTh_den_tr;
	fPh_den_tr = fSteppingAction->fPh_den_tr;

	fX_den_tf = fSteppingAction->fX_den_tf;
	fY_den_tf = fSteppingAction->fY_den_tf;
	fZ_den_tf = fSteppingAction->fZ_den_tf;
	fTh_den_tf = fSteppingAction->fTh_den_tf;
	fPh_den_tf = fSteppingAction->fPh_den_tf;

	fX_dex = fSteppingAction->fX_dex;
	fY_dex = fSteppingAction->fY_dex;
	fZ_dex = fSteppingAction->fZ_dex;
	fTh_dex = fSteppingAction->fTh_dex;
	fPh_dex = fSteppingAction->fPh_dex;

	fX_dex_tr = fSteppingAction->fX_dex_tr;
	fY_dex_tr = fSteppingAction->fY_dex_tr;
	fZ_dex_tr = fSteppingAction->fZ_dex_tr;
	fTh_dex_tr = fSteppingAction->fTh_dex_tr;
	fPh_dex_tr = fSteppingAction->fPh_dex_tr;

	fX_dex_tf = fSteppingAction->fX_dex_tf;
	fY_dex_tf = fSteppingAction->fY_dex_tf;
	fZ_dex_tf = fSteppingAction->fZ_dex_tf;
	fTh_dex_tf = fSteppingAction->fTh_dex_tf;
	fPh_dex_tf = fSteppingAction->fPh_dex_tf;

	fX_q3en = fSteppingAction->fX_q3en;
	fY_q3en = fSteppingAction->fY_q3en;
	fZ_q3en = fSteppingAction->fZ_q3en;
	fTh_q3en = fSteppingAction->fTh_q3en;
	fPh_q3en = fSteppingAction->fPh_q3en;

	fX_q3en_tr = fSteppingAction->fX_q3en_tr;
	fY_q3en_tr = fSteppingAction->fY_q3en_tr;
	fZ_q3en_tr = fSteppingAction->fZ_q3en_tr;
	fTh_q3en_tr = fSteppingAction->fTh_q3en_tr;
	fPh_q3en_tr = fSteppingAction->fPh_q3en_tr;

	fX_q3en_tf = fSteppingAction->fX_q3en_tf;
	fY_q3en_tf = fSteppingAction->fY_q3en_tf;
	fZ_q3en_tf = fSteppingAction->fZ_q3en_tf;
	fTh_q3en_tf = fSteppingAction->fTh_q3en_tf;
	fPh_q3en_tf = fSteppingAction->fPh_q3en_tf;

	fX_q3ex = fSteppingAction->fX_q3ex;
	fY_q3ex = fSteppingAction->fY_q3ex;
	fZ_q3ex = fSteppingAction->fZ_q3ex;
	fTh_q3ex = fSteppingAction->fTh_q3ex;
	fPh_q3ex = fSteppingAction->fPh_q3ex;

	fX_q3ex_tr = fSteppingAction->fX_q3ex_tr;
	fY_q3ex_tr = fSteppingAction->fY_q3ex_tr;
	fZ_q3ex_tr = fSteppingAction->fZ_q3ex_tr;
	fTh_q3ex_tr = fSteppingAction->fTh_q3ex_tr;
	fPh_q3ex_tr = fSteppingAction->fPh_q3ex_tr;

	fX_q3ex_tf = fSteppingAction->fX_q3ex_tf;
	fY_q3ex_tf = fSteppingAction->fY_q3ex_tf;
	fZ_q3ex_tf = fSteppingAction->fZ_q3ex_tf;
	fTh_q3ex_tf = fSteppingAction->fTh_q3ex_tf;
	fPh_q3ex_tf = fSteppingAction->fPh_q3ex_tf;

	fX_vdc = fSteppingAction->fX_vdc;
	fY_vdc = fSteppingAction->fY_vdc;
	fZ_vdc = fSteppingAction->fZ_vdc;
	fTh_vdc = fSteppingAction->fTh_vdc;
	fPh_vdc = fSteppingAction->fPh_vdc;

	fX_vdc_tr = fSteppingAction->fX_vdc_tr;
	fY_vdc_tr = fSteppingAction->fY_vdc_tr;
	fZ_vdc_tr = fSteppingAction->fZ_vdc_tr;
	fTh_vdc_tr = fSteppingAction->fTh_vdc_tr;
	fPh_vdc_tr = fSteppingAction->fPh_vdc_tr;

	fX_fp = fSteppingAction->fX_fp;
	fY_fp = fSteppingAction->fY_fp;
	fZ_fp = fSteppingAction->fZ_fp;
	fTh_fp = fSteppingAction->fTh_fp;
	fPh_fp = fSteppingAction->fPh_fp;

	fX_fp_tr = fSteppingAction->fX_fp_tr;
	fY_fp_tr = fSteppingAction->fY_fp_tr;
	fZ_fp_tr = fSteppingAction->fZ_fp_tr;
	fTh_fp_tr = fSteppingAction->fTh_fp_tr;
	fPh_fp_tr = fSteppingAction->fPh_fp_tr;

	fX_fp_tf = fSteppingAction->fX_fp_tf;
	fY_fp_tf = fSteppingAction->fY_fp_tf;
	fZ_fp_tf = fSteppingAction->fZ_fp_tf;
	fTh_fp_tf = fSteppingAction->fTh_fp_tf;
	fPh_fp_tf = fSteppingAction->fPh_fp_tf;


}

/*---------------------------------------------------------------------------------*/

void g4hrsIO::GrabGDMLFiles(G4String fn){
    // Reset list
    fGDMLFileNames.clear();

    g4hrsRunData *rundata = g4hrsRun::GetRun()->GetData();
    rundata->ClearGDMLFiles();

    xercesc::XMLPlatformUtils::Initialize();
    SearchGDMLforFiles(fn);
    xercesc::XMLPlatformUtils::Terminate();


    // Store filename

    unsigned int idx;

    // Copy into buffers
    for( idx = 0; idx < fGDMLFileNames.size(); idx++ ){
	G4cout << "Found GDML file " << fGDMLFileNames[idx] << G4endl;
	rundata->AddGDMLFile(fGDMLFileNames[idx]);
    }

    return;
}

void g4hrsIO::SearchGDMLforFiles(G4String fn){
    /*!  Chase down files to be included by GDML.
     *   Mainly look for file tags and perform recursively */

    struct stat thisfile;

    int ret = stat(fn.data(), &thisfile);

    if( ret != 0 ){
	G4cerr << "ERROR opening file " << fn <<  " in " << __PRETTY_FUNCTION__ << ": " << strerror(errno) << G4endl;
	exit(1);
    }

   xercesc::XercesDOMParser *xmlParser = new xercesc::XercesDOMParser();

   // Make sure file exists - otherwise freak out

   fGDMLFileNames.push_back(fn.data());

   xmlParser->parse( fn.data() );
   xercesc::DOMDocument* xmlDoc = xmlParser->getDocument();

   xercesc::DOMElement* elementRoot = xmlDoc->getDocumentElement();

   TraverseChildren( elementRoot );
   return;
}

void g4hrsIO::TraverseChildren( xercesc::DOMElement *thisel ){

   xercesc::DOMNodeList*      children = thisel->getChildNodes();
   const XMLSize_t nodeCount = children->getLength();

   for( XMLSize_t xx = 0; xx < nodeCount; ++xx ){
       xercesc::DOMNode* currentNode = children->item(xx);
       if( currentNode->getNodeType() ){   // true is not NULL

	   if( currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE ){ // is element 
	       xercesc::DOMElement* currentElement
		   = dynamic_cast< xercesc::DOMElement* >( currentNode );
	       if( xercesc::XMLString::equals(currentElement->getTagName(), xercesc::XMLString::transcode("file"))){
		   SearchGDMLforFiles(G4String(xercesc::XMLString::transcode(currentElement->getAttribute(xercesc::XMLString::transcode("name")))));
	       }

	       if( currentElement->getChildNodes()->getLength() > 0 ){
		   TraverseChildren( currentElement );
	       }
	   }
       }
   }

}









