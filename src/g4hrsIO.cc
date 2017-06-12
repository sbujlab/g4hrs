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
	fTree->Branch("x_sen", 	&fX_sen,	"x_sen/D");
	fTree->Branch("y_sen", 	&fY_sen,	"y_sen/D");
	fTree->Branch("z_sen", 	&fZ_sen,	"z_sen/D");
	fTree->Branch("p_sen", 	&fP_sen,	"p_sen/D");
	fTree->Branch("th_sen", &fTheta_sen,	"th_sen/D");
	fTree->Branch("ph_sen", &fPhi_sen,	"ph_sen/D");

	fTree->Branch("x_sm", 	&fX_sm,	"x_sm/D");
	fTree->Branch("y_sm", 	&fY_sm,	"y_sm/D");
	fTree->Branch("z_sm", 	&fZ_sm,	"z_sm/D");
	fTree->Branch("p_sm", 	&fP_sm,	"p_sm/D");
	fTree->Branch("th_sm", &fTheta_sm,	"th_sm/D");
	fTree->Branch("ph_sm", &fPhi_sm,	"ph_sm/D");

	fTree->Branch("x_sex", 	&fX_sex,	"x_sex/D");
	fTree->Branch("y_sex", 	&fY_sex,	"y_sex/D");
	fTree->Branch("z_sex", 	&fZ_sex,	"z_sex/D");
	fTree->Branch("p_sex", 	&fP_sex,	"p_sex/D");
	fTree->Branch("th_sex", &fTheta_sex,	"th_sex/D");
	fTree->Branch("ph_sex", &fPhi_sex,	"ph_sex/D");

	fTree->Branch("x_coil", 	&fX_coil,	"x_coil/D");
	fTree->Branch("y_coil", 	&fY_coil,	"y_coil/D");
	fTree->Branch("z_coil", 	&fZ_coil,	"z_coil/D");
	fTree->Branch("p_coil", 	&fP_coil,	"p_coil/D");
	fTree->Branch("th_coil", &fTheta_coil,	"th_coil/D");
	fTree->Branch("ph_coil", &fPhi_coil,	"ph_coil/D");

	fTree->Branch("x_mid", 	&fX_mid,	"x_mid/D");
	fTree->Branch("y_mid", 	&fY_mid,	"y_mid/D");
	fTree->Branch("z_mid", 	&fZ_mid,	"z_mid/D");
	fTree->Branch("p_mid", 	&fP_mid,	"p_mid/D");
	fTree->Branch("th_mid", &fTheta_mid,	"th_mid/D");
	fTree->Branch("ph_mid", &fPhi_mid,	"ph_mid/D");

	fTree->Branch("x_col", 	&fX_col,	"x_col/D");
	fTree->Branch("y_col", 	&fY_col,	"y_col/D");
	fTree->Branch("z_col", 	&fZ_col,	"z_col/D");
	fTree->Branch("p_col", 	&fP_col,	"p_col/D");
	fTree->Branch("th_col", &fTheta_col,	"th_col/D");
	fTree->Branch("ph_col", &fPhi_col,	"ph_col/D");

	fTree->Branch("x_q1en", 	&fX_q1en,	"x_q1en/D");
	fTree->Branch("y_q1en", 	&fY_q1en,	"y_q1en/D");
	fTree->Branch("z_q1en", 	&fZ_q1en,	"z_q1en/D");
	fTree->Branch("p_q1en", 	&fP_q1en,	"p_q1en/D");
	fTree->Branch("th_q1en", &fTheta_q1en,	"th_q1en/D");
	fTree->Branch("ph_q1en", &fPhi_q1en,	"ph_q1en/D");

	fTree->Branch("x_q1ex", 	&fX_q1ex,	"x_q1ex/D");
	fTree->Branch("y_q1ex", 	&fY_q1ex,	"y_q1ex/D");
	fTree->Branch("z_q1ex", 	&fZ_q1ex,	"z_q1ex/D");
	fTree->Branch("p_q1ex", 	&fP_q1ex,	"p_q1ex/D");
	fTree->Branch("th_q1ex", &fTheta_q1ex,	"th_q1ex/D");
	fTree->Branch("ph_q1ex", &fPhi_q1ex,	"ph_q1ex/D");

	fTree->Branch("x_q2en", 	&fX_q2en,	"x_q2en/D");
	fTree->Branch("y_q2en", 	&fY_q2en,	"y_q2en/D");
	fTree->Branch("z_q2en", 	&fZ_q2en,	"z_q2en/D");
	fTree->Branch("p_q2en", 	&fP_q2en,	"p_q2en/D");
	fTree->Branch("th_q2en", &fTheta_q2en,	"th_q2en/D");
	fTree->Branch("ph_q2en", &fPhi_q2en,	"ph_q2en/D");

	fTree->Branch("x_q2ex", 	&fX_q2ex,	"x_q2ex/D");
	fTree->Branch("y_q2ex", 	&fY_q2ex,	"y_q2ex/D");
	fTree->Branch("z_q2ex", 	&fZ_q2ex,	"z_q2ex/D");
	fTree->Branch("p_q2ex", 	&fP_q2ex,	"p_q2ex/D");
	fTree->Branch("th_q2ex", &fTheta_q2ex,	"th_q2ex/D");
	fTree->Branch("ph_q2ex", &fPhi_q2ex,	"ph_q2ex/D");

	fTree->Branch("x_den", 	&fX_den,	"x_den/D");
	fTree->Branch("y_den", 	&fY_den,	"y_den/D");
	fTree->Branch("z_den", 	&fZ_den,	"z_den/D");
	fTree->Branch("p_den", 	&fP_den,	"p_den/D");
	fTree->Branch("th_den", &fTheta_den,	"th_den/D");
	fTree->Branch("ph_den", &fPhi_den,	"ph_den/D");

	fTree->Branch("x_dex", 	&fX_dex,	"x_dex/D");
	fTree->Branch("y_dex", 	&fY_dex,	"y_dex/D");
	fTree->Branch("z_dex", 	&fZ_dex,	"z_dex/D");
	fTree->Branch("p_dex", 	&fP_dex,	"p_dex/D");
	fTree->Branch("th_dex", &fTheta_dex,	"th_dex/D");
	fTree->Branch("ph_dex", &fPhi_dex,	"ph_dex/D");

	fTree->Branch("x_q3en", 	&fX_q3en,	"x_q3en/D");
	fTree->Branch("y_q3en", 	&fY_q3en,	"y_q3en/D");
	fTree->Branch("z_q3en", 	&fZ_q3en,	"z_q3en/D");
	fTree->Branch("p_q3en", 	&fP_q3en,	"p_q3en/D");
	fTree->Branch("th_q3en", &fTheta_q3en,	"th_q3en/D");
	fTree->Branch("ph_q3en", &fPhi_q3en,	"ph_q3en/D");

	fTree->Branch("x_q3ex", 	&fX_q3ex,	"x_q3ex/D");
	fTree->Branch("y_q3ex", 	&fY_q3ex,	"y_q3ex/D");
	fTree->Branch("z_q3ex", 	&fZ_q3ex,	"z_q3ex/D");
	fTree->Branch("p_q3ex", 	&fP_q3ex,	"p_q3ex/D");
	fTree->Branch("th_q3ex", &fTheta_q3ex,	"th_q3ex/D");
	fTree->Branch("ph_q3ex", &fPhi_q3ex,	"ph_q3ex/D");

	fTree->Branch("x_vdc", 	&fX_vdc,	"x_vdc/D");
	fTree->Branch("y_vdc", 	&fY_vdc,	"y_vdc/D");
	fTree->Branch("z_vdc", 	&fZ_vdc,	"z_vdc/D");
	fTree->Branch("p_vdc", 	&fP_vdc,	"p_vdc/D");
	fTree->Branch("th_vdc", &fTheta_vdc,	"th_vdc/D");
	fTree->Branch("ph_vdc", &fPhi_vdc,	"ph_vdc/D");

	fTree->Branch("x_qz1", 	&fX_qz1,	"x_qz1/D");
	fTree->Branch("y_qz1", 	&fY_qz1,	"y_qz1/D");
	fTree->Branch("z_qz1", 	&fZ_qz1,	"z_qz1/D");
	fTree->Branch("p_qz1", 	&fP_qz1,	"p_qz1/D");
	fTree->Branch("th_qz1", &fTheta_qz1,	"th_qz1/D");
	fTree->Branch("ph_qz1", &fPhi_qz1,	"ph_qz1/D");

	fTree->Branch("x_qz2", 	&fX_qz2,	"x_qz2/D");
	fTree->Branch("y_qz2", 	&fY_qz2,	"y_qz2/D");
	fTree->Branch("z_qz2", 	&fZ_qz2,	"z_qz2/D");
	fTree->Branch("p_qz2", 	&fP_qz2,	"p_qz2/D");
	fTree->Branch("th_qz2", &fTheta_qz2,	"th_qz2/D");
	fTree->Branch("ph_qz2", &fPhi_qz2,	"ph_qz2/D");

	fTree->Branch("x_fp", 	&fX_fp,	"x_fp/D");
	fTree->Branch("y_fp", 	&fY_fp,	"y_fp/D");
	fTree->Branch("z_fp", 	&fZ_fp,	"z_fp/D");
	fTree->Branch("p_fp", 	&fP_fp,	"p_fp/D");
	fTree->Branch("th_fp", &fTheta_fp,	"th_fp/D");
	fTree->Branch("ph_fp", &fPhi_fp,	"ph_fp/D");

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

	fSteppingAction->fX_sen = -150.;
	fSteppingAction->fY_sen = -150.;
	fSteppingAction->fZ_sen = -150.;
	fSteppingAction->fP_sen = -150.;
	fSteppingAction->fTheta_sen = -150.;
	fSteppingAction->fPhi_sen = -150.;

}

void g4hrsIO::SetVirtualBoundaryData() {

	fX_sen = fSteppingAction->fX_sen;
	fY_sen = fSteppingAction->fY_sen;
	fZ_sen = fSteppingAction->fZ_sen;
	fP_sen = fSteppingAction->fP_sen;
	fTheta_sen = fSteppingAction->fTheta_sen;
	fPhi_sen = fSteppingAction->fPhi_sen;

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









