#include "g4hrsPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "g4hrsIO.hh"
#include "g4hrsVEventGen.hh"
#include "g4hrsEvent.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"
#include "g4hrstypes.hh"
#include "globals.hh"

#include "g4hrsGenNuclElastic.hh"
#include "g4hrsGenBeam.hh"
#include "g4hrsGenFlat.hh"

g4hrsPrimaryGeneratorAction::g4hrsPrimaryGeneratorAction() {
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);


    fDefaultEvent = new g4hrsEvent();
    fDefaultEvent->ProduceNewParticle(
        G4ThreeVector(0.*cm,0.*cm,-100.*cm),
        G4ThreeVector(0.0,0.0, gDefaultBeamE),
        "e-" );

    double kinE = sqrt(fDefaultEvent->fPartMom[0].mag()*fDefaultEvent->fPartMom[0].mag()
                       + fDefaultEvent->fPartType[0]->GetPDGMass()*fDefaultEvent->fPartType[0]->GetPDGMass() )
                  -  fDefaultEvent->fPartType[0]->GetPDGMass();

    // Default generator data
    fParticleGun->SetParticleDefinition(fDefaultEvent->fPartType[0]);
    fParticleGun->SetParticleMomentumDirection(fDefaultEvent->fPartMom[0].unit());
    fParticleGun->SetParticleEnergy( kinE  );
    fParticleGun->SetParticlePosition( fDefaultEvent->fPartPos[0] );

    fEventGen = NULL;
}

g4hrsPrimaryGeneratorAction::~g4hrsPrimaryGeneratorAction() {
    delete fParticleGun;
    delete fDefaultEvent;
}

void g4hrsPrimaryGeneratorAction::SetGenerator(G4String genname) {

    fEventGen = NULL;

    if( genname == "elastic" ) {
        fEventGen = new g4hrsGenNuclElastic();
    }else if( genname == "beam" ) {
        fEventGen = new g4hrsGenBeam();
    }else if( genname == "flat" ) {
	fEventGen = new g4hrsGenFlat();
    }

    if( !fEventGen ) {
        G4cerr << __FILE__ << " line " << __LINE__ << " - ERROR generator " << genname << " invalid" << G4endl;
        exit(1);
    } else {
        G4cout << "Setting generator to " << genname << G4endl;
    }

    g4hrsRun::GetRun()->GetData()->SetGenName(genname.data());

    return;
}

void g4hrsPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    /*  Generate event, set IO data */

    g4hrsEvent *thisev = NULL;
    if( fEventGen ) { // Specified our own generator
        thisev = fEventGen->GenerateEvent();
        for( unsigned int pidx = 0; pidx < thisev->fPartType.size(); pidx++ ) {

            double kinE = sqrt(thisev->fPartMom[pidx].mag()*thisev->fPartMom[pidx].mag() +
                               thisev->fPartType[pidx]->GetPDGMass()*thisev->fPartType[pidx]->GetPDGMass())
                          -  thisev->fPartType[pidx]->GetPDGMass();

            fParticleGun->SetParticleDefinition(thisev->fPartType[pidx]);
            fParticleGun->SetParticleMomentumDirection(thisev->fPartMom[pidx].unit());
            fParticleGun->SetParticleEnergy( kinE  );
            fParticleGun->SetParticlePosition( thisev->fPartPos[pidx] );

            fParticleGun->GeneratePrimaryVertex(anEvent);
        }

        if( thisev->fPartType.size() > 0 ) {
            fIO->SetEventData(thisev);
        }
    } else { // Use default, static single generator
        // Update this just in case things changed
        // from the command user interface
        fDefaultEvent->Reset();
        fDefaultEvent->ProduceNewParticle(
            fParticleGun->GetParticlePosition(),
            fParticleGun->GetParticleMomentumDirection()*
            fParticleGun->GetParticleMomentum(),
            fParticleGun->GetParticleDefinition()->GetParticleName() );
        fIO->SetEventData(fDefaultEvent);

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    if( thisev != NULL ) {
        delete thisev;
    }
}

G4ParticleGun* g4hrsPrimaryGeneratorAction::GetParticleGun() {
    return fParticleGun;
}

