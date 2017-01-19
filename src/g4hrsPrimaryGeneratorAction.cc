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

#include "g4hrsGenMoller.hh"
#include "g4hrsGenpElastic.hh"
#include "g4hrsGenpInelastic.hh"
#include "g4hrsGenPion.hh"
#include "g4hrsGenBeam.hh"
#include "g4hrsGenFlat.hh"
#include "g4hrsGenAl.hh"
#include "g4hrsGenLUND.hh" //Dominic Lunde adding the LUND generator command

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

    if( genname == "moller" ) {
        fEventGen = new g4hrsGenMoller();
    }else if( genname == "elastic" ) {
        fEventGen = new g4hrsGenpElastic();
    }else if( genname == "inelastic" ) {
        fEventGen = new g4hrsGenpInelastic();
    }else if( genname == "pion" ) {
        fEventGen = new g4hrsGenPion();
    }else if( genname == "beam" ) {
        fEventGen = new g4hrsGenBeam();
    }else if( genname == "flat" ) {
        fEventGen = new g4hrsGenFlat();
    }else if( genname == "inelasticAl" ) {
        fEventGen = new g4hrsGenAl(2);
    }else if( genname == "quasielasticAl" ) {
        fEventGen = new g4hrsGenAl(1);
    }else if( genname == "elasticAl" ) {
        fEventGen = new g4hrsGenAl(0);
    }else if( genname == "pion_LUND" ) { //Dominic Lunde - adding GenLUND into the generators
        fEventGen = new g4hrsGenLUND();  
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

