
#ifndef g4hrsPrimaryGeneratorAction_h
#define g4hrsPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"

class G4ParticleGun;
class G4Event;
class g4hrsIO;
class g4hrsVEventGen;
class g4hrsEvent;

class g4hrsPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    g4hrsPrimaryGeneratorAction();
    ~g4hrsPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    G4ParticleGun* GetParticleGun();
    void SetIO( g4hrsIO *io ){ fIO = io; }

    void SetGenerator( G4String );

    g4hrsVEventGen *GetGenerator(){ return fEventGen; }

  private:
    G4ParticleGun* fParticleGun;

    g4hrsVEventGen *fEventGen;
    g4hrsEvent *fDefaultEvent;
    g4hrsIO *fIO;
};

#endif


