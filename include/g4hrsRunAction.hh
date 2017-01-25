
#ifndef g4hrsRunAction_h
#define g4hrsRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;
class g4hrsIO;

class g4hrsRunAction : public G4UserRunAction
{
  public:
    g4hrsRunAction();
    ~g4hrsRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

    void SetIO( g4hrsIO *io ){ fIO = io; }

  private:
    G4Timer* timer;

    g4hrsIO *fIO;
};

#endif

