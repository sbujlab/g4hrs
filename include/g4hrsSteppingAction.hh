
#ifndef __REMOLLSTEPPINGACTION_HH
#define __REMOLLSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class g4hrsSteppingAction : public G4UserSteppingAction
{
  public:
    g4hrsSteppingAction();
    virtual ~g4hrsSteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

    void SetEnableKryptonite(G4bool k){ fEnableKryptonite = k; }

  private:
    G4bool drawFlag;

    G4bool fEnableKryptonite;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif//__REMOLLSTEPPINGACTION_HH
