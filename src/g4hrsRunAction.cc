
// Make this appear first!
#include "G4Timer.hh"

#include "g4hrsRunAction.hh"
#include "g4hrsBeamTarget.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "g4hrsIO.hh"
#include "g4hrsRun.hh"

g4hrsRunAction::g4hrsRunAction()
{
  timer = new G4Timer;
}

g4hrsRunAction::~g4hrsRunAction()
{
  delete timer;
}

void g4hrsRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //  timer->Start();
  fIO->InitializeTree();

  g4hrsRunData *rmrundata = g4hrsRun::GetRun()->GetData();

  rmrundata->SetBeamE( g4hrsBeamTarget::GetBeamTarget()->fBeamE/GeV );
  rmrundata->SetNthrown( aRun->GetNumberOfEventToBeProcessed() );

  rmrundata->Print();
}

void g4hrsRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  //       << " " << *timer << G4endl;

  fIO->WriteTree();
}

