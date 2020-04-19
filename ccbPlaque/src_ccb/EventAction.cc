#include "G4Event.hh"
#include "G4RunManager.hh"

#include "EventAction.hh"
#include "RunAction.hh"


int idontknow = 7;
extern std::vector<std::string> model_names;

EventAction::EventAction(RunAction* runAction)
    : G4UserEventAction(),
      fRunAction(runAction),
      fCollID_patient(-1),
      fEdep_0(0.),
      fEdep_1(0.),
      fEdep_2(0.),
      fEdep_3(0.),
      fEdep_4(0.),
      fEdep_5(0.),
      fEdep_6(0.),
      fEdep_7(0.),
      fEdep_8(0.),
      fEdep_9(0.),
      fEdep_10(0.),
      fEdep_11(0.),
      fEdep_12(0.),
      fEdep_13(0.),
      fEdep_14(0.),
      fEdep_15(0.),
      fEdep_16(0.),
      fEdep_17(0.),
      fEdep_18(0.),
      fEdep_19(0.),
      fEdep_20(0.),
      fEdep_21(0.),
      fEdep_22(0.),
      fEdep_23(0.),
      fEdep_24(0.),
      fEdep_25(0.),
      fEdep_26(0.),
      fEdep_27(0.),
      fEdep_28(0.),
      fEdep_29(0.),
      fEdep_30(0.),
      fEdep_31(0.),
      fEdep_32(0.),
      fEdep_33(0.),
      fEdep_34(0.)
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{


    fEdep_0 = 0.;
    fEdep_1 = 0.;
    fEdep_2 = 0.;
    fEdep_3 = 0.;
    fEdep_4 = 0.;
    fEdep_5 = 0.;
    fEdep_6 = 0.;
    fEdep_7 = 0.;
    fEdep_8 = 0.;
    fEdep_9 = 0.;
    fEdep_10 = 0.;
    fEdep_11 = 0.;
    fEdep_12 = 0.;
    fEdep_13 = 0.;
    fEdep_14 = 0.;
    fEdep_15 = 0.;
    fEdep_16 = 0.;
    fEdep_17 = 0.;
    fEdep_18 = 0.;
    fEdep_19 = 0.;
    fEdep_20 = 0.;
    fEdep_21 = 0.;
    fEdep_22 = 0.;
    fEdep_23 = 0.;
    fEdep_24 = 0.;
    fEdep_25 = 0.;
    fEdep_26 = 0.;
    fEdep_27 = 0.;
    fEdep_28 = 0.;
    fEdep_29 = 0.;
    fEdep_30 = 0.;
    fEdep_31 = 0.;
    fEdep_32 = 0.;
    fEdep_33 = 0.;
    fEdep_34 = 0.;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
    //Hits collections
    //
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    if(!HCE) return;

    // Get hits collections IDs
  //  if (fCollID_cryst < 0) {
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();
    fCollID_patient = SDMan->GetCollectionID("my_" + std::string(model_names[0]) + "_Scorer/TotalDose");
  //  }

    //Dose deposit in patient
    //
    G4double dose = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_patient));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
      ///G4int copyNb  = (itr->first);
      dose = *(itr->second);
    }
    if (dose > 0.) fRunAction->SumDose(dose);

    fRunAction->AddEdep_0(fEdep_0);
    fRunAction->AddEdep_1(fEdep_1);
    fRunAction->AddEdep_2(fEdep_2);
    fRunAction->AddEdep_3(fEdep_3);
    fRunAction->AddEdep_4(fEdep_4);
    fRunAction->AddEdep_5(fEdep_5);
    fRunAction->AddEdep_6(fEdep_6);
    fRunAction->AddEdep_7(fEdep_7);
    fRunAction->AddEdep_8(fEdep_8);
    fRunAction->AddEdep_9(fEdep_9);
    fRunAction->AddEdep_10(fEdep_10);
    fRunAction->AddEdep_11(fEdep_11);
    fRunAction->AddEdep_12(fEdep_12);
    fRunAction->AddEdep_13(fEdep_13);
    fRunAction->AddEdep_14(fEdep_14);
    fRunAction->AddEdep_15(fEdep_15);
    fRunAction->AddEdep_16(fEdep_16);
    fRunAction->AddEdep_17(fEdep_17);
    fRunAction->AddEdep_18(fEdep_18);
    fRunAction->AddEdep_19(fEdep_19);
    fRunAction->AddEdep_20(fEdep_20);
    fRunAction->AddEdep_21(fEdep_21);
    fRunAction->AddEdep_22(fEdep_22);
    fRunAction->AddEdep_23(fEdep_23);
    fRunAction->AddEdep_24(fEdep_24);
    fRunAction->AddEdep_25(fEdep_25);
    fRunAction->AddEdep_26(fEdep_26);
    fRunAction->AddEdep_27(fEdep_27);
    fRunAction->AddEdep_28(fEdep_28);
    fRunAction->AddEdep_29(fEdep_29);
    fRunAction->AddEdep_30(fEdep_30);
    fRunAction->AddEdep_31(fEdep_31);
    fRunAction->AddEdep_32(fEdep_32);
    fRunAction->AddEdep_33(fEdep_33);
    fRunAction->AddEdep_34(fEdep_34);
}
