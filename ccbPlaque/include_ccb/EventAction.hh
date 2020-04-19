#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep_0(G4double edep_0)
    {
        fEdep_0 += edep_0;
    }
    void AddEdep_1(G4double edep_1)
    {
        fEdep_1 += edep_1;
    }
    void AddEdep_2(G4double edep_2)
    {
        fEdep_2 += edep_2;
    }
    void AddEdep_3(G4double edep_3)
    {
        fEdep_3 += edep_3;
    }
    void AddEdep_4(G4double edep_4)
    {
        fEdep_4 += edep_4;
    }
    void AddEdep_5(G4double edep_5)
    {
        fEdep_5 += edep_5;
    }
    void AddEdep_6(G4double edep_6)
    {
        fEdep_6 += edep_6;
    }
    void AddEdep_7(G4double edep_7)
    {
        fEdep_7 += edep_7;
    }
    void AddEdep_8(G4double edep_8)
    {
        fEdep_8 += edep_8;
    }
    void AddEdep_9(G4double edep_9)
    {
        fEdep_9 += edep_9;
    }
    void AddEdep_10(G4double edep_10)
    {
        fEdep_10 += edep_10;
    }
    void AddEdep_11(G4double edep_11)
    {
        fEdep_11 += edep_11;
    }
    void AddEdep_12(G4double edep_12)
    {
        fEdep_12 += edep_12;
    }
    void AddEdep_13(G4double edep_13)
    {
        fEdep_13 += edep_13;
    }
    void AddEdep_14(G4double edep_14)
    {
        fEdep_14 += edep_14;
    }
    void AddEdep_15(G4double edep_15)
    {
        fEdep_15 += edep_15;
    }
    void AddEdep_16(G4double edep_16)
    {
        fEdep_16 += edep_16;
    }
    void AddEdep_17(G4double edep_17)
    {
        fEdep_17 += edep_17;
    }
    void AddEdep_18(G4double edep_18)
    {
        fEdep_18 += edep_18;
    }
    void AddEdep_19(G4double edep_19)
    {
        fEdep_19 += edep_19;
    }
    void AddEdep_20(G4double edep_20)
    {
        fEdep_20 += edep_20;
    }
    void AddEdep_21(G4double edep_21)
    {
        fEdep_21 += edep_21;
    }
    void AddEdep_22(G4double edep_22)
    {
        fEdep_22 += edep_22;
    }
    void AddEdep_23(G4double edep_23)
    {
        fEdep_23 += edep_23;
    }
    void AddEdep_24(G4double edep_24)
    {
        fEdep_24 += edep_24;
    }
    void AddEdep_25(G4double edep_25)
    {
        fEdep_25 += edep_25;
    }
    void AddEdep_26(G4double edep_26)
    {
        fEdep_26 += edep_26;
    }
    void AddEdep_27(G4double edep_27)
    {
        fEdep_27 += edep_27;
    }
    void AddEdep_28(G4double edep_28)
    {
        fEdep_28 += edep_28;
    }
    void AddEdep_29(G4double edep_29)
    {
        fEdep_29 += edep_29;
    }
    void AddEdep_30(G4double edep_30)
    {
        fEdep_30 += edep_30;
    }
    void AddEdep_31(G4double edep_31)
    {
        fEdep_31 += edep_31;
    }
    void AddEdep_32(G4double edep_32)
    {
        fEdep_32 += edep_32;
    }
    void AddEdep_33(G4double edep_33)
    {
        fEdep_33 += edep_33;
    }
    void AddEdep_34(G4double edep_34)
    {
        fEdep_34 += edep_34;
    }
private:
    RunAction* fRunAction;
    G4int fCollID_patient;  
    G4double fEdep_0,
             fEdep_1,
             fEdep_2,
             fEdep_3,
             fEdep_4,
             fEdep_5,
             fEdep_6,
             fEdep_7,
             fEdep_8,
             fEdep_9,
             fEdep_10,
             fEdep_11,
             fEdep_12,
             fEdep_13,
             fEdep_14,
             fEdep_15,
             fEdep_16,
             fEdep_17,
             fEdep_18,
             fEdep_19,
             fEdep_20,
             fEdep_21,
             fEdep_22,
             fEdep_23,
             fEdep_24,
             fEdep_25,
             fEdep_26,
             fEdep_27,
             fEdep_28,
             fEdep_29,
             fEdep_30,
             fEdep_31,
             fEdep_32,
             fEdep_33,
             fEdep_34;
};

#endif
