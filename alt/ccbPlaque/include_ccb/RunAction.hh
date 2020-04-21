#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    void SumDose(G4double dose) { fSumDose += dose; };

    void AddEdep_0 (G4double edep_0);
    void AddEdep_1 (G4double edep_1);
    void AddEdep_2 (G4double edep_2);
    void AddEdep_3 (G4double edep_3);
    void AddEdep_4 (G4double edep_4);
    void AddEdep_5 (G4double edep_5);
    void AddEdep_6 (G4double edep_6);
    void AddEdep_7 (G4double edep_7);
    void AddEdep_8 (G4double edep_8);
    void AddEdep_9 (G4double edep_9);
    void AddEdep_10 (G4double edep_10);
    void AddEdep_11 (G4double edep_11);
    void AddEdep_12 (G4double edep_12);
    void AddEdep_13 (G4double edep_13);
    void AddEdep_14 (G4double edep_14);
    void AddEdep_15 (G4double edep_15);
    void AddEdep_16 (G4double edep_16);
    void AddEdep_17 (G4double edep_17);
    void AddEdep_18 (G4double edep_18);
    void AddEdep_19 (G4double edep_19);
    void AddEdep_20 (G4double edep_20);
    void AddEdep_21 (G4double edep_21);
    void AddEdep_22 (G4double edep_22);
    void AddEdep_23 (G4double edep_23);
    void AddEdep_24 (G4double edep_24);
    void AddEdep_25 (G4double edep_25);
    void AddEdep_26 (G4double edep_26);
    void AddEdep_27 (G4double edep_27);
    void AddEdep_28 (G4double edep_28);
    void AddEdep_29 (G4double edep_29);
    void AddEdep_30 (G4double edep_30);
    void AddEdep_31 (G4double edep_31);
    void AddEdep_32 (G4double edep_32);
    void AddEdep_33 (G4double edep_33);
    void AddEdep_34 (G4double edep_34);

private:
    G4Accumulable<G4double> fSumDose; 
    G4Accumulable<G4double> fEdep_0,
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
