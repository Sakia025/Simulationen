/* GEANT4 application for calculating a depth dose profile of the CCB-type ruthenium eye plaque
 * Geometrical description from ICRU 72 Report
 * Physics list based on GEANT4 example "medical/electronScattering2/src/PhysicsList.cc"
 * Developed with GEANT4 version 10.2 patch1 http://geant4.web.cern.ch/geant4/license/LICENSE.html
 */

#ifdef G4MULTITHREADED
    #include "G4MTRunManager.hh"
#else
    #include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"

#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

//#include "/Users/smller/Simulationen/ccbPlaque/include/global_variables.hh"

int main(int argc, char** argv)
{
    G4UIExecutive* ui = 0;

    if (argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
    }
//!!!Ändern zu Set seed über runfile!!
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4int seconds =  time(NULL);
    G4Random::setTheSeed(seconds);



#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
//!!Was sind number of threads? parallelisierung
    runManager->SetNumberOfThreads(4);
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new PhysicsList);
    runManager->SetUserInitialization(new ActionInitialization());

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

//!! brauch ich das? kommt aus example B3a Activate score ntuple writer
    // The Root output type (Root) is selected in B3Analysis.hh.
    // The verbose level can be also set via UI commands
    // /score/ntuple/writerVerbose level
  //  G4TScoreNtupleWriter<G4AnalysisManager> scoreNtupleWriter;
  //  scoreNtupleWriter.SetVerboseLevel(1);


    if ( !ui )
    {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
    else
    {
        //interactive mode
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }

    delete visManager;
    delete runManager;
}
