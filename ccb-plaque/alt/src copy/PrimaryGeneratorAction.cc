
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4ChargedGeantino.hh"
#include "Randomize.hh"

#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(0)
{
    fParticleGun  = new G4ParticleGun(1);
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle
        = particleTable->FindParticle("chargedgeantino");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
    fParticleGun->SetParticleEnergy(0. * eV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();

    if (particle == G4ChargedGeantino::ChargedGeantino())
    {
        // ruthenium 106
        G4int Z_value = 44, A_value = 106;

        // rhodium 106
        // G4int Z_value = 45, A_value = 103;

        // ruthenium 103
        // G4int Z_value = 44, A_value = 103;

        // Cs 134
        // G4int Z_value = 55, A_value = 137;

        G4double ionCharge   = 0. * eplus;
        G4double excitEnergy = 0. * keV;
        G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z_value, A_value, excitEnergy);

        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(ionCharge);
    }

    G4double r_position = 12.1;
    G4double phi_position = CLHEP::twopi * G4UniformRand();
    G4double theta_position  = acos( 0.322482 * G4UniformRand() + ( 1.0 - 0.322482) );
    // cos(47.35 deg) = 0.677. We need random numbers between 1.0 and 0.677, ( 0.677 is ( 1.0 - 0.322482) )

    G4double x_position = ( r_position * std::sin(theta_position) * std::cos(phi_position) ) * mm;
    G4double y_position = ( r_position * std::sin(theta_position) * std::sin(phi_position) ) * mm;
    G4double z_position = ( r_position * std::cos(theta_position) ) * mm;
    G4ThreeVector position = (G4ThreeVector(x_position, y_position, z_position));

    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
