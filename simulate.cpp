#include <cassert>

#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "particle/particle.hpp"

void simulate() {
	std::cout << "Simulation started successfully.\n";

	double piMass = 0.13957;
  double kMass = 0.49367;
  double PMass = 0.93827;
  double KstarMass = 0.89166;
  double KstarLength = 0.05;
	
	Particle::addParticleType("pi+", piMass, 1);
	int piPlusI = 0;
	Particle::addParticleType("pi-", piMass, -1);
	int piMinI = 1;
	Particle::addParticleType("k+", kMass, 1);
	int kPlusI = 2;
	Particle::addParticleType("k-", kMass, -1);
	int kMinI = 3;
	Particle::addParticleType("P+", PMass, +1);
	int pPlusI = 4;
	Particle::addParticleType("P-", PMass, -1);
	int pMinI = 5;
	Particle::addParticleType("k*", KstarMass, 0, KstarLength);
	int kStarI = 6;

	std::cout << "Successfully added particle types.\n";

	gRandom->SetSeed();

	TH1F* pTypesH = new TH1F("pTypesH", "Particle types", 7, 0.5, 7.5);
	TH2F* anglesH = new TH2F("anglesH", "Angles", 1000, 0., 2*M_PI, 500, 0., M_PI);
	TH1F* momentumH = new TH1F("momentumH", "Momentum norm", 1000, 0., 3.);
	TH1F* tMomentumH = new TH1F("tMomentumH", "Transverse momentum", 1000, 0., 3.);
	TH1F* energyH = new TH1F("energyH", "Energy", 1000, 0., 3.);
	TH1F* invMassTotH = new TH1F("invMassTotH", "Invariant mass", 1000, 0., 3.);
	TH1F* invMassOppSgnH = new TH1F("invMassOppSgnH", "Invariant mass, opp sign", 1000, 0., 3.);
	TH1F* invMassEqSgnH = new TH1F("invMassEqSgnH", "Invariant mass, equal sign", 1000, 0., 3.);
	TH1F* invMassPKOppH = new TH1F("invMassPKOppH", "Invariant mass, opp. sign p and k", 1000, 0., 3.);
	TH1F* invMassPKEqH = new TH1F("invMassPKEqH", "Invariant mass, equal sign p and k", 1000, 0., 3.);
	TH1F* invMassKStarH = new TH1F("invMassKStarH", "Invariant mass, K* decay products", 1000, 0., 3.);

	int nEvents {static_cast<int>(1E5)};
	const int nGens {100}; // number of generated particles
	const int nParticles {nGens * 6 / 5}; // safe array size for 1% likelyhood of k* generation

	// declaring variables for the loops to use
	Particle eventParticles[nParticles];
	
	double phi {};
	double theta {};
	double pNorm {};
	double control {};

	int decayResult = {0};
	Particle dau1;
	Particle dau2;

	std::cout << "Correctly finished pre-loops setup. Starting generation.\n";
	for (int i {}; i < nEvents; ++i) {
		for (Particle& particle: eventParticles) {
			particle = {};
			assert(particle.getIndex() == -1);
		}

		for (int i {}; i < nGens; ++i) {
			phi = gRandom->Uniform(0., 2*M_PI);
			theta = gRandom->Uniform(0., M_PI);
			pNorm = gRandom->Exp(1.);
			eventParticles[i].setMomentum( { sin(phi)*cos(theta)*pNorm, sin(phi)*sin(theta)*pNorm, cos(theta)*pNorm } );
			
			control = gRandom->Uniform(0., 1.);
			if (control < 0.4) eventParticles[i].setIndex(piPlusI); 				// pi+, 40%
			else if (control < 0.8) eventParticles[i].setIndex(piMinI); 		// pi-, 40%
			else if (control < 0.85) eventParticles[i].setIndex(kPlusI); 		// k+, 5%
			else if (control < 0.9) eventParticles[i].setIndex(kMinI); 			// k-, 5%
			else if (control < 0.945) eventParticles[i].setIndex(pPlusI);		// P+, 4.5%
			else if (control < 0.99) eventParticles[i].setIndex(pMinI); 		// P-, 4.5%
			else {
				eventParticles[i].setIndex(kStarI);														// k*, 1%
				if (control < 0.995) {		// pi+ k- decay, 50%
					dau1.setIndex(piPlusI); // pi+
					dau2.setIndex(kMinI); 	// k-
				} else {									//pi- k+ decay, 50%
					dau1.setIndex(piMinI); 	// pi-
					dau2.setIndex(kPlusI); 	// k+
				}
				decayResult = eventParticles[i].decayToBody(dau1, dau2);
				assert((!decayResult));
				// std::cout << "K* generated, decay process result: " << decayResult << std::endl;
				eventParticles[tailIndex(eventParticles, nParticles, nGens)] = dau1;
				eventParticles[tailIndex(eventParticles, nParticles, nGens)] = dau2;
				invMassKStarH->Fill(invMass(dau1, dau2));
			}
			pTypesH->Fill(eventParticles[i].getIndex() + 1.);
			momentumH->Fill(pNorm);
			tMomentumH->Fill(sqrt(pNorm*pNorm*sin(theta)*sin(theta)));
			anglesH->Fill(phi, theta);
			energyH->Fill(eventParticles[i].getEnergy());
		}

		// fill histograms
		int max = tailIndex(eventParticles, nParticles, nGens);
		// std::cout << "Found tailIndex for iteration " << i << ": " << max << "\n";
		for (int i {}; i < max; ++i) {
			if (eventParticles[i].getIndex() != kStarI) {
				//momentumH->Fill(eventParticles[i].getMomentum().norm());
				//tMomentumH->Fill(eventParticles[i].getMomentum().norm2() - eventParticles[i].getPz()*eventParticles[i].getPz());
				for (int j{i+1}; j < max; ++j) {
					if (eventParticles[j].getIndex() != kStarI) {
						invMassTotH->Fill(invMass(eventParticles[i], eventParticles[j]));
						if (eventParticles[i].getCharge()*eventParticles[j].getCharge() > 0) {
							invMassEqSgnH->Fill(invMass(eventParticles[i], eventParticles[j]));
						}
						else invMassOppSgnH->Fill(invMass(eventParticles[i], eventParticles[j]));
						if (eventParticles[i].getIndex() == piPlusI) {				// pi+
							if (eventParticles[j].getIndex() == kPlusI)					// k+
								invMassPKEqH->Fill(invMass(eventParticles[i], eventParticles[j]));
							else if (eventParticles[j].getIndex() == kMinI) 		// k-
								invMassPKOppH->Fill(invMass(eventParticles[i], eventParticles[j]));
						} else if(eventParticles[i].getIndex() == piMinI) { 	// pi-
							if (eventParticles[j].getIndex() == kPlusI)					// k+
								invMassPKOppH->Fill(invMass(eventParticles[i], eventParticles[j]));
							else if (eventParticles[j].getIndex() == kMinI)			// k-
								invMassPKEqH->Fill(invMass(eventParticles[i], eventParticles[j]));
						}
					}
				}
			}
		}
	}
	std::cout << "Generation finished, saving histograms to file.\n";
	TFile* output = new TFile("output/genOutput.root", "RECREATE");
	pTypesH->Write();
	anglesH->Write();
	momentumH->Write();
	tMomentumH->Write();
	energyH->Write();
	invMassTotH->Write();
	invMassEqSgnH->Write();
	invMassOppSgnH->Write();
	invMassPKEqH->Write();
	invMassPKOppH->Write();
	invMassKStarH->Write();
	output->Close();
	delete pTypesH;
	delete anglesH;
	delete momentumH;
	delete tMomentumH;
	delete energyH;
	delete invMassTotH;
	delete invMassEqSgnH;
	delete invMassOppSgnH;
	delete invMassPKEqH;
	delete invMassPKOppH;
	delete invMassKStarH;
	delete output;
}
