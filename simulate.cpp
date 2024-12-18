#include <cassert>

#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TBenchmark.h"

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

	// std::cout << "Successfully added particle types.\n";

	gRandom->SetSeed();

	double histoExt {19.};
	TH1F* pTypesH = new TH1F("pTypesH", "Particle types", 7, 0.5, 7.5);
	TH2F* anglesH = new TH2F("anglesH", "Angles", 1000, 0., 2*M_PI, 500, 0., M_PI);
	TH1F* momentumH = new TH1F("momentumH", "Momentum norm", 1000, 0., histoExt);
	TH1F* tMomentumH = new TH1F("tMomentumH", "Transverse momentum", 1000, 0., histoExt);
	TH1F* energyH = new TH1F("energyH", "Energy", 1000, 0., histoExt);
	TH1F* invMassTotH = new TH1F("invMassTotH", "Invariant mass", 1000, 0., histoExt);
	TH1F* invMassOppSgnH = new TH1F("invMassOppSgnH", "Invariant mass, opp sign", 1000, 0., histoExt);
	invMassOppSgnH->Sumw2();
	TH1F* invMassEqSgnH = new TH1F("invMassEqSgnH", "Invariant mass, equal sign", 1000, 0., histoExt);
	invMassEqSgnH->Sumw2();
	TH1F* invMassPKOppH = new TH1F("invMassPKOppH", "Invariant mass, opp. sign p and k", 1000, 0., histoExt);
	invMassPKOppH->Sumw2();
	TH1F* invMassPKEqH = new TH1F("invMassPKEqH", "Invariant mass, equal sign p and k", 1000, 0., histoExt);
	invMassPKEqH->Sumw2();
	TH1F* invMassKStarH = new TH1F("invMassKStarH", "Invariant mass, K* decay products", 1000, 0.5, 1.2);

	int nEvents {static_cast<int>(1E5)};
	const int nGens {100}; // number of generated particles
	const int nParticles {nGens * 6 / 5}; // safe array size for 1% likelyhood of k* generation (120%)

	// declaring variables for the loops to use
	Particle eventParticles[nParticles];
	
	double phi {};
	double theta {};
	double pNorm {};
	double control {};

	int decayResult = {0};
	Particle dau1;
	Particle dau2;

	int jCurrentIndex {0};
	int kCurrentIndex {0};

	int tailIndex {100}; // index for first empty space in the eventParticles array after the 100th spot

	std::cout << "Correctly finished pre-loops setup. Starting generation.\n";
	
	gBenchmark->Start("Loop");

	for (int i {}; i < nEvents; ++i) {

		for (Particle& particle: eventParticles) {
			particle = {};
			assert(particle.getIndex() == -1);
		}

		tailIndex = 100;

		for (int j {}; j < nGens; ++j) {
			theta = gRandom->Uniform(0., 2*M_PI);
			phi = gRandom->Uniform(0., M_PI);
			pNorm = gRandom->Exp(1.);
			eventParticles[j].setMomentum( { sin(phi)*cos(theta)*pNorm, sin(phi)*sin(theta)*pNorm, cos(phi)*pNorm } );
			
			control = gRandom->Uniform(0., 1.);
			if (control < 0.4) eventParticles[j].setIndex(piPlusI); 				// pi+, 40%
			else if (control < 0.8) eventParticles[j].setIndex(piMinI); 		// pi-, 40%
			else if (control < 0.85) eventParticles[j].setIndex(kPlusI); 		// k+, 5%
			else if (control < 0.9) eventParticles[j].setIndex(kMinI); 			// k-, 5%
			else if (control < 0.945) eventParticles[j].setIndex(pPlusI);		// P+, 4.5%
			else if (control < 0.99) eventParticles[j].setIndex(pMinI); 		// P-, 4.5%
			else {
				eventParticles[j].setIndex(kStarI);														// k*, 1%
				if (control < 0.995) {		// pi+ k- decay, 50%
					dau1.setIndex(piPlusI); // pi+
					dau2.setIndex(kMinI); 	// k-
				} else {									//pi- k+ decay, 50%
					dau1.setIndex(piMinI); 	// pi-
					dau2.setIndex(kPlusI); 	// k+
				}
				decayResult = eventParticles[j].decayToBody(dau1, dau2);
				assert((!decayResult));
				// std::cout << "K* generated, decay process result: " << decayResult << std::endl;
				
				eventParticles[tailIndex] = dau1;
				++tailIndex; // increasing tailIndex
				eventParticles[tailIndex] = dau2;
				++tailIndex; // same as before
				invMassKStarH->Fill(invMass(dau1, dau2));
			}
			pTypesH->Fill(eventParticles[j].getIndex() + 1.);
			momentumH->Fill(pNorm);
			tMomentumH->Fill(sqrt(pNorm*pNorm*sin(theta)*sin(theta)));
			anglesH->Fill(theta, phi);
			energyH->Fill(eventParticles[j].getEnergy());
		}

		// fill histograms
		for (int j {}; j < tailIndex; ++j) {
			jCurrentIndex = eventParticles[j].getIndex();
			if (jCurrentIndex != kStarI) {
				//momentumH->Fill(eventParticles[i].getMomentum().norm());
				//tMomentumH->Fill(eventParticles[i].getMomentum().norm2() - eventParticles[i].getPz()*eventParticles[i].getPz());
				for (int k{j+1}; k < tailIndex; ++k) {
					kCurrentIndex = eventParticles[k].getIndex();
					if (kCurrentIndex != kStarI) {
						double invarMass = invMass(eventParticles[j], eventParticles[k]);
						invMassTotH->Fill(invarMass);
						
						if (eventParticles[j].getCharge()*eventParticles[k].getCharge() > 0) {
							invMassEqSgnH->Fill(invarMass);
						}
						else invMassOppSgnH->Fill(invarMass);
						
						if (jCurrentIndex == piPlusI) {				// pi+
							if (kCurrentIndex == kPlusI)					// k+
								invMassPKEqH->Fill(invarMass);
							else if (kCurrentIndex == kMinI) 		// k-
								invMassPKOppH->Fill(invarMass);
						} else if(jCurrentIndex == piMinI) { 	// pi-
							if (kCurrentIndex == kPlusI)					// k+
								invMassPKOppH->Fill(invarMass);
							else if (kCurrentIndex == kMinI)			// k-
								invMassPKEqH->Fill(invarMass);
						}
					}
				}
			}
		}
	}

	gBenchmark->Show("Loop");

	gBenchmark->Reset();

	std::cout << "Generation finished, saving data to file.\n";
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
	std::cout << "Histograms saved, deleting histograms from heap.\n";
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
