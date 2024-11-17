#include <cassert>

#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "particle/particle.hpp"

void simulate() {
	std::cout << "Simulation started succesfully.\n";

	double piMass = 0.13957;
  double kMass = 0.49367;
  double PMass = 0.93827;
  double KstarMass = 0.89166;
  double KstarLength = 0.05;
	
	Particle::addParticleType("pi+", piMass, 1);
	Particle::addParticleType("pi-", piMass, -1);
	Particle::addParticleType("k+", kMass, 1);
	Particle::addParticleType("k-", kMass, -1);
	Particle::addParticleType("P+", PMass, +1);
	Particle::addParticleType("P-", PMass, -1);
	Particle::addParticleType("k*", KstarMass, 0, KstarLength);
	std::cout << "Succesfully added particle types.\n";

	gRandom->SetSeed();

	TH1F* pTypesH = new TH1F("pTypesH", "Particle types", 7, 0.5, 7.5);
	TH2F* anglesH = new TH2F("anglesH", "Angles", 1000, 0., 2*M_PI, 500, 0., M_PI);
	TH1F* momentumH = new TH1F("momentumH", "Momentum norm", 1000, 0., 3.);
	TH1F* tMomentumH = new TH1F("tMomentumH", "Transverse momentum", 1000, 0., 3.);
	TH1F* energyH = new TH1F("energyH", "Energy", 1000, 0., 3.);
	TH1F* invMassTot = new TH1F("invMassTot", "Invariant mass", 1000, 0., 3.);
	TH1F* invMassOppSgn = new TH1F("invMassDisc", "Invariant mass, opp sign", 1000, 0., 3.);
	TH1F* invMassEqSgn = new TH1F("invMassConc", "Invariant mass, equal sign", 1000, 0., 3.);
	TH1F* invMassPKOpp = new TH1F("invMassPKOpp", "Invariant mass, opp. sign p and k", 1000, 0., 3.);
	TH1F* invMassPKEq = new TH1F("invMassPKEq", "Invariant mass, equal sign p and k", 1000, 0., 3.);
	TH1F* invMassKStar = new TH1F("invMassKStar", "Invariant mass, K* decay products", 1000, 0., 3.);

	int nEvents {static_cast<int>(1E5)};
	const int nGens {100}; // number of generated particles
	const int nParticles {nGens * 6 / 5}; // safe array size for 1% likelyhood of k* generation

	// declaring variables for the loops to use
	Particle eventParticles[nParticles];
	
	double phi {};
	double theta {};
	double pNorm {};
	double control {};

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
			if (control < 0.4) eventParticles[i].setIndex(0); 				// pi+, 40%
			else if (control < 0.8) eventParticles[i].setIndex(1); 		// pi-, 40%
			else if (control < 0.85) eventParticles[i].setIndex(2); 	// k+, 5%
			else if (control < 0.9) eventParticles[i].setIndex(3); 		// k-, 5%
			else if (control < 0.945) eventParticles[i].setIndex(4);	// P+, 4.5%
			else if (control < 0.99) eventParticles[i].setIndex(5); 	// P-, 4.5%
			else {
				eventParticles[i].setIndex(6);													// k*, 1%
				if (control < 0.995) {																	// pi+ k- decay, 50%
					dau1.setIndex(0); // pi+
					dau2.setIndex(3); // k-
				} else {																								//pi- k+ decay, 50%
					dau1.setIndex(1); // pi-
					dau2.setIndex(2); // k+
				}
				// std::cout << "K* generated, decay process result: " << eventParticles[i].decayToBody(dau1, dau2) << std::endl;
				eventParticles[tailIndex(eventParticles, nParticles, nGens)] = dau1;
				eventParticles[tailIndex(eventParticles, nParticles, nGens)] = dau2;
				invMassKStar->Fill(invMass(dau1, dau2));
			}
		}

		// fill histograms
		int max = tailIndex(eventParticles, nParticles, nGens);
		std::cout << "Found tailIndex for this iteration: " << max << "\n";
		for (int i {}; i <= max; ++i) {
			if (eventParticles[i].getIndex() != 6) {
				momentumH->Fill(eventParticles[i].getMomentum().norm());
				tMomentumH->Fill(eventParticles[i].getMomentum().norm2() - eventParticles[i].getPz()*eventParticles[i].getPz());
				for (int j{i+1}; j <= max; ++j) {
					if (eventParticles[i].getIndex() != 6) {
						invMassTot->Fill(invMass(eventParticles[i], eventParticles[j]));
						invMassEqSgn->Fill(invMass(eventParticles[i], eventParticles[j]));
					}
				}
			}
		}
	}
	invMassKStar->Draw();
	//pTypesH->Draw();
	//anglesH->Draw();
	//momentumH->Draw();
	//tMomentumH->Draw();
}
