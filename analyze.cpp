#include <cmath>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"

#include "particle/particle.hpp"

double pTypesProb(Double_t * x, Double_t * par) {
	if (*x < 0.5) return 0.;
	else if (*x < 2.5) return 0.4 * par[0];
	else if (*x < 4.5) return 0.05 * par[0];
	else if (*x < 6.5) return 0.045 * par[0];
	else if (*x < 7.5) return 0.01 * par[0];
	else return 0.;
}

void checkInput() {

}

void analyze() {
	checkInput();

	TFile* input = new TFile("output/genOutput.root", "READ");

	int nEvents {static_cast<int>(1E5)};
	int nParticles {100};

	TH1F* pTypesH = (TH1F*) input->Get("pTypesH");
	TH2F* anglesH = (TH2F*) input->Get("anglesH");
	TH1F* momentumH = (TH1F*) input->Get("momentumH");
	TH1F* tMomentumH = (TH1F*) input->Get("tMomentumH");
	TH1F* energyH = (TH1F*) input->Get("energyH");
	TH1F* invMassTotH = (TH1F*) input->Get("invMassTotH");
	TH1F* invMassOppSgnH = (TH1F*) input->Get("invMassOppSgnH");
	TH1F* invMassEqSgnH = (TH1F*) input->Get("invMassEqSgnH");
	TH1F* invMassPKOppH = (TH1F*) input->Get("invMassPKOppH");
	TH1F* invMassPKEqH = (TH1F*) input->Get("invMassPKEqH");
	TH1F* invMassKStarH = (TH1F*) input->Get("invMassKStarH");

	assert(pTypesH->GetEntries() == nEvents*nParticles);
	assert(anglesH->GetEntries() == nEvents*nParticles);
	assert(momentumH->GetEntries() == nEvents*nParticles);
	assert(tMomentumH->GetEntries() == nEvents*nParticles);
	assert(energyH->GetEntries() == nEvents*nParticles);
	std::cout << "Generated particle data histograms have the expected number of entries.\n";

	TF1* pTypesF = new TF1("expectedPTypes", pTypesProb, 0.5, 7.5, 1);

	pTypesH->Fit(pTypesF);

	TH1F* thetaH = (TH1F*) anglesH->ProjectionX("thetaH");
	TH1F* phiH = (TH1F*) anglesH->ProjectionY("phiH");

	TF1* thetaF = new TF1("thetaF", "[0]", 0., 2*M_PI);
	TF1* phiF = new TF1("phiF", "[0]", 0., 2*M_PI);

	thetaH->Fit(thetaF);
	phiH->Fit(phiF);

	TF1* momentumF = new TF1("momentumF", "[1]*exp([0]*x)", 0., 12.);
	momentumF->SetParameter(0, -1.);

	momentumH->Fit(momentumF);

	invMassOppSgnH->Add(invMassEqSgnH, -1.);
	invMassPKOppH->Add(invMassPKEqH, -1.);

	TFile* output = new TFile("output/processedData.root", "RECREATE");
	pTypesH->Write();
	thetaH->Write();
	phiH->Write();
	invMassOppSgnH->Write();
	invMassPKOppH->Write();
	output->Close();
}
