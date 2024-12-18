#include <cmath>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TFormula.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THistPainter.h"

#include "particle/particle.hpp"

void analyze() {

	gStyle->SetFillColor(kViolet-7);
	gStyle->SetLineColor(kGreen+3);

	TFile* input = new TFile("output/genOutput.root", "READ");

	gROOT->ForceStyle();

	int nEvents {100000};
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
	std::cout << "Generated particle data histograms have the expected number of entries.\n\n";

	const int expPTypesN = 7; // expected particle types number
	assert(expPTypesN == pTypesH->GetEntries());

	double pTypesPercentages[expPTypesN];
	double pTypesPercErrs[expPTypesN];

	for (int i{}; i < expPTypesN; ++i) {
		pTypesPercentages[i] = (pTypesH->GetBinContent(i + 1) * 100.) / pTypesH->GetEntries();
		pTypesPercErrs[i] = (pTypesH->GetBinError(i + 1) * 100.) / pTypesH->GetEntries();
	}

	const char* particleNames[expPTypesN] {
		"pi+", "pi-", "k+", "k-", "P+", "P-", "k*"
	};

	std::cout << "Generated particles %:\n";
	for(int i{0}; i < expPTypesN; ++i){
		std::cout << "Percentage for " << particleNames[i] << " particle: " << pTypesPercentages[i] << " +/- "<< pTypesPercErrs[i] << "\n";
	}

	TH1F* thetaH = (TH1F*) anglesH->ProjectionX("thetaH");
	TF1* thetaF = new TF1("thetaF", "[0]", 0., 2*M_PI);

	thetaH->SetMinimum(0.);
	thetaH->SetTitle("Occorrenze - Angolo polare");
	
	thetaH->Fit(thetaF);

	std::cout << "Theta probability density 	=	 " << thetaF->GetParameter(0) / thetaH->GetEntries()
	<< " +/- " << thetaF->GetParError(0) / thetaH->GetEntries() << std::endl;
	std::cout << "Chisquare/NDF 							=	 " << thetaF->GetChisquare() / thetaF->GetNDF() << std::endl;
	std::cout << "Fit likelyhood 							=  " << thetaF->GetProb() << std::endl;

	TH1F* phiH = (TH1F*) anglesH->ProjectionY("phiH");
	TF1* phiF = new TF1("phiF", "[0]", 0., 2*M_PI);
	
	phiH->SetMinimum(0.);
	phiH->SetTitle("Occorrenze - Angolo azimutale");

	phiH->Fit(phiF);
	std::cout << "Chisquare/NDF   =  " << phiF->GetChisquare() / thetaF->GetNDF() << std::endl;
	std::cout << "Fit likelyhood	=  " << phiF->GetProb() << "\n";
	
	TF1* momentumF = new TF1("momentumF", "[1]*exp([0]*x)", 0., 12.);
	momentumF->SetParameter(0, -1.);
	momentumF->SetNpx(1000);
	momentumH->SetTitle("Occorrenze per quantita di moto");
	momentumH->GetXaxis()->SetRangeUser(0., 9.);

	momentumH->Fit(momentumF);
	std::cout << "Chisquare/NDF     =  " << momentumF->GetChisquare() / thetaF->GetNDF() << std::endl;
	std::cout << "Fit likelyhood		=	 " << momentumF->GetProb() << std::endl;
	std::cout << "Expected mean for momentum distribution = 1.\n";

	invMassOppSgnH->Add(invMassEqSgnH, -1.);
	invMassPKOppH->Add(invMassPKEqH, -1.);

	TF1* invMassF = new TF1("invMassF", "gaus", 0.3, 1.5);
	invMassOppSgnH->SetTitle("Massa invariante - differenza totale");
	invMassOppSgnH->GetXaxis()->SetRangeUser(0.2, 3.4);
	invMassF->SetNpx(10000);
	invMassF->SetParameter(0, 5000.);
	invMassF->SetParameter(1, .89);
	invMassF->SetParameter(2, .05);
	invMassOppSgnH->Fit(invMassF, "R");

	TF1* invMassPKF = new TF1("invMassPKF", "gaus", 0.6, 1.5);
	invMassPKOppH->SetTitle("Massa invariante - differenza per K e Pi");
	invMassPKOppH->GetXaxis()->SetRangeUser(0.6, 3.4);
	invMassPKF->SetNpx(10000);
	invMassPKF->SetParameter(0, 5000.);
	invMassPKF->SetParameter(1, .89);
	invMassPKF->SetParameter(2, .05);
	invMassPKOppH->Fit(invMassPKF, "R", "", 0.75, 1.1);

	TF1* invMassKStarF = new TF1("invMassKStarF", "gaus", 0.6, 1.5);
	invMassKStarH->SetTitle("Massa invariante - decadimenti della k*");
	invMassKStarF->SetNpx(10000);
	invMassKStarH->Fit(invMassKStarF, "R");

	TFile* output = new TFile("output/processedData.root", "RECREATE");
	pTypesH->Write();
	thetaH->Write();
	phiH->Write();
	momentumH->Write();
	invMassOppSgnH->Write();
	invMassPKOppH->Write();
	invMassKStarH->Write();
	output->Close();
}
