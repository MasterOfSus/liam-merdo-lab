#include <cmath>
#include <ostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFormula.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THistPainter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "particle/particle.hpp"

void analyze() {
  TFile* input = new TFile("output/genOutput.root", "READ");

  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1111);

  int nEvents{100000};
  int nParticles{100};

  TH1F* pTypesH = (TH1F*)input->Get("pTypesH");
  TH2F* anglesH = (TH2F*)input->Get("anglesH");
  TH1F* momentumH = (TH1F*)input->Get("momentumH");
  TH1F* tMomentumH = (TH1F*)input->Get("tMomentumH");
  TH1F* energyH = (TH1F*)input->Get("energyH");
  TH1F* invMassTotH = (TH1F*)input->Get("invMassTotH");
  TH1F* invMassOppSgnH = (TH1F*)input->Get("invMassOppSgnH");
  TH1F* invMassEqSgnH = (TH1F*)input->Get("invMassEqSgnH");
  TH1F* invMassPKOppH = (TH1F*)input->Get("invMassPKOppH");
  TH1F* invMassPKEqH = (TH1F*)input->Get("invMassPKEqH");
  TH1F* invMassKStarH = (TH1F*)input->Get("invMassKStarH");

  if (pTypesH->GetEntries() == nEvents * nParticles &&
      anglesH->GetEntries() == nEvents * nParticles &&
      momentumH->GetEntries() == nEvents * nParticles &&
      tMomentumH->GetEntries() == nEvents * nParticles &&
      energyH->GetEntries() == nEvents * nParticles)
    std::cout << "Generated particle data histograms have the expected number "
                 "of entries."
              << std::endl;
  else
    std::cout << "WARNING! Generated particle data histograms do not have the "
                 "expected number of entries!";
  std::cout << std::endl;

  const int expPTypesN = 7;  // expected particle types number
  if (expPTypesN != pTypesH->GetNbinsX())
    std::cout << "WARNING! The expected number of particles is different from "
                 "the pTypesH bin number!";

  double pTypesPercentages[expPTypesN];
  double pTypesPercErrs[expPTypesN];

  for (int i{}; i < expPTypesN; ++i) {
    pTypesPercentages[i] =
        (pTypesH->GetBinContent(i + 1) * 100.) / pTypesH->GetEntries();
    pTypesPercErrs[i] =
        (pTypesH->GetBinError(i + 1) * 100.) / pTypesH->GetEntries();
  }

  const char* particleNames[expPTypesN]{"pi+", "pi-", "k+", "k-",
                                        "P+",  "P-",  "k*"};

  std::cout << "Generated particles %:" << std::endl;
  for (int i{0}; i < expPTypesN; ++i) {
    std::cout << "Percentage for " << particleNames[i]
              << " particle: " << pTypesPercentages[i] << " +/- "
              << pTypesPercErrs[i] << std::endl;
    std::cout << "Entries for " << particleNames[i]
              << " particle: " << pTypesH->GetBinContent(i + 1) << "+/-"
              << pTypesH->GetBinError(i + 1) << std::endl;
  }
  std::cout << "Total entries: " << pTypesH->GetEntries() << std::endl;

  TH1F* thetaH = (TH1F*)anglesH->ProjectionX("thetaH");
  TF1* thetaF = new TF1("thetaF", "[0]", 0., 2 * M_PI);
  thetaF->SetParName(0, "Costante");
  thetaF->SetNpx(1000);
  thetaH->SetMinimum(0.);
  thetaH->SetTitle("Occorrenze - Angolo polare");

  thetaH->Fit(thetaF, "QOR");

  std::cout << std::endl << "Theta fit results:" << std::endl;
  double probDensityTheta{thetaF->GetParameter(0) / thetaH->GetEntries() *
                          thetaH->GetNbinsX() / M_PI / 2.};
  std::cout << "Theta probability density      =  " << probDensityTheta << " +/- "
            << thetaF->GetParError(0) / thetaF->GetParameter(0) *
                   probDensityTheta
            << std::endl;
  std::cout << "Theta average bin height       =  " << thetaF->GetParameter(0)
            << " +/- " << thetaF->GetParError(0) << std::endl;
  std::cout << "Chisquare/NDF                  =  " << thetaF->GetChisquare()
            << " / " << thetaF->GetNDF() << std::endl;
  std::cout << "Fit likelyhood                 =  " << thetaF->GetProb()
            << std::endl;

  TH1F* phiH = (TH1F*)anglesH->ProjectionY("phiH");
  TF1* phiF = new TF1("phiF", "[0]", 0., 2 * M_PI);
  phiF->SetParName(0, "Costante");
  phiF->SetNpx(1000);
  phiH->SetMinimum(0.);
  phiH->SetTitle("Occorrenze - Angolo azimutale");

  phiH->Fit(phiF, "QOR");

  std::cout << std::endl << "Phi fit results:" << std::endl;
  double probDensityPhi{phiF->GetParameter(0) / phiH->GetEntries() *
                        phiH->GetNbinsX() / M_PI};
  std::cout << "Phi probability density        =  " << probDensityPhi << " +/- "
            << phiF->GetParError(0) / phiF->GetParameter(0) * probDensityPhi
            << std::endl;
  std::cout << "Phi average bin height         =  " << phiF->GetParameter(0)
            << " +/- " << phiF->GetParError(0) << std::endl;
  std::cout << "Chisquare/NDF                  =  " << phiF->GetChisquare()
            << " / " << phiF->GetNDF() << std::endl;
  std::cout << "Fit likelyhood                 =  " << phiF->GetProb()
            << std::endl;

  TF1* momentumF = new TF1("momentumF", "[1]*exp(-x/[0])", 0., 20.);
  momentumF->SetParName(0, "#tau");
  momentumF->SetParName(1, "Ampiezza");
  momentumF->SetParameter(0, 1.);
  momentumF->SetNpx(1000);
  momentumH->SetTitle("Occorrenze - Quantita di moto");

  momentumH->Fit(momentumF, "QO");

  std::cout << std::endl << "Momentum fit results:" << std::endl;
  std::cout << "Tau value                      =  "
            << momentumF->GetParameter(0) << " +/- "
            << momentumF->GetParError(0) << std::endl;
  std::cout << "Amplitude value                =  "
            << momentumF->GetParameter(1) << " +/- "
            << momentumF->GetParError(1) << std::endl;
  std::cout << "Chisquare/NDF                  =  " << momentumF->GetChisquare()
            << " / " << momentumF->GetNDF() << std::endl;
  std::cout << "Fit likelyhood                 =  " << momentumF->GetProb()
            << std::endl;
  std::cout << "Expected tau                   =  1." << std::endl;

  {
    double integral{invMassOppSgnH->GetEntries() - invMassEqSgnH->GetEntries()};
    invMassOppSgnH->Add(invMassEqSgnH, -1.);
    invMassOppSgnH->SetEntries(integral);
    integral = invMassPKOppH->GetEntries() - invMassPKEqH->GetEntries();
    invMassPKOppH->Add(invMassPKEqH, -1.);
    invMassPKOppH->SetEntries(integral);
  }

  TF1* invMassF = new TF1("invMassF", "gaus", 0., 4.5);
  invMassF->SetParName(0, "Ampiezza");
  invMassF->SetParName(1, "Valor medio");
  invMassF->SetParName(2, "Dev. standard");
  invMassOppSgnH->SetTitle("Massa invariante - Differenza totale");
  invMassF->SetNpx(10000);
  invMassF->SetParameter(1, .89);
  invMassF->SetParameter(2, .05);
  invMassF->SetParLimits(1, 0., 20.);
  invMassF->SetParLimits(0, 0., 1E5);

  invMassOppSgnH->Fit(invMassF, "QOR");

  std::cout << std::endl
            << "All particles invariant mass subtraction fit results:"
            << std::endl;
  std::cout << "Mean value (k* invariant mass) =  " << invMassF->GetParameter(1)
            << " +/- " << invMassF->GetParError(1) << std::endl;
  std::cout << "Standard deviation (k* width)  =  " << invMassF->GetParameter(2)
            << " +/- " << invMassF->GetParError(2) << std::endl;
  std::cout << "Amplitude                      =  " << invMassF->GetParameter(0)
            << " +/- " << invMassF->GetParError(0) << std::endl;
  std::cout << "Chisquare/NDF                  =  " << invMassF->GetChisquare()
            << " / " << invMassF->GetNDF() << std::endl;
  std::cout << "Fit likelyhood                 =  " << invMassF->GetProb()
            << std::endl;

  TF1* invMassPKF = new TF1("invMassPKF", "gaus", 0.3, 1.5);
  invMassPKF->SetParName(0, "Ampiezza");
  invMassPKF->SetParName(1, "Valor medio");
  invMassPKF->SetParName(2, "Dev. standard");
  invMassPKOppH->SetTitle("Massa invariante - Differenza K e Pi discordi");
  invMassPKF->SetNpx(10000);
  invMassPKF->SetParameter(0, 5000.);
  invMassPKF->SetParameter(1, .89);
  invMassPKF->SetParameter(2, .05);

  invMassPKOppH->Fit(invMassPKF, "QO");

  std::cout << std::endl
            << "Pi and K particles invariant mass subtraction fit results:"
            << std::endl;
  std::cout << "Mean value (k* invariant mass) =  "
            << invMassPKF->GetParameter(1) << " +/- "
            << invMassPKF->GetParError(1) << std::endl;
  std::cout << "Standard deviation (k* width)  =  "
            << invMassPKF->GetParameter(2) << " +/- "
            << invMassPKF->GetParError(2) << std::endl;
  std::cout << "Amplitude                      =  "
            << invMassPKF->GetParameter(0) << " +/- "
            << invMassPKF->GetParError(0) << std::endl;
  std::cout << "Chisquare/NDF                  =  "
            << invMassPKF->GetChisquare() << " / " << invMassPKF->GetNDF()
            << std::endl;
  std::cout << "Fit likelyhood                 =  " << invMassPKF->GetProb()
            << std::endl;

  TF1* invMassKStarF = new TF1("invMassKStarF", "gaus", 0., 20.);
  invMassKStarF->SetParName(0, "Ampiezza");
  invMassKStarF->SetParName(1, "Valor medio");
  invMassKStarF->SetParName(2, "Dev. standard");
  invMassKStarH->SetTitle("Massa invariante - Decadimenti della k*");
  invMassKStarF->SetNpx(10000);
  invMassKStarH->Fit(invMassKStarF, "QO");

  std::cout << std::endl << "K* invariant mass fit results:" << std::endl;
  std::cout << "Mean value (k* invariant mass) =  "
            << invMassKStarF->GetParameter(1) << " +/- "
            << invMassKStarF->GetParError(1) << std::endl;
  std::cout << "Standard deviation (k* width)  =  "
            << invMassKStarF->GetParameter(2) << " +/- "
            << invMassKStarF->GetParError(2) << std::endl;
  std::cout << "Amplitude                      =  "
            << invMassKStarF->GetParameter(0) << " +/- "
            << invMassKStarF->GetParError(0) << std::endl;
  std::cout << "Chisquare/NDF                  =  "
            << invMassKStarF->GetChisquare() << " / " << invMassKStarF->GetNDF()
            << std::endl;
  std::cout << "Fit likelyhood                 =  " << invMassKStarF->GetProb()
            << std::endl
            << std::endl;

  TFile* output = new TFile("output/processedData.root", "RECREATE");

  const char* pTypes[expPTypesN]{"#pi+", "#pi-", "K+", "K-", "P+", "P-", "K*"};
  gStyle->SetTextFont(42);
  for (int i{1}; i <= expPTypesN; ++i)
    pTypesH->GetXaxis()->SetBinLabel(i, pTypes[i - 1]);

  TCanvas* particlesCvs = new TCanvas(
      "particlesCanvas",
      "Dati su tipi e informazioni dinamiche delle particelle", 1600, 1200);

  TH1F* histograms[7]{pTypesH,        thetaH,        phiH,         momentumH,
                      invMassOppSgnH, invMassPKOppH, invMassKStarH};
  TF1* functions[6]{thetaF,   phiF,       momentumF,
                    invMassF, invMassPKF, invMassKStarF};

  {
    int i{0};
    for (TH1F* histo : histograms) {
      histo->GetXaxis()->CenterTitle(true);
      histo->GetYaxis()->CenterTitle(true);
      histo->SetFillColor(kViolet + 6);
      histo->SetLineColor(kViolet - 7);
      histo->SetName("Risultati");  // terrible but only quasi-decent way to set
                                    // the stats box title
      if (i > 0) {
        histo->GetFunction(functions[i - 1]->GetName())
            ->SetLineColor(kGreen + 2);
        histo->GetFunction(functions[i - 1]->GetName())->SetLineWidth(4);
      }
      ++i;
    }
  }

  double leftMargin{0.1375};

  particlesCvs->Divide(2, 2);
  particlesCvs->cd(1);
  gPad->SetLeftMargin(leftMargin);
  pTypesH->GetXaxis()->SetTitle("Tipo");
  pTypesH->GetXaxis()->CenterLabels(true);
  pTypesH->GetXaxis()->SetTickLength(0.01);
  pTypesH->GetYaxis()->SetTitle("Occorrenze");
  pTypesH->Draw();
  particlesCvs->cd(3);
  gPad->SetLeftMargin(leftMargin);
  thetaH->GetXaxis()->SetTitle("#theta");
  thetaH->GetXaxis()->CenterTitle(false);
  thetaH->GetYaxis()->SetTitle("Occorrenze");
  thetaH->GetYaxis()->CenterTitle(true);
  thetaH->SetMaximum(1.5 * thetaF->GetParameter(0));
  thetaH->Draw();
  particlesCvs->cd(4);
  gPad->SetLeftMargin(leftMargin);
  phiH->GetXaxis()->SetTitle("#phi");
  phiH->GetXaxis()->CenterTitle(false);
  phiH->GetYaxis()->SetTitle("Occorrenze");
  phiH->GetYaxis()->CenterTitle(true);
  phiH->SetMaximum(1.5 * phiF->GetParameter(0));
  phiH->Draw();
  particlesCvs->cd(2);
  gPad->SetLeftMargin(leftMargin);
  momentumH->GetXaxis()->SetTitle("Momento");
  momentumH->GetXaxis()->SetRangeUser(0., 9.);
  momentumH->GetYaxis()->SetTitle("Occorrenze");
  momentumH->Draw();

  particlesCvs->Write();
  particlesCvs->Print("output/particlesCvs.gif");

  TCanvas* invMassCvs = new TCanvas(
      "invMassKStarCvs",
      "Distribuzione della massa invariante delle particelle", 1800, 600);

  invMassCvs->Divide(3, 0);
  invMassCvs->cd(1);
  gPad->SetLeftMargin(leftMargin);
  invMassOppSgnH->GetXaxis()->SetTitle("Massa invariante");
  invMassOppSgnH->GetXaxis()->SetRangeUser(0.6, 1.3);
  invMassOppSgnH->GetYaxis()->SetTitle("Occorrenze");
  invMassOppSgnH->Draw();

  invMassCvs->cd(2);
  gPad->SetLeftMargin(leftMargin);
  invMassPKOppH->GetXaxis()->SetTitle("Massa invariante");
  invMassPKOppH->GetXaxis()->SetRangeUser(0.6, 1.3);
  invMassPKOppH->GetYaxis()->SetTitle("Occorrenze");
  invMassPKOppH->Draw();

  invMassCvs->cd(3);
  gPad->SetLeftMargin(leftMargin);
  invMassKStarH->GetXaxis()->SetTitle("Massa invariante");
  invMassKStarH->GetYaxis()->SetTitle("Occorrenze");
  invMassKStarH->Draw();

  invMassCvs->Write();
  invMassCvs->Print("output/invMassCvs.gif");

  pTypesH->Write();
  thetaH->Write();
  phiH->Write();
  momentumH->Write();
  invMassOppSgnH->Write();
  invMassPKOppH->Write();
  invMassKStarH->Write();

  output->Close();

  delete pTypesH;
  delete thetaH;
  delete thetaF;
  delete phiH;
  delete phiF;
  delete momentumH;
  delete momentumF;
  delete tMomentumH;
  delete invMassOppSgnH;
  delete invMassF;
  delete invMassPKOppH;
  delete invMassPKF;
  delete invMassKStarH;
  delete invMassKStarF;
  delete particlesCvs;
  delete invMassCvs;
}
