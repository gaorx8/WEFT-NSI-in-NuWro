#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>   // for Form()

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

double readXsec(const char *filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        std::exit(1);
    }

    double total = 0.0;
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int dyn;
        double events, ratio, sigma;
        if (!(iss >> dyn >> events >> ratio >> sigma)) break;
        total += sigma;
    }
    return total;
}

static std::string root2txt(const char* rootfile){
    return std::string(rootfile) + ".txt";
}

static TH1D* make_costheta_hist(const char* rootfile,
                                const char* histname,
                                const Double_t* binEdges,
                                int nbins,
                                double ngen = 100000.0)
{
    TFile *file = TFile::Open(rootfile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Unable to open ROOT file: " << rootfile << std::endl;
        std::exit(1);
    }

    TTree *tree = (TTree*)file->Get("treeout");
    if (!tree) {
        std::cerr << "Unable to find TTree 'treeout' in: " << rootfile << std::endl;
        std::exit(1);
    }

    std::vector<double> dcostheta(nbins);
    for (int i = 0; i < nbins; i++) dcostheta[i] = binEdges[i+1] - binEdges[i];

    gROOT->cd();
    if (auto old = (TH1*)gROOT->FindObject(histname)) delete old;

    TH1D* h = new TH1D(histname,
        ";cos#theta_{#mu};d#sigma/dcos#theta_{#mu} [10^{-38} cm^{2}/nucleon]",
        nbins, binEdges);
    h->SetDirectory(gROOT);
    h->Reset();

    const Long64_t nfilled = tree->Project(histname, "costheta_muon()", "select() == 1");
    if (nfilled <= 0) {
        std::cerr << "[WARN] " << rootfile
                  << " Project filled 0 entries. Check select() or costheta_muon()."
                  << std::endl;
    }

    const std::string txt = root2txt(rootfile);
    const double xsec_total = readXsec(txt.c_str());

    for (int i = 1; i <= nbins; ++i) {
        const double raw = h->GetBinContent(i);
        const double dsigma = raw * 1e38 * xsec_total / ngen / dcostheta[i - 1];
        h->SetBinContent(i, dsigma);
    }

    file->Close();
    delete file;
    return h;
}

void plotcostheta() {
    const char *filename     = "SM.root";   // SM
    const char *nsi_filename = "SM_R.root"; // SM+NSI

    const Int_t nbins = 18;
    Double_t binEdges[nbins + 1] = {
        -1, -0.85, -0.7, -0.575, -0.45, -0.325, -0.2, -0.1, 0,
         0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.725, 0.85, 0.95, 1
    };
    const double ngen = 100000.0;

    TH1D* h_sm     = make_costheta_hist(filename,     "h_sm",     binEdges, nbins, ngen);
    TH1D* h_sm_nsi = make_costheta_hist(nsi_filename, "h_sm_nsi", binEdges, nbins, ngen);

    gROOT->cd();
    if (auto old = (TH1*)gROOT->FindObject("h_delta")) delete old;

    TH1D* h_delta = (TH1D*)h_sm_nsi->Clone("h_delta");
    h_delta->SetDirectory(gROOT);
    h_delta->Reset();
    h_delta->Add(h_sm_nsi, 1.0);  
    h_delta->Add(h_sm, -1.0);    
    h_delta->Divide(h_sm);     

    h_sm->SetLineColor(kBlue);
    h_sm_nsi->SetLineColor(kRed);
    h_sm->SetLineWidth(5);
    h_sm_nsi->SetLineWidth(5);
    h_sm->SetMarkerStyle(0);
    h_sm_nsi->SetMarkerStyle(0);

    const double ymax_abs = std::max(h_sm->GetMaximum(), h_sm_nsi->GetMaximum());
    h_sm->SetMinimum(0.0);
    h_sm->SetMaximum(1.15 * ymax_abs);

    h_delta->SetLineColor(kBlack);
    h_delta->SetLineWidth(5);
    h_delta->SetMarkerStyle(0);
    h_delta->SetTitle(";cos#theta_{#mu};#Delta(#cos#theta_{#mu})");

    double maxabs = 0.0;
    for (int i = 1; i <= h_delta->GetNbinsX(); ++i) {
        const double v = h_delta->GetBinContent(i);
        if (std::isfinite(v)) maxabs = std::max(maxabs, std::abs(v));
    }
    const double pad = (maxabs > 0.0) ? 0.15 * maxabs : 1e-4;
    h_delta->SetMinimum(-maxabs - pad);
    h_delta->SetMaximum(+maxabs + pad);

    TCanvas *c = new TCanvas("c", "costheta abs + delta", 2400, 2400);
    gStyle->SetOptStat(0);

    TPad *p1 = new TPad("p1","p1", 0.0, 0.34, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(0.14);
    p1->SetRightMargin(0.04);
    p1->Draw();
    p1->cd();

    h_sm->SetTitle(";cos#theta_{#mu};d#sigma/dcos#theta_{#mu} [10^{-38} cm^{2}/nucleon]");
    h_sm->GetXaxis()->SetLabelSize(0.0);   
    h_sm->GetXaxis()->SetTitleSize(0.0);   
    h_sm->GetYaxis()->SetTitleSize(0.045);
    h_sm->Draw("HIST");
    h_sm_nsi->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.18, 0.70, 0.55, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->AddEntry(h_sm,     "SM", "l");
    leg->AddEntry(h_sm_nsi, "SM+NSI", "l");
    leg->Draw();

    c->cd();
    TPad *p2 = new TPad("p2","p2", 0.0, 0.0, 1.0, 0.34);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.32);
    p2->SetLeftMargin(0.14);
    p2->SetRightMargin(0.04);
    p2->Draw();
    p2->cd();

    h_delta->GetXaxis()->SetTitleSize(0.085);
    h_delta->GetXaxis()->SetLabelSize(0.075);
    h_delta->GetYaxis()->SetTitleSize(0.085);
    h_delta->GetYaxis()->SetLabelSize(0.075);
    h_delta->GetYaxis()->SetTitleOffset(0.5);

    h_delta->GetYaxis()->SetTitle("#delta(cos#theta_{#mu})");
    h_delta->Draw("HIST");

    TLine *l0 = new TLine(-1.0, 0.0, 1.0, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineWidth(3);
    l0->Draw("SAME");

    c->Modified();
    c->Update();
    c->Print("costheta_abs_delta.jpg");
    c->Print("costheta_abs_delta.pdf");

    delete l0;
    delete leg;
    delete p2;
    delete p1;
    delete c;
    delete h_sm;
    delete h_sm_nsi;
    delete h_delta;
}
