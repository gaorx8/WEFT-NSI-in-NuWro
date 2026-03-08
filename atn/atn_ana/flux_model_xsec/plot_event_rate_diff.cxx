// plot_fig2_diff.cxx
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <iostream>
#include <memory>
#include <algorithm>
#include <cmath>

static void set_hist_style(TH1* h, int color, int width, int style=1){
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetLineStyle(style);
    h->SetMarkerStyle(0);
}

static void make_relative_diff_inplace(TH1* h_rel, const TH1* h_den, double eps = 0.0){
    for (int ib = 1; ib <= h_rel->GetNbinsX(); ++ib){
        const double den = h_den->GetBinContent(ib);
        const double num = h_rel->GetBinContent(ib);
        const double den_use = (std::fabs(den) > eps ? den : 0.0);

        if (den_use == 0.0){
            h_rel->SetBinContent(ib, 0.0);
            h_rel->SetBinError(ib,   0.0);
        } else {
            h_rel->SetBinContent(ib, num / den_use);
            h_rel->SetBinError(ib, 0.0);
        }
    }
}

static void apply_axis_style_top(TH1* h){
    h->GetXaxis()->SetTitleSize(0.00);
    h->GetXaxis()->SetLabelSize(0.00);

    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.050);
    h->GetYaxis()->SetTitleOffset(1.10);
}

static void apply_axis_style_bottom(TH1* h){
    const double scale = 1.63;

    h->GetXaxis()->SetTitleSize(0.055 * scale);
    h->GetXaxis()->SetLabelSize(0.050 * scale);
    h->GetXaxis()->SetTitleOffset(1.2);

    h->GetYaxis()->SetTitleSize(0.055 * scale);
    h->GetYaxis()->SetLabelSize(0.050 * scale);
    h->GetYaxis()->SetTitleOffset(0.65);
}

void plot_event_rate_diff(const char* f_sm   = "atnu_mu_model_xsec_L.root",
                          const char* f_nsi  = "atnu_mu_model_xsec_R.root",
                          int iTheta = 62,
                          int jDelta = 51,
                          const char* out_pdf = "fig3.pdf")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TGaxis::SetMaxDigits(6);

    TFile Fsm(f_sm, "READ");
    TFile Fnsi(f_nsi, "READ");
    if (Fsm.IsZombie() || Fnsi.IsZombie()){
        std::cerr << "Cannot open input root files.\n";
        return;
    }

    TString hname = TString::Format("atn_%d_%d", iTheta, jDelta);
    TH1* hSM  = dynamic_cast<TH1*>(Fsm.Get(hname));
    TH1* hNSI = dynamic_cast<TH1*>(Fnsi.Get(hname));
    if (!hSM || !hNSI){
        std::cerr << "Cannot find hist: " << hname << "\n";
        return;
    }

    std::unique_ptr<TH1> h1(static_cast<TH1*>(hSM->Clone("hSM_clone")));
    std::unique_ptr<TH1> h2(static_cast<TH1*>(hNSI->Clone("hNSI_clone")));
    h1->SetDirectory(nullptr);
    h2->SetDirectory(nullptr);

    std::unique_ptr<TH1> hR(static_cast<TH1*>(h2->Clone("hRel")));
    hR->SetDirectory(nullptr);
    hR->Add(h1.get(), -1.0);
    make_relative_diff_inplace(hR.get(), h1.get(), /*eps=*/0.0);


    hR->Scale(100.0);

    const char* xTitle = "log_{10}(E_{#nu}/GeV)";

    h1->GetXaxis()->SetTitle(xTitle);
    h1->GetYaxis()->SetTitle(" "); 
    hR->GetXaxis()->SetTitle(xTitle);
    hR->GetYaxis()->SetTitle("(N_{SM+NSI}-N_{SM})/N_{SM}  [%]");

    const double x_min = std::log10(0.2);
    const double x_max = 1.0;

    int bin_min = h1->GetXaxis()->FindBin(x_min);
    int bin_max = h1->GetXaxis()->FindBin(x_max - 1e-6);

    h1->GetXaxis()->SetRange(bin_min, bin_max);
    h2->GetXaxis()->SetRange(bin_min, bin_max);
    hR->GetXaxis()->SetRange(bin_min, bin_max);

    set_hist_style(h1.get(), 1, 3, 1);
    set_hist_style(h2.get(), 2, 3, 1);
    set_hist_style(hR.get(), 4, 3, 1);

    double ymin = std::min(h1->GetMinimum(), h2->GetMinimum());
    double ymax = std::max(h1->GetMaximum(), h2->GetMaximum());
    if (ymin == ymax) { ymin *= 0.9; ymax *= 1.1; }
    h1->SetMinimum(std::max(0.0, ymin*0.9));
    h1->SetMaximum(ymax*1.25);

    double rmin = hR->GetMinimum();
    double rmax = hR->GetMaximum();
    double rabs = std::max(std::fabs(rmin), std::fabs(rmax));
    if (rabs <= 0) rabs = 1e-3;
    hR->SetMinimum(-1.25*rabs);
    hR->SetMaximum( 1.25*rabs);

    TCanvas c("c","",800,900);
    TPad p1("p1","",0,0.38,1,1.0);
    TPad p2("p2","",0,0.0, 1,0.38);

    p1.SetBottomMargin(0.06);
    p2.SetTopMargin(0.06);

    p1.SetLeftMargin(0.16);
    p1.SetRightMargin(0.05);
    p1.SetTopMargin(0.07);

    p2.SetLeftMargin(0.16);
    p2.SetRightMargin(0.05);
    p2.SetBottomMargin(0.30);

    p1.Draw(); p2.Draw();

    h1->SetTitle("");
    h2->SetTitle("");
    hR->SetTitle("");

    p1.cd();
    apply_axis_style_top(h1.get());
    apply_axis_style_top(h2.get());

    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    TLatex ylab;
    ylab.SetNDC(true);
    ylab.SetTextFont(42);
    ylab.SetTextSize(0.055);
    ylab.SetTextAngle(90); 

    ylab.SetTextSize(0.055); 

    ylab.DrawLatex(0.040, 0.10, "Oscillated atmospheric flux #times #sigma_{CCQE}");
    ylab.DrawLatex(0.075, 0.10, "per nucleon  [s^{-1} GeV^{-1}]");

    TLegend leg(0.62,0.72,0.86,0.84);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.045);
    leg.AddEntry(h1.get(), "SM", "l");
    leg.AddEntry(h2.get(), "SM+NSI", "l");
    leg.Draw();

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.045);
    lat.DrawLatex(0.62, 0.86, "NO, assuming no CPV");

    p2.cd();
    apply_axis_style_bottom(hR.get());
    hR->Draw("HIST");

    TLine z;
    z.SetLineStyle(2);
    z.DrawLine(x_min, 0.0, x_max, 0.0);

    c.SaveAs(out_pdf);
    std::cout << "Saved: " << out_pdf << "\n";

    TString out_png(out_pdf);
    if (out_png.EndsWith(".pdf")) out_png.ReplaceAll(".pdf", ".png");
    else out_png += ".png";
    c.SaveAs(out_png);
    std::cout << "Saved: " << out_png << "\n";
}
