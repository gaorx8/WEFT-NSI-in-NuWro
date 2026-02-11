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

// Safe divide: rel = (NSI - SM)/SM with protection for near-zero denominator
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

// top pad axis style
static void apply_axis_style_top(TH1* h){
    h->GetXaxis()->SetTitleSize(0.00);
    h->GetXaxis()->SetLabelSize(0.00);

    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.050);
    h->GetYaxis()->SetTitleOffset(1.10);
}

// bottom pad axis style
static void apply_axis_style_bottom(TH1* h){
    // pad height ratio: top(0.62)/bottom(0.38) ≈ 1.63
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

    // ---------- PRD-like global style ----------
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Make axis labels display enough digits (avoid "all zeros" when numbers are tiny)
    // This suppresses the "×10^n" compact notation by allowing more digits on the axis.
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

    // rel = (NSI - SM) / SM
    std::unique_ptr<TH1> hR(static_cast<TH1*>(h2->Clone("hRel")));
    hR->SetDirectory(nullptr);
    hR->Add(h1.get(), -1.0);
    make_relative_diff_inplace(hR.get(), h1.get(), /*eps=*/0.0);

    // percent
    hR->Scale(100.0);

    // --------- Axis labels ----------
    const char* xTitle = "log_{10}(E_{#nu}/GeV)";

    h1->GetXaxis()->SetTitle(xTitle);
    // h1->GetYaxis()->SetTitle("Oscillated atmospheric flux\n#times CCQE cross section (arb. units)");
    h1->GetYaxis()->SetTitle(" ");  // 占位，避免ROOT自己画复杂title

    hR->GetXaxis()->SetTitle(xTitle);
    hR->GetYaxis()->SetTitle("(N_{SM+NSI}-N_{SM})/N_{SM}  [%]");

    // IMPORTANT: do NOT disable exponent here; instead, show enough significant digits.
    // Otherwise tiny numbers (~1e-33) will format to 0.00.


    // --------- Crop to 0–10 GeV => log10(E/GeV) in [-1, 1] ----------
    // const double x_min = -1.0; // log10(0.1 GeV)
    const double x_min = std::log10(0.2);
    const double x_max = 1.0;   // log10(10 GeV)

    // 找到 bin index
    int bin_min = h1->GetXaxis()->FindBin(x_min);
    int bin_max = h1->GetXaxis()->FindBin(x_max - 1e-6); // 防止踩到下一 bin

    h1->GetXaxis()->SetRange(bin_min, bin_max);
    h2->GetXaxis()->SetRange(bin_min, bin_max);
    hR->GetXaxis()->SetRange(bin_min, bin_max);

    // Styles
    set_hist_style(h1.get(), 1, 3, 1);
    set_hist_style(h2.get(), 2, 3, 1);
    set_hist_style(hR.get(), 4, 3, 1);

    // Y-range for top pad
    double ymin = std::min(h1->GetMinimum(), h2->GetMinimum());
    double ymax = std::max(h1->GetMaximum(), h2->GetMaximum());
    if (ymin == ymax) { ymin *= 0.9; ymax *= 1.1; }
    h1->SetMinimum(std::max(0.0, ymin*0.9));
    h1->SetMaximum(ymax*1.25);

    // Relative diff range (auto, symmetric)
    double rmin = hR->GetMinimum();
    double rmax = hR->GetMaximum();
    double rabs = std::max(std::fabs(rmin), std::fabs(rmax));
    if (rabs <= 0) rabs = 1e-3;
    hR->SetMinimum(-1.25*rabs);
    hR->SetMaximum( 1.25*rabs);

    // Canvas/pads
    TCanvas c("c","",800,900);
    TPad p1("p1","",0,0.38,1,1.0);
    TPad p2("p2","",0,0.0, 1,0.38);

    p1.SetBottomMargin(0.06);
    p2.SetTopMargin(0.06);

    // give more left margin to avoid y-title clipping
    p1.SetLeftMargin(0.16);
    p1.SetRightMargin(0.05);
    p1.SetTopMargin(0.07);

    p2.SetLeftMargin(0.16);
    p2.SetRightMargin(0.05);
    p2.SetBottomMargin(0.30);

    p1.Draw(); p2.Draw();

    // remove histogram titles
    h1->SetTitle("");
    h2->SetTitle("");
    hR->SetTitle("");

    // ---------------- TOP ----------------
    p1.cd();
    apply_axis_style_top(h1.get());
    apply_axis_style_top(h2.get());

    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    // ---- two-line y-label using TLatex (robust, vertical) ----
    TLatex ylab;
    ylab.SetNDC(true);
    ylab.SetTextFont(42);
    ylab.SetTextSize(0.055);
    ylab.SetTextAngle(90);   // vertical text

    ylab.SetTextSize(0.055);   // 放大一档（PRD友好）

    // 竖排时：x 控制上下，y 控制左右
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

    // ---------------- BOTTOM ----------------
    p2.cd();
    apply_axis_style_bottom(hR.get());
    hR->Draw("HIST");

    TLine z;
    z.SetLineStyle(2);
    z.DrawLine(x_min, 0.0, x_max, 0.0);

    // Save
    c.SaveAs(out_pdf);
    std::cout << "Saved: " << out_pdf << "\n";

    TString out_png(out_pdf);
    if (out_png.EndsWith(".pdf")) out_png.ReplaceAll(".pdf", ".png");
    else out_png += ".png";
    c.SaveAs(out_png);
    std::cout << "Saved: " << out_png << "\n";
}
