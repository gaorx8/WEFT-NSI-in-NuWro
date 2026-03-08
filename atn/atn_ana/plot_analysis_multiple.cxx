#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH2.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TList.h"
#include "TPad.h"
#include "TMarker.h"

#include <vector>
#include <iostream>
#include <memory>
#include <cmath>

static std::vector<TGraph*> ExtractContourGraphs(TH2* h_dchi2, double level) {
    std::vector<TGraph*> out;
    if (!h_dchi2) return out;

    TPad* pad_save = (TPad*)gPad;

    if (auto specials = gROOT->GetListOfSpecials()) {
        if (auto old = specials->FindObject("contours")) specials->Remove(old);
    }

    auto hc = (TH2*)h_dchi2->Clone(Form("tmp_for_contours_%p_%.6f", (void*)h_dchi2, level));
    hc->SetDirectory(nullptr);
    hc->SetContour(1);
    hc->SetContourLevel(0, level);

    TCanvas ctmp("ctmp_contour","ctmp_contour",600,600);
    ctmp.cd();
    hc->Draw("CONT Z LIST");
    ctmp.Update();

    auto conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    if (conts && conts->GetSize() > 0) {
        auto llist = (TList*)conts->At(0);
        if (llist && llist->GetSize() > 0) {
            for (int i = 0; i < llist->GetSize(); ++i) {
                auto gr = (TGraph*)llist->At(i);
                if (gr && gr->GetN() > 2) {
                    out.push_back((TGraph*)gr->Clone(Form("contour_level%.6f_seg%d", level, i)));
                }
            }
        }
    }

    delete hc;
    if (pad_save) pad_save->cd();
    return out;
}

static TH2* MakeDeltaChi2Hist(const TH2* raw) {
    if (!raw) return nullptr;

    int ixmin, iymin, izmin;
    ((TH2*)raw)->GetBinXYZ(((TH2*)raw)->GetMinimumBin(), ixmin, iymin, izmin);
    const double chi2min = raw->GetBinContent(ixmin, iymin, izmin);

    auto* hd = (TH2*)raw->Clone(Form("%s_dchi2", raw->GetName()));
    hd->SetDirectory(nullptr);
    for (int ix=1; ix<=hd->GetNbinsX(); ++ix)
        for (int iy=1; iy<=hd->GetNbinsY(); ++iy)
            hd->SetBinContent(ix, iy, hd->GetBinContent(ix, iy) - chi2min);
    return hd;
}

static void GetTruthPointDeg(const TString& order, double& theta23_deg, double& deltaCP_deg) {
    const bool no = (order != "_io");
    const double sin2_theta23 = no ? 0.558 : 0.553;
    const double deltaCP      = no ? 180.0 : 270.0;

    const double theta23_rad = std::asin(std::sqrt(sin2_theta23));
    theta23_deg = theta23_rad * 180.0 / TMath::Pi();
    deltaCP_deg = deltaCP;
}

static void draw_overlay_new(TH2* raw_chi2_fit,          // SM-only fit surface
                             TCanvas* c,
                             const TString& order,
                             const TString& title_text,
                             const TString& eps_text,
                             const TString& unc_text,
                             const TString& out_png,
                             int dof = 2)
{
    if (!raw_chi2_fit) { ::Error("draw_overlay_new","raw_chi2_fit is null"); return; }

    c->cd();
    c->Clear();

    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.09);

    int ixmin, iymin, izmin;
    raw_chi2_fit->GetBinXYZ(raw_chi2_fit->GetMinimumBin(), ixmin, iymin, izmin);
    const double chi2min = raw_chi2_fit->GetBinContent(ixmin, iymin, izmin);
    const double theta_hat   = raw_chi2_fit->GetXaxis()->GetBinCenter(ixmin);
    const double deltacp_hat = raw_chi2_fit->GetYaxis()->GetBinCenter(iymin);

    std::cout << "[BESTFIT] (Fit: SM only) chi2min=" << chi2min
              << "  at (theta23, deltaCP)=(" << theta_hat << ", " << deltacp_hat << ")\n";

    double theta_true, delta_true;
    GetTruthPointDeg(order, theta_true, delta_true);
    std::cout << "[TRUTH] (Data: SM+NSI) (theta23, deltaCP)=(" << theta_true << ", " << delta_true << ")\n";

    int ix_t = raw_chi2_fit->GetXaxis()->FindBin(theta_true);
    int iy_t = raw_chi2_fit->GetYaxis()->FindBin(delta_true);

    const double chi2_truth = raw_chi2_fit->GetBinContent(ix_t, iy_t);
    const double dchi2_bias = chi2_truth - chi2min;

    std::cout << "[BIAS] Δχ²_bias = χ²(truth) − χ²_min = "
            << dchi2_bias << std::endl;

    std::unique_ptr<TH2> hd( MakeDeltaChi2Hist(raw_chi2_fit) );

    const double X_MIN_FORCE = 40.0;
    const double X_MAX_FORCE = 52.0;
    const double Y_MIN_FORCE = 0.0;
    const double Y_MAX_FORCE = 360.0;
 
    TH2F* hframe = new TH2F("hframe_overlay_new","",
                            1, X_MIN_FORCE, X_MAX_FORCE,
                            1, Y_MIN_FORCE, Y_MAX_FORCE);
    hframe->SetDirectory(nullptr);
    hframe->SetStats(0);
    hframe->GetXaxis()->SetTitle("#theta_{23} (deg)");
    hframe->GetYaxis()->SetTitle("#delta_{CP} (deg)");
    hframe->GetXaxis()->SetNdivisions(510);
    hframe->GetYaxis()->SetNdivisions(510);
    hframe->GetXaxis()->SetTitleSize(0.045);
    hframe->GetYaxis()->SetTitleSize(0.045);
    hframe->GetXaxis()->SetLabelSize(0.040);
    hframe->GetYaxis()->SetLabelSize(0.040);
    hframe->Draw();

    const double lev68 = TMath::ChisquareQuantile(0.6827, dof);
    const double lev90 = TMath::ChisquareQuantile(0.90,   dof);
    const double lev95 = TMath::ChisquareQuantile(0.95,   dof);

    auto seg68 = ExtractContourGraphs(hd.get(), lev68);
    auto seg90 = ExtractContourGraphs(hd.get(), lev90);
    auto seg95 = ExtractContourGraphs(hd.get(), lev95);

    std::cout << "[DIAG] Fit contour seg counts: 68=" << seg68.size()
              << " 90=" << seg90.size()
              << " 95=" << seg95.size() << "\n";

    c->cd(); gPad->cd();

    for (auto* g : seg68) { if (!g) continue; g->SetLineColor(kRed);     g->SetLineWidth(3); g->SetLineStyle(1); g->Draw("L SAME"); }
    for (auto* g : seg90) { if (!g) continue; g->SetLineColor(kBlack);   g->SetLineWidth(3); g->SetLineStyle(2); g->Draw("L SAME"); }
    for (auto* g : seg95) { if (!g) continue; g->SetLineColor(kGreen+2); g->SetLineWidth(3); g->SetLineStyle(3); g->Draw("L SAME"); }

    auto m_truth = new TMarker(theta_true, delta_true, kFullCircle);
    m_truth->SetMarkerColor(kBlue+2);
    m_truth->SetMarkerSize(2.0);
    m_truth->Draw("SAME");

    auto m_best = new TMarker(theta_hat, deltacp_hat, kFullStar);
    m_best->SetMarkerColor(kBlack);
    m_best->SetMarkerSize(1.5);
    m_best->Draw("SAME");

    {
        auto ttl = new TLatex(0.50, 0.96, title_text);
        ttl->SetNDC();
        ttl->SetTextAlign(22);
        ttl->SetTextSize(0.045);
        ttl->Draw();
    }

    const double xN = 0.20;
    const double yN = 0.88;
    const double dy = 0.04;
    const int    kF = 42;

    auto t0 = new TLatex(xN, yN, Form("Normal Ordering, assuming no CPV"));
    t0->SetNDC(); t0->SetTextFont(kF); t0->SetTextSize(0.026); t0->SetTextAlign(13);
    t0->Draw();

    auto tU = new TLatex(xN, yN - 1*dy, unc_text);
    tU->SetNDC(); tU->SetTextFont(kF); tU->SetTextSize(0.026); tU->SetTextAlign(13);
    tU->Draw();

    auto tE = new TLatex(xN, yN - 2*dy, "10^{6} events");
    tE->SetNDC(); tE->SetTextFont(kF); tE->SetTextSize(0.026); tE->SetTextAlign(13);
    tE->Draw();

    auto leg = new TLegend(0.19, 0.16, 0.42, 0.40);
    leg->SetTextSize(0.022);
    leg->SetMargin(0.18);

    {
        auto d68 = new TGraph(); d68->SetLineColor(kRed);     d68->SetLineStyle(1); d68->SetLineWidth(3);
        auto d90 = new TGraph(); d90->SetLineColor(kBlack);   d90->SetLineStyle(2); d90->SetLineWidth(3);
        auto d95 = new TGraph(); d95->SetLineColor(kGreen+2); d95->SetLineStyle(3); d95->SetLineWidth(3);

        leg->AddEntry(d68, "68% C.L.", "L");
        leg->AddEntry(d90, "90% C.L.", "L");
        leg->AddEntry(d95, "95% C.L.", "L");
    }
    leg->AddEntry(m_truth, "Truth (SM + NSI)", "P");
    leg->AddEntry(m_best,  "Best fit (SM only)",   "P");
    leg->Draw();

    gPad->RedrawAxis();
    hframe->Draw("AXIS SAME");
    gPad->Modified();
    gPad->Update();

    c->SaveAs(out_png);

    delete hframe;
}

auto plot_analysis_multiple(TString order) -> void {
    gStyle->Reset("Plain");
    gStyle->SetLabelColor(kBlack, "XY");
    gStyle->SetTitleColor(kBlack, "XY");
    gStyle->SetAxisColor(kBlack, "XY");
    gStyle->SetLabelSize(0.04,  "XY");
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gROOT->ForceStyle();

    const TString fit_file = "analysis-0.010000-0.010000.root"; 
    TFile* f_fit = TFile::Open(fit_file);
    if (!f_fit || f_fit->IsZombie()) {
        std::cerr << "[ERR] cannot open " << fit_file << "\n";
        return;
    }
    TH2* chi2_raw_fit = f_fit->Get<TH2>("raw_chi2");
    if (!chi2_raw_fit) {
        std::cerr << "[ERR] cannot find TH2 'raw_chi2' in " << fit_file << "\n";
        delete f_fit;
        return;
    }

    TString title_text = "";
    TString eps_text   = "#epsilon_{R} = 4 #times 10^{-2}";
    TString unc_text   = "1% flux and cross-section uncertainties";
    TString out_png    = Form("p_bias_R_0.04%s.pdf", order.Data());

    auto c = new TCanvas("atn_ana","atn_ana", (int)(1400/0.98), 1400);
    c->cd();
    draw_overlay_new(chi2_raw_fit, c, order, title_text, eps_text, unc_text, out_png, 2);

    delete f_fit;
}
