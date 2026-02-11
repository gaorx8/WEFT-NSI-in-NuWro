#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH2.h"
#include "TMarker.h"

#include <iostream>

auto plot_analysis_single() -> void {
    const auto ll{TFile::Open("analysis-0.020000-0.020000.root")->Get<TH2>("log_likelihood")};

    int i, j, k;
    ll->GetBinXYZ(ll->GetMaximumBin(), i, j, k);
    const auto theta_hat{ll->GetXaxis()->GetBinCenter(i)};
    const auto delta_m2_hat{ll->GetYaxis()->GetBinCenter(j)};
    std::cout << theta_hat << '\n'
              << delta_m2_hat << '\n';

    const auto c(new TCanvas{"atn_ana", "atn_ana", static_cast<int>(800 / 0.95), 800});
    c->SetLeftMargin(0.15);

    ll->SetTitle("Atmospheric neutrino mixing parameters;#theta_{23} (deg);#delta_{CP} (deg)");
    ll->SetStats(false);
    TGaxis::SetMaxDigits(3);

    ll->SetContour(3);
    ll->SetLineColor(kBlue);
    ll->SetContourLevel(2, -0.5);
    ll->SetContourLevel(1, -1.0);
    ll->SetContourLevel(0, -1.5);

    ll->Draw("CONT2");

    const auto best_fit(new TMarker{theta_hat, delta_m2_hat, kFullStar});
    best_fit->SetMarkerColor(kBlue);
    best_fit->Draw("SAME");
}
