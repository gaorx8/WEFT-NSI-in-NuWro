#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

auto analysis(double model_unc, double xsec_unc, TString order) -> void {
    const double theta_min{30};
    const double theta_max{60};
    const int n_theta{100};
    const double delta_cp_min{0};
    const double delta_cp_max{360};
    const int n_delta_cp{100};

    // -------------------------
    // Open input files
    // -------------------------
    TFile model_file{"flux_model_xsec_res/atnu_mu_all" + order + "_model_xsec_res.root"};
    if (model_file.IsZombie()) {
        std::cerr << "[ERR] cannot open model file\n";
        return;
    }

    TFile data_file{"asimov_dataset.root"};
    if (data_file.IsZombie()) {
        std::cerr << "[ERR] cannot open data file\n";
        return;
    }

    // NOTE: use fixed-format filename to avoid mismatched names
    const TString out_name = Form("analysis-%0.6f-%0.6f.root", model_unc, xsec_unc);
    TFile output_file{out_name, "RECREATE", "",
                      ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose};
    if (output_file.IsZombie()) {
        std::cerr << "[ERR] cannot create output file " << out_name << "\n";
        return;
    }

    // -------------------------
    // Load data histogram (do NOT own it via unique_ptr)
    // -------------------------
    TH1* data = data_file.Get<TH1>("data");
    if (!data) {
        std::cerr << "[ERR] cannot find TH1 'data' in asimov_dataset.root\n";
        return;
    }

    const auto Cut = [&](int /*bin*/) {
        // return data->GetXaxis()->GetBinUpEdge(bin) > 0;
        return true;
    };

    const double n_event = data->Integral();
    double n_event_used{};

    for (int k{1}; k <= data->GetNcells(); ++k) {
        if (Cut(k)) n_event_used += data->GetBinContent(k);
    }

    // -------------------------
    // Create output histograms
    // IMPORTANT: lock output directory, then detach hist from directory to avoid gDirectory issues
    // -------------------------
    output_file.cd();

    TH2D log_likelihood{"log_likelihood", "log_likelihood",
                        n_theta, theta_min, theta_max,
                        n_delta_cp, delta_cp_min, delta_cp_max};

    TH2D raw_chi2{"raw_chi2", "raw_chi2",
                  n_theta, theta_min, theta_max,
                  n_delta_cp, delta_cp_min, delta_cp_max};

    // Detach from any TFile directory ownership (robust against gDirectory switching)
    log_likelihood.SetDirectory(nullptr);
    raw_chi2.SetDirectory(nullptr);

    // -------------------------
    // Scan grid
    // -------------------------
    double best_chi2 = std::numeric_limits<double>::infinity();
    int best_i = -1, best_j = -1;

    for (int i = 1; i <= n_theta; ++i) {
        for (int j = 1; j <= n_delta_cp; ++j) {

            // --- Get model hist pointer from model_file (owned by model_file) ---
            TH1* htmp = model_file.Get<TH1>(TString::Format("atn_%d_%d", i, j));
            if (!htmp) {
                std::cerr << "[ERR] cannot find model hist atn_" << i << "_" << j
                          << " in model file\n";
                return;
            }

            // --- Clone to local owned hist so we can scale safely and independent of file ---
            std::unique_ptr<TH1> model{ (TH1*)htmp->Clone(Form("model_%d_%d_tmp", i, j)) };
            model->SetDirectory(nullptr);

            if (model->GetNcells() != data->GetNcells()) {
                std::cerr << "[ERR] model->GetNcells() != data->GetNcells()\n";
                return;
            }

            const double model_int = model->Integral();
            if (model_int <= 0) {
                // Avoid divide-by-zero; set huge chi2
                const auto paramIndex = log_likelihood.GetBin(i, j);
                log_likelihood.SetBinContent(paramIndex, -1e100);
                log_likelihood.SetBinError(paramIndex, 0);
                raw_chi2.SetBinContent(paramIndex, 1e100);
                raw_chi2.SetBinError(paramIndex, 0);
                continue;
            }

            // Normalize model to data total events
            model->Scale(n_event / model_int);

            // --- compute binned Gaussian chi2 with combined (stat + sys) variance ---
            double ll{};
            for (int k = 1; k <= data->GetNcells(); ++k) {
                if (!Cut(k)) continue;

                const double n  = data->GetBinContent(k);
                const double mu = model->GetBinContent(k);

                const double sigma_n  = data->GetBinError(k);
                const double sigma_mu = mu * TMath::Hypot(xsec_unc, model_unc);

                const double denom = TMath::Sq(sigma_n) + TMath::Sq(sigma_mu);
                if (denom <= 0) continue;

                const double numer = TMath::Sq(n - mu);
                ll -= (numer / denom) / 2.0; // ll = -0.5*chi2
            }

            const double chi2_here = -2.0 * ll;

            if (chi2_here < best_chi2) {
                best_chi2 = chi2_here;
                best_i = i;
                best_j = j;
            }

            const auto paramIndex = log_likelihood.GetBin(i, j);
            log_likelihood.SetBinContent(paramIndex, ll);
            log_likelihood.SetBinError(paramIndex, 0);

            raw_chi2.SetBinContent(paramIndex, chi2_here);
            raw_chi2.SetBinError(paramIndex, 0);
        }
    }

    // -------------------------
    // Print best fit + save best-fit spectrum (optional)
    // -------------------------
    if (best_i > 0 && best_j > 0) {
        const double theta_hat   = log_likelihood.GetXaxis()->GetBinCenter(best_i);
        const double delta_hat   = log_likelihood.GetYaxis()->GetBinCenter(best_j);

        int ixmin, iymin, izmin;
        raw_chi2.GetBinXYZ(raw_chi2.GetMinimumBin(), ixmin, iymin, izmin);
        const double chi2min_hist   = raw_chi2.GetBinContent(ixmin, iymin, izmin);
        const double theta_hat_hist = raw_chi2.GetXaxis()->GetBinCenter(ixmin);
        const double delta_hat_hist = raw_chi2.GetYaxis()->GetBinCenter(iymin);

        std::cout << "[BESTFIT] scan-min: chi2min=" << best_chi2
                  << "  (best_i,best_j)=(" << best_i << "," << best_j << ")"
                  << "  (theta23,deltaCP)=(" << theta_hat << "," << delta_hat << ")\n";

        std::cout << "[BESTFIT] hist-min: chi2min=" << chi2min_hist
                  << "  (ixmin,iymin)=(" << ixmin << "," << iymin << ")"
                  << "  (theta23,deltaCP)=(" << theta_hat_hist << "," << delta_hat_hist << ")\n";

        // Save best-fit model spectrum PNG (still optional)
        TH1* hbest_src = model_file.Get<TH1>(Form("atn_%d_%d", best_i, best_j));
        if (hbest_src) {
            std::unique_ptr<TH1> hbest{ (TH1*)hbest_src->Clone("best_model") };
            hbest->SetDirectory(nullptr);

            TCanvas c("bestfit","bestfit",1200,800);
            hbest->SetLineWidth(2);
            hbest->Draw("HIST");
            // c.SaveAs(Form("best_model_%d_%d.png", best_i, best_j));
        }
    } else {
        std::cout << "[WARN] Best-fit not found (indices invalid).\n";
    }

    // -------------------------
    // Write outputs (critical: cd to output_file, overwrite, flush, close)
    // -------------------------
    output_file.cd();
    log_likelihood.Write("log_likelihood", TObject::kOverwrite);
    raw_chi2.Write("raw_chi2", TObject::kOverwrite);

    output_file.Write();
    output_file.Close();

    std::cout << "Output saved to: " << out_name << "\n";
    std::cout << "N event of interest: " << n_event_used << std::endl;
}
