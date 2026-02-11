#include "TFile.h"
#include "TH1D.h"

#include <algorithm>
#include <iostream>
#include <unordered_set>

auto smear_model(TString flavor, TString order) -> void {
    const int nTheta{100};
    const int nDeltaCP{100};

    TFile model_file{"../flux_model_xsec/atnu_" + flavor + order + "_model_xsec.root"};
    TFile resolution_file{"resolution_model_nu_" + flavor + ".root"};
    TFile output_file{"atnu_" + flavor + order + "_model_xsec_res.root", "RECREATE",
                      "", ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose};

    const auto model_template{model_file.Get<TH1>("atn_1_1")};
    const auto n_lgE{model_template->GetNbinsX()};
    const auto lgE_min{model_template->GetXaxis()->GetBinLowEdge(1)};
    const auto lgE_max{model_template->GetXaxis()->GetBinUpEdge(model_template->GetNbinsX())};
    // calculate bin range where free of marginal effect
    std::unordered_set<int> underflowed_lgE_bin;
    std::unordered_set<int> overflowed_lgE_bin;
    for (auto k{1}; k <= n_lgE; ++k) {
        std::unique_ptr<TH1> res_ker{resolution_file.Get<TH1>(TString::Format("res_ker_%d", k))};
        const auto n_ker_lgE{res_ker->GetNbinsX()};
        const auto m_c{res_ker->FindBin(0)};
        for (auto m{1}; m <= n_ker_lgE; ++m) {
            const auto k_prime{k + (m - m_c)};
            if (k_prime < 1) { underflowed_lgE_bin.emplace(k); }
            if (k_prime > n_lgE) { overflowed_lgE_bin.emplace(k); }
        }
    }
    const auto last_underflowed_lgE_bin{*std::ranges::max_element(underflowed_lgE_bin)};
    const auto first_overflowed_lgE_bin{*std::ranges::min_element(overflowed_lgE_bin)};
    const auto new_n_lgE{first_overflowed_lgE_bin - last_underflowed_lgE_bin - 1};
    const auto new_lgE_min{model_template->GetXaxis()->GetBinUpEdge(last_underflowed_lgE_bin)};
    const auto new_lgE_max{model_template->GetXaxis()->GetBinLowEdge(first_overflowed_lgE_bin)};
    delete model_template;

    for (auto i{1}; i <= nTheta; i++) {
        // std::cout << i << ": " << std::flush;
        for (auto j{1}; j <= nDeltaCP; j++) {
            std::unique_ptr<TH1> original_model{model_file.Get<TH1>(TString::Format("atn_%d_%d", i, j))};
            TH1D smeard_model{"htemp", "htemp", n_lgE, lgE_min, lgE_max};
            for (auto k{1}; k <= n_lgE; ++k) {
                const auto original_flux{original_model->GetBinContent(k)};
                std::unique_ptr<TH1> res_ker{resolution_file.Get<TH1>(TString::Format("res_ker_%d", k))};
                const auto n_ker_lgE{res_ker->GetNbinsX()};
                const auto m_c{res_ker->FindBin(0.)};
                for (auto m{1}; m <= n_ker_lgE; ++m) {
                    const auto bin{k + (m - m_c)};
                    const auto updated_flux{smeard_model.GetBinContent(bin) + original_flux * res_ker->GetBinContent(m)};
                    smeard_model.SetBinContent(bin, updated_flux);
                }
            }
            
            TH1D truncated_smeard_model{original_model->GetName(), original_model->GetTitle(),
                                        new_n_lgE, new_lgE_min, new_lgE_max};
            for (auto k{1}; k <= new_n_lgE; ++k) {
                const auto flux{smeard_model.GetBinContent(last_underflowed_lgE_bin + k)};
                truncated_smeard_model.SetBinContent(k, flux);
            }
            truncated_smeard_model.Write();
            // std::cout << j << ' ' << std::flush;
        }
        // std::cout << std::endl;
    }
}
