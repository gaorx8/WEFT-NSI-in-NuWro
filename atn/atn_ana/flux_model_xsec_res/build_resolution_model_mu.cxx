#include "Math/PdfFunc.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"

#include <cmath>
#include <iostream>

constexpr auto pi{3.14159265358979323846};
constexpr auto deg{pi / 180};
constexpr auto ln10{2.30258509299404568402};

auto sigma_E(double E) -> double {
    // return 0.03 * std::sqrt(E * 1000) / 1000;
    return 0.043 * std::sqrt(E) + 0.011; // nu_mu
    // return 0.049 * std::sqrt(E) + 0.006; // nu_mu_bar
}

auto resolution_kernel_pdf(const double var[1], const double param[1]) -> double {
    const auto delta_lgE{var[0]};
    const auto lgE0{param[0]};

    const auto E0{std::pow(10, lgE0)};
    const auto E{E0 * std::pow(10, delta_lgE)};

    return ln10 * E * ROOT::Math::gaussian_pdf(E - E0, sigma_E(E0));
}

auto build_resolution_model_mu() -> void {
    TFile model_file{"../flux_model_xsec/atnu_mu_model_xsec.root"};
    std::unique_ptr<TH1> model{model_file.Get<TH1>("atn_1_1")};
    const auto n_lgE{model->GetNbinsX()};
    const auto lgE_axis{model->GetXaxis()};
    const auto lgE_min{lgE_axis->GetBinLowEdge(1)};
    const auto lgE_max{lgE_axis->GetBinUpEdge(model->GetNbinsX())};
    const auto width_lgE{lgE_axis->GetBinWidth(1)};

    TFile output_file{"resolution_model_nu_mu.root", "RECREATE",
                      "", ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose};

    const auto delta_lgE_bound{std::max(std::abs(lgE_min), std::abs(lgE_max))};
    TF1 res_ker_pdf{"resolution_kernel_pdf", resolution_kernel_pdf,
                    -delta_lgE_bound, delta_lgE_bound, 1};

    constexpr auto kernel_range_factor{4}; // n sigma range
 
    for (auto i{1}; i <= n_lgE; ++i) {
        // bin center
        const auto lgE0{lgE_axis->GetBinCenter(i)};
        res_ker_pdf.SetParameter(0, lgE0);
        // calculate lgE kernel range
        const auto E0{std::pow(10, lgE0)};
        const auto sigma_lgE0{sigma_E(E0) / (E0 * ln10)};
        const auto ker_lgE_index_range{std::max(2, static_cast<int>(kernel_range_factor * sigma_lgE0 / width_lgE) + 1)};
        const auto n_ker_lgE{2 * ker_lgE_index_range + 1};
        const auto ker_lgE_range{(ker_lgE_index_range + 0.5) * width_lgE};

        TH1D kernel{TString::Format("res_ker_%d", i),
                    "Resolution kernel;#Delta log_{10}(E)",
                    n_ker_lgE, -ker_lgE_range, ker_lgE_range};
        for (auto k{1}; k <= n_ker_lgE; ++k) {
            const auto delta_lgE_low{kernel.GetXaxis()->GetBinLowEdge(k)};
            const auto delta_lgE_up{kernel.GetXaxis()->GetBinUpEdge(k)};
            const auto p{res_ker_pdf.Integral(delta_lgE_low, delta_lgE_up)};
            kernel.SetBinContent(k, p);
            kernel.SetBinError(k, 0);
        }
        kernel.Scale(1 / kernel.Integral());
        kernel.Write();

        std::cout << i << ' ' << std::flush;
    }
    std::cout << std::endl;
}
