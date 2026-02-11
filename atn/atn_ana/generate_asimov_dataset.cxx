#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"

constexpr auto pi{3.14159265358979323846};
constexpr auto deg{pi / 180};

auto generate_asimov_dataset(double n_event, TString order) -> void {
    const auto no{order != "_io"};
  
    const double sin2_theta{no ? 0.558 : 0.553};
    const double delta_cp{no ? 180 * deg : 270 * deg};
 
    const double theta_min{30 * deg};
    const double theta_max{60 * deg};
    const int n_theta{100};
    const double delta_cp_min{0 * deg};
    const double delta_cp_max{360 * deg};
    const int n_delta_cp{100};
    TH2D param_space{"param_space", "param_space",
                     n_theta, theta_min, theta_max,
                     n_delta_cp, delta_cp_min, delta_cp_max};
    const auto i{param_space.GetXaxis()->FindBin(std::asin(std::sqrt(sin2_theta)))};
    const auto j{param_space.GetYaxis()->FindBin(delta_cp)};

    TFile model_file{"flux_model_xsec_res/atnu_mu_all" + order + "_model_xsec_res_R_0.04.root"};
    TFile output_file{"asimov_dataset.root", "RECREATE",
                      "", ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose};

    std::unique_ptr<TH1> model{model_file.Get<TH1>(TString::Format("atn_%d_%d", i, j))};
    model->SetName("data");
    model->Scale(n_event / model->Integral());
    
    for (int k{1}; k <= model->GetNcells(); ++k) {
        model->SetBinError(k, std::sqrt(model->GetBinContent(k)));
    }
    model->Write();   
} 
 