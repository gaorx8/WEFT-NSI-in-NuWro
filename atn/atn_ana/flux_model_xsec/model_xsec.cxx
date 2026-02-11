#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "TROOT.h"
#include "TError.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cctype>
#include <cmath>
#include <memory>
#include <numeric>

static inline std::string trim(const std::string& s){
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a==std::string::npos) return "";
    return s.substr(a, b-a+1);
}
static std::vector<std::string> split_ws(const std::string& line){
    std::vector<std::string> v; std::stringstream ss(line); std::string t;
    while (ss >> t) v.push_back(t);
    return v;
}
static inline bool is_unit_token(const std::string& s){ 
    return !s.empty() && s.front()=='[' && s.back()==']';
}
static inline std::string col_p1(int case_id){ 
    return "sigma_case" + std::to_string(case_id) + "_p1";
}

struct XSecTable {
    std::vector<double> E_GeV;                                
    std::map<std::string, std::vector<double>> cols;       

    double interp(const std::string& col, double E) const {
        auto it = cols.find(col);
        if (it == cols.end()) return 0.0;
        const auto& Y = it->second;
        if (E <= E_GeV.front()) return Y.front();
        if (E >= E_GeV.back())  return Y.back();
        auto hi = std::upper_bound(E_GeV.begin(), E_GeV.end(), E);
        size_t i = std::max<size_t>(1, hi - E_GeV.begin()) - 1;
        double x0 = E_GeV[i], x1 = E_GeV[i+1];
        double y0 = Y[i],     y1 = Y[i+1];
        double t = (E - x0) / std::max(1e-300, x1 - x0);
        return y0 + t*(y1 - y0);
    }
};

static bool load_xsec_txt_p1(const char* path, XSecTable& tab){
    std::ifstream fin(path);
    if (!fin.is_open()){
        std::cerr << "❌ cannot open " << path << "\n";
        return false;
    }
    std::string header;
    if (!std::getline(fin, header)){
        std::cerr << "❌ empty file\n";
        return false;
    }
    auto head = split_ws(header);
    if (head.size() < 2){
        std::cerr << "❌ header needs >=2 columns (E + at least one xsec)\n";
        return false;
    }

    bool energy_is_MeV = false;
    {
        std::string H = header; for (auto &c : H) c = std::toupper(c);
        if (H.find("MEV") != std::string::npos) energy_is_MeV = true;
    }

    std::vector<std::string> colnames;
    for (size_t i=1; i<head.size(); ++i){
        if (is_unit_token(head[i])) continue; 
        tab.cols[head[i]] = {};
        colnames.push_back(head[i]);
    }

    std::string line;
    while (std::getline(fin, line)){
        line = trim(line);
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.size() < 1 + colnames.size()) continue; 

        double E = std::stod(f[0]);
        if (energy_is_MeV) E /= 1000.0; 
        tab.E_GeV.push_back(E);

        for (size_t j=0; j<colnames.size(); ++j){
            double v = std::stod(f[1+j]);
            tab.cols[colnames[j]].push_back(v);
        }
    }
    if (tab.E_GeV.size() < 2){
        std::cerr << "❌ not enough rows\n";
        return false;
    }

    std::vector<size_t> idx(tab.E_GeV.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
        [&](size_t a,size_t b){ return tab.E_GeV[a] < tab.E_GeV[b]; });

    bool sorted = true;
    for (size_t i=0; i<idx.size(); ++i) if (idx[i] != i) { sorted = false; break; }
    if (!sorted){
        std::vector<double> E2; E2.reserve(idx.size());
        for (auto id: idx) E2.push_back(tab.E_GeV[id]);
        tab.E_GeV.swap(E2);
        for (auto &kv : tab.cols){
            const auto &Y = kv.second;
            std::vector<double> Y2; Y2.reserve(idx.size());
            for (auto id: idx) Y2.push_back(Y[id]);
            kv.second.swap(Y2);
        }
    }
    return true;
}

enum CombinedCase { R_case=0, S0_case=1, S1_case=2, P_case=3, T_case=4 };

static bool check_required_columns(const XSecTable& tab, CombinedCase cc, std::string& miss){
    auto has = [&](const std::string& c){ return tab.cols.count(c)>0; };
    if (!has(col_p1(0))) { miss = col_p1(0); return false; } // LL
    switch(cc){
        case R_case:  if (!has(col_p1(2))){ miss=col_p1(2); return false; } // LR
                      if (!has(col_p1(1))){ miss=col_p1(1); return false; } // RR
                      break;
        case S0_case: if (!has(col_p1(4))){ miss=col_p1(4); return false; } // LS0
                      if (!has(col_p1(3))){ miss=col_p1(3); return false; } // SS0
                      break;
        case S1_case: if (!has(col_p1(6))){ miss=col_p1(6); return false; } // LS1
                      if (!has(col_p1(5))){ miss=col_p1(5); return false; } // SS1
                      break;
        case P_case:  if (!has(col_p1(8))){ miss=col_p1(8); return false; } // LP
                      if (!has(col_p1(7))){ miss=col_p1(7); return false; } // PP
                      break;
        case T_case:  if (!has(col_p1(10))){ miss=col_p1(10); return false; } // LT
                      if (!has(col_p1(9))){  miss=col_p1(9);  return false; } // TT
                      break;
        default: break;
    }
    return true;
}

static double combined_sigma_p1(const XSecTable& tab, CombinedCase cc, double eps, double E, bool allow_negative){
    const double LL = tab.interp(col_p1(0), E);
    double s = LL;
    switch (cc){
        case R_case:  s = LL + eps*tab.interp(col_p1(2),E) + eps*eps*tab.interp(col_p1(1),E);  break;
        case S0_case: s = LL + eps*tab.interp(col_p1(4),E) + eps*eps*tab.interp(col_p1(3),E);  break;
        case S1_case: s = LL + eps*tab.interp(col_p1(6),E) + eps*eps*tab.interp(col_p1(5),E);  break;
        case P_case:  s = LL + eps*tab.interp(col_p1(8),E) + eps*eps*tab.interp(col_p1(7),E);  break;
        case T_case:  s = LL + eps*tab.interp(col_p1(10),E)+ eps*eps*tab.interp(col_p1(9),E);  break;
        default: break;
    }
    if (!std::isfinite(s)) s = 0.0;
    if (!allow_negative) s = std::max(0.0, s);
    return s;
}

static void apply_one(const std::string& input_root,
                      const std::string& output_root,
                      const std::string& xsec_txt,
                      int combined_case,
                      double epsilon,
                      int nTheta,
                      int nDeltaCP,
                      bool allow_negative)
{
    XSecTable tab;
    if (!load_xsec_txt_p1(xsec_txt.c_str(), tab)) {
        std::cerr << "❌ Failed to read txt: " << xsec_txt << "\n";
        return;
    }

    CombinedCase CC = static_cast<CombinedCase>(combined_case);
    std::string miss;
    if (!check_required_columns(tab, CC, miss)){
        std::cerr << "❌ Missing column: " << miss << "  in " << xsec_txt << "\n";
        return;
    }

    const double logEmin = std::log10(0.2);
    const double logEmax = std::log10(10.0);

    TF1 fXsec("fXsec",
        [&](const double* x, const double*) {
            double E = std::pow(10.0, x[0]);          // GeV
            E = std::min(std::max(E, tab.E_GeV.front()), tab.E_GeV.back());

            double sig_cm2_nucleon = combined_sigma_p1(tab, CC, epsilon, E, allow_negative);
            double sig_m2_O  = sig_cm2_nucleon * 1e-4;
            return sig_m2_O;
        },
        logEmin, logEmax, 0
    );

    TFile model_file(input_root.c_str(), "READ");
    if (model_file.IsZombie()){
        std::cerr << "❌ cannot open input ROOT: " << input_root << "\n";
        return;
    }

    TFile output_file(output_root.c_str(),
                      "RECREATE", "",
                      ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);

    auto multiply_with_fun = [&](TH1* h){
        for (int ib=1; ib<=h->GetNbinsX(); ++ib){
            double x  = h->GetXaxis()->GetBinCenter(ib);
            double w  = fXsec.Eval(x);
            double yc = h->GetBinContent(ib);
            double ye = h->GetBinError(ib);
            h->SetBinContent(ib, yc * w);
            h->SetBinError(ib,   ye * std::fabs(w));
        }
    };

    for (int i=1; i<=nTheta; ++i){
        for (int j=1; j<=nDeltaCP; ++j){
            TString hname = TString::Format("atn_%d_%d", i, j);
            TH2* h2 = model_file.Get<TH2>(hname);
            if (!h2) continue;
            std::unique_ptr<TH1> projX{ h2->ProjectionX(hname) };
            const double dcos = h2->GetYaxis()->GetBinWidth(1);
            projX->Scale(2.0 * TMath::Pi() * dcos);
            multiply_with_fun(projX.get());
            projX->Write();
        }
    }

    output_file.Close();
    std::cout << "✅ Written: " << output_root << "\n";
}

void model_xsec_from_txt_p1_run(int combined_case = 0,      // 0=R,1=S0,2=S1,3=P,4=T
                                double epsilon   = 0.,
                                int nTheta       = 100,
                                int nDeltaCP     = 100,
                                bool allow_negative = true,
                                const char* in_nu_root   = "atnu_mu_model.root",
                                const char* in_nubar_root= "atnu_mu_bar_model.root",
                                const char* xsec_nu_txt  = "weft_sigma_mu.txt",
                                const char* xsec_nubar_txt="weft_sigma_mu_bar.txt")
{
    // νμ
    apply_one(in_nu_root,
              "atnu_mu_model_xsec.root",
              xsec_nu_txt,
              combined_case, epsilon, nTheta, nDeltaCP, allow_negative);

    // ν̄μ
    apply_one(in_nubar_root,
              "atnu_mu_bar_model_xsec.root",
              xsec_nubar_txt,
              combined_case, epsilon, nTheta, nDeltaCP, allow_negative);
}

void model_xsec()
{
    model_xsec_from_txt_p1_run(/*combined_case=*/0,
                              /*epsilon=*/0,
                              /*nTheta=*/100,
                              /*nDeltaCP=*/100,
                              /*allow_negative=*/true);
} 
 