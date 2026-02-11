#include <fstream>
#include <random>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "vect.h"
#include "vec.h"
#include "ff.h"

#include "event1.h"
#include "kinematics.h"
#include "nucleus.h"
#include "params.h"
#include "sf/CSFOptions.h"
#include "sf/CSpectralFunc.h"
#include "sf/GConstants.h"
#include "sfevent_weft.h"
#include "sfevent.h"
#include <iomanip> 
  
static inline vect on_shell_4vec(particle& N, double M) {
    const double p2 = N.p() * N.p();
    const double E  = std::sqrt(M * M + p2);
    return vect(E, N.p()); 
}

static inline double safe_sqrt(double x) {
    return std::sqrt(std::max(0.0, x));
}

static inline double pow2(double x) { return x * x; }

double V_ud = 0.97154;
double vac = 246220;
double D_0 = 2 * pow(V_ud, 2) / pow(vac, 4);

double m_n = 939.56542052;
double m_p = 938.27208816;
double M_N = (m_n + m_p) / 2;

double m_pi = 139.57039;
double M_2 = pow(M_N, 2);
double MV2_0 = 840 * 840;
double mu_p_0 = 2.79284;
double mu_n_0 = -1.91304;

double sfevent_weft(params &par, event &e, nucleus &t) {
    particle &l0 = e.in[0];  
    particle &N0 = e.in[1]; 
    double E_v = l0.E();

    const bool is_anti = l0.pdg < 0;               
    const bool is_on_n = N0.pdg == pdg_neutron;   
    const bool is_cc_possible = is_anti xor is_on_n;
    
    if (e.flag.cc and not is_cc_possible) return 0;

    particle l1;  
    particle N1; 

    N1.r = N0.r;

    is_on_n xor e.flag.cc ? N1.set_neutron() : N1.set_proton();

    l1.pdg = l0.pdg;
    if (e.flag.cc) l1.pdg += is_anti ? 1 : -1;

    const double m_l = mass(l1.pdg);
    const double M_i = N0.mass();
    const double M_f = N1.mass();
    l1.set_mass(m_l);

    double p = 0;  
    double E = 0;  

    CSFOptions options(par, e.flag.cc, !is_on_n, is_anti); 
    CSpectralFunc *sf = options.get_SF();     

    if (par.weft_param != 0) {             
        p = sf->MomDist()->generate();  
        E = get_E(sf, p);  
    } 

    N0.set_momentum(rand_dir() * p);
    N0.t = N0.mass() - E;
    vect p_i_eff = N0.p4();

    vect s = l0.p4() + N0.p4();

    const double s2 = s * s;

    if (s2 < pow2(M_f + m_l)) return 0;  

    const vec v = s.v(); 

    const double mom_cms = sqrt(0.25 * pow2(s2 + pow(m_l, 2) - pow(M_f, 2)) / s2 - pow(m_l, 2));
    const vec dir_cms = rand_dir();

    l1.set_momentum(mom_cms * dir_cms);
    N1.set_momentum(-l1.p());

    N0.boost2(v);

    double N0_p_length = N0.p().length();
    double Q2_a = pow(N0_p_length - mom_cms, 2) - pow(N1.E() - N0.E(), 2);
    double Q2_b = pow(N0_p_length + mom_cms, 2) - pow(N1.E() - N0.E(), 2);

    N0.boost(v);

    l1.boost(v);
    N1.boost(v);

    if (par.pauli_blocking) {
        if (par.sf_pb == 0 and N1.momentum() < sf->get_pBlock())
        return 0;
        else if (par.sf_pb == 1 and N1.momentum() < t.localkf(N1))
        return 0;
        else if (par.sf_pb == 2 and frandom() < sf->MomDist()->Tot(N1.momentum()) / sf->MomDist()->Tot(0))
        return 0;
    }
  
    vect p_v = l0.p4();
    vect p_l = l1.p4();
    vect p_i = on_shell_4vec(N0, M_i);  // initial nucleon 4-vector
    vect p_f = on_shell_4vec(N1, M_f);  // final   nucleon 4-vector
    // vect p_n = N0.p4();
    // vect p_p = N1.p4();

    vect q = p_f - p_i;
    double q2 = q * q;
    double Q2 = - q2;

    const double rel = 5e-4; 

    double s_num = pow(M_N, 2) + 2 * M_N * E_v;
    double Q2_min = (1 / (2 * E_v + M_N )) * (2 * M_N * pow(E_v, 2) - pow(m_l, 2) * (M_N + E_v) 
                - E_v * sqrt(- 2 * pow(M_N, 2) * (s_num + pow(m_l, 2)) + pow(s_num - pow(m_l, 2), 2) + pow(M_N, 4)));
    double Q2_max = (1 / (2 * E_v + M_N )) * (2 * M_N * pow(E_v, 2) - pow(m_l, 2) * (M_N + E_v) 
                + E_v * sqrt(- 2 * pow(M_N, 2) * (s_num + pow(m_l, 2)) + pow(s_num - pow(m_l, 2), 2) + pow(M_N, 4)));

    double a2 = (1.0 + Q2 / MV2_0) * (1.0 + Q2 / MV2_0); 
    double G_D = 1 / a2;
    double tau = Q2 / (4 * M_2);

    double GEp = 1.0 / a2;
    double GEn = 0;
    double GMp = mu_p_0 / a2;
    double GMn = mu_n_0 / a2;

    double Ge = GEp - GEn;
    double Gm = GMp - GMn;

    double F1 = (Ge + tau * Gm) / (1 + tau);
    double F2 = (Gm - Ge) / (1 + tau);

    double g_A = 1.2728; 
    double m_A = 961;
    double m_q = 3.410;

    double G_A = g_A / pow(1 + Q2 / pow(m_A, 2), 2);
    double G_t_P = - 4 * pow(M_N, 2) * G_A / (Q2 + pow(m_pi, 2));
    double G_P = M_N * G_A / m_q + Q2 * G_t_P / (4 * M_N * m_q);

    // LQCD
    if (par.weft_param == 2) {
        double tc = 9 * m_pi * m_pi;
        double t0 = -tc;
        double sqrt_tc_plus_Q2 = std::sqrt(tc + Q2);
        double sqrt_tc_minus_t0 = std::sqrt(tc - t0);
        double z = (sqrt_tc_plus_Q2 - sqrt_tc_minus_t0) / (sqrt_tc_plus_Q2 + sqrt_tc_minus_t0);
        
        std::vector<double> a_a = {1.009, -1.756, -1.059, 1.621, 3.919, -5.739, 2.005}; 
        double sum_a = 0.0;
        for (size_t k = 0; k < a_a.size(); ++k) {
            sum_a += a_a[k] * std::pow(z, k);
        }
        G_A = sum_a;

        std::vector<double> a_t_p = {1.008, -1.831, -1.713, 4.994, -1.522, -1.984, 1.047};   

        double sum_t_p = 0.0;
        for (size_t k = 0; k < a_t_p.size(); ++k) {
            sum_t_p += a_t_p[k] * std::pow(z, k);
        }
        G_t_P = -4 * M_N * M_N * sum_t_p / (Q2 + m_pi * m_pi);
        
        std::vector<double> a_p = {1.066, -1.461, -1.053, -2.504, 12.446, -12.260, 3.766};   
        double sum_p = 0.0;
        for (size_t k = 0; k < a_p.size(); ++k) {
            sum_p += a_p[k] * std::pow(z, k);
        }
        G_P = M_N * pow(m_pi, 2) * sum_p / (m_q * (Q2 + pow(m_pi, 2)));
    }

    // D2
    if (par.weft_param == 3) {
        double GA, GtP;
        list(GA, GtP)=fap(q2, 0);
        G_A = -GA; 
        G_t_P = - 4.0 * M_N * M_N * G_A / (Q2 + m_pi * m_pi);
        G_P = M_N * G_A / m_q + Q2 * G_t_P / (4 * M_N * m_q);
    }

    // S, P, T
    double G_V = F1;
    double delta_M_N_QCD = 2.580;
    double delta_m_q = 2.527;
    double G_S = - delta_M_N_QCD * G_V / delta_m_q;
    double G_t_S = 3 * M_N * G_V * (delta_M_N_QCD * g_A / (2 * M_N) - 3 * delta_m_q / (2 * m_p)) / m_p;

    double F1_u_0 = 0.784;
    double F1_d_0 = -0.204;
    double B_pi_u_0 = 0.195;
    double m_b1 = 1229;
    double m_h1 = 1166;

    double F2_p = (mu_p_0 - 1) * G_D / (1.0 + tau);
    double F2_n = mu_n_0 * G_D / (1.0 + tau);  

    double D_h1 = pow(m_h1, 2) / (pow(m_h1, 2) + Q2);
    double D_b1 = pow(m_b1, 2) / (pow(m_b1, 2) + Q2);

    double F1_u = 0.5 * (F1_u_0 - F1_d_0) * D_b1 + 0.5 * (F1_u_0 + F1_d_0) * D_h1;
    double F1_d = - 0.5 * (F1_u_0 - F1_d_0) * D_b1 + 0.5 * (F1_u_0 + F1_d_0) * D_h1;
    double F2_u = - M_N * B_pi_u_0 * (2 * GMp + GMn) / (2 * m_pi) + 2 * pow(M_N, 2) * F1_u_0 * D_b1 / pow(m_b1, 2);
    double F2_d = - M_N * B_pi_u_0 * (GMp + 2 * GMn) / (2 * m_pi) + 2 * pow(M_N, 2) * F1_d_0 * D_b1 / pow(m_b1, 2);
    double F3_u = M_N * B_pi_u_0 * (2 * F2_p + F2_n) / (4 * m_pi) - pow(M_N, 2) * F1_u_0 * D_b1 / pow(m_b1, 2);
    double F3_d = M_N * B_pi_u_0 * (F2_p + 2 * F2_n) / (4 * m_pi) - pow(M_N, 2) * F1_d_0 * D_b1 / pow(m_b1, 2);

    double G_T = F1_u - F1_d;
    double G_T_1 = F2_u - F2_d;
    double G_T_2 = F3_u - F3_d;

    double C0 = 4 * pow(V_ud, 2) * pow(M_N, 4) / pow(vac, 4);
    const double s_u = (p_v + p_i)*(p_v + p_i) - (p_v - p_f)*(p_v - p_f);

    const double MN2 = M_N * M_N;
    const double a   = Q2 / (4.0 * MN2);

    const double termA = (1.0 + a) * (G_A * G_A);

    const double termB = (1.0 - a) * ( (F1 * F1) - a * (F2 * F2) );

    const double termC = (Q2 / MN2) * (F1 * F2);

    const double termD = (m_l * m_l) / (4.0 * MN2) *
        ( (F1 + F2) * (F1 + F2)
        + (G_A - G_t_P) * (G_A - G_t_P)
        - (1.0 + a) * (G_t_P * G_t_P) );

    const double B_LL = ((m_l * m_l + Q2) / MN2) * ( termA - termB + termC - termD );
    
    double C_LL = (Q2 / MN2) * ( G_A * (F1 + F2) );    
    
    if (is_anti) C_LL = -C_LL;                 
    
    double D_LL = 0.25 * ( (G_A * G_A) + (F1 * F1) + a * (F2 * F2) ); 

    double A2_LL = C0 * (B_LL + C_LL * s_u / pow(M_N, 2) + D_LL * pow(s_u, 2) / pow(M_N, 4));


    // double A2_LL = D_0 * ((((p_l * p_n) * (p_v * p_p)) + ((p_l * p_p) * (p_v * p_n)))
    //                             * (pow((2 * F1 + F2), 2) + pow((G_t_P - 2 * G_A), 2) + (p_n * p_p) * (pow(F2, 2) - pow(G_t_P, 2)) / pow(M_N, 2))
    //                         - (((p_l * p_n) * (p_v * p_n)) + ((p_l * p_p) * (p_v * p_p))) 
    //                             * (4 * F1 * F2 + 3 * pow(F2, 2) + G_t_P * (G_t_P - 4 * G_A) - (p_n * p_p) * (pow(F2, 2) + pow(G_t_P, 2)) / pow(M_N, 2))
    //                         - 8 * (((p_l * p_n) * (p_v * p_p))  - ((p_l * p_p) * (p_v * p_n)))* (G_A * (F1 + F2))
    //                         - pow(M_N, 2) * (p_l * p_v) * ((1 - (p_n * p_p) / pow(M_N, 2)) * 
    //                                                     (4 * F1 * F2 + pow(F2, 2) + G_t_P * (4 * G_A - G_t_P) - (p_n * p_p) * (pow(F2, 2) - pow(G_t_P, 2)) / pow(M_N, 2)) + 4 * (pow(F1, 2) - pow(G_A, 2))));

    double B_RR = B_LL;     
    double C_RR = -C_LL;
    if (is_anti) C_RR = -C_RR;       
    double D_RR = D_LL;       
    double A2_RR = C0 * (B_RR + C_RR * s_u / pow(M_N, 2) + D_RR * pow(s_u, 2) / pow(M_N, 4));

    // double A2_RR = D_0 * ((((p_l * p_n) * (p_v * p_p)) + ((p_l * p_p) * (p_v * p_n)))
    //                             * (pow((2 * F1 + F2), 2) + pow((G_t_P - 2 * G_A), 2) + (p_n * p_p) * (pow(F2, 2) - pow(G_t_P, 2)) / pow(M_N, 2))
    //                         - (((p_l * p_n) * (p_v * p_n)) + ((p_l * p_p) * (p_v * p_p))) 
    //                             * (4 * F1 * F2 + 3 * pow(F2, 2) + G_t_P * (G_t_P - 4 * G_A) - (p_n * p_p) * (pow(F2, 2) + pow(G_t_P, 2)) / pow(M_N, 2))
    //                         + 8 * (((p_l * p_n) * (p_v * p_p))  - ((p_l * p_p) * (p_v * p_n)))* (G_A * (F1 + F2))
    //                         - pow(M_N, 2) * (p_l * p_v) * ((1 - (p_n * p_p) / pow(M_N, 2)) * 
    //                                                     (4 * F1 * F2 + pow(F2, 2) + G_t_P * (4 * G_A - G_t_P) - (p_n * p_p) * (pow(F2, 2) - pow(G_t_P, 2)) / pow(M_N, 2)) + 4 * (pow(F1, 2) - pow(G_A, 2))));

    double B_LR = ((pow(m_l, 2) + Q2) / pow(M_N, 2)) * (- (1 + Q2 / (4 * pow(M_N, 2))) * pow(G_A, 2) 
                - (1 - Q2 / (4 * pow(M_N, 2))) * (pow(F1, 2) - Q2 * pow(F2, 2) / (4 * pow(M_N, 2))) 
                + Q2 * F1 * F2 / pow(M_N, 2)
                - (pow(m_l, 2) / (4 * pow(M_N, 2))) * (pow(F1 + F2, 2) - pow(G_A - G_t_P, 2) + (1 + Q2 / (4 * pow(M_N, 2))) * pow(G_t_P, 2)));
    double C_LR = 0;
    double D_LR = (- pow(G_A, 2) + pow(F1, 2) + Q2 * pow(F2, 2) / (4 * pow(M_N, 2))) / 4;
    double A2_LR = C0 * (B_LR + C_LR * s_u / MN2 + D_LR * (s_u*s_u) / (MN2*MN2));

    // double A2_LR = D_0 * ((((p_l * p_n) * (p_v * p_p)) + ((p_l * p_p) * (p_v * p_n)))
    //                     * (pow((2 * F1 + F2), 2) - pow((G_t_P - 2 * G_A), 2) + (p_n * p_p) * (pow(F2, 2) + pow(G_t_P, 2)) / pow(M_N, 2))
    //                 - (((p_l * p_n) * (p_v * p_n)) + ((p_l * p_p) * (p_v * p_p))) 
    //                     * (4 * F1 * F2 + 3 * pow(F2, 2) - G_t_P * (G_t_P - 4 * G_A) - (p_n * p_p) * (pow(F2, 2) - pow(G_t_P, 2)) / pow(M_N, 2))
    //                 + pow(M_N, 2) * (p_l * p_v) * ((1 - (p_n * p_p) / pow(M_N, 2)) * 
    //                                             (4 * F1 * F2 + pow(F2, 2) - G_t_P * (4 * G_A - G_t_P) - (p_n * p_p) * (pow(F2, 2) + pow(G_t_P, 2)) / pow(M_N, 2)) + 4 * (pow(F1, 2) + pow(G_A, 2))));
    
    double B_SS0 = ((m_l * m_l + Q2) / MN2) * (1.0 + a) * (G_S * G_S); 
    double C_SS0 = 0.0;                                                
    double D_SS0 = 0.0;                                                
    double A2_SS0 = C0 * (B_SS0 + C_SS0 * s_u / pow(M_N, 2) + D_SS0 * pow(s_u, 2) / pow(M_N, 4));

    // double A2_SS0 = 2 * D_0 * pow(G_S, 2) * (p_l * p_v) * ((p_n * p_p) + pow(M_N, 2));

    double B_LS0 = 0;
    double C_LS0 = (m_l / (2.0 * M_N)) * (F1 - Q2 * F2 / (4 * pow(M_N, 2))) * G_S;
    if (is_anti) C_LS0 = -C_LS0;   
    double D_LS0 = 0;
    double A2_LS0 = C0 * (B_LS0 + C_LS0 * s_u / pow(M_N, 2) + D_LS0 * pow(s_u, 2) / pow(M_N, 4));

    // double A2_LS0 = D_0 * m_l * M_N * G_S * ((p_v * p_n) + (p_v * p_p)) * (2 * F1 + F2 * (1 - (p_n * p_p) / pow(M_N, 2)));

    G_S = - delta_M_N_QCD * G_V / delta_m_q + Q2 * G_t_S / (2 * M_N * delta_m_q);
    double B_SS1 = ((m_l * m_l + Q2) / MN2) * (1.0 + a) * (G_S * G_S); 
    double C_SS1 = 0.0;                                                
    double D_SS1 = 0.0;                                                
    double A2_SS1 = C0 * (B_SS1 + C_SS1 * s_u / pow(M_N, 2) + D_SS1 * pow(s_u, 2) / pow(M_N, 4));

    // double A2_SS1 = 2 * D_0 * pow(G_S, 2) * (p_l * p_v) * ((p_n * p_p) + pow(M_N, 2));

    // G_S = - delta_M_N_QCD * G_V / delta_m_q + Q2 * G_t_S / (2 * M_N * delta_m_q);
    double B_LS1 = 0;
    double C_LS1 = (m_l / (2.0 * M_N)) * (F1 - Q2 * F2 / (4 * pow(M_N, 2))) * G_S;
    if (is_anti) C_LS1 = -C_LS1;   
    double D_LS1 = 0;
    double A2_LS1 = C0 * (B_LS1 + C_LS1 * s_u / pow(M_N, 2) + D_LS1 * pow(s_u, 2) / pow(M_N, 4));

    // double A2_LS1 = D_0 * m_l * M_N * G_S * ((p_v * p_n) + (p_v * p_p)) * (2 * F1 + F2 * (1 - (p_n * p_p) / pow(M_N, 2)));

    double B_PP = ((m_l * m_l + Q2) / MN2) * (Q2 / (4.0 * MN2)) * (G_P * G_P);
    double C_PP = 0.0;
    double D_PP = 0.0;
    double A2_PP = C0 * (B_PP + C_PP * s_u / pow(M_N, 2) + D_PP * pow(s_u, 2) / pow(M_N, 4));

    // double A2_PP = 2 * D_0 * pow(G_P, 2) * (p_l * p_v) * ((p_n * p_p) - pow(M_N, 2));

    double B_LP = - (m_l / (2 * M_N)) * ((m_l * m_l + Q2) / (M_N * M_N)) * (G_A + Q2 * G_t_P / (4 * M_N * M_N)) * G_P;
    if (is_anti) B_LP = -B_LP;   
    double C_LP = 0;
    double D_LP = 0;
    double A2_LP = C0 * (B_LP + C_LP * s_u / pow(M_N, 2) + D_LP * pow(s_u, 2) / pow(M_N, 4)); 

    // double A2_LP = - D_0 * m_l * M_N * G_P * ((p_v * p_n) - (p_v * p_p)) * (2 * G_A - G_t_P * (1 - (p_n * p_p) / pow(M_N, 2)));    

    // double B_TT = - ((pow(m_l, 2) + Q2) / pow(M_N, 2)) * (pow(G_T, 2) 
    //             + (Q2 / (4 * pow(M_N, 2))) * (4 * G_T * G_T_1 - pow(m_l, 2) * (2 * G_T_2 * (G_T - 4 * G_T_2) + pow(G_T_1, 2) - 4 * G_T_1 * G_T_2) / pow(M_N, 2))
    //             + (pow(m_l, 2) / (2 * pow(M_N, 2))) * (pow(G_T, 2) - 4 * G_T * (G_T_1 + G_T_2) + 2 * pow(G_T_1 + 2 * G_T_2, 2))
    //             + 4 * pow(Q2 / (4 * pow(M_N, 2)), 2) * (pow(G_T_1, 2) + pow(m_l, 2) * pow(G_T_2, 2) / pow(M_N, 2)));
    // B = - ((pow(m_l, 2) + Q2) / pow(M_N, 2)) * (pow(G_T, 2) * (m_l * m_l / (2 * M_N * M_N) + 1) 
    //         + G_T_1 * G_T_1 * (Q2 * Q2 + m_l * m_l * (4 * M_N * M_N - Q2)) / (4 * pow(M_N, 4))
    //         + G_T_2 * G_T_2 * (4 * M_N * M_N + Q2) * (-2 * Q2 * Q2 * Q2 + pow(m_l, 4) * (4 * M_N * M_N + Q2) + m_l * m_l * Q2 * (4 * M_N * M_N + Q2)) / (4 * pow(M_N, 6) * (m_l * m_l + Q2))
    //         + G_T * G_T_1 * (Q2 - 2 * m_l * m_l) / (M_N * M_N)
    //         - G_T * G_T_2 * (-2 * Q2 * Q2 * Q2 + pow(m_l, 4) * (4 * M_N * M_N + Q2) + m_l * m_l * Q2 * (4 * M_N * M_N + Q2)) / (2 * pow(M_N, 4) * (m_l * m_l + Q2))
    //         + G_T_1 * G_T_2 * (-2 * Q2 * Q2 * Q2 + pow(m_l, 4) * (4 * M_N * M_N + Q2) + m_l * m_l * Q2 * (4 * M_N * M_N + Q2)) / (pow(M_N, 4) * (m_l * m_l + Q2)));
    const double ml2 = m_l * m_l;

    const double GT  = G_T;
    const double GT1 = G_T_1;
    const double GT2 = G_T_2;

    const double inner1 = GT*GT + a * ( 4.0*GT*GT1 - (ml2/MN2) * ( 2.0*GT2*(GT - 4.0*GT2) + GT1*GT1 - 4.0*GT1*GT2 ) );

    const double inner2 = (ml2 / (2.0 * MN2)) * ( GT*GT - 4.0*GT * ( (GT1 + GT2) + 2.0 * (GT1 + 2.0*GT2) * (GT1 + 2.0*GT2) ) );

    const double inner3 = 4.0 * a * a * ( GT1*GT1 + (ml2/MN2) * (GT2*GT2) );

    double B_TT = - ((ml2 + Q2) / MN2) * ( inner1 + inner2 + inner3 );
    double C_TT = 0;
    double D_TT = pow(G_T, 2) / 2 + (Q2 / (4 * pow(M_N, 2))) * (-2 * G_T * G_T_2 + pow(G_T_1 + 2 * G_T_2, 2))
                    + 4 * pow(G_T_2, 2) * pow(Q2 / (4 * pow(M_N, 2)), 2);          
    double A2_TT = C0 * (B_TT + C_TT * s_u / pow(M_N, 2) + D_TT * pow(s_u, 2) / pow(M_N, 4)); 
    
    // A2 = 4 * D_0 * ((((p_l * p_n) * (p_v * p_n)) + ((p_l * p_p) * (p_v * p_p))) * (2 * (G_T_1 + G_T_2) * (G_T - G_T_1 - G_T_2) + 
    //                 (p_n * p_p) * (pow(G_T_1, 2) - 2 * pow(G_T_2, 2)) / pow(M_N, 2))
    //             + (((p_l * p_n) * (p_v * p_p)) + ((p_l * p_p) * (p_v * p_n))) * (2 * (p_n * p_p) * G_T_2 
    //                 * (-G_T + 2 * G_T_1 + G_T_2 * (1 + (p_n * p_p) / pow(M_N, 2))) / pow(M_N, 2) + pow(G_T - G_T_1, 2)) 
    //             + pow(M_N, 2) * (p_l * p_v) * (pow(G_T, 2) * (p_n * p_p) / pow(M_N, 2) + 2 * (1 - (p_n * p_p) / pow(M_N, 2)) 
    //                 * (G_T_1 + G_T_2 * (1 + (p_n * p_p) / pow(M_N, 2))) * (G_T - G_T_1 - G_T_2 * (1 + (p_n * p_p) / pow(M_N, 2)))));
    // double A2_TT = 8 * D_0 * ((((p_l * p_n) * (p_v * p_n)) + ((p_l * p_p) * (p_v * p_p))) * (2 * (G_T_1 + G_T_2) * (G_T - G_T_1 - G_T_2) + 
    //                 (p_n * p_p) * (pow(G_T_1, 2) - 2 * pow(G_T_2, 2)) / pow(M_N, 2))
    //             + (((p_l * p_n) * (p_v * p_p)) + ((p_l * p_p) * (p_v * p_n))) * (2 * (p_n * p_p) * G_T_2 
    //                 * (-G_T + 2 * G_T_1 + G_T_2 * (1 + (p_n * p_p) / pow(M_N, 2))) / pow(M_N, 2) + pow(G_T - G_T_1, 2)) 
    //             - pow(M_N, 2) * (p_l * p_v) * (pow(G_T, 2) * (p_n * p_p) / pow(M_N, 2) + 2 * (1 - (p_n * p_p) / pow(M_N, 2)) 
    //                 * (G_T_1 + G_T_2 * (1 + (p_n * p_p) / pow(M_N, 2))) * (G_T - G_T_1 - G_T_2 * (1 + (p_n * p_p) / pow(M_N, 2)))));
       
    // double B_LT = - (m_l * m_l) * (m_l * m_l + Q2) / (M_N * M_N * M_N) * (F1 * (3 * G_T - 2 * G_T_1 - 4 * G_T_2) + 2 * F1 * G_T + (Q2 / (4 * M_N * M_N)) * 
    //         (4 * F1 * (G_T_1 - G_T_2) + F2 * (-G_T + 6 * G_T_1 + 4 * G_T_2)) + 4 * F2 * G_T_2 * pow(Q2 / (4 * M_N * M_N), 2));

    double B_LT = - ((m_l * m_l + Q2) / MN2) * (m_l / M_N) *
    (
        F1 * (3.0 * G_T - 2.0 * G_T_1 - 4.0 * G_T_2)
    + 2.0 * F2 * G_T
    + a * ( 4.0 * F1 * (G_T_1 - G_T_2) + F2 * (-G_T + 6.0 * G_T_1 + 4.0 * G_T_2) )
    + 4.0 * F2 * G_T_2 * (a * a)
    );

    double C_LT = -(m_l / M_N) * (3 * G_A * G_T + (G_t_P * G_T + 4 * G_A * G_T_1) * Q2 / (4 * M_N * M_N));
    if (is_anti) C_LT = -C_LT;   
    double D_LT = 0;
    double A2_LT = 0.5 * C0 * (B_LT + C_LT * s_u / pow(M_N, 2) + D_LT * pow(s_u, 2) / pow(M_N, 4)); 

    // B = - m_l * (m_l * m_l + Q2) / (M_N * M_N * M_N) * (F1 * (3 * G_T - 2 * G_T_1 - 4 * G_T_2) + 2 * F2 * G_T + (Q2 / (4 * M_N * M_N)) * 
    //         (4 * F1 * (G_T_1 - G_T_2) + F2 * (-G_T + 6 * G_T_1 + 4 * G_T_2)) + 4 * F2 * G_T_2 * pow(Q2 / (4 * M_N * M_N), 2));
    // C = -(m_l / M_N) * (3 * G_A * G_T + (G_t_P * G_T + 4 * G_A * G_T_1) * Q2 / (4 * M_N * M_N));
    // A2 = 0.5 * C0 * (B + C * s_u / pow(M_N, 2) + D * pow(s_u, 2) / pow(M_N, 4)); 
    
    // double A2_LT = - D_0 * m_l * M_N * (((p_v * p_n) - (p_v * p_p)) * (2 * F1 * (3 * G_T - 4 * G_T_1 - 2 * G_T_2)
    //             + F2 * (5 * G_T - 6 * G_T_1 - 2 * G_T_2)
    //             + (p_n * p_p) * (4 * F1 * (G_T_1 - G_T_2) - F2 * (G_T - 6 * G_T_1) + 2 * F2 * G_T_2 * (p_n * p_p) /  pow(M_N, 2)) / pow(M_N, 2))
    //             + ((p_v * p_n) + (p_v * p_p)) * (2 * G_A * (3 * G_T - 2 * G_T_1) - G_t_P * G_T + (4 * G_A * G_T_1 + G_T * G_t_P) * (p_n * p_p) / pow(M_N, 2)));     

    // double A2_LT = - D_0 * m_l * M_N * (((p_v * p_i) - (p_v * p_f)) * (2 * F1 * (3 * G_T - 4 * G_T_1 - 2 * G_T_2)
    //             + F2 * (5 * G_T - 6 * G_T_1 - 2 * G_T_2)
    //             + (p_i * p_f) * (4 * F1 * (G_T_1 - G_T_2) - F2 * (G_T - 6 * G_T_1) + 2 * F2 * G_T_2 * (p_i * p_f) /  pow(M_N, 2)) / pow(M_N, 2))
    //             + ((p_v * p_i) + (p_v * p_f)) * (2 * G_A * (3 * G_T - 2 * G_T_1) - G_t_P * G_T + (4 * G_A * G_T_1 + G_T * G_t_P) * (p_i * p_f) / pow(M_N, 2)));     

    double A2 = 0;
    if (par.weft_individual_or_combined == 0) {
        switch (par.weft_individual_case) {
            case 0: // LL
                A2 = A2_LL;
                break;
            
            case 1: // RR                 
                A2 = A2_RR;
                break;

            case 2: // LR
                A2 = A2_LR;
                break;

            case 3: // SS0
                A2 = A2_SS0;
                break;

            case 4: // LS0
                A2 = A2_LS0;
                break;

            case 5: // SS1
                A2 = A2_SS1;
                break;

            case 6: // LS1                      
                A2 = A2_LS1;
                break;

            case 7: // PP
                A2 = A2_PP;  
                break;

            case 8:  // LP 
                A2 = A2_LP; 
                break;

            case 9: // TT
                A2 = A2_TT;
                break;

            case 10: // LT
                A2 = A2_LT; 
                break;

            default:
                std::cout << "Please provide weft_individual_case!\n";
        }
    } else {
        double epi = par.weft_epi;

        switch (par.weft_combined_case) {
            case 0: // R
                A2 = A2_LL + epi * A2_LR + epi * epi * A2_RR;
                break;

            case 1: // S0
                A2 = A2_LL + epi * A2_LS0 + epi * epi * A2_SS0;
                break;

            case 2: // S1
                A2 = A2_LL + epi * A2_LS1 + epi * epi * A2_SS1;
                break;

            case 3: // P
                A2 = A2_LL + epi * A2_LP + epi * epi * A2_PP;
                break;
                
            case 4: // T
                A2 = A2_LL + epi * A2_LT + epi * epi * A2_TT;
                break;

            default:
                std::cout << "Please provide weft_combined_case!\n";
        }
    }

    double dsigma_dQ2 = (1 / (64 * M_PI)) * (1 / pow(p_v * p_i_eff, 2)) * A2;

    double result = dsigma_dQ2 * (std::min(Q2_b, Q2_max) - std::max(Q2_a, Q2_min));

    if (std::isnan(result) || std::isinf(result)) {
        std::cerr << "f(Q2, E_v) returned invalid value! Q2: " << Q2 << " E_v: " << E_v << std::endl;
        return 0;
    }
    
    double val = result / cm2;

    double q0_shift = 0.0; 

    if (par.sf_coulomb and e.flag.cc) q0_shift += coulomb_correction(is_anti, par.nucleus_p, par.nucleus_n);

    // modify lepton kinetic energy or xsec = 0 if not possible
    if (l1.Ek() > q0_shift)
        l1.set_energy(l1.E() - q0_shift);
    else
        return 0;

    if (N1.E() > N1.mass() - q0_shift)
        N1.set_energy(N1.E() + q0_shift);
    else
        return 0;

    e.weight = val;

    e.in[1] = N0;
    e.out.push_back(l1);
    e.out.push_back(N1);

    // selection of events for electron scattering using acceptance information
    if(l0.pdg==11)
    {
        double kosine=l1.z/l1.momentum();  
        if ( kosine < (par.el_costh_lab-par.el_costh_del) || kosine > (par.el_costh_lab+par.el_costh_del) )
        {
        e.weight=0;
        return 0;
        }
        else
        {
        e.weight /= 2*par.el_costh_del;
        val      /= 2*par.el_costh_del;
        }
        // KN: the output should be the differential in costh!
        //     (so the user doesn't have to remember what the width was)
    }

    return result;
}
