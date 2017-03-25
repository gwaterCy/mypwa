#include "TROOT.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "DPFPWAPdf.h"
#include "DPFPWAPoint.h"
//#include "RooStringVar.h"
//#include "RooAbsPdf.h"
//#include "RooAbsData.h"
//#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPrintable.h"
#include "TTree.h"
//#include "RooRealProxy.h"
//#include "RooListProxy.h"
//#include "RooComplex.h"
//#include "RooFitResult.h"
//#include "RooArgSet.h"
//#include "RooRealConstant.h"

#include "fitproxy.h"
#include "../codes_frac/phikk_structure.h"
#include "../codes_frac/phipipi_structure.h"
#include <omp.h>
#include <iomanip>

#if defined(TEST)
const bool m_test = true;
#else
const bool m_test = false;
#endif

#if defined(FIXPHI)
const bool m_fitphi = true;
#else
const bool m_fitphi = false;
#endif

using namespace RooFit;
using namespace std;

const Double_t mpsip=3.686,mka=0.493677,mpi=0.13957;
const Double_t pi=3.14159265358979312;
const Double_t pi2=2 * pi;
const Double_t high=mpsip;
const Double_t low=0-high;
const Double_t Nreal=50000;
//double high=mjpsi/2;
void Info(RooRealVar *bb) {
    double ll = bb->getMax() - bb->getMin();
    if (bb->getAttribute("Constant")) {
    cout << "*" << bb->GetName() << " = " << setw(8) << bb->getVal() << " In(" << setw(8) << bb->getMin() << ", " << setw(8) << bb->getMax() << ")"<< " L" << setw(6) << (bb->getVal() - bb->getMin()) / ll << " R" << setw(6) << (bb->getMax() - bb->getVal()) / ll << endl;
    } else {
    cout << " " << bb->GetName() << " = " << setw(8) << bb->getVal() << " In(" << setw(8) << bb->getMin() << ", " << setw(8) << bb->getMax() << ")"<< " L" << setw(6) << (bb->getVal() - bb->getMin()) / ll << " R" << setw(6) << (bb->getMax() - bb->getVal()) / ll << endl;
    }
}
void fitproxy::createlist_allparas() {
    for(int i = 0; i < vpars.size(); i++) {
        for(int j = 0; j < vpars[i].size(); j++) {
            allparas.add(vpars[i][j]);
        }
    }
}
void fitproxy::add_fitparas(TString rn) {
    for(int i = 0; i < vn.size(); i++) {
        if (rn == vn[i]) {
            for(int j = 0; j < vpars[i].size(); j++) {
                fitparas.add(vpars[i][j]);
            }
            break;
        }
    }
}
void fitproxy::set_all_constant(TString rn) {
    for(int i = 0; i < vn.size(); i++) {
        if (rn == vn[i]) {
            for(int j = 0; j < vpars[i].size(); j++) {
                vpars[i][j].setConstant();
            }
            break;
        }
    }
}
void fitproxy::set_all_constant() {
    for(int i = 0; i < vn.size(); i++) {
        for(int j = 0; j < vpars[i].size(); j++) {
            vpars[i][j].setConstant();
        }
    }
}
void fitproxy::add_res980_list(TString rn, double mss_min, double mss_max) {
    vector<RooRealVar> tmp;
    if (m_fitphi) {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 2));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2, mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_g10_", rn + "a_g10_", 0.165, 0.1, 1.0));
        tmp.push_back(RooRealVar(rn + "a_g20_", rn + "a_g20_", 0.695, 0.1, 1.5));
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 23.6, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 23.6, -100, 100));
        // f0(980) 0+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 1));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "a_phi1", rn + "p_phi1", 0, -pi2, pi2));
        // f0(980) 0+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 2));
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi2", rn + "a_phi2", 0, -pi2, pi2));
    } else {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 2));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2, mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_g10_", rn + "a_g10_", 0.165, 0.1, 1.0));
        tmp.push_back(RooRealVar(rn + "a_g20_", rn + "a_g20_", 0.695, 0.1, 1.5));
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 23.6, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 23.6, -100, 100));
        // f0(980) 0+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 1));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "p_phi1", rn + "p_phi1", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi1", rn + "k_phi1", 0, -pi2, pi2));
        // f0(980) 0+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 2));
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi2", rn + "p_phi2", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi2", rn + "k_phi2", 0, -pi2, pi2));
    }
    vpars.push_back(tmp);
    vn.push_back(rn);
}
void fitproxy::act_res980(TString rn, bool act_pp, bool act_kk) {
    if (m_fitphi) {
        if (act_pp) {
            pdfphipp->addResonance980(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance980(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance980(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance980(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
        }

    } else {
        if (act_pp) {
            pdfphipp->addResonance980(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "p_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance980(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "p_phi2"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance980(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "k_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance980(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_g10_"), *gp(rn + "a_g20_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "k_phi2"), *gp(rn + "a_typ_"));
        }
    }
    add_fitparas(rn);
}
void fitproxy::add_res0_list(TString rn, double mss_min, double mss_max, double wdt_min, double wdt_max) {
    vector<RooRealVar> tmp;
    if (m_fitphi) {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 1));                                         // 0
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2, mss_min, mss_max)); // 1
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2, wdt_min, wdt_max)); // 2
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 23.6, -100, 100));                           // 3
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 23.6, -100, 100));                           // 4
        // f0 0+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 1));                                         // 5
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));                                         // 6
        tmp.push_back(RooRealVar(rn + "a_phi1", rn + "a_phi1", 0, -pi2, pi2));                                // 7
        // f0 0+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 2));                                         // 9
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));                              // 10
        tmp.push_back(RooRealVar(rn + "a_phi2", rn + "a_phi2", 0, -pi2, pi2));                                // 11
    } else {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 1));                                         // 0
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2, mss_min, mss_max)); // 1
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2, wdt_min, wdt_max)); // 2
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 23.6, -100, 100));                           // 3
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 23.6, -100, 100));                           // 4
        // f0 0+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 1));                                         // 5
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));                                         // 6
        tmp.push_back(RooRealVar(rn + "p_phi1", rn + "p_phi1", 0, -pi2, pi2));                                // 7
        tmp.push_back(RooRealVar(rn + "k_phi1", rn + "k_phi1", 0, -pi2, pi2));                              // 8
        // f0 0+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 2));                                         // 9
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));                              // 10
        tmp.push_back(RooRealVar(rn + "p_phi2", rn + "p_phi2", 0, -pi2, pi2));                                // 11
        tmp.push_back(RooRealVar(rn + "k_phi2", rn + "k_phi2", 0, -pi2, pi2));                              // 12
    }
    vpars.push_back(tmp);
    vn.push_back(rn);
}
RooRealVar* fitproxy::gp(TString tn) {
    cout << allparas.find(tn)->GetName() << endl;
    return (RooRealVar*)allparas.find(tn);
}
void fitproxy::act_res0(TString rn, bool act_pp, bool act_kk) {
    if (m_fitphi) {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
        }

    } else {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "p_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "p_phi2"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "k_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "k_phi2"), *gp(rn + "a_typ_"));
        }
    }
    add_fitparas(rn);
}
void fitproxy::act_res2(TString rn, bool act_pp, bool act_kk) {
    if (m_fitphi) {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_3", rn + "p_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc3"), *gp(rn + "a_phi3"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_4", rn + "p_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc4"), *gp(rn + "a_phi4"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_5", rn + "p_5", *gp(rn + "a_spn5"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc5"), *gp(rn + "a_phi5"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_3", rn + "k_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc3"), *gp(rn + "a_phi3"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_4", rn + "k_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc4"), *gp(rn + "a_phi4"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_5", rn + "k_5", *gp(rn + "a_spn5"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc5"), *gp(rn + "a_phi5"), *gp(rn + "a_typ_"));
        }

    } else {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "p_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "p_phi2"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_3", rn + "p_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc3"), *gp(rn + "p_phi3"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_4", rn + "p_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc4"), *gp(rn + "p_phi4"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_5", rn + "p_5", *gp(rn + "a_spn5"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc5"), *gp(rn + "p_phi5"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "k_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "k_phi2"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_3", rn + "k_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc3"), *gp(rn + "k_phi3"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_4", rn + "k_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc4"), *gp(rn + "k_phi4"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_5", rn + "k_5", *gp(rn + "a_spn5"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc5"), *gp(rn + "k_phi5"), *gp(rn + "a_typ_"));
        }
    }
    add_fitparas(rn);
}
void fitproxy::add_res2_list(TString rn, double mss_min, double mss_max, double wdt_min, double wdt_max) {
    vector<RooRealVar> tmp;
    if (m_fitphi) {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 6));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2,  mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2,  wdt_min, wdt_max));
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 23.6, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 23.6, -100, 100));
        // f2 2+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 21));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "a_phi1", rn + "a_phi1", 0, -pi2, pi2));
        // f2 2+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 22));
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi2", rn + "a_phi2", 0, -pi2, pi2));
        // f2 2+ 3
        tmp.push_back(RooRealVar(rn + "a_spn3", rn + "a_spn3", 23));
        tmp.push_back(RooRealVar(rn + "a_frc3", rn + "a_frc3", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi3", rn + "a_phi3", 0, -pi2, pi2));
        // f2 2+ 4
        tmp.push_back(RooRealVar(rn + "a_spn4", rn + "a_spn4", 24));
        tmp.push_back(RooRealVar(rn + "a_frc4", rn + "a_frc4", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi4", rn + "a_phi4", 0, -pi2, pi2));
        // f2 2+ 5
        tmp.push_back(RooRealVar(rn + "a_spn5", rn + "a_spn5", 25));
        tmp.push_back(RooRealVar(rn + "a_frc5", rn + "a_frc5", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi5", rn + "a_phi5", 0, -pi2, pi2));
    } else {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 6));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2,  mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2,  wdt_min, wdt_max));
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 23.6, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 23.6, -100, 100));
        // f2 2+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 21));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "p_phi1", rn + "p_phi1", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi1", rn + "k_phi1", 0, -pi2, pi2));
        // f2 2+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 22));
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi2", rn + "p_phi2", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi2", rn + "k_phi2", 0, -pi2, pi2));
        // f2 2+ 3
        tmp.push_back(RooRealVar(rn + "a_spn3", rn + "a_spn3", 23));
        tmp.push_back(RooRealVar(rn + "a_frc3", rn + "a_frc3", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi3", rn + "p_phi3", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi3", rn + "k_phi3", 0, -pi2, pi2));
        // f2 2+ 4
        tmp.push_back(RooRealVar(rn + "a_spn4", rn + "a_spn4", 24));
        tmp.push_back(RooRealVar(rn + "a_frc4", rn + "a_frc4", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi4", rn + "p_phi4", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi4", rn + "k_phi4", 0, -pi2, pi2));
        // f2 2+ 5
        tmp.push_back(RooRealVar(rn + "a_spn5", rn + "a_spn5", 25));
        tmp.push_back(RooRealVar(rn + "a_frc5", rn + "a_frc5", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi5", rn + "p_phi5", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi5", rn + "k_phi5", 0, -pi2, pi2));
    }
    vpars.push_back(tmp);
    vn.push_back(rn);
}


void fitproxy::add_res1m_list(TString rn, double mss_min, double mss_max, double wdt_min, double wdt_max) {
    vector<RooRealVar> tmp;
    if (m_fitphi) {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 4));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2, mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2, wdt_min, wdt_max));
        tmp.push_back(RooRealVar(rn + "a_rho_", rn + "a_rho_", 0, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 0, -100, 100));
        // z 1- 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 111));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "a_phi1", rn + "a_phi1", 0, -pi2, pi2));
    } else {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 4));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2, mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2, wdt_min, wdt_max));
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 0, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 0, -100, 100));
        // z 1- 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 111));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "p_phi1", rn + "p_phi1", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi1", rn + "k_phi1", 0, -pi2, pi2));
    }
    vpars.push_back(tmp);
    vn.push_back(rn);
}

void fitproxy::act_res1m(TString rn, bool act_pp, bool act_kk) {
    if (m_fitphi) {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
        }
    } else {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "p_phi1"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "k_phi1"), *gp(rn + "a_typ_"));
        }
    }
    add_fitparas(rn);
}

void fitproxy::add_res1p_list(TString rn, double mss_min, double mss_max, double wdt_min, double wdt_max) {
    vector<RooRealVar> tmp;
    if (m_fitphi) {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 4));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2,  mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2,  wdt_min, wdt_max));
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 0, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 0, -100, 100));
        // z 1+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 11));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "a_phi1", rn + "a_phi1", 0, -pi2, pi2));
        // z 1+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 12));
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi2", rn + "a_phi2", 0, -pi2, pi2));
        // z 1+ 3
        tmp.push_back(RooRealVar(rn + "a_spn3", rn + "a_spn3", 13));
        tmp.push_back(RooRealVar(rn + "a_frc3", rn + "a_frc3", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi3", rn + "a_phi3", 0, -pi2, pi2));
        // z 1+ 4
        tmp.push_back(RooRealVar(rn + "a_spn4", rn + "a_spn4", 14));
        tmp.push_back(RooRealVar(rn + "a_frc4", rn + "a_frc4", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "a_phi4", rn + "a_phi4", 0, -pi2, pi2));
    } else {
        tmp.push_back(RooRealVar(rn + "a_typ_", rn + "a_typ_", 4));
        tmp.push_back(RooRealVar(rn + "a_mss_", rn + "a_mss_", (mss_min + mss_max) / 2,  mss_min, mss_max));
        tmp.push_back(RooRealVar(rn + "a_wdt_", rn + "a_wdt_", (wdt_min + wdt_max) / 2,  wdt_min, wdt_max));
        tmp.push_back(RooRealVar(rn + "p_rho_", rn + "p_rho_", 0, -100, 100));
        tmp.push_back(RooRealVar(rn + "k_rho_", rn + "k_rho_", 0, -100, 100));
        // z 1+ 1
        tmp.push_back(RooRealVar(rn + "a_spn1", rn + "a_spn1", 11));
        tmp.push_back(RooRealVar(rn + "a_frc1", rn + "a_frc1", 0));
        tmp.push_back(RooRealVar(rn + "p_phi1", rn + "p_phi1", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi1", rn + "k_phi1", 0, -pi2, pi2));
        // z 1+ 2
        tmp.push_back(RooRealVar(rn + "a_spn2", rn + "a_spn2", 12));
        tmp.push_back(RooRealVar(rn + "a_frc2", rn + "a_frc2", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi2", rn + "p_phi2", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi2", rn + "k_phi2", 0, -pi2, pi2));
        // z 1+ 3
        tmp.push_back(RooRealVar(rn + "a_spn3", rn + "a_spn3", 13));
        tmp.push_back(RooRealVar(rn + "a_frc3", rn + "a_frc3", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi3", rn + "p_phi3", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi3", rn + "k_phi3", 0, -pi2, pi2));
        // z 1+ 4
        tmp.push_back(RooRealVar(rn + "a_spn4", rn + "a_spn4", 14));
        tmp.push_back(RooRealVar(rn + "a_frc4", rn + "a_frc4", 0, -4.6, 4.6));
        tmp.push_back(RooRealVar(rn + "p_phi4", rn + "p_phi4", 0, -pi2, pi2));
        tmp.push_back(RooRealVar(rn + "k_phi4", rn + "k_phi4", 0, -pi2, pi2));
    }
    vpars.push_back(tmp);
    vn.push_back(rn);
}

void fitproxy::act_res1p(TString rn, bool act_pp, bool act_kk) {
    if (m_fitphi) {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_3", rn + "p_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc3"), *gp(rn + "a_phi3"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_4", rn + "p_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc4"), *gp(rn + "a_phi4"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "a_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "a_phi2"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_3", rn + "k_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc3"), *gp(rn + "a_phi3"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_4", rn + "k_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc4"), *gp(rn + "a_phi4"), *gp(rn + "a_typ_"));
        }
    } else {
        if (act_pp) {
            pdfphipp->addResonance(rn + "p_1", rn + "p_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc1"), *gp(rn + "p_phi1"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_2", rn + "p_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc2"), *gp(rn + "p_phi2"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_3", rn + "p_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc3"), *gp(rn + "p_phi3"), *gp(rn + "a_typ_"));
            pdfphipp->addResonance(rn + "p_4", rn + "p_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "p_rho_"), *gp(rn + "a_frc4"), *gp(rn + "p_phi4"), *gp(rn + "a_typ_"));
        }
        if (act_kk) {
            pdfphikk->addResonance(rn + "k_1", rn + "k_1", *gp(rn + "a_spn1"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc1"), *gp(rn + "k_phi1"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_2", rn + "k_2", *gp(rn + "a_spn2"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc2"), *gp(rn + "k_phi2"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_3", rn + "k_3", *gp(rn + "a_spn3"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc3"), *gp(rn + "k_phi3"), *gp(rn + "a_typ_"));
            pdfphikk->addResonance(rn + "k_4", rn + "k_4", *gp(rn + "a_spn4"), *gp(rn + "a_mss_"), *gp(rn + "a_wdt_"), *gp(rn + "k_rho_"), *gp(rn + "a_frc4"), *gp(rn + "k_phi4"), *gp(rn + "a_typ_"));
        }
    }
    add_fitparas(rn);
}

void fitproxy::init_all_files_name(const PWA_CTRL & pwa_ctrl) {
    outf_phipp = pwa_ctrl.outf_phipp;
    outf_phikk = pwa_ctrl.outf_phikk;
    proj_phipp = pwa_ctrl.proj_phipp;
    proj_phikk = pwa_ctrl.proj_phikk;
//    if(m_test){
//      phsp_phipp = "../phipipi/phsp_pwa_pipi_small_sample.dat";
//      phsp_phikk = "../phikk/phsp_pwa_kk_small_sample.dat";
//      data_phipp = "../phipipi/data_pwa_pipi_weight_small_sample.dat";
//      idx_pp = "../phipipi/idp_pp.dat";
//      data_phikk = "../phikk/data_pwa_kk_weight_small_sample.dat";
//      idx_kk = "../phikk/idp_kk.dat";
//      cout << "idx_pp = " << idx_pp << endl;
//    }
//    else{
        phsp_phipp = pwa_ctrl.phsp_phipp;
        phsp_phikk = pwa_ctrl.phsp_phikk;
        data_phipp = pwa_ctrl.data_phipp;
        data_phikk = pwa_ctrl.data_phikk;
        idx_pp = pwa_ctrl.idx_pp;
        idx_kk = pwa_ctrl.idx_kk;
      //phsp_phipp = "../phipipi/phsp_pwa_pipi.dat";
      //phsp_phikk = "../phikk/phsp_pwa_kk.dat";
      //data_phipp = "../phipipi/data_pwa_pipi_weight.dat";
      //data_phikk = "../phikk/data_pwa_kk_weight.dat";
      //idx_pp = "../phipipi/idp_pp.dat";
      //idx_kk = "../phikk/idp_kk.dat";
//    }
}
void fitproxy::init_input_argset() {
    vv.push_back(RooRealVar("idp", "idp", 10));
    vv.push_back(RooRealVar("weight", "weight", low, high));
//    for (int i = 1; i <= 5; i++) {
//        for (int j = 1; j <= 4; j++) {
//            TString vname = "v" + Int2Str(i) + Int2Str(j);
//            TString vdesp = Int2Str(i) + Int2Str(j);
//            cout << vname << "   " << vdesp << endl;
//            vv.push_back(RooRealVar(vname, vdesp, low, high));
//        }
//    }

    //for (vector<RooRealVar>::iterator iter = vv.begin(); iter != vv.end(); iter++) {
    //    theSet.add(*iter);
    //}

    //idp = RooRealVar("idp", "index", 1);
    theSet.add(vv[0]);
    theSet.add(vv[1]);
    cout << "haha: " << __LINE__ << endl;
    cout << vv[0].getVal() << endl;

//    theSet.writeToFile("para_v.txt");
//    theSet.readFromFile("para_v.txt");
//    theSet.writeToStream(cout, false);

}
void fitproxy::read_data() {
    cout<<"-----------------"<<endl;
    //RooDataSet *data11 = RooDataSet::read(data_phipp,theSet);
    RooDataSet *data11 = RooDataSet::read(idx_pp,theSet);
    data11->Print();
    cout << "haha: " << __LINE__ << endl;
    datapp = new RooDataSet(data11->GetName(),data11->GetTitle(),data11,*data11->get(),0,(theSet.find("weight"))->GetName());
    cout << "haha: " << __LINE__ << endl;
    datapp->Print();

    cout<<"-----------------"<<endl;
    //RooDataSet *data22 = RooDataSet::read(data_phikk,theSet);
    RooDataSet *data22 = RooDataSet::read(idx_kk,theSet);
    data22->Print();
    datakk = new RooDataSet(data22->GetName(),data22->GetTitle(),data22,*data22->get(),0,(theSet.find("weight"))->GetName());

    datakk->Print();
    cout<<"RooDataSet is initialized!!"<<endl;
}
void fitproxy::init_pdf(const PWA_CTRL & pwa_ctrl) {
    dphipp = new DPFPWAPoint(mka, mka, mpi, mpi, mpsip, phsp_phipp, data_phipp);
    dphikk = new DPFPWAPoint(mka, mka, mka, mka, mpsip, phsp_phikk, data_phikk);
    cout<<"DPFPWAPoint is intialized!!"<<endl;
//    pdfphipp = new DPFPWAPdf("pdfphipp", "pdfphipp", vv[0], vv[1], vv[2], vv[3], vv[4], vv[5], vv[6], vv[7], vv[8], vv[9], vv[10], vv[11], vv[12], vv[13], vv[14], vv[15], vv[16], vv[17], vv[18], vv[19], dphipp);
//    pdfphikk = new DPFPWAPdf("pdfphikk", "pdfphikk", vv[0], vv[1], vv[2], vv[3], vv[4], vv[5], vv[6], vv[7], vv[8], vv[9], vv[10], vv[11], vv[12], vv[13], vv[14], vv[15], vv[16], vv[17], vv[18], vv[19], dphikk);
    //cout << "vv[0]=" << vv[0] << endl;
    pdfphipp = new DPFPWAPdf("pdfphipp", "pdfphipp", vv[0], dphipp);
    pdfphikk = new DPFPWAPdf("pdfphikk", "pdfphikk", vv[0], dphikk);
    pdfphipp->work_path = pwa_ctrl.workPath;
    pdfphikk->work_path = pwa_ctrl.workPath;
    cout<<"DPFPWAPdf is intialized!!"<<endl;
}
void fitproxy::setup_resonances() {
    add_res0_list("f01000", 1.2, 1.5, 0, 20000);
    add_res2_list("f21000", 1.2, 1.5, 0, 20000);

    // baseline 0++
    add_res980_list("f00980", 0.93, 1.05);
    add_res0_list("f01370", 1.2, 1.5, 0.2, 0.5);                    // weak, pipi seen, kkbar seen
    add_res0_list("f01500", 1.486, 1.522, 0.088, 0.130);            // definite, p+p-(pi0pi0) 34.9%, kkbar 8.6%
    add_res0_list("f01710", 1.705, 1.741, 0.115, 0.163);            // weak, pipi seen, kkbar seen

    // PDG left 0++
    add_res0_list("f02020", 1.944, 2.04, 0.262, 0.622);             // weak, pi0pi0 seen
    add_res0_list("f02100", 2.079, 2.18, 0.155, 0.63);              // weak
    add_res0_list("f02200", 2.15, 2.228, 0.088, 0.388);             // weak
    add_res0_list("f02330", 2.239, 2.389, 0.084, 0.316);            // weak
    add_res0_list("f02201", 1.944, 4.0, 0.084, 0.63);             // homemade f02020, f02100 f02200 f02330

    // baseline 2++
    //add_res2_list("f21270", 1.2731, 1.2779, 0.1779, 0.1933);        // definite, pp 84.2%, kkbar 4.6%
    add_res2_list("f21270", 1.1731, 1.2779, 0.1779, 0.1933);        // definite, pp 84.2%, kkbar 4.6%
    add_res2_list("f21525", 1.51, 1.54, 0.058, 0.091);              // definite, pp 0.82%, kkbar 88.7%

    // PDG left 2++
    //add_res2_list("f21430", 1.33, 1.63, 0, 0.3);                    // weak
    //add_res2_list("f21565", 1.523, 1.629, 0.11, 0.29);              // weak, pp seen
    //add_res2_list("f21640", 1.572, 1.668, 0, 0.32);                 // weak, kkbar seen
    add_res2_list("f21501", 1.33, 1.668, 0, 0.32);                  // homemade f21430, f21565, f21640

    add_res2_list("f21810", 1.71, 1.962, 0.131, 0.34);              // weak
    add_res2_list("f21910", 1.86, 1.982, 0.021, 0.346);             // weak, k+k- seen
    add_res2_list("f21950", 1.79, 2.13, 0.418, 0.7);                // weak, p+p- seen, kkbar seen
    add_res2_list("f22150", 2.121, 2.239, 0.062, 0.42);             // weak, kkbar seen
    add_res2_list("f22300", 2.213, 2.381, 0.029, 0.289);            // weak, kkbar seen
    add_res2_list("f22201", 2.0, 4.0, 0.021, 0.418);             // homemade f21810, f21910, f21950, f22150, f22300


    // PDG left 1--
    add_res1m_list("1m1800", 1.5, 2.0, 0, 0.3);

    // PDG left 1+-
    add_res1p_list("1p1800", 1.5, 2.0, 0, 0.3);

}
void fitproxy::print_all_paras() {
    for(int i = 0; i < vpars.size(); i++) {
        for(int j = 0; j < vpars[i].size(); j++) {
            Info(&vpars[i][j]);
        }
    }
}
void fitproxy::print_fit_paras() {
    for(int i = 0; i < vpars.size(); i++) {
        for(int j = 0; j < vpars[i].size(); j++) {
            //cout <<  fitparas.find(vpars[i][j].GetName()) << endl;
            if (fitparas.find(vpars[i][j].GetName()) != NULL) {
                Info(&vpars[i][j]);
            }
        }
    }
}
void fitproxy::store_all_paras(TString fname) {
    TFile *ff2 = new TFile(fname, "RECREATE");
    TIterator *it = allparas.createIterator();
    for(int i = 0; i < allparas.getSize(); i++) {
        RooRealVar *bb = (RooRealVar*)it->Next();
        RooRealVar aa = *bb;
        //bb->Write(bb->GetName());
        aa.Write();
    }
    ff2->Close();
    cout << "all parameters are stored in " << fname << endl;
}
void fitproxy::store_fit_paras(TString fname) {
    TFile *ff2 = new TFile(fname, "RECREATE");
    for(int i = 0; i < vpars.size(); i++) {
        for(int j = 0; j < vpars[i].size(); j++) {
            //cout <<  fitparas.find(vpars[i][j].GetName()) << endl;
            if (fitparas.find(vpars[i][j].GetName()) != NULL) {
                RooRealVar aa = vpars[i][j];
                cout << "WRITE vpars" << endl;
                aa.Write();
                //Info(&vpars[i][j]);
            }
        }
    }
    //allparas.Write("all_paras");
    //fitparas.Write("fit_paras");
    //TIterator *it = fitparas.createIterator();
    //cout << "fitparas size is " << fitparas.getSize() << endl;
    //for(int i = 0; i < fitparas.getSize(); i++) {
    //    RooRealVar *bb = (RooRealVar*)it->Next();
    //    cout << "test roopdf" << " i = " << i << endl;
    //    //bb->Write(bb->GetName());
    //    cout << bb->GetName() << endl;
    //    cout << "bb value is " << bb->getValV() << endl;
    //    bb->Write();
    //    cout << "end roopdf" << endl;
    //}
    ff2->Close();
    cout << "fit parameters are stored in " << fname << endl;
}
void fitproxy::reload_paras(TString fname) {
    TFile *ff3 = new TFile(fname);
//    allparas.Read("all_paras");
//    fitparas.Read("fit_paras");
//    TIterator *it = allparas.createIterator();
//    for(int i = 0; i < allparas.getSize(); i++) {
    TIterator *it = fitparas.createIterator();
    cout <<"@@@@@@@ fitparas size = " << fitparas.getSize() << endl;
    for(int i = 0; i < fitparas.getSize(); i++) {
        RooRealVar *bb = (RooRealVar*)it->Next();
        //RooRealVar bb_;
        //bb_.Read(bb->GetName());
        //cout << "bb_ name is " << bb_.GetName() << endl;
        //bb->Copy(bb_);
        bb->Read(bb->GetName());
        cout << "+++++++++++++++++++++++++++" << endl;
        cout << bb->GetName() << "the point is " << bb->getValV() << endl;
    }
    ff3->Close();
    cout << "these parameters are reloaded from " << fname << endl;
}
void fitproxy::FIT() {
    double startTime_Total=omp_get_wtime();
//    RooRealVar *bb = (RooRealVar*)allparas.find("f01000a_mss_");
//    cout << bb->GetName() << "---------->" << bb->getVal() << endl;
    RooCategory kapi("kapi","kapi");
    kapi.defineType("phipp");
    kapi.defineType("phikk");

    RooDataSet combdata("combdata","combined data",theSet, Index(kapi),Import("phipp",*datapp),Import("phikk",*datakk));

    RooSimultaneous simPdf("simPdf","simultaneous pdf",kapi);
    simPdf.addPdf(*pdfphipp,"phipp");
    simPdf.addPdf(*pdfphikk,"phikk");

    cout << "Begin fitTo" << endl;
    RooFitResult* res = simPdf.fitTo(combdata,Save());
    //   RooFitResult* res = pdfphikk->fitTo(*datakk,Save());

    double endTime_Total = omp_get_wtime();
    cout<<"Fit time : "<< endTime_Total - startTime_Total <<"\n";
    res->Print("V");

    cout << "The minNll of res is " << res->minNll() << endl;
    cout << "The number of floating parameters is " << (res->floatParsFinal()).getSize() << endl;

//    cout <<((RooRealVar*)fitparas.find("f00980a_mss_"))->getVal() << endl;
////    ((RooRealVar*)fitparas.find("f00980a_mss_"))->setVal(2);
//    pdfphikk->change_value("f00980a_mss_", 0.4);
////    pdfphikk->change_value("f00980k_rho_", 10);
////    gp("f00980a_mss_")->setVal(1.2);
//    cout <<((RooRealVar*)fitparas.find("f00980a_mss_"))->getVal() << endl;
//    cout << "The minNll of res is " << res->minNll() << endl;
//    cout << "The number of floating parameters is " << (res->floatParsFinal()).getSize() << endl;
}
fitproxy::fitproxy () {
//    init_all_files_name();
//    init_input_argset();
//    setup_resonances();
//    createlist_allparas();
//    read_data();
//    init_pdf();
//
////    TFile *ff = new TFile("fit_results.root", "RECREATE");
////    allparas.Write("init_all_paras");
////    ff->Close();
//    print_all_paras();
//
//    act_res980("f00980");
////    act_res2("f21270");
////    FIT();
//
//    print_fit_paras();
//    gp("f00980a_mss_")->setVal(1.05);
//    cout << "End Write" << endl;
//    gp("f00980a_mss_")->setVal(0.93);
//    print_fit_paras();
//
//    cout << "Store fit results!!!" << endl;
////    pdfphipp->writeToFile("pipi_paras_results.root");
////    pdfphikk->writeToFile("kk_paras_results.root");
//
//    cout << "Fit is done!!!!" << endl;
//    delete datapp;
//    delete datakk;
}



