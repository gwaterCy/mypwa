#ifndef PHIKK_STR_CC
#define PHIKK_STR_CC

#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <sstream>
#include "phikk_structure.h"

TString MEM_ANA_KK = "Mphi/F:MKp1Km1:MKp2Km2:MKp1Km2:MKp2Km1:MphiKp2:MphiKm2:M2phiKp2:M2phiKm2:Chi2_4c:QKp1:QKm1:QKp2:QKm2:QKp1Km1:QKp2Km2:MomKp1:MomKm1:MomKp2:MomKm2:MKp1Km1Kp2:MKp1Km1Km2:MKp1Kp2Km2:MKm1Kp2Km2:weight:itopo/I";

TString MEM_PWA_KK = "Kp1X/D:Kp1Y:Kp1Z:Kp1E:Km1X:Km1Y:Km1Z:Km1E:Kp2X:Kp2Y:Kp2Z:Kp2E:Km2X:Km2Y:Km2Z:Km2E:weight";

double value_kk(TString var, const DATA_ANA_KK& ss) {
    if (var == "Mphi") {
        return ss.Mphi;
    } else if (var == "MKp1Km1") {
        return ss.MKp1Km1;
    } else if (var == "MKp2Km2") {
        return ss.MKp2Km2;
    } else if (var == "MKp1Km2") {
        return ss.MKp1Km2;
    } else if (var == "MKp2Km1") {
        return ss.MKp2Km1;
    } else if (var == "MphiKp2") {
        return ss.MphiKp2;
    } else if (var == "MphiKm2") {
        return ss.MphiKm2;
    } else if (var == "M2phiKp2") {
        return ss.M2phiKp2;
    } else if (var == "M2phiKm2") {
        return ss.M2phiKm2;
    } else if (var == "QKp1") {
        return ss.QKp1;
    } else if (var == "QKm1") {
        return ss.QKm1;
    } else if (var == "QKp2") {
        return ss.QKp2;
    } else if (var == "QKm2") {
        return ss.QKm2;
    } else if (var == "MomKp1") {
        return ss.MomKp1;
    } else if (var == "MomKm1") {
        return ss.MomKm1;
    } else if (var == "MomKp2") {
        return ss.MomKp2;
    } else if (var == "MomKm2") {
        return ss.MomKm2;
    } else if (var == "QKp1Km1") {
        return ss.QKp1Km1;
    } else if (var == "QKp2Km2") {
        return ss.QKp2Km2;
    } else if (var == "MKp1Km1Kp2") {
        return ss.MKp1Km1Kp2;
//        if ((ss.PKm2).E() > (ss.PKm1).E()) {
//            return (ss.PKp1 + ss.PKm1 + ss.PKp2).M();
//        } else {
//            return (ss.PKp1 + ss.PKp2 + ss.PKm2).M();
//        }
    } else if (var == "MKp1Km1Km2") {
        return ss.MKp1Km1Km2;
//        if ((ss.PKp2).E() > (ss.PKp1).E()) {
//            return (ss.PKp1 + ss.PKm1 + ss.PKm2).M();
//        } else {
//            return (ss.PKm1 + ss.PKp2 + ss.PKm2).M();
//        }
    } else if (var == "MKp1Kp2Km2") {
        return ss.MKp1Kp2Km2;
    } else if (var == "MKm1Kp2Km2") {
        return ss.MKm1Kp2Km2;
    } else {
        return -1000;
    }
}

double value_kk(TString var, const DATA_ORIG_KK& ss) {
    if (var == "Mphi") {
        return (ss.PKp1 + ss.PKm1).M();
    } else if (var == "MKp1Km1") {
        return (ss.PKp1 + ss.PKm1).M();
    } else if (var == "MKp2Km2") {
        return (ss.PKp2 + ss.PKm2).M();
    } else if (var == "MKp1Km2") {
        return (ss.PKp1 + ss.PKm2).M();
    } else if (var == "MKp2Km1") {
        return (ss.PKp2 + ss.PKm1).M();
    } else if (var == "MphiKp2") {
        return (ss.PKp1 + ss.PKm1 + ss.PKp2).M();
    } else if (var == "MphiKm2") {
        return (ss.PKp1 + ss.PKm1 + ss.PKm2).M();
    } else if (var == "M2phiKp2") {
        return (ss.PKp1 + ss.PKm1 + ss.PKp2).M2();
    } else if (var == "M2phiKm2") {
        return (ss.PKp1 + ss.PKm1 + ss.PKm2).M2();
    } else if (var == "QKp1") {
        return ss.PKp1.CosTheta();
    } else if (var == "QKm1") {
        return ss.PKm1.CosTheta();
    } else if (var == "QKp2") {
        return ss.PKp2.CosTheta();
    } else if (var == "QKm2") {
        return ss.PKm2.CosTheta();
    } else if (var == "MomKp1") {
        return ss.PKp1.P();
    } else if (var == "MomKm1") {
        return ss.PKm1.P();
    } else if (var == "MomKp2") {
        return ss.PKp2.P();
    } else if (var == "MomKm2") {
        return ss.PKm2.P();
    } else if (var == "QKp1Km1") {
        return (ss.PKp1 + ss.PKm1).CosTheta();
    } else if (var == "QKp2Km2") {
        return (ss.PKp2 + ss.PKm2).CosTheta();
    } else if (var == "QphiKp2") {
        return (ss.PKp1 + ss.PKm1 + ss.PKp2).CosTheta();
    } else if (var == "QphiKm2") {
        return (ss.PKp1 + ss.PKm1 + ss.PKm2).CosTheta();
    } else if (var == "MKp1Km1Kp2") {
        if ((ss.PKm2).E() > (ss.PKm1).E()) {
            return (ss.PKp1 + ss.PKm1 + ss.PKp2).M();
        } else {
            return (ss.PKp1 + ss.PKp2 + ss.PKm2).M();
        }
    } else if (var == "MKp1Km1Km2") {
        if ((ss.PKp2).E() > (ss.PKp1).E()) {
            return (ss.PKp1 + ss.PKm1 + ss.PKm2).M();
        } else {
            return (ss.PKm1 + ss.PKp2 + ss.PKm2).M();
        }
    } else if (var == "MKp1Kp2Km2") {
        return (ss.PKp1 + ss.PKp2 + ss.PKm2).M();
    } else if (var == "MKm1Kp2Km2") {
        return (ss.PKm1 + ss.PKp2 + ss.PKm2).M();
    } else {
        return -1000;
    }
}

TString TexName_KK(TString var) {
    if (var == "Mphi") {
        return "M(#phi)";
    } else if (var == "MKp1Km1") {
        return "M(K_{#phi}^{+}K_{#phi}^{-})";
    } else if (var == "MKp2Km2") {
        return "M(K^{+}K^{-})";
    } else if (var == "MKp1Km2") {
        return "M(K_{#phi}^{+}K^{-})";
    } else if (var == "MKp2Km1") {
        return "M(K_{#phi}^{-}K^{+})";
    } else if (var == "MphiKp2") {
        return "M(#phi K^{+})";
    } else if (var == "MphiKm2") {
        return "M(#phi K^{-})";
    } else if (var == "M2phiKp2") {
        return "M^{2}(#phi K^{+})";
    } else if (var == "M2phiKm2") {
        return "M^{2}(#phi K^{-})";
    } else if (var == "Chi2_4c") {
        return "#chi^{2}_{4c}";
    } else if (var == "QKp1") {
        return "#theta(K_{#phi}^{+})";
    } else if (var == "QKm1") {
        return "#theta(K_{#phi}^{-})";
    } else if (var == "QKp2") {
        return "#theta(K^{+})";
    } else if (var == "QKm2") {
        return "#theta(K^{-})";
    } else if (var == "MomKp1") {
        return "p(K_{#phi}^{+})";
    } else if (var == "MomKm1") {
        return "p(K_{#phi}^{-})";
    } else if (var == "MomKp2") {
        return "p(K^{+})";
    } else if (var == "MomKm2") {
        return "p(K^{-})";
    } else if (var == "QKp1Km1") {
        return "#theta(K_{#phi}^{+}K_{#phi}^{-})";
    } else if (var == "QKp2Km2") {
        return "#theta(K^{+}K^{-})";
    } else if (var == "QphiKp2") {
        return "#theta(#phi K^{+})";
    } else if (var == "QphiKm2") {
        return "#theta(#phi K^{-})";
    } else if (var == "MKp1Km1Kp2") {
        return "M(K_{#phi}^{+}K_{#phi}^{-}K^{+})";
    } else if (var == "MKp1Km1Km2") {
        return "M(K_{#phi}^{+}K_{#phi}^{-}K^{-})";
    } else if (var == "MKp1Kp2Km2") {
        return "M(K_#phi^{+}K^{+}K^{-})";
    } else if (var == "MKm1Kp2Km2") {
        return "M(K_#phi^{-}K^{+}K^{-})";
    } else {
        return "UNKOWN";
    }
}

void orig_to_ana(const DATA_ORIG_KK& ss, DATA_ANA_KK& bb) {
    bb.Mphi = (ss.Pphi).M();
    bb.MKp1Km1 = (ss.PKp1 + ss.PKm1).M();
    bb.MKp2Km2 = (ss.PKp2 + ss.PKm2).M();
    bb.MKp1Km2 = (ss.PKp1 + ss.PKm2).M();
    bb.MKp2Km1 = (ss.PKp2 + ss.PKm1).M();
    bb.MphiKp2 = (ss.Pphi + ss.PKp2).M();
    bb.MphiKm2 = (ss.Pphi + ss.PKm2).M();
    bb.M2phiKp2 = (ss.Pphi + ss.PKp2).M2();
    bb.M2phiKm2 = (ss.Pphi + ss.PKm2).M2();
    bb.Chi2_4c = ss.Chi2_4c;
    bb.QKp2 = (ss.PKp2).CosTheta();
    bb.QKm2 = (ss.PKm2).CosTheta();
    bb.MomKp2 = (ss.PKp2).P();
    bb.MomKm2 = (ss.PKm2).P();
    bb.QKp1Km1 = (ss.PKp1 + ss.PKm1).CosTheta();
    bb.QKp2Km2 = (ss.PKp2 + ss.PKm2).CosTheta();
    bb.MKp1Km1Kp2 = (ss.PKp1 + ss.PKm1 + ss.PKp2).M();
    bb.MKp1Km1Km2 = (ss.PKp1 + ss.PKm1 + ss.PKm2).M();
    bb.MKp1Kp2Km2 = (ss.PKp1 + ss.PKp2 + ss.PKm2).M();
    bb.MKm1Kp2Km2 = (ss.PKm1 + ss.PKp2 + ss.PKm2).M();
    bb.itopo = ss.itopo;
    TVector3 b = (ss.Pphi).BoostVector();
    TLorentzVector bPKp1 = ss.PKp1;
    TLorentzVector bPKm1 = ss.PKm1;
    bPKp1.Boost(-b);
    bPKm1.Boost(-b);
    bb.QKp1 = bPKp1.CosTheta();
    bb.QKm1 = bPKm1.CosTheta();
    bb.MomKp1 = bPKp1.P();
    bb.MomKm1 = bPKm1.P();
    bb.weight = ss.weight;
}

void orig_to_pwa(const DATA_ORIG_KK& ss, DATA_PWA_KK& bb) {
    bb.Kp1X = ss.PKp1.X(); bb.Kp1Y = ss.PKp1.Y(); bb.Kp1Z = ss.PKp1.Z(); bb.Kp1E = ss.PKp1.E();
    bb.Km1X = ss.PKm1.X(); bb.Km1Y = ss.PKm1.Y(); bb.Km1Z = ss.PKm1.Z(); bb.Km1E = ss.PKm1.E();
    bb.Kp2X = ss.PKp2.X(); bb.Kp2Y = ss.PKp2.Y(); bb.Kp2Z = ss.PKp2.Z(); bb.Kp2E = ss.PKp2.E();
    bb.Km2X = ss.PKm2.X(); bb.Km2Y = ss.PKm2.Y(); bb.Km2Z = ss.PKm2.Z(); bb.Km2E = ss.PKm2.E();
    bb.weight = ss.weight;
}

void get_orig(DATA_ORIG_KK& ss, Double_t m4xyz[4][4], Double_t chi_4c, Int_t itopo) {
    V2LV(ss.PKp1, m4xyz[0]);
    V2LV(ss.PKm1, m4xyz[1]);
    V2LV(ss.PKp2, m4xyz[2]);
    V2LV(ss.PKm2, m4xyz[3]);
    ss.Pphi = ss.PKp1 + ss.PKm1;
    ss.Chi2_4c = chi_4c;
    ss.itopo = itopo;
    ss.weight = 1;
}

void pwa_to_orig(const DATA_PWA_KK& bb, DATA_ORIG_KK& ss) {
    ss.PKp1 = TLorentzVector(bb.Kp1X, bb.Kp1Y, bb.Kp1Z, bb.Kp1E);
    ss.PKm1 = TLorentzVector(bb.Km1X, bb.Km1Y, bb.Km1Z, bb.Km1E);
    ss.PKp2 = TLorentzVector(bb.Kp2X, bb.Kp2Y, bb.Kp2Z, bb.Kp2E);
    ss.PKm2 = TLorentzVector(bb.Km2X, bb.Km2Y, bb.Km2Z, bb.Km2E);
    ss.weight = bb.weight;
    ss.Pphi = ss.PKp1 + ss.PKm1;
}

bool good_event(const DATA_ORIG_KK& ss) {
  if (OutRange((ss.PKp1 + ss.PKm1).M(), 1.006, 1.032)) return false;
//  if ((ss.PKp2 + ss.PKm2).M() < 1.2) return false;
    return true;
}

bool good_event_sideband(const DATA_ORIG_KK& ss) {
    if (OutRange((ss.PKp1 + ss.PKm1).M(), 1.04, 1.06)) return false;
//  if ((ss.PKp2 + ss.PKm2).M() < 1.2) return false;
    return true;
}

//void copy_ana(const DATA_ANA& rr, DATA_ANA& aa) {
//    aa.Mphi == rr.Mphi;
//    aa.MKp1Km1 == rr.MKp1Km1;
//    aa.MKp2Km2 == rr.MKp2Km2;
//    aa.MKp1Km2 == rr.MKp1Km2;
//    aa.MKp2Km1 == rr.MKp2Km1;
//    aa.MphiKp2 == rr.MphiKp2;
//    aa.MphiKm2 == rr.MphiKm2;
//    aa.M2phiKp2 == rr.M2phiKp2;
//    aa.M2phiKm2 == rr.M2phiKm2;
//    aa.Chi2_4c == rr.Chi2_4c;
//    aa.QKp1 = rr.QKp1;
//    aa.QKm1 = rr.QKm1;
//    aa.QKp2 = rr.QKp2;
//    aa.QKm2 = rr.QKm2;
//    aa.QKp1Km1 = rr.QKp1Km1;
//    aa.QKp2Km2 = rr.QKp2Km2;
//    aa.MomKp1 = rr.MomKp1;
//    aa.MomKm1 = rr.MomKm1;
//    aa.MomKp2 = rr.MomKp2;
//    aa.MomKm2 = rr.MomKm2;
//    aa.itopo = rr.itopo;
//}


#endif
