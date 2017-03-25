#ifndef PHIPIPI_STR_CC
#define PHIPIPI_STR_CC

#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <sstream>
#include "phipipi_structure.h"

TString MEM_ANA_PIPI = "Mphi/F:MKppim:MKmpip:M2phipip:M2phipim:Mpippim:Chi2_4c:MKpKm:Mphipip:Mphipim:QKp:QKm:Qpip:Qpim:QKpKm:Qpippim:weight:itopo/I";

TString MEM_PWA_PIPI = "pipX/D:pipY:pipZ:pipE:pimX:pimY:pimZ:pimE:KpX:KpY:KpZ:KpE:KmX:KmY:KmZ:KmE:weight";

double value_pipi(TString var, const DATA_ANA_PIPI& ss) {
    if (var == "Mphi") {
        return ss.Mphi;
    } else if (var == "MKppim") {
        return ss.MKppim;
    } else if (var == "MKmpip") {
        return ss.MKmpip;
    } else if (var == "M2phipip") {
        return ss.M2phipip;
    } else if (var == "M2phipim") {
        return ss.M2phipim;
    } else if (var == "Mpippim") {
        return ss.Mpippim;
    } else if (var == "MKpKm") {
        return ss.MKpKm;
    } else if (var == "Mphipip") {
        return ss.Mphipip;
	} else if (var == "Mphipim") {
        return ss.Mphipim;
    } else if (var == "QKp") {
        return ss.QKp;
    } else if (var == "QKm") {
        return ss.QKm;
    } else if (var == "Qpip") {
        return ss.Qpip;
    } else if (var == "Qpim") {
        return ss.Qpim;
    } else if (var == "QKpKm") {
        return ss.QKpKm;
    } else if (var == "Qpippim") {
        return ss.Qpippim;
    } else {
        return -10000;
    }

}
double value_pipi(TString var, const DATA_ORIG_PIPI& ss) {
    if (var == "Mphi") {
        return (ss.PKp + ss.PKm).M();
    } else if (var == "MKppim") {
        return (ss.PKp + ss.Ppim).M();
    } else if (var == "MKmpip") {
        return (ss.PKm + ss.Ppip).M();
    } else if (var == "M2phipip") {
        return (ss.PKp + ss.PKm + ss.Ppip).M2();
    } else if (var == "M2phipim") {
        return (ss.PKp + ss.PKm + ss.Ppim).M2();
    } else if (var == "Mpippim") {
        return (ss.Ppip + ss.Ppim).M();
    } else if (var == "MKpKm") {
        return (ss.PKp + ss.PKm).M();
    } else if (var == "Mphipip") {
        return (ss.PKp + ss.PKm + ss.Ppip).M();
	} else if (var == "Mphipim") {
		return (ss.PKp + ss.PKm + ss.Ppim).M();
    } else if (var == "QKp") {
        return ss.PKp.CosTheta();
    } else if (var == "QKm") {
        return ss.PKm.CosTheta();
    } else if (var == "Qpip") {
        return ss.Ppip.CosTheta();
    } else if (var == "Qpim") {
        return ss.Ppim.CosTheta();
    } else if (var == "QKpKm") {
        return (ss.PKp + ss.PKm).CosTheta();
    } else if (var == "Qpippim") {
        return (ss.Ppip + ss.Ppim).CosTheta();
    } else {
        return -10000;
    }

}

TString TexName_PIPI(TString var) {
    if (var == "Mphi") {
        return "M(#phi)";
    } else if (var == "MKppim") {
        return "M(K^{+}#pi^{-})";
    } else if (var == "MKmpip") {
        return "M(K^{-}#pi^{+})";
    } else if (var == "M2phipip") {
        return "M^{2}(#phi#pi^{+})";
    } else if (var == "M2phipim") {
        return "M^{2}(#phi#pi^{-})";
    } else if (var == "Mpippim") {
        return "M(#pi^{+}#pi^{-})";
    } else if (var == "Chi2_4c") {
        return "#chi^{2}_{4c}";
    } else if (var == "MKpKm") {
        return "M(K^{+}K^{-})";
    } else if (var == "Mphipip") {
        return "M(#phi#pi^{+})";
	} else if (var == "Mphipim") {
		return "M(#phi#pi^{-})";
    } else if (var == "QKp") {
        return "#theta(K^{+})";
    } else if (var == "QKm") {
        return "#theta(K^{-})";
    } else if (var == "Qpip") {
        return "#theta(#pi^{+})";
    } else if (var == "Qpim") {
        return "#theta(#pi^{-})";
    } else if (var == "QKpKm") {
        return "#theta(K^{+}K^{-})";
    } else if (var == "Qpippim") {
        return "#theta(#pi^{+}#pi^{-})";
    } else {
        return "UNKOWN";
    }
}


bool good_event(const DATA_ORIG_PIPI& ss) {
    if (OutRange((ss.PKp + ss.PKm).M(), 1.006, 1.032)) return false;
    if (InRange((ss.PKp + ss.Ppim).M(), 0.85, 0.95)) return false;
    if (InRange((ss.PKm + ss.Ppip).M(), 0.85, 0.95)) return false;
//  if ((ss.Ppip + ss.Ppim).M() < 1.2) return false;
    return true;
}
bool good_event_sideband(const DATA_ORIG_PIPI& ss) {
    if (OutRange((ss.PKp + ss.PKm).M(), 1.04, 1.06)) return false;
    if (InRange((ss.PKp + ss.Ppim).M(), 0.85, 0.95)) return false;
    if (InRange((ss.PKm + ss.Ppip).M(), 0.85, 0.95)) return false;
//  if ((ss.Ppip + ss.Ppim).M() < 1.2) return false;
    return true;
}

void orig_to_ana(const DATA_ORIG_PIPI& ss, DATA_ANA_PIPI& bb) {
    bb.Mphi = (ss.Pphi).M();
    bb.MKppim = (ss.PKp + ss.Ppim).M();
    bb.MKmpip = (ss.PKm + ss.Ppip).M();
    bb.M2phipip = (ss.Pphi + ss.Ppip).M2();
    bb.M2phipim = (ss.Pphi + ss.Ppim).M2();
    bb.Mpippim = (ss.Ppip + ss.Ppim).M();
    bb.Chi2_4c = ss.Chi2_4c;
    bb.MKpKm = (ss.PKp + ss.PKm).M();
    bb.Mphipip = (ss.Pphi + ss.Ppip).M();
    bb.Mphipim = (ss.Pphi + ss.Ppim).M();
    bb.QKp = (ss.PKp).CosTheta();
    bb.QKm = (ss.PKm).CosTheta();
    bb.Qpip = (ss.Ppip).CosTheta();
    bb.Qpim = (ss.Ppim).CosTheta();
    bb.QKpKm = (ss.PKp + ss.PKm).CosTheta();
    bb.Qpippim = (ss.Ppip + ss.Ppim).CosTheta();
    bb.weight = ss.weight;
    bb.itopo = ss.itopo;
}
void pwa_to_orig(const DATA_PWA_PIPI& bb, DATA_ORIG_PIPI& ss) {
    ss.Ppip = TLorentzVector(bb.pipX, bb.pipY, bb.pipZ, bb.pipE);
    ss.Ppim = TLorentzVector(bb.pimX, bb.pimY, bb.pimZ, bb.pimE);
    ss.PKp = TLorentzVector(bb.KpX, bb.KpY, bb.KpZ, bb.KpE);
    ss.PKm = TLorentzVector(bb.KmX, bb.KmY, bb.KmZ, bb.KmE);
    ss.Pphi = ss.PKp + ss.PKm;
    ss.weight = bb.weight;
}

void orig_to_pwa(const DATA_ORIG_PIPI& ss, DATA_PWA_PIPI& bb) {
    bb.pipX = ss.Ppip.X(); bb.pipY = ss.Ppip.Y(); bb.pipZ = ss.Ppip.Z(); bb.pipE = ss.Ppip.E();
    bb.pimX = ss.Ppim.X(); bb.pimY = ss.Ppim.Y(); bb.pimZ = ss.Ppim.Z(); bb.pimE = ss.Ppim.E();
    bb.KpX = ss.PKp.X(); bb.KpY = ss.PKp.Y(); bb.KpZ = ss.PKp.Z(); bb.KpE = ss.PKp.E();
    bb.KmX = ss.PKm.X(); bb.KmY = ss.PKm.Y(); bb.KmZ = ss.PKm.Z(); bb.KmE = ss.PKm.E();
    bb.weight = ss.weight;
}


void get_orig(DATA_ORIG_PIPI& ss, Double_t m4xyz[4][4], Double_t chi_4c, Int_t itopo) {
    V2LV(ss.Ppip, m4xyz[0]);
    V2LV(ss.Ppim, m4xyz[1]);
    V2LV(ss.PKp, m4xyz[2]);
    V2LV(ss.PKm, m4xyz[3]);
    ss.Pphi = ss.PKp + ss.PKm;
    ss.Chi2_4c = chi_4c;
    ss.itopo = itopo;
    ss.weight = 1;
}


#endif
