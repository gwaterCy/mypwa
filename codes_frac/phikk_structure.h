#ifndef PHIKK_STR_H
#define PHIKK_STR_H

#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <sstream>
#include "common_tools.h"

struct DATA_ANA_KK {
    Float_t Mphi;
    Float_t MKp1Km1;
    Float_t MKp2Km2;
    Float_t MKp1Km2;
    Float_t MKp2Km1;
    Float_t MphiKp2;
    Float_t MphiKm2;
    Float_t M2phiKp2;
    Float_t M2phiKm2;
    Float_t Chi2_4c;
    Float_t QKp1;
    Float_t QKm1;
    Float_t QKp2;
    Float_t QKm2;
    Float_t QKp1Km1;
    Float_t QKp2Km2;
    Float_t MomKp1;
    Float_t MomKm1;
    Float_t MomKp2;
    Float_t MomKm2;
    Float_t MKp1Km1Kp2;
    Float_t MKp1Km1Km2;
    Float_t MKp1Kp2Km2;
    Float_t MKm1Kp2Km2;
    Float_t weight;
    Int_t   itopo;
};

struct DATA_PWA_KK {
    Double_t Kp1X, Kp1Y, Kp1Z, Kp1E;
    Double_t Km1X, Km1Y, Km1Z, Km1E;
    Double_t Kp2X, Kp2Y, Kp2Z, Kp2E;
    Double_t Km2X, Km2Y, Km2Z, Km2E;
    Double_t weight;
};

struct DATA_ORIG_KK {
    TLorentzVector PKp1, PKm1, PKp2, PKm2, Pphi;
    Double_t Chi2_4c;
    Double_t weight;
    Int_t itopo;
};

extern TString MEM_ANA_KK;

extern TString MEM_PWA_KK;

double value_kk(TString var, const DATA_ORIG_KK& ss);
double value_kk(TString var, const DATA_ANA_KK& ss);

TString TexName_KK(TString var);

void orig_to_ana(const DATA_ORIG_KK& ss, DATA_ANA_KK& bb);

void orig_to_pwa(const DATA_ORIG_KK& ss, DATA_PWA_KK& bb);

void get_orig(DATA_ORIG_KK& ss, Double_t m4xyz[4][4], Double_t chi_4c, Int_t itopo);

void pwa_to_orig(const DATA_PWA_KK& bb, DATA_ORIG_KK& ss);

bool good_event(const DATA_ORIG_KK& ss);

bool good_event_sideband(const DATA_ORIG_KK& ss);

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
