#ifndef PHIPIPI_STR_H
#define PHIPIPI_STR_H

#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <sstream>
#include "common_tools.h"

struct DATA_ANA_PIPI {
    Float_t Mphi;
    Float_t MKppim;
    Float_t MKmpip;
    Float_t M2phipip;
    Float_t M2phipim;
    Float_t Mpippim;
    Float_t Chi2_4c;
    Float_t MKpKm;
    Float_t Mphipip;
    Float_t Mphipim;
    Float_t QKp;
    Float_t QKm;
    Float_t Qpip;
    Float_t Qpim;
    Float_t QKpKm;
    Float_t Qpippim;
    Float_t weight;
    Int_t   itopo;
};

struct DATA_PWA_PIPI {
    Double_t pipX, pipY, pipZ, pipE;
    Double_t pimX, pimY, pimZ, pimE;
    Double_t KpX, KpY, KpZ, KpE;
    Double_t KmX, KmY, KmZ, KmE;
    Double_t weight;
};

struct DATA_ORIG_PIPI {
    TLorentzVector Ppip, Ppim, PKp, PKm, Pphi;
    Double_t Chi2_4c;
    Double_t weight;
    Int_t itopo;
};

extern TString MEM_ANA_PIPI;

extern TString MEM_PWA_PIPI;

double value_pipi(TString var, const DATA_ORIG_PIPI& ss);
double value_pipi(TString var, const DATA_ANA_PIPI& ss);

TString TexName_PIPI(TString var);

bool good_event(const DATA_ORIG_PIPI& ss);

bool good_event_sideband(const DATA_ORIG_PIPI& ss);

void orig_to_ana(const DATA_ORIG_PIPI& ss, DATA_ANA_PIPI& bb);

void pwa_to_orig(const DATA_PWA_PIPI& bb, DATA_ORIG_PIPI& ss);

void orig_to_pwa(const DATA_ORIG_PIPI& ss, DATA_PWA_PIPI& bb);

void get_orig(DATA_ORIG_PIPI& ss, Double_t m4xyz[4][4], Double_t chi_4c, Int_t itopo);

#endif






