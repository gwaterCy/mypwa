#include "../phikk_structure.h"
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void get_para_kk(TString fname, TString tname, TString bname, TString var, double& sum_weight, double& dmin, double& dmax) {
    TFile *ff = new TFile(fname);
    TTree *tt = (TTree*)ff->Get(tname);
    DATA_PWA_KK dp;
    TBranch *bb;   //!
    tt->SetMakeClass(1);
    tt->SetBranchAddress(bname, &dp, &bb);

//    cout << "Set Adress" << endl;
    sum_weight = 0;
    dmin = 1000;
    dmax = -1000;

    Long64_t nentries = tt->GetEntriesFast();
    Long64_t nb = 0; 

    DATA_ORIG_KK ss;
//    cout << "Begin scan" << endl;
    cout << nentries << endl;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        nb = tt->GetEntry(jentry);
        pwa_to_orig(dp, ss); 
//        cout << "Done converter" << endl;
//        if (jentry % 1000 == 0) cout << jentry << endl;
//        if (value_kk("Mphi", ss) > 1.032) continue;
//        if (value_kk("Mphi", ss) < 1.006) continue;
        if (!good_event(ss)) continue; 
        double tmp = value_kk(var, ss);
        dmin = dmin < tmp ? dmin : tmp;
        dmax = dmax > tmp ? dmax : tmp;
        sum_weight += dp.weight;
    }
    cout << dmin << "   " << dmax << endl;
    cout << "Done parameters!" << endl;
} 

TH1F* hist_kk(TString fname, TString tname, TString bname, TString var, double dmin, double dmax) {
    TH1F* hh = new TH1F(var, TexName_KK(var), 100, dmin, dmax);
    TFile *ff = new TFile(fname);
    TTree *tt = (TTree*)ff->Get(tname);
    DATA_PWA_KK dp;
    TBranch *bb;   //!
    tt->SetMakeClass(1);
    tt->SetBranchAddress(bname, &dp, &bb);

    Long64_t nentries = tt->GetEntriesFast();
    Long64_t nb = 0; 

    DATA_ORIG_KK ss;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        nb = tt->GetEntry(jentry);
        pwa_to_orig(dp, ss); 
        if (!good_event(ss)) continue; 
//        if (value_kk("Mphi", ss) > 1.032) continue;
//        if (value_kk("Mphi", ss) < 1.006) continue;
        hh->Fill(value_kk(var, ss), dp.weight);
    }
    return hh;
}

void draw_kk(TString var) {
    double data_sum_weight, data_dmin, data_dmax;
    double phsp_sum_weight, phsp_dmin, phsp_dmax;
    get_para_kk("data_pwa_kk.root", "pwa_tr", "data_pwa_kk", var, data_sum_weight, data_dmin, data_dmax);
    get_para_kk("phsp_pwa_kk_weight.root", "pwa_tr", "kk_weight", var, phsp_sum_weight, phsp_dmin, phsp_dmax);
    double dmin = data_dmin < phsp_dmin ? data_dmin : phsp_dmin;
    double dmax = data_dmax > phsp_dmax ? data_dmax : phsp_dmax;
    TH1F* data_h = hist_kk("data_pwa_kk.root", "pwa_tr", "data_pwa_kk", var, dmin, dmax);
    TH1F* phsp_h = hist_kk("phsp_pwa_kk_weight.root", "pwa_tr", "kk_weight", var, dmin, dmax);
    phsp_h->Scale(abs(data_sum_weight) / phsp_sum_weight);
    TCanvas *c1 = new TCanvas();
    phsp_h->Draw("h");
    data_h->Draw("esame");
    c1->SaveAs("EPS/" + var + ".eps");
    c1->SaveAs("PNG/" + var + ".png");
}

void UpdataFigure() {
    draw_kk("Mphi");
    draw_kk("MKp1Km1");
    draw_kk("MKp2Km2");
    draw_kk("MKp1Km2");
    draw_kk("MKp2Km1");
    draw_kk("MphiKp2");
    draw_kk("MphiKm2");
    draw_kk("M2phiKp2");
    draw_kk("M2phiKm2");
    draw_kk("QKp1");
    draw_kk("QKm1");
    draw_kk("QKp2");
    draw_kk("QKm2");
    draw_kk("MomKp1");
    draw_kk("MomKm1");
    draw_kk("MomKp2");
    draw_kk("MomKm2");
    draw_kk("QKp1Km1");
    draw_kk("QKp2Km2");
}


