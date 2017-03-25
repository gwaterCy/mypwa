#include "../phikk/phikk_structure.h"
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void get_para_kk(TString fname, TString tname, TString bname, TString var1, TString var2, double& sum_weight, double& dmin, double& dmax) {
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
        if (!good_event(ss)) continue;
//  if ((ss.PKp2 + ss.PKm2).M() < 1.2) continue;
//        cout << "Done converter" << endl;
//        if (jentry % 1000 == 0) cout << jentry << endl;
//        if (value_kk("Mphi", ss) > 1.032) continue;
//        if (value_kk("Mphi", ss) < 1.006) continue;
        double tmp = value_kk(var1, ss);
        dmin = dmin < tmp ? dmin : tmp;
        dmax = dmax > tmp ? dmax : tmp;
        tmp = value_kk(var2, ss);
        dmin = dmin < tmp ? dmin : tmp;
        dmax = dmax > tmp ? dmax : tmp;
        sum_weight += dp.weight;
    }
    cout << "Done parameters!" << endl;
} 

TH1F* hist_kk(TString fname, TString tname, TString bname, TString var1, TString var2, double dmin, double dmax) {
    TH1F* hh = new TH1F(var1 + "_" + var2, TexName_KK(var1) + " and " + TexName_KK(var2), 100, dmin, dmax);
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
        if (ss.Chi2_4c > 50) continue;
//  if ((ss.PKp2 + ss.PKm2).M() < 1.2) continue;
//        if (value_kk("Mphi", ss) > 1.032) continue;
//        if (value_kk("Mphi", ss) < 1.006) continue;
        hh->Fill(value_kk(var1, ss), dp.weight);
        hh->Fill(value_kk(var2, ss), dp.weight);
    }
    return hh;
}

void draw_kk(TString var1, TString var2, TString path) {
    double data_sum_weight, data_dmin, data_dmax;
    double phsp_sum_weight, phsp_dmin, phsp_dmax;
    double base_sum_weight, base_dmin, base_dmax;
    get_para_kk("../phikk/data_pwa_kk.root", "pwa_tr", "data_pwa_kk", var1, var2, data_sum_weight, data_dmin, data_dmax);
    get_para_kk("../" + path + "/phsp_pwa_kk_weight.root", "pwa_tr", "kk_weight", var1, var2, phsp_sum_weight, phsp_dmin, phsp_dmax);
    get_para_kk("../all01/phsp_pwa_kk_weight.root", "pwa_tr", "kk_weight", var1, var2, base_sum_weight, base_dmin, base_dmax);
//    data_dmin = phsp_dmin = base_dmin = 1.5;
//    data_dmax = phsp_dmax = base_dmax = 2.2;
    TH1F* data_h = hist_kk("../phikk/data_pwa_kk.root", "pwa_tr", "data_pwa_kk", var1, var2, data_dmin, data_dmax);
    TH1F* phsp_h = hist_kk("../" + path + "/phsp_pwa_kk_weight.root", "pwa_tr", "kk_weight", var1, var2, phsp_dmin, phsp_dmax);
    TH1F* base_h = hist_kk("../all01/phsp_pwa_kk_weight.root", "pwa_tr", "kk_weight", var1, var2, base_dmin, base_dmax);
    phsp_h->Scale(abs(data_sum_weight) / phsp_sum_weight);
    base_h->Scale(abs(data_sum_weight) / base_sum_weight);
    TCanvas *c1 = new TCanvas();
    double Ymax = phsp_h->GetMaximum() > data_h->GetMaximum() ? phsp_h->GetMaximum() : data_h->GetMaximum();
    phsp_h->GetYaxis()->SetRangeUser(0, 1.1 * Ymax);
    data_h->GetYaxis()->SetRangeUser(0, 1.1 * Ymax);
    base_h->GetYaxis()->SetRangeUser(0, 1.1 * Ymax);
    phsp_h->SetLineColor(8);
    phsp_h->SetLineWidth(2.5);
    phsp_h->Draw("h");
    base_h->SetLineColor(2);
    base_h->Draw("hsame");
    data_h->Draw("esame");
    c1->SaveAs("../" + path + "/kk_EPS/" + var1 + "_" + var2 + ".eps");
    c1->SaveAs("../" + path + "/kk_PNG/" + var1 + "_" + var2 + ".png");
}

//void UpdataFigure(TString path = ".") {
//    draw_kk("Mphi", path);
//    draw_kk("MKp1Km1", path);
//    draw_kk("MKp2Km2", path);
//    draw_kk("MKp1Km2", path);
//    draw_kk("MKp2Km1", path);
//    draw_kk("MphiKp2", path);
//    draw_kk("MphiKm2", path);
//    draw_kk("M2phiKp2", path);
//    draw_kk("M2phiKm2", path);
//    draw_kk("QKp1", path);
//    draw_kk("QKm1", path);
//    draw_kk("QKp2", path);
//    draw_kk("QKm2", path);
//    draw_kk("MomKp1", path);
//    draw_kk("MomKm1", path);
//    draw_kk("MomKp2", path);
//    draw_kk("MomKm2", path);
//    draw_kk("QKp1Km1", path);
//    draw_kk("QKp2Km2", path);
//}
void test(TString var, TString mycut, Int_t bin_num = 100, Float_t dmin = -100, Float_t dmax = -100)
{
    TFile* f1 = TFile::Open("../phikk/data_ana_kk_signal.root");
    TTree* t1 = (TTree*)f1->Get("signal_tr");
    TFile* f2 = TFile::Open("../phikk/data_ana_kk_sideband.root");
    TTree* t2 = (TTree*)f2->Get("sideband_tr");
    if (dmax < -10) {
        dmax = (t1->GetMaximum(var) > t2->GetMaximum(var))? t1->GetMaximum(var) : t2->GetMaximum(var);
    }
    if (dmin < -10) {
        dmin = (t1->GetMinimum(var) < t2->GetMinimum(var))? t1->GetMinimum(var) : t2->GetMinimum(var);
    }
    if (mycut == "") {
        mycut = var + ">" + Float2Str(dmin) + "&" + var + "<" + Float2Str(dmax);
    } else {
        mycut += "&" + var + ">" + Float2Str(dmin) + "&" + var + "<" + Float2Str(dmax);
    }
    cout << mycut << endl;
    TH1F *h1 = new TH1F("h1", TexName_KK(var), bin_num, dmin, dmax);
    h1->GetXaxis()->SetTitle(TexName_KK(var));
    h1->GetYaxis()->SetTitle("Entries");
    TH1F *h2 = new TH1F("h2", TexName_KK(var), bin_num, dmin, dmax);
    h2->GetXaxis()->SetTitle(TexName_KK(var));
    h2->GetYaxis()->SetTitle("Entries");
    t1->Project("h1", var, mycut);
    t2->Project("h2", var, mycut);
    cout << h1->GetEntries() << endl;
    cout << h2->GetEntries() << endl;
    Float_t scale_factor = 1.0 * h1->GetEntries() / h2->GetEntries();
    h2->Scale(scale_factor);
    h2->SetLineColor(kRed);
    h2->SetLineWidth(1.5);
    h1->SetLineWidth(2);
    h1->Sumw2();
    TCanvas *c1 = new TCanvas();
    if (h1->GetMaximum() > h2->GetMaximum()) {
        h1->Draw("e");
        h2->Draw("esame");
    } else {
        h2->Draw("e");
        h1->Draw("esame");
    }
    c1->SaveAs("kk_signal_sideband_PNG/" + var + ".png");
}


void UpdataFigure_sideband() {
    test("Mphi", "");
    test("MKp1Km1", "");
    test("MKp2Km2", "");
    test("MKp1Km2", "");
    test("MKp2Km1", "");
    test("MphiKp2", "");
    test("MphiKm2", "");
    test("M2phiKp2", "");
    test("M2phiKm2", "");
    test("QKp1", "");
    test("QKm1", "");
    test("QKp2", "");
    test("QKm2", "");
    test("MomKp1", "");
    test("MomKm1", "");
    test("MomKp2", "");
    test("MomKm2", "");
    test("QKp1Km1", "");
    test("QKp2Km2", "");
}
