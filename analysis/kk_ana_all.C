#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <iostream>
#include <algorithm>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>

#include "../codes_frac/phikk_structure.cc"
#include "../codes_frac/PWA_CTRL.C"
#include "../codes_frac/common_tools.cc"

#include <iomanip>

using namespace std;

void get_para_kk(TString fname, TString tname, TString bname, TString var, double& sum_weight, double& dmin, double& dmax) {
    cout << "fname = " << fname << endl;
    cout << "var = " << var << endl;
    TFile *ff = new TFile(fname);
    TTree *tt = (TTree*)ff->Get(tname);
    DATA_ANA_KK dp;
    TBranch *bb;   //!
    tt->SetMakeClass(1);
    tt->SetBranchAddress(bname, &dp, &bb);

//    cout << "Set Adress" << endl;
    sum_weight = 0;
    dmin = 1000;
    dmax = -1000;

    Long64_t nentries = tt->GetEntriesFast();
    Long64_t nb = 0;

//    cout << "Begin scan" << endl;
    cout << nentries << endl;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        nb = tt->GetEntry(jentry);
//        if (!good_event(ss)) continue;
//        cout << "Done converter" << endl;
//        if (jentry % 1000 == 0) cout << jentry << endl;
        //if (value_kk("Mphi", dp) > 1.032) continue;
        //if (value_kk("Mphi", dp) < 1.006) continue;
        double tmp = value_kk(var, dp);
        dmin = dmin < tmp ? dmin : tmp;
        dmax = dmax > tmp ? dmax : tmp;
        sum_weight += dp.weight;
    }
    cout << "Done parameters!" << endl;
}

TH1F* hist_kk(TString fname, TString tname, TString bname, TString resName, TString var, double dmin, double dmax) {
    TH1F* hh = new TH1F(var, TexName_KK(var) + "_" +resName, 100, dmin, dmax);
    TFile *ff = new TFile(fname);
    TTree *tt = (TTree*)ff->Get(tname);
    DATA_ANA_KK dp;
    TBranch *bb;   //!
    tt->SetMakeClass(1);
    tt->SetBranchAddress(bname, &dp, &bb);

    Long64_t nentries = tt->GetEntriesFast();
    Long64_t nb = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        nb = tt->GetEntry(jentry);
        //if (!good_event(ss)) continue;
        //if (value_kk("Mphi", dp) > 1.032) continue;
        //if (value_kk("Mphi", dp) < 1.006) continue;
        hh->Fill(value_kk(var, dp), dp.weight);
    }
    return hh;
    ff->Close();
}

void draw_kk(TString var, string cfgFile) {
    vector<double> weight, dmin, dmax;
    double _weight, _dmin, _dmax;
    string path;
    readConfigFile(cfgFile, "work_path", path);
    cout << "path = " << path << endl;
    //TString path = "/public/users/caihao1/tmp/baseline_small/";
    get_para_kk(path + "data_ana_kk_signal.root", "signal_tr", "signal_ana_kk", var, _weight, _dmin, _dmax);
    weight.push_back(_weight);
    dmin.push_back(_dmin);
    dmax.push_back(_dmax);
    get_para_kk(path + "phsp_pwa_kk_weight_all.root", "ana_tr", "kk", var, _weight, _dmin, _dmax);
    weight.push_back(_weight);
    dmin.push_back(_dmin);
    dmax.push_back(_dmax);
    string value, value1;
    readConfigFile(cfgFile, "active_resonances", value);
    readConfigFile(cfgFile, "active_resonances_kk", value1);
    if (value1 != "NONE") {
        value = value + "," + value1;
    }
    vector<string> resNameList;
    string_to_vector(value, resNameList);
    for(vector<string>::iterator it = resNameList.begin(); it != resNameList.end(); it++) {
        get_para_kk(path + "phsp_pwa_kk_weight_" + *it + ".root", "ana_tr", "kk", var, _weight, _dmin, _dmax);
        weight.push_back(_weight);
        dmin.push_back(_dmin);
        dmax.push_back(_dmax);
        TString ss = *it;
        cout << ss << "------" << _weight << endl;
    }
    _dmin = *min_element(dmin.begin(), dmin.end());
    _dmax = *max_element(dmax.begin(), dmax.end());

    THStack *hs = new THStack("hs_" + var, "");
    double scaleFactor = weight[0] / weight[1];
    int tcolor = 1;
    TH1F *signal_h = hist_kk(path + "data_ana_kk_signal.root", "signal_tr", "signal_ana_kk", "signal", var, _dmin, _dmax);
    signal_h->SetLineColor(tcolor++);
    signal_h->SetLineWidth(3);
    signal_h->GetXaxis()->SetTitle(TexName_KK(var));

    double sp = valid_digits((_dmax - _dmin) / 100, 2);
    signal_h->GetYaxis()->SetTitle("Entries / " + Double2Str(sp));

    TH1F *all_h = hist_kk(path + "phsp_pwa_kk_weight_all.root", "ana_tr", "kk", "all", var, _dmin, _dmax);
    all_h->SetLineColor(tcolor++);
    all_h->SetLineWidth(2);
    all_h->GetXaxis()->SetTitle(TexName_KK(var));
    //double scaleFactor = 1.0 * signal_h->GetEntries() / all_h->GetEntries();
    std::cout << "scaleFactor = " << scaleFactor << endl;
    all_h->Scale(scaleFactor);
    TCanvas *tc = new TCanvas();
    gStyle->SetOptStat(0);
    TLegend *leg = new TLegend(0.7, 0.4, 0.9, 0.9);
    leg->SetHeader("PWA");
    leg->AddEntry(signal_h, "Data", "le");
    leg->AddEntry(all_h, "Total", "l");
    signal_h->Draw("e");
    all_h->Draw("same");

    hs->Add(signal_h);
    hs->Add(all_h);
    double weight_all = 0;
    for(unsigned int i = 2; i < weight.size(); i++) {
        weight_all += weight[i];
        cout << resNameList[i - 2] << " has " << weight[i] << endl;
        cout << "weight_all = " << weight_all << endl;
    }
    for(unsigned int i = 2; i < weight.size(); i++) {
        cout << "fraction of " << resNameList[i - 2] << " = " <<  weight[i] / weight[1] << endl;
    }
    int iw = 2;
    for(vector<string>::iterator it = resNameList.begin(); it != resNameList.end(); it++) {
        TH1F *res_h = hist_kk(path +"phsp_pwa_kk_weight_" + *it + ".root", "ana_tr", "kk", *it, var,  _dmin, _dmax);
        res_h->Scale(weight[0] / weight[1]);
        //res_h->Scale(weight[0] * weight[iw] / (weight[1] * weight[1]));
        //res_h->Scale(scaleFactor * weight[iw] / weight_all);
        //res_h->Scale(scaleFactor);
        if (tcolor == 10) tcolor++;
        cout << *it << "--color--" << tcolor << endl;
        res_h->SetLineColor(tcolor++);
        res_h->SetLineWidth(2);
        res_h->Draw("same");
        hs->Add(res_h);
        TString ss = *it;
        leg->AddEntry(res_h, ss + " " + Double2Str(valid_digits(weight[iw] / weight[1], 3) * 100) + "%", "l");
        iw++;
    }
    leg->SetFillStyle(kFDotted1);
    leg->Draw();
    //hs->Draw();
}


void UpdateFigure(TString var, string cfgFile) {
    draw_kk(var, cfgFile);
//    draw_kk(var, "/public/users/caihao1/tmp/baseline_small/baseline_small_read_conf");
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
}
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
