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

#include "../codes_frac/phipipi_structure.cc"
#include "../codes_frac/PWA_CTRL.C"
#include "../codes_frac/common_tools.cc"

#include <iomanip>

using namespace std;

void get_para_pipi(TString fname, TString tname, TString bname, TString var, double& sum_weight, double& dmin, double& dmax) {
    cout << "fname = " << fname << endl;
    cout << "var = " << var << endl;
    TFile *ff = new TFile(fname);
    TTree *tt = (TTree*)ff->Get(tname);
    DATA_ANA_PIPI dp;
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
//        if (value_pipi("Mphi", ss) > 1.032) continue;
//        if (value_pipi("Mphi", ss) < 1.006) continue;
        double tmp = value_pipi(var, dp);
        dmin = dmin < tmp ? dmin : tmp;
        dmax = dmax > tmp ? dmax : tmp;
        sum_weight += dp.weight;
    }
    cout << "Done parameters!" << endl;
    ff->Close();
}

TH1F* hist_pipi(TString fname, TString tname, TString bname, TString resName, TString var, double dmin, double dmax) {
    TH1F* hh = new TH1F(var, TexName_PIPI(var) + "_" +resName, 100, dmin, dmax);
    TFile *ff = new TFile(fname);
    TTree *tt = (TTree*)ff->Get(tname);
    DATA_ANA_PIPI dp;
    TBranch *bb;   //!
    tt->SetMakeClass(1);
    tt->SetBranchAddress(bname, &dp, &bb);

    Long64_t nentries = tt->GetEntriesFast();
    Long64_t nb = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        nb = tt->GetEntry(jentry);
//        if (!good_event(ss)) continue;
//        if (value_pipi("Mphi", ss) > 1.032) continue;
//        if (value_pipi("Mphi", ss) < 1.006) continue;
        hh->Fill(value_pipi(var, dp), dp.weight);
    }
    return hh;
    ff->Close();
}

void draw_pipi(TString var, string cfgFile) {
    vector<float> weight, dmin, dmax;
    double _weight, _dmin, _dmax;
    string path;
    readConfigFile(cfgFile, "work_path", path);
    cout << "path = " << path << endl;
    get_para_pipi(path + "data_ana_pipi_signal.root", "signal_tr", "signal_ana_pipi", var, _weight, _dmin, _dmax);
    weight.push_back(_weight);
    dmin.push_back(_dmin);
    dmax.push_back(_dmax);
    get_para_pipi(path + "phsp_pwa_pipi_weight_all.root", "ana_tr", "pp", var, _weight, _dmin, _dmax);
    weight.push_back(_weight);
    dmin.push_back(_dmin);
    dmax.push_back(_dmax);
    string value, value1;
    readConfigFile(cfgFile, "active_resonances", value);
    readConfigFile(cfgFile, "active_resonances_pp", value1);
    if (value1 != "NONE") {
        value = value + "," + value1;
    }
    vector<string> resNameList;
    string_to_vector(value, resNameList);
    for(vector<string>::iterator it = resNameList.begin(); it != resNameList.end(); it++) {
        get_para_pipi(path + "phsp_pwa_pipi_weight_" + *it + ".root", "ana_tr", "pp", var, _weight, _dmin, _dmax);
        weight.push_back(_weight);
        dmin.push_back(_dmin);
        dmax.push_back(_dmax);
    }
    _dmin = *min_element(dmin.begin(), dmin.end());
    _dmax = *max_element(dmax.begin(), dmax.end());

    double scaleFactor = weight[0] / weight[1];
    int tcolor = 1;
    TH1F *signal_h = hist_pipi(path + "data_ana_pipi_signal.root", "signal_tr", "signal_ana_pipi", "signal", var, _dmin, _dmax);
    signal_h->SetLineColor(tcolor++);
    signal_h->SetLineWidth(3);
    signal_h->GetXaxis()->SetTitle(TexName_PIPI(var));

    double sp = valid_digits((_dmax - _dmin) / 100, 2);
    signal_h->GetYaxis()->SetTitle("Entries / " + Double2Str(sp));

    TH1F *all_h = hist_pipi(path + "phsp_pwa_pipi_weight_all.root", "ana_tr", "pp", "all", var, _dmin, _dmax);
    all_h->SetLineColor(tcolor++);
    all_h->SetLineWidth(2);
    all_h->GetXaxis()->SetTitle(TexName_PIPI(var));
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

    double weight_all = 0;
    for(unsigned int i = 2; i < weight.size(); i++) {
        weight_all += weight[i];
        cout << "weight[i] = " << weight[i] << endl;
        cout << "weight_all = " << weight_all << endl;
    }
    for(unsigned int i = 2; i < weight.size(); i++) {
        cout << resNameList[i - 2] << " has " << weight[i] / weight_all << endl;
    }
    int iw = 2;
    for(vector<string>::iterator it = resNameList.begin(); it != resNameList.end(); it++) {
        TH1F *res_h = hist_pipi(path +"phsp_pwa_pipi_weight_" + *it + ".root", "ana_tr", "pp", *it, var,  _dmin, _dmax);
        res_h->Scale(weight[0] / weight[1]);
        //res_h->Scale(weight[0] * weight[iw] / (weight[1] * weight[1]));
        //iw++;
        //res_h->Scale(scaleFactor * weight[iw] / weight_all);
        //res_h->Scale(scaleFactor);
        if (tcolor == 10) tcolor++;
        res_h->SetLineColor(tcolor++);
        res_h->SetLineWidth(2);
        res_h->Draw("same");
        TString ss = *it;
        leg->AddEntry(res_h, ss + " " + Double2Str(valid_digits(weight[iw] / weight[1], 3) * 100) + "%", "l");
        iw++;
    }
    leg->SetFillStyle(kFDotted1);
    leg->Draw();
    //hs->Draw();

}


void UpdateFigure(TString var, string cfgFile) {
    draw_pipi(var, cfgFile);
    //draw_pipi(var, "/public/users/caihao1/tmp/baseline_small/baseline_small_read_conf");
//    TFile* ff = new TFile("../"+ path + "/pipi_hists.root", "RECREATE");
//    TH1F* th = new TH1F("test", "test", 100, -1, 1);
//    th->Write();
//    ff->Close();
//    draw_pipi("MKppim", path);
//    draw_pipi("MKmpip", path);
//    draw_pipi("Mpippim", path);
//    draw_pipi("MKpKm", path);
//    draw_pipi("Mphipip", path);
//    draw_pipi("Mphipim", path);
//    draw_pipi("QKp", path);
//    draw_pipi("QKm", path);
//    draw_pipi("Qpip", path);
//    draw_pipi("Qpim", path);
//    draw_pipi("QKpKm", path);
//    draw_pipi("Qpippim", path);
}
void test(TString var, TString mycut, Int_t bin_num = 100, Float_t dmin = -100, Float_t dmax = -100)
{
    TFile* f1 = TFile::Open("../phipipi/data_ana_pipi_signal.root");
    TTree* t1 = (TTree*)f1->Get("signal_tr");
    TFile* f2 = TFile::Open("../phipipi/data_ana_pipi_sideband.root");
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
    TH1F *h1 = new TH1F("h1", TexName_PIPI(var), bin_num, dmin, dmax);
    h1->GetXaxis()->SetTitle(TexName_PIPI(var));
    h1->GetYaxis()->SetTitle("Entries");
    TH1F *h2 = new TH1F("h2", TexName_PIPI(var), bin_num, dmin, dmax);
    h2->GetXaxis()->SetTitle(TexName_PIPI(var));
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
    c1->SaveAs("pipi_signal_sideband_PNG/" + var + ".png");
}

void UpdataFigure_sideband() {
    test("Mpippim", "");
    test("MKppim", "");
    test("MKmpip", "");
    test("Mpippim", "");
    test("MKpKm", "");
    test("Mphipip", "");
    test("Mphipim", "");
    test("QKp", "");
    test("QKm", "");
    test("Qpip", "");
    test("Qpim", "");
    test("QKpKm", "");
    test("Qpippim", "");
}

