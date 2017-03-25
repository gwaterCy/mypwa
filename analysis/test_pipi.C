#include "phipipi_structure.h"

TString con1 = "Chi2_4c<60&(MKmpip<0.744|MKmpip>1.048)&(MKppim<0.744|MKppim>1.048)";
TString con2 = "Mphi>1.006&Mphi<1.032&" + con1;
TString con3 = "Mphi>1.006&Mphi<1.032&(MKmpip<0.744|MKmpip>1.048)&(MKppim<0.744|MKppim>1.048)";


void test(TString var, TString mycut, Int_t bin_num = 100, Float_t dmin = -100, Float_t dmax = -100)
{
    TFile* f1 = TFile::Open("data_ana_pipi.root");
    TTree* t1 = (TTree*)f1->Get("data_tr");
    TFile* f2 = TFile::Open("topo_ana_pipi.root");
    TTree* t2 = (TTree*)f2->Get("topo_tr");
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
    TH1F *h1 = new TH1F("h1", var + " as " + mycut, bin_num, dmin, dmax);
    h1->GetXaxis()->SetTitle(TexName_PIPI(var));
    h1->GetYaxis()->SetTitle("Entries");
    TH1F *h2 = new TH1F("h2", var + " as " + mycut, bin_num, dmin, dmax);
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
    if (h1->GetMaximum() > h2->GetMaximum()) {
        h1->Draw("e");
        h2->Draw("same");
    } else {
        h2->Draw("");
        h1->Draw("esame");
    }
}
void test2(TString mycut, Int_t bin_num = 100, Float_t dmin = -1, Float_t dmax = -1)
{
    TFile* f1 = TFile::Open("data_ana_pipi.root");
    TTree* t1 = (TTree*)f1->Get("data_tr");
    TFile* f2 = TFile::Open("topo_ana_pipi.root");
    TTree* t2 = (TTree*)f2->Get("topo_tr");
//    if (mycut == "") {
//        mycut = "Chi2_4c<60&(MKmpip<0.744|MKmpip>1.048)&(MKppim<0.744|MKppim>1.048)";
//    } else {
//        mycut += "&Chi2_4c<60&(MKmpip<0.744|MKmpip>1.048)&(MKppim<0.744|MKppim>1.048)";
//    }
    cout << mycut << endl;
    new TCanvas;
    t1->Draw("M2phipip:M2phipim", mycut);
    new TCanvas;
    t2->Draw("M2phipip:M2phipim", mycut);

}
void test3(TString mycut, Int_t bin_num = 100, Float_t dmin = -1, Float_t dmax = -1)
{
    TFile* f1 = TFile::Open("data_ana_pipi.root");
    TTree* t1 = (TTree*)f1->Get("data_tr");
    TFile* f2 = TFile::Open("topo_ana_pipi.root");
    TTree* t2 = (TTree*)f2->Get("topo_tr");
    if (mycut == "") {
        mycut = "Chi2_4c<60&(MKmpip<0.744|MKmpip>1.048)&(MKppim<0.744|MKppim>1.048)";
    } else {
        mycut += "&Chi2_4c<60&(MKmpip<0.744|MKmpip>1.048)&(MKppim<0.744|MKppim>1.048)";
    }
    new TCanvas;
    t1->Draw("MKpKm:Mpippim", mycut);
    new TCanvas;
    t2->Draw("MKpKm:Mpippim", mycut);

}

void draw_Mphi() {
    c1 = new TCanvas();
    test("Mphi", con1, 100, 0.98, 1.09);
    c1->SaveAs("PNG/Mphi.png");
    c1->SaveAs("EPS/Mphi.eps");
}
void draw_MKpKm() {
    c1 = new TCanvas();
    test("MKpKm", con1, 100, 0.98, 1.09);
    c1->SaveAs("PNG/MKpKm.png");
    c1->SaveAs("EPS/MKpKm.eps");
}
void draw_MKmpip() {
    c1 = new TCanvas();
    test("MKmpip", con2);
    c1->SaveAs("PNG/MKmpip.png");
    c1->SaveAs("EPS/MKmpip.eps");
}
void draw_MKppim() {
    c1 = new TCanvas();
    test("MKppim", con2);
    c1->SaveAs("PNG/MKppim.png");
    c1->SaveAs("EPS/MKppim.eps");
}
void draw_M2phipip() {
    c1 = new TCanvas();
    test("M2phipip", con2);
    c1->SaveAs("PNG/M2phipip.png");
    c1->SaveAs("EPS/M2phipip.eps");
}
void draw_M2phipim() {
    c1 = new TCanvas();
    test("M2phipim", con2);
    c1->SaveAs("PNG/M2phipim.png");
    c1->SaveAs("EPS/M2phipim.eps");
}
void draw_Mpippim() {
    c1 = new TCanvas();
    test("Mpippim", con2);
    c1->SaveAs("PNG/Mpippim.png");
    c1->SaveAs("EPS/Mpippim.eps");
}
void draw_Mpippim_f0() {
    c1 = new TCanvas();
    test("Mpippim", con2, 50, 1.5, 2.5);
    c1->SaveAs("PNG/Mpippim_f0.png");
    c1->SaveAs("EPS/Mpippim_f0.eps");
}
void draw_Chi2_4c() {
    c1 = new TCanvas();
    test("Chi2_4c", con3, 100, 0);
    c1->SaveAs("PNG/Chi2_4c.png");
    c1->SaveAs("EPS/Chi2_4c.eps");
}
void draw_Mphipip() {
    c1 = new TCanvas();
    test("Mphipip", con2);
    c1->SaveAs("PNG/Mphipip.png");
    c1->SaveAs("EPS/Mphipip.eps");
}
void draw_Mphipim() {
    c1 = new TCanvas();
    test("Mphipim", con2);
    c1->SaveAs("PNG/Mphipim.png");
    c1->SaveAs("EPS/Mphipim.eps");
}
void draw_QKp() {
    c1 = new TCanvas();
    test("QKp", con2, 100, -1, 1);
    c1->SaveAs("PNG/QKp.png");
    c1->SaveAs("EPS/QKp.eps");
}
void draw_QKm() {
    c1 = new TCanvas();
    test("QKm", con2, 100, -1, 1);
    c1->SaveAs("PNG/QKm.png");
    c1->SaveAs("EPS/QKm.eps");
}
void draw_Qpip() {
    c1 = new TCanvas();
    test("Qpip", con2);
    c1->SaveAs("PNG/Qpip.png");
    c1->SaveAs("EPS/Qpip.eps");
}
void draw_Qpim() {
    c1 = new TCanvas();
    test("Qpim", con2);
    c1->SaveAs("PNG/Qpim.png");
    c1->SaveAs("EPS/Qpim.eps");
}
void draw_QKpKm() {
    c1 = new TCanvas();
    test("QKpKm", con2);
    c1->SaveAs("PNG/QKpKm.png");
    c1->SaveAs("EPS/QKpKm.eps");
}
void draw_Qpippim() {
    c1 = new TCanvas();
    test("Qpippim", con2);
    c1->SaveAs("PNG/Qpippim.png");
    c1->SaveAs("EPS/Qpippim.eps");
}

void UpdateAllFigures() {
	draw_Mphi();
	draw_MKpKm();
	draw_MKmpip();
	draw_MKppim();
	draw_M2phipip();
	draw_M2phipim();
	draw_Mpippim();
	draw_Mpippim_f0();
	draw_Chi2_4c();
	draw_Mphipip();
	draw_Mphipim();
	draw_QKp();
	draw_QKm();
	draw_Qpip();
	draw_Qpim();
	draw_QKpKm();
	draw_Qpippim();

}

