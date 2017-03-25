#ifndef FITPROXY_H
#define FITPROXY_H

#include "DPFPWAPdf.h"
#include "DPFPWAPoint.h"
#include "TString.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "PWA_CTRL.H"


//void Info(RooRealVar *bb);

class fitproxy {
public:
    TString outf_phipp, outf_phikk;
    TString data_phipp, data_phikk;
    TString proj_phipp, proj_phikk;
    TString phsp_phipp, phsp_phikk;
    TString idx_pp, idx_kk;
    RooArgSet theSet;

    vector<RooRealVar> vv;
    vector< vector<RooRealVar> > vpars;
    vector<TString> vn;

    RooRealVar idp;

    RooArgSet allparas, fitparas;
    RooDataSet *datapp, *datakk;

    DPFPWAPoint *dphipp, *dphikk;
    DPFPWAPdf *pdfphipp, *pdfphikk;

    fitproxy();
    ~fitproxy() {
        cout << "*****" << endl;
        cout << "delete dphipp" << endl;
        delete dphipp;
        delete dphikk;
        cout << "delete pdfphipp" << endl;
        delete pdfphipp;
        delete pdfphikk;
        cout << "delete datapp" << endl;
        datapp->Delete();
        delete datapp;
        datakk->Delete();
        delete datakk;
        cout << "End delete" << endl;
        vpars.resize(0);
        vv.resize(0);
        vn.resize(0);
        cout << "end vpars resize" << endl;
        theSet.Delete();
        allparas.Delete();
        fitparas.Delete();
        vv.resize(0);
        vpars.resize(0);
        vn.resize(0);
    };
    void init_all_files_name(const PWA_CTRL&);
    void init_input_argset();
    void read_data();
    void init_pdf(const PWA_CTRL&);
    void setup_resonances();
    void print_all_paras();
    void print_fit_paras();
    void add_res0_list(TString, double, double, double, double);
    void add_res2_list(TString, double, double, double, double);
    void add_res980_list(TString, double, double);
    void add_res1m_list(TString, double, double, double, double);
    void add_res1p_list(TString, double, double, double, double);
    void createlist_allparas();
    void add_fitparas(TString rn);

    void act_res0(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res2(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res980(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res1m(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res1p(TString rn, bool act_pp = true, bool act_kk = true);
    RooRealVar* gp(TString tn);

    void set_all_constant(TString rn);
    void set_all_constant();
    void store_fit_paras(TString fname);
    void store_all_paras(TString fname);
    void reload_paras(TString fname);

    void FIT();

};

#endif


