#include "TROOT.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "DPFPWAPdf.h"
#include "DPFPWAPoint.h"
//#include "RooStringVar.h"
//#include "RooAbsPdf.h"
//#include "RooAbsData.h"
//#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPrintable.h"
#include "TTree.h"
//#include "RooRealProxy.h"
//#include "RooListProxy.h"
//#include "RooComplex.h"
//#include "RooFitResult.h"
//#include "RooArgSet.h"
//#include "RooRealConstant.h"

#include <omp.h>

void fillpipi(RooFitResult *, DPFPWAPdf &, DPFPWAPdf &, TString, TString);
void fillkk(RooFitResult *, DPFPWAPdf &, DPFPWAPdf &, TString, TString);

int main()
{
//    gSystem->Load("libRooFit");
//    gSystem->Load("libPhysics");
//    gSystem->Load("../tmp/libDPF.so");
    using namespace RooFit;
    //   gROOT->ProcessLine(".L PWAPdf.cxx+");
    //   gROOT->ProcessLine(".L bkgpdf.cxx+");
    const Double_t mpsip=3.686,mka=0.493677,mpi=0.13957;
    const Double_t pi=3.14159265358979312;
    // double high=mjpsi/2;
    const Double_t high=mpsip;
    const Double_t low=0-high;
    const Double_t Nreal=50000;

    TString indata,bkgdata;
    TString outfitpipi="outpipi.rep";
    TString outfit="outkk.rep";
    TString indatapipi="../phipipi/data_pwa_pipi_weight.dat";
    TString indatakaka="../phikk/data_pwa_kk_weight.dat";
    TString projectfilepipi="resultpipi.root";
    TString projectfilekaka="resultkaka.root";
    bkgdata="kstarkp.dat";

    RooRealVar v11( "v11", "11", low, high);
    RooRealVar v12( "v12", "12", low, high);
    RooRealVar v13( "v13", "13", low, high);
    RooRealVar v14( "v14", "14", low, high);
    RooRealVar v21( "v21", "21", low, high);
    RooRealVar v22( "v22", "22", low, high);
    RooRealVar v23( "v23", "23", low, high);
    RooRealVar v24( "v24", "24", low, high);
    RooRealVar v31( "v31", "31", low, high);
    RooRealVar v32( "v32", "32", low, high);
    RooRealVar v33( "v33", "33", low, high);
    RooRealVar v34( "v34", "34", low, high);
    RooRealVar v41( "v41", "41", low, high);
    RooRealVar v42( "v42", "42", low, high);
    RooRealVar v43( "v43", "43", low, high);
    RooRealVar v44( "v44", "44", low, high);
    RooRealVar v51( "v51", "51", low, high);
    RooRealVar v52( "v52", "52", low, high);
    RooRealVar v53( "v53", "53", low, high);
    RooRealVar v54( "v54", "54", low, high);
    RooRealVar weight("weight","weight",low,high);

    //  RooArgSet theSet;
    //  theSet.add(RooArgSet(v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54,weight));
    RooArgSet theSet1,theSet2,theSet3;
    theSet1.add(RooArgSet(v11,v12,v13,v14,v21,v22,v23,v24));
    theSet2.add(RooArgSet(v31,v32,v33,v34,v41,v42,v43,v44));
    theSet3.add(RooArgSet(v51,v52,v53,v54,weight));
    RooArgSet theSet4(theSet1,theSet2,"");
    RooArgSet theSet(theSet4,theSet3,"");

    cout<<"-----------------"<<endl;
    //  RooDataSet *data = RooDataSet::read(indata,theSet,"Q");
    RooDataSet *data11 = RooDataSet::read(indatapipi,theSet);
    //  RooDataSet *data2 = RooDataSet::read(bkgdata,theSet);
    //  PWAPdf pdf("pdf","pdf", v11, v12, v13, v14, v21, v22, v23, v24);
    //  bkgpdf bkg("bkgpdf","bkgpdf", v11, v12, v13, v14, v21, v22, v23, v24);
    data11->Print();
    RooDataSet *datapipi= new RooDataSet(data11->GetName(),data11->GetTitle(),data11,*data11->get(),0,weight.GetName());
    datapipi->Print();

    cout<<"-----------------"<<endl;
    //  RooDataSet *data = RooDataSet::read(indata,theSet,"Q");
    RooDataSet *data22 = RooDataSet::read(indatakaka,theSet);
    //  RooDataSet *data2 = RooDataSet::read(bkgdata,theSet);
    //  PWAPdf pdf("pdf","pdf", v11, v12, v13, v14, v21, v22, v23, v24);
    //  bkgpdf bkg("bkgpdf","bkgpdf", v11, v12, v13, v14, v21, v22, v23, v24);
    data22->Print();
    RooDataSet *datakaka= new RooDataSet(data22->GetName(),data22->GetTitle(),data22,*data22->get(),0,weight.GetName());
    datakaka->Print();
    cout<<"RooDataSet is initialized!!"<<endl;

    DPFPWAPoint dphipipi(mka,mka, mpi, mpi, mpsip, "../phipipi/phsp_pwa_pipi.dat");
    DPFPWAPoint dphikaka(mka,mka, mka, mka, mpsip, "../phikk/phsp_pwa_kk.dat");
    cout<<"DPFPWAPoint is intialized!!"<<endl;
    DPFPWAPdf pdfpipi("pdfpipi","pdfpipi",v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54,&dphipipi);
    DPFPWAPdf pdfkaka("pdfkaka","pdfkaka",v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54,&dphikaka);
    cout<<"DPFPWAPdf is intialized!!"<<endl;


//  construct_resonances(pdfpipi, pdfkaka);
//  EM radiactive
    RooRealVar f01000a_spn_("f01000a_spn_", "f01000a_spn_",       1);
    RooRealVar f01000a_mss_("f01000a_mss_", "f01000a_mss_",   1.350,   1.200,   1.500);
    RooRealVar f01000a_wdt_("f01000a_wdt_", "f01000a_wdt_",   10000,       0,   20000);
    RooRealVar f01000p_rho_("f01000p_rho_", "f01000p_rho_", 65.8879,    -100,     100);
    RooRealVar f01000p_phi_("f01000p_phi_", "f01000p_phi_",  0.4468,     -pi,      pi);
    RooRealVar f01000a_typ_("f01000a_typ_", "f01000a_typ_",       1);
    //         f01000a_spn_.setConstant();
    //         f01000a_mss_.setConstant();
    //         f01000a_wdt_.setConstant();
               f01000p_rho_.setConstant();
               f01000p_phi_.setConstant();
    //         f01000a_typ_.setConstant();

//  f0(980) 0+ 1 pp
    RooRealVar f00980a_spn1("f00980a_spn1", "f00980a_spn1",       1);
    RooRealVar f00980a_mss_("f00980a_mss_", "f00980a_mss_",   0.965,     0.9,       1);
    RooRealVar f00980a_g10_("f00980a_g10_", "f00980a_g10_",   0.165,     0.1,     0.5);
    RooRealVar f00980a_g20_("f00980a_g20_", "f00980a_g20_",   0.695,     0.1,     1.5);
    RooRealVar f00980p_rho1("f00980p_rho1", "f00980p_rho1", 25.4288,    -100,     100);
    RooRealVar f00980p_phi1("f00980p_phi1", "f00980p_phi1",  4.6788,   -2*pi,    2*pi);
    RooRealVar f00980a_typ_("f00980a_typ_", "f00980a_typ_",       2);
               f00980a_spn1.setConstant();
    //         f00980a_mss_.setConstant();
    //         f00980a_g10_.setConstant();
    //         f00980a_g20_.setConstant();
    //         f00980p_rho1.setConstant();
    //         f00980p_phi1.setConstant();
               f00980a_typ_.setConstant();
//  f0(980) 0+ 2 pp
    RooRealVar f00980a_spn2("f00980a_spn2", "f00980a_spn2",       2);
    RooRealVar f00980p_rho2("f00980p_rho2", "f00980p_rho2",  8.0024,     -100,     100);
    RooRealVar f00980p_phi2("f00980p_phi2", "f00980p_phi2",  4.3268,    -2*pi,    2*pi);
               f00980a_spn2.setConstant();
    //         f00980p_rho2.setConstant();
    //         f00980p_phi2.setConstant();

//  f0(1370) 0+ 1 pp
    RooRealVar f01370a_spn1("f01370a_spn1", "f01370a_spn1",       1);
    RooRealVar f01370a_mss_("f01370a_mss_", "f01370a_mss_",    1.35,      1.2,     1.5);
    RooRealVar f01370a_wdt_("f01370a_wdt_", "f01370a_wdt_",   0.265,      0.2,     0.5);
    RooRealVar f01370p_rho1("f01370p_rho1", "f01370p_rho1", 65.8879,     -100,     100);
    RooRealVar f01370p_phi1("f01370p_phi1", "f01370p_phi1",  0.4468,      -pi,      pi);
    RooRealVar f01370a_typ_("f01370a_typ_", "f01370a_typ_",       1);
               f01370a_spn1.setConstant();
               f01370a_mss_.setConstant();
               f01370a_wdt_.setConstant();
               f01370p_rho1.setConstant();
               f01370p_phi1.setConstant();
               f01370a_typ_.setConstant();
//  f0(1370) 0+ 2 pp
    RooRealVar f01370a_spn2("f01370a_spn2", "f01370a_spn2",       2);
    RooRealVar f01370p_rho2("f01370p_rho2", "f01370p_rho2", 16.1689,     -100,     100);
    RooRealVar f01370p_phi2("f01370p_phi2", "f01370p_phi2",  0.7596,      -pi,      pi);
               f01370a_spn2.setConstant();
               f01370p_rho2.setConstant();
               f01370p_phi2.setConstant();

//  f0(1750) 0+ 1 pp
    RooRealVar f01750a_spn1("f01750a_spn1", "f01750a_spn1",       1);
    RooRealVar f01750a_mss_("f01750a_mss_", "f01750a_mss_",    1.79,      1.7,     1.9);
    RooRealVar f01750a_wdt_("f01750a_wdt_", "f01750a_wdt_",    0.27,     0.15,    0.35);
    RooRealVar f01750p_rho1("f01750p_rho1", "f01750p_rho1", 69.9798,     -100,     100);
    RooRealVar f01750p_phi1("f01750p_phi1", "f01750p_phi1",   0.543,      -pi,      pi);
    RooRealVar f01750a_typ_("f01750a_typ_", "f01750a_typ_",       1);
    //         f01750a_mss_.setConstant();
    //         f01750a_wdt_.setConstant();
               f01750a_spn1.setConstant();
    //         f01750p_rho1.setConstant();
    //         f01750p_phi1.setConstant();
               f01750a_typ_.setConstant();
//  f0(1750) 0+ 2 pp
    RooRealVar f01750a_spn2("f01750a_spn2", "f01750a_spn2",       2);
    RooRealVar f01750p_rho2("f01750p_rho2", "f01750p_rho2", 39.7695,     -100,   100.0);
    RooRealVar f01750p_phi2("f01750p_phi2", "f01750p_phi2",  -0.377,      -pi,      pi);
               f01750a_spn2.setConstant();
    //         f01750p_rho2.setConstant();
    //         f01750p_phi2.setConstant();

//  f0(sigma600) 0+ 1 pp
    RooRealVar f0s600a_spn1("f0s600a_spn1", "f0s600a_spn1",       1);
    RooRealVar f0s600a_mss_("f0s600a_mss_", "f0s600a_mss_",  0.9264,      0.8,       1);
    RooRealVar f0s600a_b1__("f0s600a_b1__", "f0s600a_b1__",  0.5843,      0.5,     0.6);
    RooRealVar f0s600a_b2__("f0s600a_b2__", "f0s600a_b2__",  1.6663,      1.5,     1.7);
    RooRealVar f0s600a_b3__("f0s600a_b3__", "f0s600a_b3__",   1.082,      0.8,     1.5);
    RooRealVar f0s600a_b4__("f0s600a_b4__", "f0s600a_b4__",  0.0024,    0.001,   0.005);
    RooRealVar f0s600a_b5__("f0s600a_b5__", "f0s600a_b5__",       1,      0.5,     1.5);
    RooRealVar f0s600p_rho1("f0s600p_rho1", "f0s600p_rho1", 34.8889,     -100,     100);
    RooRealVar f0s600p_phi1("f0s600p_phi1", "f0s600p_phi1",  1.8491,      -pi,      pi);
    RooRealVar f0s600a_typ_("f0s600a_typ_", "f0s600a_typ_",       3);
               f0s600a_spn1.setConstant();
               f0s600a_mss_.setConstant();
               f0s600a_b1__.setConstant();
               f0s600a_b2__.setConstant();
               f0s600a_b3__.setConstant();
               f0s600a_b4__.setConstant();
               f0s600a_b5__.setConstant();
               f0s600p_rho1.setConstant();
               f0s600p_phi1.setConstant();
               f0s600a_typ_.setConstant();
//  f0(sigma600) 0+ 2 pp
    RooRealVar f0s600a_spn2("f0s600a_spn2", "f0s600a_spn2",       2);
    RooRealVar f0s600p_rho2("f0s600p_rho2", "f0s600p_rho2", 13.1566,     -100,     100);
    RooRealVar f0s600p_phi2("f0s600p_phi2", "f0s600p_phi2",  1.9269,      -pi,      pi);
               f0s600a_spn2.setConstant();
               f0s600p_rho2.setConstant();
               f0s600p_phi2.setConstant();

//  f0(1500) 0+ 1 pp
    RooRealVar f01500a_spn1("f01500a_spn1", "f01500a_spn1",       1);
    RooRealVar f01500a_mss_("f01500a_mss_", "f01500a_mss_",  1.4950,    1.00,     1.90);
    RooRealVar f01500a_wdt_("f01500a_wdt_", "f01500a_wdt_",  0.1220,   0.100,    0.150);
    RooRealVar f01500p_rho1("f01500p_rho1", "f01500p_rho1", 53.6250,  -100.0,    100.0);
    RooRealVar f01500p_phi1("f01500p_phi1", "f01500p_phi1",  2.7549,     -pi,       pi);
    RooRealVar f01500a_typ_("f01500a_typ_", "f01500a_typ_",       1);
               f01500a_spn1.setConstant();
    //         f01500a_mss_.setConstant();
    //         f01500a_wdt_.setConstant();
    //         f01500p_rho1.setConstant();
    //         f01500p_phi1.setConstant();
               f01500a_typ_.setConstant();
//  f0(1500) 0+ 2 pp
    RooRealVar f01500a_spn2("f01500a_spn2", "f01500a_spn2",       2);
    RooRealVar f01500p_rho2("f01500p_rho2", "f01500p_rho2", 23.4287,  -100.0,    100.0);
    RooRealVar f01500p_phi2("f01500p_phi2", "f01500p_phi2",  2.5820,     -pi,       pi);
               f01500a_spn2.setConstant();
    //         f01500p_rho2.setConstant();
    //         f01500p_phi2.setConstant();

//  f0(1270) 2+ 1 pp
    RooRealVar f01270a_spn1("f01270a_spn1", "f01270a_spn1",      21);
    RooRealVar f01270a_mss_("f01270a_mss_", "f01270a_mss_",  1.2750,     1.0,    1.800);
    RooRealVar f01270a_wdt_("f01270a_wdt_", "f01270a_wdt_",  0.1900,    0.10,     0.90);
    RooRealVar f01270p_rho1("f01270p_rho1", "f01270p_rho1", 23.4619,  -100.0,      100);
    RooRealVar f01270p_phi1("f01270p_phi1", "f01270p_phi1",  0.6866,     -pi,       pi);
    RooRealVar f01270a_typ_("f01270a_typ_", "f01270a_typ_",       6);
               f01270a_mss_.setConstant();
               f01270a_wdt_.setConstant();
               f01270a_spn1.setConstant();
               f01270p_rho1.setConstant();
               f01270p_phi1.setConstant();
               f01270a_typ_.setConstant();
//  f0(1270) 2+ 2 pp
    RooRealVar f01270a_spn2("f01270a_spn2", "f01270a_spn2",      22);
    RooRealVar f01270p_rho2("f01270p_rho2", "f01270p_rho2",  7.5520,    -100,      100);
    RooRealVar f01270p_phi2("f01270p_phi2", "f01270p_phi2", -0.1037,     -pi,       pi);
               f01270a_spn2.setConstant();
               f01270p_rho2.setConstant();
               f01270p_phi2.setConstant();
//  f0(1270) 2+ 3 pp
    RooRealVar f01270a_spn3("f01270a_spn3", "f01270a_spn3",      23);
    RooRealVar f01270p_rho3("f01270p_rho3", "f01270p_rho3",  5.3346,    -100,      100);
    RooRealVar f01270p_phi3("f01270p_phi3", "f01270p_phi3", -0.0750,     -pi,       pi);
               f01270a_spn3.setConstant();
               f01270p_rho3.setConstant();
               f01270p_phi3.setConstant();
//  f0(1270) 2+ 4 pp
    RooRealVar f01270a_spn4("f01270a_spn4", "f01270a_spn4",      24);
    RooRealVar f01270p_rho4("f01270p_rho4", "f01270p_rho4",  4.3929,    -100,      100);
    RooRealVar f01270p_phi4("f01270p_phi4", "f01270p_phi4",  3.2811,   -2*pi,     2*pi);
               f01270a_spn4.setConstant();
               f01270p_rho4.setConstant();
               f01270p_phi4.setConstant();
//  f0(1270) 2+ 5 pp
    RooRealVar f01270a_spn5("f01270a_spn5", "f01270a_spn5",      25);
    RooRealVar f01270p_rho5("f01270p_rho5", "f01270p_rho5",  4.3929,    -100,     100);
    RooRealVar f01270p_phi5("f01270p_phi5", "f01270p_phi5",  3.2811,   -2*pi,    2*pi);
               f01270a_spn5.setConstant();
               f01270p_rho5.setConstant();
               f01270p_phi5.setConstant();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  EM radiactive
    RooRealVar f01000k_rho_("f01000k_rho_", "f01000k_rho_", 65.8879,    -100,     100);
    RooRealVar f01000k_phi_("f01000k_phi_", "f01000k_phi_",  0.4468,     -pi,      pi);
    //         f01000k_rho_.setConstant();
    //         f01000k_phi_.setConstant();

//  f0(980) 0+ 1 kk 
    RooRealVar f00980k_rho1("f00980k_rho1", "f00980k_rho1",  0.7628,  -100.0,    100.0);
    RooRealVar f00980k_phi1("f00980k_phi1", "f00980k_phi1",  1.2817,     -pi,       pi);
    //         f00980k_rho1.setConstant();
    //         f00980k_phi1.setConstant();
//  f0(980) 0+ 2 kk 
    RooRealVar f00980k_rho2("f00980k_rho2", "f00980k_rho2",  0.2400,  -100.0,    100.0);
    RooRealVar f00980k_phi2("f00980k_phi2", "f00980k_phi2",  1.0694,     -pi,       pi);
    //         f00980k_rho2.setConstant();
    //         f00980k_phi2.setConstant();

//  f0(1370) 0+ 1 kk
    RooRealVar f01370k_rho1("f01370k_rho1", "f01370k_rho1",  0.3360,  -100.0,    100.0);
    RooRealVar f01370k_phi1("f01370k_phi1", "f01370k_phi1",  3.3365,   -2*pi,     2*pi);
               f01370k_rho1.setConstant();
               f01370k_phi1.setConstant();
//  f0(1370) 0+ 2 kk
    RooRealVar f01370k_rho2("f01370k_rho2", "f01370k_rho2", 0.08246,  -100.0,    100.0);
    RooRealVar f01370k_phi2("f01370k_phi2", "f01370k_phi2", -2.7921,     -pi,       pi);
               f01370k_rho2.setConstant();
               f01370k_phi2.setConstant();

//  f0(1750) 0+ 1 kk 
    RooRealVar f01750k_rho1("f01750k_rho1", "f01750k_rho1",  0.5598,  -100.0,    100.0);
    RooRealVar f01750k_phi1("f01750k_phi1", "f01750k_phi1",  5.6406,   -2*pi,     2*pi);
    //         f01750k_rho1.setConstant();
    //         f01750k_phi1.setConstant();
//  f0(1750) 0+ 2 kk
    RooRealVar f01750k_rho2("f01750k_rho2", "f01750k_rho2",  0.3181,  -100.0,    100.0);
    RooRealVar f01750k_phi2("f01750k_phi2", "f01750k_phi2",  4.2865,   -2*pi,     2*pi);
    //         f01750k_rho2.setConstant();
    //         f01750k_phi2.setConstant();

//  f0(sigma) 0+ 1 kk
    RooRealVar f0s600k_rho1("f0s600k_rho1", "f0s600k_rho1",  0.5477,  -100.0,    100.0);
    RooRealVar f0s600k_phi1("f0s600k_phi1", "f0s600k_phi1",  3.0705,   -2*pi,     2*pi);
               f0s600k_rho1.setConstant();
               f0s600k_phi1.setConstant();
//  f0(sigma) 0+ 2 kk
    RooRealVar f0s600k_rho2("f0s600k_rho2", "f0s600k_rho2",  0.2065,  -100.0,    100.0);
    RooRealVar f0s600k_phi2("f0s600k_phi2", "f0s600k_phi2",  2.1623,     -pi,       pi);
               f0s600k_rho2.setConstant();
               f0s600k_phi2.setConstant();

//  f0(1500) 0+ 1 kk 
    RooRealVar f01500k_rho1("f01500k_rho1", "f01500k_rho1",  0.6166,  -100.0,    100.0);
    RooRealVar f01500k_phi1("f01500k_phi1", "f01500k_phi1",  3.9934,   -2*pi,     2*pi);
    //         f01500k_rho1.setConstant();
    //         f01500k_phi1.setConstant();
//  f0(1500) 0+ 2 kk
    RooRealVar f01500k_rho2("f01500k_rho2", "f01500k_rho2",  0.2694,  -100.0,    100.0);
    RooRealVar f01500k_phi2("f01500k_phi2", "f01500k_phi2",  4.3124,   -2*pi,     2*pi);
    //         f01500k_rho2.setConstant();
    //         f01500k_phi2.setConstant();

//  f0(1710) 0+ 1 kk 
    RooRealVar f01710a_spn1("f01710a_spn1", "f01710a_spn1",       1);
    RooRealVar f01710a_mss_("f01710a_mss_", "f01710a_mss_",  1.7070,   1.650,    1.800);
    RooRealVar f01710a_wdt_("f01710a_wdt_", "f01710a_wdt_",  0.1250,    0.10,   0.1500);
    RooRealVar f01710k_rho1("f01710k_rho1", "f01710k_rho1",  1.0527,  -100.0,    100.0);
    RooRealVar f01710k_phi1("f01710k_phi1", "f01710k_phi1",  6.8807,   -3*pi,     3*pi);
    RooRealVar f01710a_typ_("f01710a_typ_", "f01710a_typ_",       1);
               f01710a_spn1.setConstant();
    //         f01710a_mss_.setConstant();
    //         f01710a_wdt_.setConstant();
    //         f01710k_rho1.setConstant();
    //         f01710k_phi1.setConstant();
               f01710a_typ_.setConstant();
//  f0(1710) 0+ 2 kk 
    RooRealVar f01710a_spn2("f01710a_spn2", "f01710a_spn2",       2);
    RooRealVar f01710k_rho2("f01710k_rho2", "f01710k_rho2",  0.2545,  -100.0,    100.0);
    RooRealVar f01710k_phi2("f01710k_phi2", "f01710k_phi2", -0.0052,     -pi,       pi);
               f01710a_spn2.setConstant();
    //         f01710k_rho2.setConstant();
    //         f01710k_phi2.setConstant();

//  f0(1270) 2+ 1 kk
    RooRealVar f01270k_rho1("f01270k_rho1", "f01270k_rho1",  0.1000,  -100.0,    100.0);
    RooRealVar f01270k_phi1("f01270k_phi1", "f01270k_phi1",  2.8428,   -2*pi,     2*pi);
               f01270k_rho1.setConstant();
               f01270k_phi1.setConstant();
//  f0(1270) 2+ 2 kk
    RooRealVar f01270k_rho2("f01270k_rho2", "f01270k_rho2",  0.0322,  -100.0,    100.0);
    RooRealVar f01270k_phi2("f01270k_phi2", "f01270k_phi2", -1.9102,     -pi,       pi);
               f01270k_rho2.setConstant();
               f01270k_phi2.setConstant();
//  f0(1270) 2+ 3 kk
    RooRealVar f01270k_rho3("f01270k_rho3", "f01270k_rho3",  0.0226,  -100.0,    100.0);
    RooRealVar f01270k_phi3("f01270k_phi3", "f01270k_phi3",  0.3628,     -pi,       pi);
               f01270k_rho3.setConstant();
               f01270k_phi3.setConstant();
//  f0(1270) 2+ 4 kk
    RooRealVar f01270k_rho4("f01270k_rho4", "f01270k_rho4",  0.0187,  -100.0,    100.0);
    RooRealVar f01270k_phi4("f01270k_phi4", "f01270k_phi4",  1.4555,     -pi,       pi);
               f01270k_rho4.setConstant();
               f01270k_phi4.setConstant();
//  f0(1270) 2+ 5 kk
    RooRealVar f01270k_rho5("f01270k_rho5", "f01270k_rho5",  0.0000,  -100.0,    100.0);
    RooRealVar f01270k_phi5("f01270k_phi5", "f01270k_phi5",  0.0000,     -pi,       pi);
               f01270k_rho5.setConstant();
               f01270k_phi5.setConstant();

//  f0(1525) 2+ 1 kk
    RooRealVar f01525a_spn1("f01525a_spn1", "f01525a_spn1",      21);
    RooRealVar f01525a_mss_("f01525a_mss_", "f01525a_mss_",  1.5210,   0.500,     2.50);
    RooRealVar f01525a_wdt_("f01525a_wdt_", "f01525a_wdt_",  0.0770,  0.0010,    0.800);
    RooRealVar f01525k_rho1("f01525k_rho1", "f01525k_rho1",  1.0000,  -100.0,    100.0);
    RooRealVar f01525k_phi1("f01525k_phi1", "f01525k_phi1",  0.2000,     -pi,       pi);
    RooRealVar f01525a_typ_("f01525a_typ_", "f01525a_typ_",       1);
    //         f01525a_mss_.setConstant();
    //         f01525a_wdt_.setConstant();
               f01525a_spn1.setConstant();
    //         f01525k_rho1.setConstant();
    //         f01525k_phi1.setConstant();
               f01525a_typ_.setConstant();
//  f0(1525) 2+ 2 kk
    RooRealVar f01525a_spn2("f01525a_spn2", "f01525a_spn2",      22);
    RooRealVar f01525k_rho2("f01525k_rho2", "f01525k_rho2",  0.5878,  -100.0,    100.0);
    RooRealVar f01525k_phi2("f01525k_phi2", "f01525k_phi2",  5.9124,   -3*pi,     3*pi);
               f01525a_spn2.setConstant();
    //         f01525k_rho2.setConstant();
    //         f01525k_phi2.setConstant();
//  f0(1525) 2+ 3 kk
    RooRealVar f01525a_spn3("f01525a_spn3", "f01525a_spn3",      23);
    RooRealVar f01525k_rho3("f01525k_rho3", "f01525k_rho3",  0.3493,  -100.0,    100.0);
    RooRealVar f01525k_phi3("f01525k_phi3", "f01525k_phi3",  0.1663,     -pi,       pi);
               f01525a_spn3.setConstant();
    //         f01525k_rho3.setConstant();
    //         f01525k_phi3.setConstant();
//  f0(1525) 2+ 4 kk
    RooRealVar f01525a_spn4("f01525a_spn4", "f01525a_spn4",      24);
    RooRealVar f01525k_rho4("f01525k_rho4", "f01525k_rho4",  0.0000,  -100.0,    100.0);
    RooRealVar f01525k_phi4("f01525k_phi4", "f01525k_phi4",  0.0000,     -pi,       pi);
               f01525a_spn4.setConstant();
    //         f01525k_rho4.setConstant();
    //         f01525k_phi4.setConstant();
//  f0(1525) 2+ 5 kk
    RooRealVar f01525a_spn5("f01525a_spn5", "f01525a_spn5",      25);
    RooRealVar f01525k_rho5("f01525k_rho5", "f01525k_rho5",  0.0000,  -100.0,    100.0);
    RooRealVar f01525k_phi5("f01525k_phi5", "f01525k_phi5",  0.0000,     -pi,       pi);
               f01525a_spn5.setConstant();
    //         f01525k_rho5.setConstant();
    //         f01525k_phi5.setConstant();

//  phi(1680) 1- 1 kk
    RooRealVar ph1680a_spn1("ph1680a_spn1", "ph1680a_spn1",     191);
    RooRealVar ph1680a_mss_("ph1680a_mss_", "ph1680a_mss_",  1.6150,   0.500,    2.500);
    RooRealVar ph1680a_m98_("ph1680a_m98_", "ph1680a_m98_",  0.9650,   0.600,    1.200);
    RooRealVar ph1680a_wdt_("ph1680a_wdt_", "ph1680a_wdt_",  0.1300,    0.01,     0.50);
    RooRealVar ph1680a_g1__("ph1680a_g1__", "ph1680a_g1__",  0.1650,    0.01,    0.500);
    RooRealVar ph1680a_g2__("ph1680a_g2__", "ph1680a_g2__",  0.6950,     0.1,      1.0);
    RooRealVar ph1680k_rho1("ph1680k_rho1", "ph1680k_rho1",  0.0192,  -100.0,    100.0);
    RooRealVar ph1680k_phi1("ph1680k_phi1", "ph1680k_phi1", -0.5701,     -pi,       pi);
    RooRealVar ph1680a_typ_("ph1680a_typ_", "ph1680a_typ_",       5);
               ph1680a_spn1.setConstant();
               ph1680a_mss_.setConstant();
               ph1680a_m98_.setConstant();
               ph1680a_wdt_.setConstant();
               ph1680a_g1__.setConstant();
               ph1680a_g2__.setConstant();
               ph1680k_rho1.setConstant();
               ph1680k_phi1.setConstant();
               ph1680a_typ_.setConstant();
//  phi(1680) 1- 2 kk
    RooRealVar ph1680a_spn2("ph1680a_spn2", "ph1680a_spn2",     192);
    RooRealVar ph1680k_rho2("ph1680k_rho2", "ph1680k_rho2",  0.0297,  -100.0,    100.0);
    RooRealVar ph1680k_phi2("ph1680k_phi2", "ph1680k_phi2",  0.1496,     -pi,       pi);
               ph1680a_spn2.setConstant();
               ph1680k_rho2.setConstant();
               ph1680k_phi2.setConstant();


    pdfpipi.addResonance   ("f01000p"  , "f01000p"  , f01000a_spn_, f01000a_mss_, f01000a_wdt_, f01000p_rho_, f01000p_phi_, f01000a_typ_);
    pdfpipi.addResonance600("f0s600p_1", "f0s600p_1", f0s600a_spn1, f0s600a_mss_, f0s600a_b1__, f0s600a_b2__, f0s600a_b3__, f0s600a_b4__, f0s600a_b5__, f0s600p_rho1, f0s600p_phi1, f0s600a_typ_);
    pdfpipi.addResonance600("f0s600p_2", "f0s600p_2", f0s600a_spn2, f0s600a_mss_, f0s600a_b1__, f0s600a_b2__, f0s600a_b3__, f0s600a_b4__, f0s600a_b5__, f0s600p_rho2, f0s600p_phi2, f0s600a_typ_);
    pdfpipi.addResonance980("f00980p_1", "f00980p_1", f00980a_spn1, f00980a_mss_, f00980a_g10_, f00980a_g20_, f00980p_rho1, f00980p_phi1, f00980a_typ_);
    pdfpipi.addResonance980("f00980p_2", "f00980p_2", f00980a_spn2, f00980a_mss_, f00980a_g10_, f00980a_g20_, f00980p_rho2, f00980p_phi2, f00980a_typ_);
    pdfpipi.addResonance   ("f01370p_1", "f01370p_1", f01370a_spn1, f01370a_mss_, f01370a_wdt_, f01370p_rho1, f01370p_phi1, f01370a_typ_);
    pdfpipi.addResonance   ("f01370p_2", "f01370p_2", f01370a_spn2, f01370a_mss_, f01370a_wdt_, f01370p_rho2, f01370p_phi2, f01370a_typ_);
    pdfpipi.addResonance   ("f01750p_1", "f01750p_1", f01750a_spn1, f01750a_mss_, f01750a_wdt_, f01750p_rho1, f01750p_phi1, f01750a_typ_);
    pdfpipi.addResonance   ("f01750p_2", "f01750p_2", f01750a_spn2, f01750a_mss_, f01750a_wdt_, f01750p_rho2, f01750p_phi2, f01750a_typ_);
    pdfpipi.addResonance   ("f01500p_1", "f01500p_1", f01500a_spn1, f01500a_mss_, f01500a_wdt_, f01500p_rho1, f01500p_phi1, f01500a_typ_);
    pdfpipi.addResonance   ("f01500p_2", "f01500p_2", f01500a_spn2, f01500a_mss_, f01500a_wdt_, f01500p_rho2, f01500p_phi2, f01500a_typ_);
    pdfpipi.addResonance   ("f01270p_1", "f01270p_1", f01270a_spn1, f01270a_mss_, f01270a_wdt_, f01270p_rho1, f01270p_phi1, f01270a_typ_);
    pdfpipi.addResonance   ("f01270p_2", "f01270p_2", f01270a_spn2, f01270a_mss_, f01270a_wdt_, f01270p_rho2, f01270p_phi2, f01270a_typ_);
    pdfpipi.addResonance   ("f01270p_3", "f01270p_3", f01270a_spn3, f01270a_mss_, f01270a_wdt_, f01270p_rho3, f01270p_phi3, f01270a_typ_);
    pdfpipi.addResonance   ("f01270p_4", "f01270p_4", f01270a_spn4, f01270a_mss_, f01270a_wdt_, f01270p_rho4, f01270p_phi4, f01270a_typ_);
    pdfpipi.addResonance   ("f01270p_5", "f01270p_5", f01270a_spn5, f01270a_mss_, f01270a_wdt_, f01270p_rho5, f01270p_phi5, f01270a_typ_);


    pdfkaka.addResonance   ("f01000k"  , "f01000k"  , f01000a_spn_, f01000a_mss_, f01000a_wdt_, f01000k_rho_, f01000k_phi_, f01000a_typ_);
    pdfkaka.addResonance600("f0s600k_1", "f0s600k_1", f0s600a_spn1, f0s600a_mss_, f0s600a_b1__, f0s600a_b2__, f0s600a_b3__, f0s600a_b4__, f0s600a_b5__, f0s600k_rho1, f0s600k_phi1, f0s600a_typ_);
    pdfkaka.addResonance600("f0s600k_2", "f0s600k_2", f0s600a_spn2, f0s600a_mss_, f0s600a_b1__, f0s600a_b2__, f0s600a_b3__, f0s600a_b4__, f0s600a_b5__, f0s600k_rho2, f0s600k_phi2, f0s600a_typ_);
    pdfkaka.addResonance980("f00980k_1", "f00980k_1", f00980a_spn1, f00980a_mss_, f00980a_g10_, f00980a_g20_, f00980k_rho1, f00980k_phi1, f00980a_typ_);
    pdfkaka.addResonance980("f00980k_2", "f00980k_2", f00980a_spn2, f00980a_mss_, f00980a_g10_, f00980a_g20_, f00980k_rho2, f00980k_phi2, f00980a_typ_);
    pdfkaka.addResonance   ("f01370k_1", "f01370k_1", f01370a_spn1, f01370a_mss_, f01370a_wdt_, f01370k_rho1, f01370k_phi1, f01370a_typ_);
    pdfkaka.addResonance   ("f01370k_2", "f01370k_2", f01370a_spn2, f01370a_mss_, f01370a_wdt_, f01370k_rho2, f01370k_phi2, f01370a_typ_);
    pdfkaka.addResonance   ("f01750k_1", "f01750k_1", f01750a_spn1, f01750a_mss_, f01750a_wdt_, f01750k_rho1, f01750k_phi1, f01750a_typ_);
    pdfkaka.addResonance   ("f01750k_2", "f01750k_2", f01750a_spn2, f01750a_mss_, f01750a_wdt_, f01750k_rho2, f01750k_phi2, f01750a_typ_);
    pdfkaka.addResonance   ("f01500k_1", "f01500k_1", f01500a_spn1, f01500a_mss_, f01500a_wdt_, f01500k_rho1, f01500k_phi1, f01500a_typ_);
    pdfkaka.addResonance   ("f01500k_2", "f01500k_2", f01500a_spn2, f01500a_mss_, f01500a_wdt_, f01500k_rho2, f01500k_phi2, f01500a_typ_);
    pdfkaka.addResonance   ("f01710k_1", "f01710k_1", f01710a_spn1, f01710a_mss_, f01710a_wdt_, f01710k_rho1, f01710k_phi1, f01710a_typ_);
    pdfkaka.addResonance   ("f01710k_2", "f01710k_2", f01710a_spn2, f01710a_mss_, f01710a_wdt_, f01710k_rho2, f01710k_phi2, f01710a_typ_);
    pdfkaka.addResonance   ("f01270k_1", "f01270k_1", f01270a_spn1, f01270a_mss_, f01270a_wdt_, f01270k_rho1, f01270k_phi1, f01270a_typ_);
    pdfkaka.addResonance   ("f01270k_2", "f01270k_2", f01270a_spn2, f01270a_mss_, f01270a_wdt_, f01270k_rho2, f01270k_phi2, f01270a_typ_);
    pdfkaka.addResonance   ("f01270k_3", "f01270k_3", f01270a_spn3, f01270a_mss_, f01270a_wdt_, f01270k_rho3, f01270k_phi3, f01270a_typ_);
    pdfkaka.addResonance   ("f01270k_4", "f01270k_4", f01270a_spn4, f01270a_mss_, f01270a_wdt_, f01270k_rho4, f01270k_phi4, f01270a_typ_);
    pdfkaka.addResonance   ("f01270k_5", "f01270k_5", f01270a_spn5, f01270a_mss_, f01270a_wdt_, f01270k_rho5, f01270k_phi5, f01270a_typ_);
    pdfkaka.addResonance   ("f01525k_1", "f01525k_1", f01525a_spn1, f01525a_mss_, f01525a_wdt_, f01525k_rho1, f01525k_phi1, f01525a_typ_);
    pdfkaka.addResonance   ("f01525k_2", "f01525k_2", f01525a_spn2, f01525a_mss_, f01525a_wdt_, f01525k_rho2, f01525k_phi2, f01525a_typ_);
    pdfkaka.addResonance   ("f01525k_3", "f01525k_3", f01525a_spn3, f01525a_mss_, f01525a_wdt_, f01525k_rho3, f01525k_phi3, f01525a_typ_);
    pdfkaka.addResonance   ("f01525k_4", "f01525k_4", f01525a_spn4, f01525a_mss_, f01525a_wdt_, f01525k_rho4, f01525k_phi4, f01525a_typ_);
    pdfkaka.addResonance   ("f01525k_5", "f01525k_5", f01525a_spn5, f01525a_mss_, f01525a_wdt_, f01525k_rho5, f01525k_phi5, f01525a_typ_);
   pdfkaka.addResonance1680("ph1680k_1", "ph1680k_1", ph1680a_spn1, ph1680a_mss_, ph1680a_m98_, ph1680a_wdt_, ph1680a_g1__, ph1680a_g2__, ph1680k_rho1, ph1680k_phi1, ph1680a_typ_);
   pdfkaka.addResonance1680("ph1680k_2", "ph1680k_2", ph1680a_spn2, ph1680a_mss_, ph1680a_m98_, ph1680a_wdt_, ph1680a_g1__, ph1680a_g2__, ph1680k_rho2, ph1680k_phi2, ph1680a_typ_);

    //    pdf.addResonance("R1","R1",22,m1,w1,rho1,phi1,1);
    //    pdf.addResonance980("R21","R2",01,m3,g1,g2,rho3,phi4,2);
    //    pdf.addResonance980("R22","R2",02,m3,g1,g2,rho4,phi4,2);
    //    pdf.addResonance600("R31","R3",01,m4,b1,b2,b3,b4,b5,rho4,phi5,3);
    //    pdf.addResonance600("R32","R3",02,m4,b1,b2,b3,b4,b5,rho5,phi5,3);
    //  pdf.addResonance21("R21","R2",m2,g2,rho21,phi21);
    //  pdf.addResonance22("R22","R2",rho22,phi21);
    //  pdf.addResonance23("R23","R2",rho23,phi21);

    //  bkg.addResonance0("R0","R0",m0,g0,rho0,phi0);
    //  bkg.addResonance0("R1","R1",m1,g1,rho1,phi1);
    //  bkg.addResonance21("R21","R2",m2,g2,rho21,phi21);
    //  bkg.addResonance22("R22","R2",rho22,phi21);
    //  bkg.addResonance23("R23","R2",rho23,phi21);

    //
    //  RooCategory d("d","d");
    //  d.defineType("data");
    //  d.defineType("bkg");
    //  d.setLabel("data");data->addColumn(d);
    //  d.setLabel("bkg");data2->addColumn(d);
    //  RooDataSet *dataTotal=data->Clone("dataTotal");
    //  dataTotal->append(*data2);
    //  RooDataSet dataTotal("dataTotal","total",theSet,Index(d),Import("data",*data),Import("bkg",*data2));
    //  RooSimultaneous sim_model("sim_model","",d);
    //  sim_model.addPdf(pdf,"data");
    //  sim_model.addPdf(bkg,"bkg");

    /*
       TStopwatch timer;
       timer.Start();
       cout<<"hello2"<<endl;
    //  RooNLLVar  nll("nll","nll",sim_model,*dataTotal) ;
    RooNLLVar  nllkaka("nllkaka","nllkaka",pdfkaka,*datakaka) ;
    cout<<"hello3"<<endl;
    RooMinuit m(nllkaka) ;
    cout<<"hello4"<<endl;
    */

    /*
       RooPlot* frame1 = g1.frame(Range(0.1,0.2),Title("-log(L) scan vs width")) ;
       RooPlot* frame2 = m1.frame(Range(1.6,1.8),Title("-log(L) scan vs mass")) ;
    //  nll.plotOn(frame1,PrintEvalErrors(0),ShiftToZero(),LineColor(kRed),Precision(1e-1)) ;
    nll.plotOn(frame1,LineColor(kRed),Precision(1e-1)) ;
    nll.plotOn(frame2,LineColor(kRed),Precision(1e-2)) ;
    TCanvas* c = new TCanvas("nll","nllshow",1000,400) ;
    c->Divide(2) ;
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
    c->Print("nll.eps");
    */

    /*	  
    // Activate verbose logging of MINUIT parameter space stepping
    m.setVerbose(kTRUE) ;
    cout<<"hello5"<<endl;

    // Call MIGRAD to minimize the likelihood
    m.migrad() ;
    cout<<"hello6"<<endl;
    // Print values of all paramaters, that reflect values (and error estimates)
    // that are back propagated from MINUIT
    //  pdf.getParameters(x)->Print("s") ;

    // Disable verbose logging
    m.setVerbose(kFALSE) ;

    // Run HESSE to calculate errors from d2L/dp2
    m.hesse() ;

    RooFitResult* res = m.save();
    cout<<"hello:"<<endl;
    //  RooFitResult* res = pdf.fitTo(*data,"mHrt");
    timer.Stop();
    cout<<"Fit time:"<<timer.CpuTime()<<endl;
    // Print the fit result snapshot
    res->Print("v") ;
    cout<<"hello2:"<<endl;
    */

    cout<<"Resonances are added!!"<<endl;
    double startTime_Total=omp_get_wtime();

    RooCategory kapi("kapi","kapi");
    kapi.defineType("phipipi");
    kapi.defineType("phikaka");

    RooDataSet combdata("combdata","combined data",theSet, Index(kapi),Import("phipipi",*datapipi),Import("phikaka",*datakaka));

    RooSimultaneous simPdf("simPdf","simultaneous pdf",kapi);
    simPdf.addPdf(pdfpipi,"phipipi");
    simPdf.addPdf(pdfkaka,"phikaka");

//   RooFitResult* res = simPdf.fitTo(combdata,Save(),RooFit::NumCPU(8));
    RooFitResult* res = simPdf.fitTo(combdata,Save());

    double endTime_Total = omp_get_wtime();   
    cout<<"Fit time : "<< endTime_Total - startTime_Total <<"\n";
    res->Print("V");

    cout << "The minNll of res is " << res->minNll() << endl;
    cout << "The number of floating parameters is " << (res->floatParsFinal()).getSize() << endl;

    pdfpipi.projectpipi();
    pdfkaka.projectkk();
    fillpipi(res, pdfpipi, pdfkaka, projectfilepipi, projectfilekaka);
    fillkk(res, pdfpipi, pdfkaka, projectfilepipi, projectfilekaka);

    cout<<"123456789000000000000"<<endl;
#ifdef PSIP

#endif
    return 0;
}


void fillkk(RooFitResult *res, DPFPWAPdf &pdfpipi, DPFPWAPdf &pdfkaka, TString projectfilepipi, TString projectfilekaka) {
    TString outfit="outkk.rep";
    RooArgList cpar = res->constPars();
    RooArgList fpar = res->floatParsFinal();

    pdfpipi.project(cpar,fpar,projectfilepipi);
    pdfkaka.project(cpar,fpar,projectfilekaka);

    TIterator* parIter =fpar.createIterator() ;
    RooRealVar *par(0);
    parIter->Reset();
    ofstream paramout("parakaka.rep");
    while(0 != (par= (RooRealVar*)parIter->Next())) {
        paramout<<par->getVal()<<"  "<<par->getError()<<endl;
    }
    paramout.close();

    //
    ofstream fout(outfit);
    res->printStream(fout,4,RooPrintable::StyleOption(3),"ostream");
    double n1,n2,n3,n4,n5,n6,n7,n8,n9;
    RooArgSet *s1 = pdfkaka.fitFractions(fpar,kTRUE,fout);
    TIterator *i1 = s1->createIterator();
    RooRealVar *v1 = (RooRealVar*)i1->Next();
    n1=v1->getVal();
    cout<<"n1_ini="<<n1<<endl;
    RooRealVar *v2 = (RooRealVar*)i1->Next();
    n2=v2->getVal();
    cout<<"n2_ini="<<n2<<endl;
    RooRealVar *v3 = (RooRealVar*)i1->Next();
    n3=v3->getVal();
    cout<<"n3_ini="<<n3<<endl;
    RooRealVar *v4 = (RooRealVar*)i1->Next();
    n4=v4->getVal();
    cout<<"n4_ini="<<n4<<endl;
    RooRealVar *v5 = (RooRealVar*)i1->Next();
    n5=v5->getVal();
    cout<<"n5_ini="<<n5<<endl;
    RooRealVar *v6 = (RooRealVar*)i1->Next();
    n6=v6->getVal();
    cout<<"n6_ini="<<n6<<endl;
    RooRealVar *v7 = (RooRealVar*)i1->Next();
    n7=v7->getVal();
    cout<<"n7_ini="<<n7<<endl;
    RooRealVar *v8 = (RooRealVar*)i1->Next();
    n8=v8->getVal();
    cout<<"n8_ini="<<n8<<endl;
    RooRealVar *v9 = (RooRealVar*)i1->Next();
    n9=v9->getVal();
    cout<<"n9_ini="<<n9<<endl;

    TFile *fo=new TFile("fractionkaka.root","recreate");
    TTree *tree=new TTree("fraction","fraction");
    tree->Branch("NumR1",&n1,"n1/D");
    tree->Branch("NumR2",&n2,"n2/D");
    tree->Branch("NumR3",&n3,"n3/D");
    tree->Branch("NumR4",&n4,"n4/D");
    tree->Branch("NumR5",&n5,"n5/D");
    tree->Branch("NumR6",&n6,"n6/D");
    tree->Branch("NumR7",&n7,"n7/D");
    tree->Branch("NumR8",&n8,"n8/D");
    tree->Branch("NumR9",&n9,"n9/D");
    for (int itoy=0; itoy<300; itoy++) {
        RooArgList check = res->randomizePars();
        RooArgSet *mf = pdfkaka.fitFractions(check,kTRUE,fout);
        TIterator *iF = mf->createIterator();
        RooRealVar *fb1 = (RooRealVar*)iF->Next();
        n1=fb1->getVal();
        RooRealVar *fb2 = (RooRealVar*)iF->Next();
        n2=fb2->getVal();
        RooRealVar *fb3 = (RooRealVar*)iF->Next();
        n3=fb3->getVal();
        RooRealVar *fb4 = (RooRealVar*)iF->Next();
        n4=fb4->getVal();
        RooRealVar *fb5 = (RooRealVar*)iF->Next();
        n5=fb5->getVal();
        RooRealVar *fb6 = (RooRealVar*)iF->Next();
        n6=fb6->getVal();
        RooRealVar *fb7 = (RooRealVar*)iF->Next();
        n7=fb7->getVal();
        RooRealVar *fb8 = (RooRealVar*)iF->Next();
        n8=fb8->getVal();
        RooRealVar *fb9 = (RooRealVar*)iF->Next();
        n9=fb9->getVal();
        tree->Fill(); 
    }
    tree->Write();
    fo->Close();

    fout.close();

    cout << "Done fillkk! " << endl;

}

void fillpipi(RooFitResult *res, DPFPWAPdf &pdfpipi, DPFPWAPdf &pdfkaka, TString projectfilepipi, TString projectfilekaka) {
    TString outfitpipi="outpipi.rep";
    RooArgList cpar = res->constPars();
    RooArgList fpar = res->floatParsFinal();

    pdfpipi.project(cpar,fpar,projectfilepipi);
    pdfkaka.project(cpar,fpar,projectfilekaka);

    TIterator* parIterpipi =fpar.createIterator() ;
    RooRealVar *parpipi(0);
    parIterpipi->Reset();
    ofstream paramoutpipi("para1.rep");
    while(0 != (parpipi= (RooRealVar*)parIterpipi->Next())) {
        paramoutpipi<<parpipi->getVal()<<"  "<<parpipi->getError()<<endl;
    }
    paramoutpipi.close();

    //
    ofstream foutpipi(outfitpipi);
    res->printStream(foutpipi,4,RooPrintable::StyleOption(3),"ostream");
    double n1,n2,n3,n4,n5,n6;
    TFile *fopipi=new TFile("fractionpipi.root","recreate");
    TTree *treepipi=new TTree("fraction","fraction");
    treepipi->Branch("NumR1",&n1,"n1/D");
    treepipi->Branch("NumR2",&n2,"n2/D");
    treepipi->Branch("NumR3",&n3,"n3/D");
    treepipi->Branch("NumR4",&n4,"n4/D");
    treepipi->Branch("NumR5",&n5,"n5/D");
    treepipi->Branch("NumR6",&n6,"n6/D");
    //tree->Branch("NumR7",&n7,"n7/D");
    //tree->Branch("NumR8",&n8,"n8/D");
    //tree->Branch("NumR9",&n9,"n9/D");
    //tree->Branch("NumR10",&n10,"n10/D");
    //tree->Branch("NumR11",&n11,"n11/D");
    //tree->Branch("NumR12",&n12,"n12/D");
    //tree->Branch("NumR13",&n13,"n13/D");
    //tree->Branch("NumR14",&n14,"n14/D");
    RooArgSet *s1 = pdfpipi.fitFractions(fpar,kTRUE,foutpipi);
    cout<<"-----------------"<<endl;
    fpar.Print();
    cout<<"-----------------"<<endl;
    TIterator *i1 = s1->createIterator();
    RooRealVar *v1 = (RooRealVar*)i1->Next();
    n1=v1->getVal();
    cout<<"n1_ini="<<n1<<endl;
    RooRealVar *v2 = (RooRealVar*)i1->Next();
    n2=v2->getVal();
    cout<<"n2_ini="<<n2<<endl;
    RooRealVar *v3 = (RooRealVar*)i1->Next();
    n3=v3->getVal();
    cout<<"n3_ini="<<n3<<endl;
    RooRealVar *v4 = (RooRealVar*)i1->Next();
    n4=v4->getVal();
    cout<<"n4_ini="<<n4<<endl;
    RooRealVar *v5 = (RooRealVar*)i1->Next();
    n5=v5->getVal();
    cout<<"n5_ini="<<n5<<endl;
    RooRealVar *v6 = (RooRealVar*)i1->Next();
    n6=v6->getVal();
    cout<<"n6_ini="<<n6<<endl;
    //RooRealVar *v7 = (RooRealVar*)i1->Next();
    //double n7=v7->getVal();
    //cout<<"n7_ini="<<n7<<endl;
    //RooRealVar *v8 = (RooRealVar*)i1->Next();
    //double n8=v8->getVal();
    //cout<<"n8_ini="<<n8<<endl;
    //RooRealVar *v9 = (RooRealVar*)i1->Next();
    //double n9=v9->getVal();
    //cout<<"n9_ini="<<n9<<endl;
    //RooRealVar *v10 = (RooRealVar*)i1->Next();
    //double n10=v10->getVal();
    //cout<<"n10_ini="<<n10<<endl;
    //RooRealVar *v11 = (RooRealVar*)i1->Next();
    //double n11=v11->getVal();
    //cout<<"n11_ini="<<n11<<endl;
    //RooRealVar *v12 = (RooRealVar*)i1->Next();
    //double n12=v12->getVal();
    //cout<<"n12_ini="<<n12<<endl;
    //RooRealVar *v13 = (RooRealVar*)i1->Next();
    //double n13=v13->getVal();
    //cout<<"n13_ini="<<n13<<endl;
    //RooRealVar *v14 = (RooRealVar*)i1->Next();
    //double n14=v14->getVal();
    //cout<<"n14_ini="<<n14<<endl;
    for (int itoy=0; itoy<300; itoy++) {
        RooArgList check = res->randomizePars();
        cout<<"-----------------"<<endl;
        check.Print();
        cout<<"-----------------"<<endl;
        RooArgSet *mf = pdfpipi.fitFractions(check,kTRUE,foutpipi);
        TIterator *iF = mf->createIterator();
        RooRealVar *fb1 = (RooRealVar*)iF->Next();
        n1=fb1->getVal();
        RooRealVar *fb2 = (RooRealVar*)iF->Next();
        n2=fb2->getVal();
        RooRealVar *fb3 = (RooRealVar*)iF->Next();
        n3=fb3->getVal();
        RooRealVar *fb4 = (RooRealVar*)iF->Next();
        n4=fb4->getVal();
        RooRealVar *fb5 = (RooRealVar*)iF->Next();
        n5=fb5->getVal();
        RooRealVar *fb6 = (RooRealVar*)iF->Next();
        n6=fb6->getVal();
        treepipi->Fill(); 
    }
    treepipi->Write();
    fopipi->Close();
    foutpipi.close();
    cout << "Done fillpipi! " << endl;
}
