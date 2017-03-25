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
    Double_t mpsip=3.686,mka=0.493677,mpi=0.13957;
    Double_t pi=3.14159265358979312;
    // double high=mjpsi/2;
    Double_t high=mpsip;
    Double_t low=0-high;
    Double_t Nreal=50000;

    TString indata,bkgdata;
    TString outfitpipi="outpipi.rep";
    TString outfit="outkk.rep";
    TString indatapipi="../phipipi/data_pwa_pipi.dat";
    TString indatakaka="../phikk/data_pwa_kk.dat";
    TString projectfilepipi="resultpipi.root";
    TString projectfilekaka="resultkaka.root";
    bkgdata="kstarkp.dat";

    DPFPWAPoint dphipipi(mka,mka, mpi, mpi, mpsip, "../phipipi/phsp_pwa_pipi.dat");
    DPFPWAPoint dphikaka(mka,mka, mka, mka, mpsip, "../phikk/phsp_pwa_kk.dat");
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
    ////  f0(980) 0+ 1 pp
    RooRealVar spin0("spin0","spin0", 1);
    RooRealVar m980("m980","m980", 0.9650, 0.900, 1.00);
    RooRealVar g10("g10","g10", 0.1650,0.100,0.500);
    RooRealVar g20("g20","g20", 0.6950,0.100,1.500);
    RooRealVar rho980pp1("rho980pp1","rho980pp1", 25.4288,-100.,100.0);
    RooRealVar phi980pp1("phi980pp1","phi980pp1", 4.6788,-2*pi,2*pi);
    RooRealVar propType0("propType0","propType0", 2);
    spin0.setConstant();
    //  m980.setConstant();
    //  g10.setConstant();
    //  g20.setConstant();
    //  rho980pp1.setConstant();
    //  phi980pp1.setConstant();
    propType0.setConstant();
    ////  f0(980) 0+ 2 pp
    RooRealVar spin1("spin1","spin1", 2);
    //  RooRealVar m1("m1","m1", 0.9650, 0.600, 1.9700);
    //  RooRealVar g11("g11","g11", 0.1650,0.100,0.700);
    //  RooRealVar g21("g21","g21", 0.6950,0.600,0.8000);
    RooRealVar rho980pp2("rho980pp2","rho980pp2", 8.0024,-100.,100.0);
    RooRealVar phi980pp2("phi980pp2","phi980pp2", 4.3268,-2*pi,2*pi);
    RooRealVar propType1("propType1","propType1", 2);
    spin1.setConstant();
    //  m1.setConstant();
    //  g11.setConstant();
    //  g21.setConstant();
    //  rho980pp2.setConstant();
    //  phi980pp2.setConstant();
    propType1.setConstant();
    ////****************  f0(1370) 0+ 1 pp
    RooRealVar spin2("spin2","spin2", 1);
    RooRealVar m1370("m1370","m1370", 1.350, 1.200, 1.500);
    RooRealVar w1370("w1370","w1370", 0.265,0.2,0.5);
    RooRealVar rho1370pp1("rho1370pp1","rho1370pp1", 65.8879,-100.0,100.0);
    RooRealVar phi1370pp1("phi1370pp1","phi1370pp1", 0.4468,-pi,pi);
    RooRealVar propType2("propType2","propType2", 1);
    m1370.setConstant();
    w1370.setConstant();
    spin2.setConstant();
    rho1370pp1.setConstant();
    phi1370pp1.setConstant();
    propType2.setConstant();
    ////**************  f0(1370) 0+ 2 pp
    RooRealVar spin3("spin3","spin3", 2);
    //  RooRealVar m3("m3","m3", 1.350, 1.300, 1.400);
    //  RooRealVar w3("w3","w3", 0.265,0.01,0.9);
    RooRealVar rho1370pp2("rho1370pp2","rho1370pp2", 16.1689,-100.0,100.0);
    RooRealVar phi1370pp2("phi1370pp2","phi1370pp2", 0.7596,-pi,pi);
    RooRealVar propType3("propType3","propType3", 1);
    //  m3.setConstant();
    //  w3.setConstant();
    spin3.setConstant();
    rho1370pp2.setConstant();
    phi1370pp2.setConstant();
    propType3.setConstant();
    ////  f0(1750) 0+ 1 pp
    RooRealVar spin4("spin4","spin4", 1);
    RooRealVar m1750("m1750","m1750", 1.7900, 1.70, 1.90);
    RooRealVar w1750("w1750","w1750", 0.2700,0.150,0.350);
    RooRealVar rho1750pp1("rho1750pp1","rho1750pp1", 69.9798,-100.0,100.0);
    RooRealVar phi1750pp1("phi1750pp1","phi1750pp1", 0.5430,-pi,pi);
    RooRealVar propType4("propType4","propType4", 1);
    //  m1750.setConstant();
    //  w1750.setConstant();
    spin4.setConstant();
    //  rho1750pp1.setConstant();
    //  phi1750pp1.setConstant();
    propType4.setConstant();
    ////  f0(1750) 0+ 2 pp
    RooRealVar spin5("spin5","spin5", 2);
    //  RooRealVar m5("m5","m5", 1.7900, 1.70, 1.90);
    //  RooRealVar w5("w5","w5", 0.2700,0.200,0.30);
    RooRealVar rho1750pp2("rho1750pp2","rho1750pp2", 39.7695,-100.0,100.0);
    RooRealVar phi1750pp2("phi1750pp2","phi1750pp2", -0.3770,-pi,pi);
    RooRealVar propType5("propType5","propType5", 1);
    //  m5.setConstant();
    //  w5.setConstant();
    spin5.setConstant();
    //  rho1750pp2.setConstant();
    //  phi1750pp2.setConstant();
    propType5.setConstant();
    ////  f0(sigma600) 0+ 1 pp
    RooRealVar spin6("spin6","spin6", 1);
    RooRealVar m6("m6","m6", 0.9264, 0.800, 1.00);
    RooRealVar b16("b16","b16", 0.5843,0.500,0.600);
    RooRealVar b26("b26","b26", 1.6663,1.500,1.7000);
    RooRealVar b36("b36","b36", 1.0820,0.800,1.500);
    RooRealVar b46("b46","b46", 0.0024,0.0010,0.0050);
    RooRealVar b56("b56","b56", 1.0000,0.5000,1.500);
    RooRealVar rho6("rho6","rho6", 34.8889,-100.0,100.0);
    RooRealVar phi6("phi6","phi6", 1.8491,-pi,pi);
    RooRealVar propType6("propType6","propType6", 3);
    m6.setConstant();
    b16.setConstant();
    b26.setConstant();
    b36.setConstant();
    b46.setConstant();
    b56.setConstant();
    spin6.setConstant();
    rho6.setConstant();
    phi6.setConstant();
    propType6.setConstant();
    ////  f0(sigma600) 0+ 2 pp
    RooRealVar spin7("spin7","spin7", 2);
    RooRealVar m7("m7","m7", 0.9264, 0.800, 1.500);
    RooRealVar b17("b17","b17", 0.5843,0.200,0.900);
    RooRealVar b27("b27","b27", 1.6663,1.00,1.9000);
    RooRealVar b37("b37","b37", 1.0820,0.800,1.500);
    RooRealVar b47("b47","b47", 0.0024,0.0010,0.0050);
    RooRealVar b57("b57","b57", 1.0000,0.5000,1.500);
    RooRealVar rho7("rho7","rho7", 13.1566,-100.0,100.0);
    RooRealVar phi7("phi7","phi7", 1.9269,-pi,pi);
    RooRealVar propType7("propType7","propType7", 3);
    m7.setConstant();
    b17.setConstant();
    b27.setConstant();
    b37.setConstant();
    b47.setConstant();
    b57.setConstant();
    spin7.setConstant();
    rho7.setConstant();
    phi7.setConstant();
    propType7.setConstant();
    ////  f0(1500) 0+ 1 pp
    RooRealVar spin8("spin8","spin8", 1);
    RooRealVar m8("m8","m8", 1.4950, 1.00, 1.90);
    RooRealVar w8("w8","w8", 0.1220,0.100,0.150);
    RooRealVar rho8("rho8","rho8", 53.6250,-100.0,100.0);
    RooRealVar phi8("phi8","phi8", 2.7549,-pi,pi);
    RooRealVar propType8("propType8","propType8", 1);
    m8.setConstant();
    w8.setConstant();
    spin8.setConstant();
    rho8.setConstant();
    phi8.setConstant();
    propType8.setConstant();
    ////  f0(1500) 0+ 2 pp
    RooRealVar spin9("spin9","spin9", 2);
    RooRealVar m9("m9","m9", 1.4950, 1.00, 1.90);
    RooRealVar w9("w9","w9", 0.1220,0.100,0.50);
    RooRealVar rho9("rho9","rho9", 23.4287,-100.0,100.0);
    RooRealVar phi9("phi9","phi9", 2.5820,-pi,pi);
    RooRealVar propType9("propType9","propType9", 1);
    m9.setConstant();
    w9.setConstant();
    spin9.setConstant();
    rho9.setConstant();
    phi9.setConstant();
    propType9.setConstant();
    ////  f0(1270) 2+ 1 pp
    RooRealVar spin10("spin10","spin10", 21);
    RooRealVar m10("m10","m10", 1.2750, 1.0, 1.800);
    RooRealVar w10("w10","w10", 0.1900,0.10,0.90);
    RooRealVar rho10("rho10","rho10", 23.4619,-100.0,100.0);
    RooRealVar phi10("phi10","phi10", 0.6866,-pi,pi);
    RooRealVar propType10("propType10","propType10", 6);
    m10.setConstant();
    w10.setConstant();
    spin10.setConstant();
    rho10.setConstant();
    phi10.setConstant();
    propType10.setConstant();
    ////  f0(1270) 2+ 2 pp
    RooRealVar spin11("spin11","spin11", 22);
    RooRealVar m11("m11","m11", 1.2750, 1.0, 1.800);
    RooRealVar w11("w11","w11", 0.1900,0.10,0.950);
    RooRealVar rho11("rho11","rho11", 7.5520,-100.,100.0);
    RooRealVar phi11("phi11","phi11", -0.1037,-pi,pi);
    RooRealVar propType11("propType11","propType11", 6);
    m11.setConstant();
    w11.setConstant();
    spin11.setConstant();
    rho11.setConstant();
    phi11.setConstant();
    propType11.setConstant();
    ////  f0(1270) 2+ 3 pp
    RooRealVar spin12("spin12","spin12", 23);
    RooRealVar m12("m12","m12", 1.2750, 1.00, 1.800);
    RooRealVar w12("w12","w12", 0.1900,0.10,0.950);
    RooRealVar rho12("rho12","rho12", 5.3346,-100.,100.0);
    RooRealVar phi12("phi12","phi12", -0.0750,-pi,pi);
    RooRealVar propType12("propType12","propType12", 6);
    m12.setConstant();
    w12.setConstant();
    spin12.setConstant();
    rho12.setConstant();
    phi12.setConstant();
    propType12.setConstant();
    ////  f0(1270) 2+ 4 pp
    RooRealVar spin13("spin13","spin13", 24);
    RooRealVar m13("m13","m13", 1.2750, 1.00, 1.800);
    RooRealVar w13("w13","w13", 0.1900,0.10,0.950);
    RooRealVar rho13("rho13","rho13", 4.3929,-100.,100.0);
    RooRealVar phi13("phi13","phi13", 3.2811,-2*pi,2*pi);
    RooRealVar propType13("propType13","propType13", 6);
    m13.setConstant();
    w13.setConstant();
    spin13.setConstant();
    rho13.setConstant();
    phi13.setConstant();
    propType13.setConstant();
    ////  f0(1270) 2+ 5 pp
    RooRealVar spin14("spin14","spin14", 25);
    RooRealVar m14("m14","m14", 1.2750, 1.00, 1.800);
    RooRealVar w14("w14","w14", 0.1900,0.10,0.950);
    RooRealVar rho14("rho14","rho14", 4.3929,-100.,100.0);
    RooRealVar phi14("phi14","phi14", 3.2811,-2*pi,2*pi);
    RooRealVar propType14("propType14","propType14", 6);
    m14.setConstant();
    w14.setConstant();
    spin14.setConstant();
    rho14.setConstant();
    phi14.setConstant();
    propType14.setConstant();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////  f0(980) 0+ 1 kk 
    RooRealVar spin50("spin50","spin50", 1);
    //  RooRealVar m50("m50","m50", 0.9650, 0.900, 1.200);
    //  RooRealVar g150("g150","g150", 0.1650,0.1,0.5);
    //  RooRealVar g250("g250","g250", 0.6950,0.10,0.900);
    RooRealVar rho980kk1("rho980kk1","rho980kk1", 0.7628,-100.0,100.0);
    RooRealVar phi980kk1("phi980kk1","phi980kk1", 1.2817,-pi,pi);
    RooRealVar propType50("propType50","propType50", 2);
    //  m50.setConstant();
    //  g150.setConstant();
    spin50.setConstant();
    //  g250.setConstant();
    //  rho980kk1.setConstant();
    //  phi980kk1.setConstant();
    propType50.setConstant();
    ////  f0(980) 0+ 2 kk 
    RooRealVar spin51("spin51","spin51", 2);
    //  RooRealVar m51("m51","m51", 0.9650, 0.600, 1.20);
    //  RooRealVar g151("g151","g151", 0.1650,0.10,0.70);
    //  RooRealVar g251("g251","g251", 0.6950,0.100,0.900);
    RooRealVar rho980kk2("rho980kk2","rho980kk2", 0.2400,-100.0,100.0);
    RooRealVar phi980kk2("phi980kk2","phi980kk2", 1.0694,-pi,pi);
    RooRealVar propType51("propType51","propType51", 2);
    //  m51.setConstant();
    //  g151.setConstant();
    spin51.setConstant();
    //  rho980kk2.setConstant();
    //  g251.setConstant();
    //  phi980kk2.setConstant();
    propType51.setConstant();
    ////**************  f0(1370) 0+ 1 kk
    RooRealVar spin52("spin52","spin52", 1);
    //  RooRealVar m52("m52","m52", 1.350, 1.00, 1.500);
    //  RooRealVar w52("w52","w52", 0.265,0.01,0.9);
    RooRealVar rho1370kk1("rho1370kk1","rho1370kk1", 0.3360,-100.0,100.0);
    RooRealVar phi1370kk1("phi1370kk1","phi1370kk1", 3.3365,-2*pi,2*pi);
    RooRealVar propType52("propType52","propType52", 1);
    //  m52.setConstant();
    //  w52.setConstant();
    spin52.setConstant();
    rho1370kk1.setConstant();
    phi1370kk1.setConstant();
    propType52.setConstant();
    ////***************  f0(1370) 0+ 2 kk
    RooRealVar spin53("spin53","spin53", 2);
    //  RooRealVar m53("m53","m53", 1.350, 1.0, 1.50);
    //  RooRealVar w53("w53","w53", 0.265,0.01,0.9);
    RooRealVar rho1370kk2("rho1370kk2","rho1370kk2", 0.08246,-100.0,100.0);
    RooRealVar phi1370kk2("phi1370kk2","phi1370kk2", -2.7921,-pi,pi);
    RooRealVar propType53("propType53","propType53", 1);
    //  m53.setConstant();
    //  w53.setConstant();
    spin53.setConstant();
    rho1370kk2.setConstant();
    phi1370kk2.setConstant();
    propType53.setConstant();
    ////  f0(1750) 0+ 1 kk 
    RooRealVar spin54("spin54","spin54", 1);
    //  RooRealVar m54("m54","m54", 1.7900, 1.0, 2.50);
    //  RooRealVar w54("w54","w54", 0.2700,0.10,0.80);
    RooRealVar rho1750kk1("rho1750kk1","rho1750kk1", 0.5598,-100.0,100.0);
    RooRealVar phi1750kk1("phi1750kk1","phi1750kk1", 5.6406,-2*pi,2*pi);
    RooRealVar propType54("propType54","propType54", 1);
    //  m54.setConstant();
    //  w54.setConstant();
    spin54.setConstant();
    //  rho1750kk1.setConstant();
    //  phi1750kk1.setConstant();
    propType54.setConstant();
    ////  f0(1750) 0+ 2 kk
    RooRealVar spin55("spin55","spin55", 2);
    //  RooRealVar m55("m55","m55", 1.7900, 1.0, 2.50);
    //  RooRealVar w55("w55","w55", 0.2700,0.10,0.80);
    RooRealVar rho1750kk2("rho1750kk2","rho1750kk2", 0.3181,-100.0,100.0);
    RooRealVar phi1750kk2("phi1750kk2","phi1750kk2", 4.2865,-2*pi,2*pi);
    RooRealVar propType55("propType55","propType55", 1);
    //  m55.setConstant();
    //  w55.setConstant();
    spin55.setConstant();
    //  rho1750kk2.setConstant();
    //  phi1750kk2.setConstant();
    propType55.setConstant();
    ////  f0(sigma) 0+ 1 kk
    RooRealVar spin56("spin56","spin56", 1);
    RooRealVar m56("m56","m56", 0.9264, 0.20, 1.500);
    RooRealVar b156("b156","b156", 0.5843,0.10,0.90);
    RooRealVar b256("b256","b256", 1.6663,1.0,2.50);
    RooRealVar b356("b356","b356", 1.0820,0.500,1.500);
    RooRealVar b456("b456","b456", 0.0024,0.0010,0.0080);
    RooRealVar b556("b556","b556", 1.0000,0.500,1.500);
    RooRealVar rho56("rho56","rho56", 0.5477,-100.0,100.0);
    RooRealVar phi56("phi56","phi56", 3.0705,-2*pi,2*pi);
    RooRealVar propType56("propType56","propType56", 3);
    spin56.setConstant();
    m56.setConstant();
    b156.setConstant();
    b256.setConstant();
    b356.setConstant();
    b456.setConstant();
    b556.setConstant();
    rho56.setConstant();
    phi56.setConstant();
    propType56.setConstant();
    ////  f0(sigma) 0+ 2 kk
    RooRealVar spin57("spin57","spin57", 2);
    RooRealVar m57("m57","m57", 0.9264, 0.200, 1.500);
    RooRealVar b157("b157","b157", 0.5843,0.100,1.500);
    RooRealVar b257("b257","b257", 1.6663,0.500,2.7000);
    RooRealVar b357("b357","b357", 1.0820,0.500,1.500);
    RooRealVar b457("b457","b457", 0.0024,0.001,0.0080);
    RooRealVar b557("b557","b557", 1.0000,0.5000,1.500);
    RooRealVar rho57("rho57","rho57", 0.2065,-100.0,100.0);
    RooRealVar phi57("phi57","phi57", 2.1623,-pi,pi);
    RooRealVar propType57("propType57","propType57", 3);
    m57.setConstant();
    b157.setConstant();
    b257.setConstant();
    b357.setConstant();
    b457.setConstant();
    b557.setConstant();
    spin57.setConstant();
    rho57.setConstant();
    phi57.setConstant();
    propType57.setConstant();
    ////  f0(1500) 0+ 1 kk 
    RooRealVar spin58("spin58","spin58", 1);
    RooRealVar m58("m58","m58", 1.4950, 0.500, 2.50);
    RooRealVar w58("w58","w58", 0.1220,0.01,0.50);
    RooRealVar rho58("rho58","rho58", 0.6166,-100.0,100.0);
    RooRealVar phi58("phi58","phi58", 3.9934,-2*pi,2*pi);
    RooRealVar propType58("propType58","propType58", 1);
    m58.setConstant();
    w58.setConstant();
    spin58.setConstant();
    rho58.setConstant();
    phi58.setConstant();
    propType58.setConstant();
    ////  f0(1500) 0+ 2 kk
    RooRealVar spin59("spin59","spin59", 2);
    RooRealVar m59("m59","m59", 1.4950, 0.500, 2.50);
    RooRealVar w59("w59","w59", 0.1220,0.01,0.50);
    RooRealVar rho59("rho59","rho59", 0.2694,-100.0,100.0);
    RooRealVar phi59("phi59","phi59", 4.3124,-2*pi,2*pi);
    RooRealVar propType59("propType59","propType59", 1);
    m59.setConstant();
    w59.setConstant();
    spin59.setConstant();
    rho59.setConstant();
    phi59.setConstant();
    propType59.setConstant();
    ////  f0(1710) 0+ 1 kk 
    RooRealVar spin60("spin60","spin60", 1);
    RooRealVar m1710("m1710","m1710", 1.7070, 1.650, 1.800);
    RooRealVar w1710("w1710","w1710", 0.1250,0.10,0.1500);
    RooRealVar rho1710kk1("rho1710kk1","rho1710kk1", 1.0527,-100.0,100.0);
    RooRealVar phi1710kk1("phi1710kk1","phi1710kk1", 6.8807,-3*pi,3*pi);
    RooRealVar propType60("propType60","propType60", 1);
    //  m1710.setConstant();
    //  w1710.setConstant();
    spin60.setConstant();
    //  rho1710kk1.setConstant();
    //  phi1710kk1.setConstant();
    propType60.setConstant();
    ////  f0(1710) 0+ 2 kk 
    RooRealVar spin61("spin61","spin61", 2);
    //  RooRealVar m61("m61","m61", 1.7070, 0.500, 2.700);
    //  RooRealVar w61("w61","w61", 0.1250,0.0100,0.500);
    RooRealVar rho1710kk2("rho1710kk2","rho1710kk2", 0.2545,-100.0,100.0);
    RooRealVar phi1710kk2("phi1710kk2","phi1710kk2", -0.0052,-pi,pi);
    RooRealVar propType61("propType61","propType61", 1);
    //  m61.setConstant();
    //  w61.setConstant();
    spin61.setConstant();
    //  rho1710kk2.setConstant();
    //  phi1710kk2.setConstant();
    propType61.setConstant();
    ////  f0(1270) 2+ 1 kk
    RooRealVar spin62("spin62","spin62", 21);
    RooRealVar m62("m62","m62", 1.2750, 0.500, 1.800);
    RooRealVar w62("w62","w62", 0.1900,0.010,0.90);
    RooRealVar rho62("rho62","rho62", 0.1000,-100.0,100.0);
    RooRealVar phi62("phi62","phi62", 2.8428,-2*pi,2*pi);
    RooRealVar propType62("propType62","propType62", 6);
    m62.setConstant();
    w62.setConstant();
    spin62.setConstant();
    rho62.setConstant();
    phi62.setConstant();
    propType62.setConstant();
    ////  f0(1270) 2+ 2 kk
    RooRealVar spin63("spin63","spin63", 22);
    RooRealVar m63("m63","m63", 1.2750, 0.500, 2.500);
    RooRealVar w63("w63","w63", 0.1900,0.010,0.90);
    RooRealVar rho63("rho63","rho63", 0.0322,-100.0,100.0);
    RooRealVar phi63("phi63","phi63", -1.9102,-pi,pi);
    RooRealVar propType63("propType63","propType63", 6);
    m63.setConstant();
    w63.setConstant();
    spin63.setConstant();
    rho63.setConstant();
    phi63.setConstant();
    propType63.setConstant();
    ////  f0(1270) 2+ 3 kk
    RooRealVar spin64("spin64","spin64", 23);
    RooRealVar m64("m64","m64", 1.2750, 0.500, 1.800);
    RooRealVar w64("w64","w64", 0.1900,0.010,0.90);
    RooRealVar rho64("rho64","rho64", 0.0226,-100.0,100.0);
    RooRealVar phi64("phi64","phi64", 0.3628,-pi,pi);
    RooRealVar propType64("propType64","propType64", 6);
    m64.setConstant();
    w64.setConstant();
    spin64.setConstant();
    rho64.setConstant();
    phi64.setConstant();
    propType64.setConstant();
    ////  f0(1270) 2+ 4 kk
    RooRealVar spin65("spin65","spin65", 24);
    RooRealVar m65("m65","m65", 1.2750, 0.500, 2.00);
    RooRealVar w65("w65","w65", 0.1900,0.010,0.90);
    RooRealVar rho65("rho65","rho65", 0.0187,-100.0,100.0);
    RooRealVar phi65("phi65","phi65", 1.4555,-pi,pi);
    RooRealVar propType65("propType65","propType65", 6);
    m65.setConstant();
    w65.setConstant();
    spin65.setConstant();
    rho65.setConstant();
    phi65.setConstant();
    propType65.setConstant();
    ////  f0(1270) 2+ 5 kk
    RooRealVar spin66("spin66","spin66", 25);
    RooRealVar m66("m66","m66", 1.2750, 0.500, 1.800);
    RooRealVar w66("w66","w66", 0.1900,0.010,0.90);
    RooRealVar rho66("rho66","rho66", 0.0000,-100.0,100.0);
    RooRealVar phi66("phi66","phi66", 0.0000,-pi,pi);
    RooRealVar propType66("propType66","propType66", 6);
    m66.setConstant();
    w66.setConstant();
    spin66.setConstant();
    rho66.setConstant();
    phi66.setConstant();
    propType66.setConstant();
    ////  f0(1525) 2+ 1 kk
    RooRealVar spin67("spin67","spin67", 21);
    RooRealVar m67("m67","m67", 1.5210, 0.500, 2.50);
    RooRealVar w67("w67","w67", 0.0770,0.0010,0.800);
    RooRealVar rho67("rho67","rho67", 1.0000,-100.0,100.0);
    RooRealVar phi67("phi67","phi67", 0.2000,-pi,pi);
    RooRealVar propType67("propType67","propType67", 1);
    m67.setConstant();
    w67.setConstant();
    spin67.setConstant();
    rho67.setConstant();
    phi67.setConstant();
    propType67.setConstant();
    ////  f0(1525) 2+ 2 kk
    RooRealVar spin68("spin68","spin68", 22);
    RooRealVar m68("m68","m68", 1.5210, 0.500, 2.50);
    RooRealVar w68("w68","w68", 0.0770,0.0010,0.800);
    RooRealVar rho68("rho68","rho68", 0.5878,-100.0,100.0);
    RooRealVar phi68("phi68","phi68", 5.9124,-3*pi,3*pi);
    RooRealVar propType68("propType68","propType68", 1);
    m68.setConstant();
    w68.setConstant();
    spin68.setConstant();
    rho68.setConstant();
    phi68.setConstant();
    propType68.setConstant();
    ////  f0(1525) 2+ 3 kk
    RooRealVar spin69("spin69","spin69", 23);
    RooRealVar m69("m69","m69", 1.5210, 0.500, 2.50);
    RooRealVar w69("w69","w69", 0.0770,0.0010,0.800);
    RooRealVar rho69("rho69","rho69", 0.3493,-100.0,100.0);
    RooRealVar phi69("phi69","phi69", 0.1663,-pi,pi);
    RooRealVar propType69("propType69","propType69", 1);
    m69.setConstant();
    w69.setConstant();
    spin69.setConstant();
    rho69.setConstant();
    phi69.setConstant();
    propType69.setConstant();
    ////  f0(1525) 2+ 4 kk
    RooRealVar spin70("spin70","spin70", 24);
    RooRealVar m70("m70","m70", 1.5210, 0.500, 2.50);
    RooRealVar w70("w70","w70", 0.0770,0.0010,0.800);
    RooRealVar rho70("rho70","rho70", 0.0000,-100.0,100.0);
    RooRealVar phi70("phi70","phi70", 0.0000,-pi,pi);
    RooRealVar propType70("propType70","propType70", 1);
    m70.setConstant();
    w70.setConstant();
    spin70.setConstant();
    rho70.setConstant();
    phi70.setConstant();
    propType70.setConstant();
    ////  f0(1525) 2+ 5 kk
    RooRealVar spin71("spin71","spin71", 25);
    RooRealVar m71("m71","m71", 1.5210, 0.500, 2.50);
    RooRealVar w71("w71","w71", 0.0770,0.0010,0.800);
    RooRealVar rho71("rho71","rho71", 0.0000,-100.0,100.0);
    RooRealVar phi71("phi71","phi71", 0.0000,-pi,pi);
    RooRealVar propType71("propType71","propType71", 1);
    m71.setConstant();
    w71.setConstant();
    spin71.setConstant();
    rho71.setConstant();
    phi71.setConstant();
    propType71.setConstant();
    ////  phi(1680) 1- 1 kk
    RooRealVar spin72("spin72","spin72", 191);
    RooRealVar m72("m72","m72", 1.6150, 0.500, 2.500);
    RooRealVar m72980("m72980","m72980", 0.9650, 0.600, 1.200);
    RooRealVar w72("w72","w72", 0.1300,0.01,0.50);
    RooRealVar g172("g172","g172", 0.1650,0.01,0.500);
    RooRealVar g272("g272","g272", 0.6950,0.1,1.0);
    RooRealVar rho72("rho72","rho72", 0.0192,-100.0,100.0);
    RooRealVar phi72("phi72","phi72", -0.5701,-pi,pi);
    RooRealVar propType72("propType72","propType72", 5);
    m72.setConstant();
    m72980.setConstant();
    w72.setConstant();
    g172.setConstant();
    g272.setConstant();
    spin72.setConstant();
    rho72.setConstant();
    phi72.setConstant();
    propType72.setConstant();
    ////  phi(1680) 1- 2 kk
    RooRealVar spin73("spin73","spin73", 192);
    RooRealVar m73("m73","m73", 1.6150, 0.500, 2.500);
    RooRealVar m73980("m73980","m73980", 0.9650, 0.600, 1.200);
    RooRealVar w73("w73","w73", 0.1300,0.010,0.500);
    RooRealVar g173("g173","g173", 0.1650,0.010,0.500);
    RooRealVar g273("g273","g273", 0.6950,0.1,1.0);
    RooRealVar rho73("rho73","rho73", 0.0297,-100.0,100.0);
    RooRealVar phi73("phi73","phi73", 0.1496,-pi,pi);
    RooRealVar propType73("propType73","propType73", 5);
    m73.setConstant();
    m73980.setConstant();
    w73.setConstant();
    g173.setConstant();
    g273.setConstant();
    spin73.setConstant();
    rho73.setConstant();
    phi73.setConstant();
    propType73.setConstant();
    //////  f0(400) 0+ 1 kk 
    //  RooRealVar spin56("spin56","spin56", 1);
    //  RooRealVar m56("m56","m56", 0.4930, 0.4900, 0.4950);
    //  RooRealVar w56("w56","w56", 0.5380,0.5350,0.5400);
    //  RooRealVar rho56("rho56","rho56", 0.0157,-100.,100.0);
    //  RooRealVar phi56("phi56","phi56", 3.0721,-pi,pi);
    //  RooRealVar propType56("propType56","propType56", 3);
    //  m56.setConstant();
    //  w56.setConstant();
    //  spin56.setConstant();
    //  rho56.setConstant();
    //  phi56.setConstant();
    //  propType56.setConstant();
    //////  f0(400) 0+ 2 kk
    //  RooRealVar spin57("spin57","spin57", 2);
    //  RooRealVar m57("m57","m57", 0.4930, 0.4900, 0.4950);
    //  RooRealVar w57("w57","w57", 0.5380,0.5350,0.5400);
    //  RooRealVar rho57("rho57","rho57", 0.0000,-100.,100.0);
    //  RooRealVar phi57("phi57","phi57", 2.1614,-pi,pi);
    //  RooRealVar propType57("propType57","propType57", 3);
    //  m57.setConstant();
    //  w57.setConstant();
    //  spin57.setConstant();
    //  rho57.setConstant();
    //  phi57.setConstant();
    //  propType57.setConstant();

    // RooRealVar m1("m1","m1", 1.79, 1.7, 1.9);
    // RooRealVar g1("g1","g1", 0.27,0.20,0.35);
    // RooRealVar rho1("rho1","rho1", 35.,-100.,100.0);
    // RooRealVar phi1("phi1","phi1", -0.5,-pi,pi);

    //  RooRealVar m2("m2","m2", 1.5,1.4,1.6);
    //  RooRealVar g2("g2","g2", 0.109,0.,0.2);
    //  RooRealVar rho21("rho21","rho21", -16.,-100.,100.);
    //  RooRealVar phi21("phi21","phi21", 0.5,-pi,pi);
    //
    //  RooRealVar rho22("rho22","rho22", -1., -100.,100.);
    //  RooRealVar rho23("rho23","rho22", -3., -100.,100.);
    //  g1.setConstant();
    //  g2.setConstant();
    //  m1.setConstant();
    //  m2.setConstant();
    //  rho21.setConstant();
    //  rho22.setConstant();
    //  rho23.setConstant();
    //  phi21.setConstant();

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
    cout<<"hello1"<<endl;
    DPFPWAPdf pdfpipi("pdfpipi","pdfpipi",v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54,&dphipipi);
    pdfpipi.addResonance980("R9801","R980pipi",spin0,m980,g10,g20,rho980pp1,phi980pp1,propType0);
    //    pdfpipi.addResonance980("R9802","R980pipi",spin1,m1,g11,g21,rho1,phi1,propType1);
    pdfpipi.addResonance980("R9802","R980pipi",spin1,m980,g10,g20,rho980pp2,phi980pp2,propType1);
    pdfpipi.addResonance("R13701","R1370pipi",spin2,m1370,w1370,rho1370pp1,phi1370pp1,propType2);
    //    pdfpipi.addResonance("R13702","R1370pipi",spin3,m3,w3,rho3,phi3,propType3);
    pdfpipi.addResonance("R13702","R1370pipi",spin3,m1370,w1370,rho1370pp2,phi1370pp2,propType3);
    pdfpipi.addResonance("R17501","R1750pipi",spin4,m1750,w1750,rho1750pp1,phi1750pp1,propType4);
    //    pdfpipi.addResonance("R17502","R1750pipi",spin5,m5,w5,rho5,phi5,propType5);
    pdfpipi.addResonance("R17502","R1750pipi",spin5,m1750,w1750,rho1750pp2,phi1750pp2,propType5);
    pdfpipi.addResonance600("sigma01","sigma",spin6, m6, b16, b26, b36, b46, b56, rho6, phi6, propType6);
    pdfpipi.addResonance600("sigma02","sigma",spin7, m7, b17, b27, b37, b47, b57, rho7, phi7, propType7);
    pdfpipi.addResonance("R15001","R1500pipi",spin8,m8,w8,rho8,phi8,propType8);
    pdfpipi.addResonance("R15002","R1500pipi",spin9,m9,w9,rho9,phi9,propType9);
    pdfpipi.addResonance("R12701","R1270pipi",spin10,m10,w10,rho10,phi10,propType10);
    pdfpipi.addResonance("R12702","R1270pipi",spin11,m11,w11,rho11,phi11,propType11);
    pdfpipi.addResonance("R12703","R1270pipi",spin12,m12,w12,rho12,phi12,propType12);
    pdfpipi.addResonance("R12704","R1270pipi",spin13,m13,w13,rho13,phi13,propType13);
    //    pdfpipi.addResonance("R12705","R12705",spin14,m14,w14,rho14,phi14,propType14);
    // pdf.addResonance("R1","R1",22,m1,w1,rho1,phi1,1);
    // pdf.addResonance980("R21","R2",01,m3,g1,g2,rho3,phi4,2);
    // pdf.addResonance980("R22","R2",02,m3,g1,g2,rho4,phi4,2);
    // pdf.addResonance600("R31","R3",01,m4,b1,b2,b3,b4,b5,rho4,phi5,3);
    // pdf.addResonance600("R32","R3",02,m4,b1,b2,b3,b4,b5,rho5,phi5,3);
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
    RooNLLVar  nllpipi("nllpipi","nllpipi",pdfpipi,*datapipi) ;
    cout<<"hello3"<<endl;
    RooMinuit m(nllpipi) ;
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

    RooFitResult* respipi = m.save();
    cout<<"hello:"<<endl;
    //  RooFitResult* res = pdf.fitTo(*data,"mHrt");
    timer.Stop();
    cout<<"Fit time:"<<timer.CpuTime()<<endl;
    // Print the fit result snapshot
    respipi->Print("v") ;
    cout<<"hello2:pipi"<<endl;

    RooArgList cparpipi = respipi->constPars();
    RooArgList fparpipi = respipi->floatParsFinal();
    //  pdfpipi.project(cparpipi,fparpipi,projectfilepipi);
    pdfpipi.project(cparpipi,fparpipi,projectfilepipi);

    TIterator* parIter =fparpipi.createIterator() ;
    RooRealVar *par(0);
    parIter->Reset();
    ofstream paramout("para1.rep");
    while(0 != (par= (RooRealVar*)parIter->Next())) {
    paramout<<par->getVal()<<"  "<<par->getError()<<endl;
    }
    paramout.close();

    //
    ofstream foutpipi(outfitpipi);
    respipi->printStream(foutpipi,4,3,"ostream");
    double n1,n2,n3,n4,n5,n6;
    TFile *fo=new TFile("fractionpipi.root","recreate");
    TTree *tree=new TTree("fraction","fraction");
    tree->Branch("NumR1",&n1,"n1/D");
    tree->Branch("NumR2",&n2,"n2/D");
    tree->Branch("NumR3",&n3,"n3/D");
    tree->Branch("NumR4",&n4,"n4/D");
    tree->Branch("NumR5",&n5,"n5/D");
    tree->Branch("NumR6",&n6,"n6/D");
    //tree->Branch("NumR7",&n7,"n7/D");
    //tree->Branch("NumR8",&n8,"n8/D");
    //tree->Branch("NumR9",&n9,"n9/D");
    //tree->Branch("NumR10",&n10,"n10/D");
    //tree->Branch("NumR11",&n11,"n11/D");
    //tree->Branch("NumR12",&n12,"n12/D");
    //tree->Branch("NumR13",&n13,"n13/D");
    //tree->Branch("NumR14",&n14,"n14/D");
    RooArgSet *s1 = pdfpipi.fitFractions(fparpipi,kTRUE,foutpipi);
    cout<<"-----------------"<<endl;
    fparpipi->Print();
    cout<<"-----------------"<<endl;
    TIterator *i1 = s1->createIterator();
    RooRealVar *v1 = (RooRealVar*)i1->Next();
    double n1=v1->getVal();
    cout<<"n1_ini="<<n1<<endl;
    RooRealVar *v2 = (RooRealVar*)i1->Next();
    double n2=v2->getVal();
    cout<<"n2_ini="<<n2<<endl;
    RooRealVar *v3 = (RooRealVar*)i1->Next();
    double n3=v3->getVal();
    cout<<"n3_ini="<<n3<<endl;
    RooRealVar *v4 = (RooRealVar*)i1->Next();
    double n4=v4->getVal();
    cout<<"n4_ini="<<n4<<endl;
    RooRealVar *v5 = (RooRealVar*)i1->Next();
    double n5=v5->getVal();
    cout<<"n5_ini="<<n5<<endl;
    RooRealVar *v6 = (RooRealVar*)i1->Next();
    double n6=v6->getVal();
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
    for (int itoy=0; itoy<1; itoy++) {
        RooArgList check = respipi->randomizePars();
        cout<<"-----------------"<<endl;
        check->Print();
        cout<<"-----------------"<<endl;
        RooArgSet *mf = pdfpipi.fitFractions(check,kTRUE,foutpipi);
        TIterator *iF = mf->createIterator();
        RooRealVar *fb1 = (RooRealVar*)iF->Next();
        double n1=fb1->getVal();
        RooRealVar *fb2 = (RooRealVar*)iF->Next();
        double n2=fb2->getVal();
        RooRealVar *fb3 = (RooRealVar*)iF->Next();
        double n3=fb3->getVal();
        RooRealVar *fb4 = (RooRealVar*)iF->Next();
        double n4=fb4->getVal();
        RooRealVar *fb5 = (RooRealVar*)iF->Next();
        double n5=fb5->getVal();
        RooRealVar *fb6 = (RooRealVar*)iF->Next();
        double n6=fb6->getVal();
        tree->Fill(); 
    }
    tree->Write();
    fo->Close();
    foutpipi.close();
    */
        /*
           RooRealVar mcv11("mcv11","11",low,high);
           RooRealVar mcv12("mcv12","12",low,high);
           RooRealVar mcv13("mcv13","13",low,high);
           RooRealVar mcv14("mcv14","14",low,high);
           RooRealVar mcv21("mcv21","21",low,high);
           RooRealVar mcv22("mcv22","22",low,high);
           RooRealVar mcv23("mcv23","23",low,high);
           RooRealVar mcv24("mcv24","24",low,high);
           RooArgSet mcSet;
           mcSet.add(RooArgSet(mcv11,mcv12,mcv13,mcv14,mcv21,mcv22,mcv23,mcv24));
           RooDataSet *mcdata = RooDataSet::read("phspdata.dat",theSet);
        //  RooFormulaVar mean("mean","mcv11+mcv22",RooArgList(mcv11,mcv22));
        //  RooPlot *meanframe=mean.frame();
        RooPlot *v11frame=v11.frame();
        mcdata->plotOn(v11frame);
        pdf.plotOn(v11frame,ProjWData(*mcdata));
        v11frame->Draw();
        */

        /*
           RooPlot* frame1 = g1.frame(Range(0.11,0.125),Title("-log(L) scan vs width")) ;
           RooPlot* frame2 = m1.frame(Range(1.715,1.73),Title("-log(L) scan vs mass")) ;
        //  nll.plotOn(frame1,PrintEvalErrors(0),ShiftToZero(),LineColor(kRed),Precision(1e-1)) ;
        nll.plotOn(frame1,ShiftToZero(),LineColor(kRed),Precision(1e-1)) ;
        nll.plotOn(frame2,ShiftToZero(),LineColor(kRed),Precision(1e-2)) ;
        TCanvas* c = new TCanvas("nll","nllshow",1000,400) ;
        c->Divide(2) ;
        c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
        c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
        c->Print("nll.eps");
        */





        cout<<"-----------------"<<endl;
    //  RooDataSet *data = RooDataSet::read(indata,theSet,"Q");
    RooDataSet *data22 = RooDataSet::read(indatakaka,theSet);
    //  RooDataSet *data2 = RooDataSet::read(bkgdata,theSet);
    //  PWAPdf pdf("pdf","pdf", v11, v12, v13, v14, v21, v22, v23, v24);
    //  bkgpdf bkg("bkgpdf","bkgpdf", v11, v12, v13, v14, v21, v22, v23, v24);
    data22->Print();
    RooDataSet *datakaka= new RooDataSet(data22->GetName(),data22->GetTitle(),data22,*data22->get(),0,weight.GetName());
    datakaka->Print();
    cout<<"hello1"<<endl;
    DPFPWAPdf pdfkaka("pdfkaka","pdfkaka",v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54,&dphikaka);
    //    pdfkaka.addResonance980("R9801","R980kk",spin50,m50,g150,g250,rho50,phi50,propType50);
    pdfkaka.addResonance980("R9801","R980kk",spin50,m980,g10,g20,rho980kk1,phi980kk1,propType50);
    //    pdfkaka.addResonance980("R9802","R980kk",spin51,m51,g151,g251,rho51,phi51,propType51);
    pdfkaka.addResonance980("R9802","R980kk",spin51,m980,g10,g20,rho980kk2,phi980kk2,propType51);
    //    pdfkaka.addResonance("R13701","R1370kk",spin52,m52,w52,rho52,phi52,propType52);
    pdfkaka.addResonance("R13701","R1370kk",spin52,m1370,w1370,rho1370kk1,phi1370kk1,propType52);
    //    pdfkaka.addResonance("R13702","R1370kk",spin53,m53,w53,rho53,phi53,propType53);
    pdfkaka.addResonance("R13702","R1370kk",spin53,m1370,w1370,rho1370kk2,phi1370kk2,propType53);
    //    pdfkaka.addResonance("R17501","R1750kk",spin54,m54,w54,rho54,phi54,propType54);
    pdfkaka.addResonance("R17501","R1750kk",spin54,m1750,w1750,rho1750kk1,phi1750kk1,propType54);
    //    pdfkaka.addResonance("R17502","R1750kk",spin55,m55,w55,rho55,phi55,propType55);
    pdfkaka.addResonance("R17502","R1750kk",spin55,m1750,w1750,rho1750kk2,phi1750kk2,propType55);
    pdfkaka.addResonance600("R6001","R600kk",spin56, m56, b156, b256, b356, b456, b556, rho56,phi56,propType56);
    pdfkaka.addResonance600("R6002","R600kk",spin57, m57, b157, b257, b357, b457, b557, rho57,phi57,propType57);
    pdfkaka.addResonance("R15001","R1500kk",spin58,m58,w58,rho58,phi58,propType58);
    pdfkaka.addResonance("R15002","R1500kk",spin59,m59,w59,rho59,phi59,propType59);
    pdfkaka.addResonance("R17101","R1710kk",spin60,m1710,w1710,rho1710kk1,phi1710kk1,propType60);
    //    pdfkaka.addResonance("R17102","R1710kk",spin61,m61,w61,rho61,phi61,propType61);
    pdfkaka.addResonance("R17102","R1710kk",spin61,m1710,w1710,rho1710kk2,phi1710kk2,propType61);
    pdfkaka.addResonance("R12701","R1270kk",spin62,m62,w62,rho62,phi62,propType62);
    pdfkaka.addResonance("R12702","R1270kk",spin63,m63,w63,rho63,phi63,propType63);
    pdfkaka.addResonance("R12703","R1270kk",spin64,m64,w64,rho64,phi64,propType64);
    pdfkaka.addResonance("R12704","R1270kk",spin65,m65,w65,rho65,phi65,propType65);
    pdfkaka.addResonance("R12705","R1270kk",spin66,m66,w66,rho66,phi66,propType66);
    pdfkaka.addResonance("R15251","R1525kk",spin67,m67,w67,rho67,phi67,propType67);
    pdfkaka.addResonance("R15252","R1525kk",spin68,m68,w68,rho68,phi68,propType68);
    pdfkaka.addResonance("R15253","R1525kk",spin69,m69,w69,rho69,phi69,propType69);
    pdfkaka.addResonance("R15254","R1525kk",spin70,m70,w70,rho70,phi70,propType70);
    pdfkaka.addResonance("R15255","R1525kk",spin71,m71,w71,rho71,phi71,propType71);
    pdfkaka.addResonance1680("R16801","R1680kk",spin72, m72, m72980, w72, g172, g272,rho72,phi72,propType72);
    pdfkaka.addResonance1680("R16802","R1680kk",spin73, m73, m73980, w73, g173, g273,rho73,phi73,propType73);
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
