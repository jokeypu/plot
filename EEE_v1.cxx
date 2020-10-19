#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooNovosibirsk.h"
#include "TFile.h"
#include "TTree.h"
#include "RooFitResult.h"
#include <sstream>
#include <string>

#include "bes3plotstyle.C"
using namespace std;
using namespace RooFit;
Double_t Novosibirsk(Double_t x,Double_t peak=0.,Double_t width=0.,Double_t tail=0.)
{
  if (TMath::Abs(tail) < 1.e-7) {
    return TMath::Exp( -0.5 * TMath::Power( ( (x - peak) / width ), 2 ));
  }

  Double_t arg = 1.0 - ( x - peak ) * tail / width;

  if (arg < 1.e-7) {
    //Argument of logarithm negative. Real continuation -> function equals zero
    return 0.0;
  }

  Double_t log = TMath::Log(arg);
  static const Double_t xi = 2.3548200450309494; // 2 Sqrt( Ln(4) )

  Double_t width_zero = ( 2.0 / xi ) * TMath::ASinH( tail * xi * 0.5 );
  Double_t width_zero2 = width_zero * width_zero;
  Double_t exponent = ( -0.5 / (width_zero2) * log * log ) - ( width_zero2 * 0.5 );

  return TMath::Exp(exponent);
}
void EEE_v1(){
    gROOT->SetStyle("Plain");
    gSystem->Load("libRooFit");
    
    Double_t y;
    int min = 320000;
    string str;
    
    
    TString namet(""), namex("Energy"), namey("Entries");
    
    int bin1(500),bin2(200);
    float tx(1200),ty(900);
    double xmin(0.85),xmax(1.2);
    //double xmin(0.3),xmax(0.7);
    string out1_name("doc/out1.txt"), out2_name("doc/out2.txt"), out3_name("doc/out3.txt");
    Double_t x1_min(xmin),x1_max(xmax);
    
    TCanvas* c=new TCanvas("PANDA1","c1",tx,ty);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(1);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    TH1D* h1D1 = new TH1D("Hist1_1","h1_1", bin1, xmin, xmax);
    h1D1->SetLineColor(kBlue);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("Energy");
    h1D1->GetYaxis()->SetTitle("Entries");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    TH1D* h1D2 = new TH1D("Hist1_2","h1_2", bin1, xmin, xmax);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);
    h1D2->GetXaxis()->SetTitle("Energy");
    h1D2->GetYaxis()->SetTitle("Entries");
    h1D2->GetXaxis()->CenterTitle();
    h1D2->GetYaxis()->CenterTitle();
    
    TH1D* h1D3 = new TH1D("Hist1_3","h1_3", bin1, xmin, xmax);
    h1D3->SetLineColor(kGreen);
    h1D3->SetLineWidth(2);
    h1D3->GetXaxis()->SetTitle("Energy");
    h1D3->GetYaxis()->SetTitle("Entries");
    h1D3->GetXaxis()->CenterTitle();
    h1D3->GetYaxis()->CenterTitle();
    
    ifstream out1, out2, out3;
    out1.open(out1_name, ios::in);
    out2.open(out2_name, ios::in);
    out3.open(out3_name, ios::in);
    int cunt;
    cunt = 0;
    for (int i= 0;i < min; i++) {
        getline(out1,str);
        double value= atof(str.c_str());
        h1D1->Fill(value);
    }
    for (int i= 0;i < min; i++) {
        getline(out2,str);
        double value= atof(str.c_str());
        h1D2->Fill(value);
    }
    for (int i= 0;i < min; i++) {
        getline(out3,str);
        double value= atof(str.c_str());
        h1D3->Fill(value);
    }
    out1.close();
    out2.close();
    out3.close();
    TF1* f1=new TF1("f1","[0]*Novosibirsk(x,[1],[2],[3])",0.85,1.2);
    f1->SetParameters(3000,1.,.02,.3);
    f1->SetParLimits(0, 0, 10000);
    f1->SetParLimits(1, 0.9, 1.3);
    f1->SetParLimits(2, 0.001, 0.1);
    f1->SetParLimits(3, 0, 10);
    //f1->SetParLimits(4, 0, 1000);
    h1D3->Fit(f1,"R");
    return;
}
