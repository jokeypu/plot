#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "TFile.h"
#include "TTree.h"
#include "RooFitResult.h"
#include <sstream>
#include <string>

#include "bes3plotstyle.C"
using namespace std;
using namespace RooFit;
void EEE(){
    gROOT->SetStyle("Plain");
    gSystem->Load("libRooFit");
    
    Double_t y;
    int min = 320000;
    string str;
    
    
    TString namet(""), namex("Energy"), namey("Entries");
    
    int bin1(50),bin2(200);
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
    
    //TCanvas* c = new TCanvas("bes3fit","BESIII Fit", 1200,800);
    
    RooRealVar x("x",namex,x1_min,x1_max);
    //RooRealVar mean("mean","mean",1.0,0.7,1.3);
    //RooRealVar sigma("sigma","sigma",0.026,0.00,0.1);
    //RooRealVar sigmap("sigmap","sigmap",0.01,0.00,0.18);
    RooRealVar mean("mean","mean",1.0,0.7,1.3);
    RooRealVar sigma("sigma","sigma",0.026,0.00,1);
    RooRealVar sigmap("sigmap","sigmap",0.1,-50,50);
    RooRealVar nn("nn","nn",1,-100,100);
    RooVoigtian sig("sig","signal p.d.f.",x,mean,sigma,sigmap);
    RooVoigtian sig1("sig1","signal p.d.f.",x,mean,sigma,sigmap);
    RooCBShape sig2("sig2","signal p.d.f.2",x,mean,sigma,sigmap,nn);
    
    RooRealVar c0("c0","coefficient #0", 1.0,-1.5,1.5);
    RooRealVar c1("c1","coefficient #1", 0.1,-1.5,1.5);
    RooRealVar c2("c2","coefficient #2",-0.1,-1.5,1.5);
    RooChebychev bkg("bkg","background p.d.f.",x,RooArgList(c0,c1,c2));
    
    RooRealVar nsig("nsig","signal fraction",2150,0.,40000.);
    RooRealVar nbkg("nbkg","background fraction",700,0.,50000.);
    RooAddPdf model("model","model",RooArgList(sig),RooArgList(nsig));
    RooAddPdf model1("model1","model1",RooArgList(sig),RooArgList(nsig));
    
    RooPlot* frame = x.frame(Title(namet));
    RooDataHist data("data","dataset with x",x,h1D1);
    
    model.fitTo(data,Extended());
    
    data.plotOn(frame, MarkerStyle(21), MarkerSize(1),MarkerColor(kBlue));
    
    model.plotOn(frame,LineColor(kBlue),LineWidth(3));
    //model.plotOn(frame, Components(sig),LineStyle(kDashed));
    
    frame->GetXaxis()->CenterTitle( kTRUE );
    frame->GetYaxis()->CenterTitle( kTRUE );
    frame->GetXaxis()->SetLabelFont( 42 );
    frame->GetYaxis()->SetLabelFont( 42 );
    frame->SetXTitle( namex );
    frame->SetYTitle( namey );
    frame->Draw();
    
    RooPlot* frame1 = x.frame(Title(namet));
    RooDataHist data1("data","dataset with x1",x,h1D2);
    model1.fitTo(data1,Extended());
    data1.plotOn(frame1, MarkerStyle(22), MarkerSize(1),MarkerColor(kRed));
    model1.plotOn(frame1,LineColor(kRed),LineWidth(3));
    frame1->GetXaxis()->CenterTitle( kTRUE );
    frame1->GetYaxis()->CenterTitle( kTRUE );
    frame1->GetXaxis()->SetLabelFont( 42 );
    frame1->GetYaxis()->SetLabelFont( 42 );
    frame1->SetXTitle( namex );
    frame1->SetYTitle( namey );
    frame1->Draw("SAME");
    //c->Print("Picture/"+tree1+"-"+branch1+"-F.png");
    return;
}
