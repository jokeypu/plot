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
void EEE_v2(){
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
    //Double_t x1_min(xmin),x1_max(xmax);
    
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
    
    TF1* f1=new TF1("f1","[0]*Novosibirsk(x,[1],[2],[3])",0.85,1.2);
    f1->SetParameters(3000,1.,.02,.3);
    f1->SetParLimits(0, 0, 10000);
    f1->SetParLimits(1, 0.9, 1.3);
    f1->SetParLimits(2, 0.001, 0.1);
    f1->SetParLimits(3, 0, 10);
    //f1->SetParLimits(4, 0, 1000);
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/Gamma_tow_1G_last_o/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TClonesArray* fBumpArray = (TClonesArray*) ioman->GetObject("EmcBump");
    if (!fBumpArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    int cunt1(0);
    
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt);
        int nbump = fBumpArray->GetEntriesFast();
        for (int n = 0; n < nbump ; n++){
            PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(n);
            Double_t bump_E = bump->energy();
            //if (bump_E<1.5)
            h1D1->Fill(bump_E);
            cunt1++;
        }
    }
      
    int cunt2(0);
    
    TFile* f = new TFile("../data/Gamma_tow_1G_last/evtcomplete_digi.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fBumpArray_fix = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray_fix);
    if (!fBumpArray_fix) return -1;
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        //if (cunt2>=cunt1) continue;
        t->GetEntry(ievt);
        int nbump = fBumpArray_fix->GetEntriesFast();
        for (int n = 0; n < nbump ; n++){
            PndEmcBump* bump = (PndEmcBump*)fBumpArray_fix->At(n);
            Double_t bump_E = bump->energy();
            //if (bump_E<1.5)
            h1D2->Fill(bump_E);
            cunt2++;
            //if (cunt2>=cunt1) break;
        }
    }
    cout << "cunt1:" << cunt1 << ", cunt2:" << cunt2 << endl;
    
    h1D1->Fit(f1,"R");
    return;
}
