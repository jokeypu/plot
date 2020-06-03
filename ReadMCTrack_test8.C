#include "TrackFinder.h"

void hist(TClonesArray* fDigiArray, std::map<int, std::vector<int> > &digilist) {
    int bin1(40),bin2(40);
    double xmin(0.9),xmax(1.8),ymin(-0.6),ymax(0.6);
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy2=new TH2D("hvx0vy02","vx vs vy2",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy3=new TH2D("hvx0vy03","vx vs vy3",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy4=new TH2D("hvx0vy04","vx vs vy4",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy5=new TH2D("hvx0vy05","vx vs vy5",bin1,xmin,xmax,bin2,ymin,ymax);
    
    for (int j = 0; j < digilist.size(); j++) {
    for (int i = 0; i < digilist[j].size(); i++) {
            PndEmcDigi* digi = (PndEmcDigi* ) fDigiArray->At(digilist[j].at(i));
            double theta = digi->GetTheta();
            double phi = digi->GetPhi();
            double EE = 8000*pow((digi->GetEnergy()),0.4);
            //cout << EE << endl;
            for (int k = 0; k < EE; k++) {
                switch (j) {
                    case 0:
                        histxy->Fill(theta,phi);
                        break;
                    case 1:
                        histxy1->Fill(theta,phi);
                        break;
                    case 2:
                        histxy2->Fill(theta,phi);
                        break;
                    case 3:
                        histxy3->Fill(theta,phi);
                        break;
                    case 4:
                        histxy4->Fill(theta,phi);
                        break;
                    case 5:
                        histxy5->Fill(theta,phi);
                        break;
                    case -1:
                        break;
                    default:
                        cout << "-E-" << endl;
                        break;
                }
            }
    }
    }
    TCanvas* c1=new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->SetFillColor(46);
    histxy->Draw("BOX");
    histxy1->SetFillColor(13);
    histxy1->Draw("SAMEBOX");
    histxy2->SetFillColor(7);
    histxy2->Draw("SAMEBOX");
    histxy3->SetFillColor(9);
    histxy3->Draw("SAMEBOX");
    histxy4->SetFillColor(36);
    histxy4->Draw("SAMEBOX");
    histxy5->SetFillColor(6);
    histxy5->Draw("SAMEBOX");
}

int ReadMCTrack_test8()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/pi_phi_0/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/pi_phi_0/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    //TH1D* h1=new TH1D("h1","htemp1",200,0.,1.2);
   // TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",28,0.75,1.75,28,0,1);
    //TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",28,0.75,1.75,28,0,1);
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    Int_t nlnks;
    std::map<int, std::vector<int> > table;
    //cout << "maxEvtNo " << maxEvtNo << endl;
    Int_t ievt = 0;
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nclusters = fClusterArray->GetEntriesFast();
        cout << "\n------------ EvtNo: " << ievt+1 << " ------------ ClustersNo: "<< nclusters << " ------------" << endl;
        for (int i = 0; i < nclusters; i++) {
            //cout << "\n***cluster*** " << i+1 <<endl;
            PndEmcCluster* cluster = (PndEmcCluster* ) fClusterArray->At(i);
            //double theta = cluster->theta();
            //double phi = cluster->phi();
            std::vector<int> digilist = cluster->DigiList();
            int nlist = digilist.size();
            cout << nlist <<endl;
            table[i] = digilist;
            
        }
        hist(fDigiArray, table);
    
   /* h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Energy(GeV)");
    h1->GetYaxis()->SetTitle("Entries");
    h1->SetLineColor(kRed);
    h1->Draw();*/
    return 0;
}

