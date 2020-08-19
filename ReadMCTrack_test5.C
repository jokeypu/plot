#include "TrackFinder.h"
int ReadMCTrack_test5()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/Gamma/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/Gamma/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    //TH1D* h1=new TH1D("h1","htemp1",200,0.,1.2);
    int bin1(40),bin2(40);
    double xmin(0.9),xmax(1.8),ymin(-0.6),ymax(0.6);
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",bin1,xmin,xmax,bin2,ymin,ymax);
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    Int_t nlnks;
    //cout << "maxEvtNo " << maxEvtNo << endl;
    for (Int_t ievt = 0; ievt < 1; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nclusters = fClusterArray->GetEntriesFast();
        cout << "Number of Clusters: " << nclusters << endl;
        //cout << "\n------------ EvtNo: " << ievt+1 << " ------------ ClustersNo: "<< nclusters << " ------------" << endl;
        for (int i = 0; i < nclusters; i++) {
            cout << "\n***cluster*** " << i+1 <<endl;
            PndEmcCluster* cluster = (PndEmcCluster* ) fClusterArray->At(i);
            //double theta = cluster->theta();
            //double phi = cluster->phi();
            std::vector<Int_t> digilist = cluster->DigiList();
            int nlist = digilist.size();
            //cout << nlist <<endl;
            for (int j = 0; j < nlist; j++) {
                PndEmcDigi* digi = (PndEmcDigi* ) fDigiArray->At(digilist[j]);
                double theta = digi->GetTheta();
                double phi = digi->GetPhi();
                double EE = 8000 * pow(digi->GetEnergy(), 0.4);
                //cout << digi->GetEnergy() << endl;
                for (int k = 0; k < EE; k++) {
                    if ( i == 0 ) histxy->Fill(theta,phi);
                    if ( i == 1 ) histxy1->Fill(theta,phi);
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
    histxy1->SetFillColor(27);
    histxy1->Draw("SAMEBOX");
   /* h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Energy(GeV)");
    h1->GetYaxis()->SetTitle("Entries");
    h1->SetLineColor(kRed);
    h1->Draw();*/
    return 0;
}

