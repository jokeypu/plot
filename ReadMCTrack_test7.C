#include "TrackFinder.h"
int ReadMCTrack_test7()
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
    TH2D* histxy = new TH2D("hvx0vy0","vx vs vy",500,1.32,1.34,500,0,0.1);
    TGraph* gr = new TGraph();
    long nn(0);

    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    Int_t nlnks;
    //cout << "maxEvtNo " << maxEvtNo << endl;
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nclusters = fClusterArray->GetEntriesFast();
        //cout << "\n------------ EvtNo: " << ievt+1 << " ------------ ClustersNo: "<< nclusters << " ------------" << endl;
        for (int i = 0; i < nclusters; i++) {
            //cout << "\n***cluster*** " << i+1 <<endl;
            PndEmcCluster* cluster = (PndEmcCluster* ) fClusterArray->At(i);
            TVector3 pos = cluster->position();
            double Ptheta = cluster->theta();
            double Pphi = cluster->phi();
            double Eclus = cluster->energy();
            if (Eclus < 0.2) continue;
            //if ( theta > 1.2 && theta < 1.45 && phi > -0.1 && phi < 0.1 ){
            gr->SetPoint(nn,Ptheta,Pphi);
            nn++;
            //}
            //cout<<cluster->GetModule()<<endl;
            //cout<<"position:"<<pos.Mag()<<endl;
            /*std::vector<FairLink> links = (cluster->GetTrackEntering()).GetSortedMCTracks();
            std::vector<FairLink> linksout = (cluster->GetTrackExiting()).GetSortedMCTracks();
            //cout << (cluster->GetTrackEntering()).GetNLinks() << endl;
            nlnks = links.size();
            for (int j = 0; j < nlnks; j++){
                int id = links[j].GetIndex();
                PndMCTrack* track = (PndMCTrack*)fMCtrackArray->At(id);
                TLorentzVector s4momentum = track->Get4Momentum();
                Double_t EE = s4momentum.E();
                TVector3 svertex = track->GetStartVertex();
                double phi = svertex.Phi();
                double theta = svertex.Theta();
                double r = svertex.Mag();
                //double r = sqrt(svertex.x()*svertex.x()+svertex.y()*svertex.y());
                //if (EE > 0.2) continue;
                if (r < 66 ||  r > 70 ) continue;
                histxy->Fill(theta,phi);
		if ( theta > 1.2 && theta < 1.45 && phi > -0.1 && phi < 0.1 ){
		gr->SetPoint(nn,theta,phi);
		nn++;
		}
            }*/
            
        }
    }
    TCanvas* c1=new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
//    histxy->Draw("SCAT");
    gr->GetXaxis()->SetTitle("#theta");
    gr->GetYaxis()->SetTitle("#phi");
    gr->SetMarkerStyle(6);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(1);
//    gr->ComputeRange(xmin, ymin, xmax, ymax);
    gr->Draw("ap");
//    histxy->Draw("box");
   /* h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Energy(GeV)");
    h1->GetYaxis()->SetTitle("Entries");
    h1->SetLineColor(kRed);
    h1->Draw();*/
    return 0;
}

