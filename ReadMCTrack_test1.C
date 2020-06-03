#include "TrackFinder.h"
int ReadMCTrack_test1()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    TH1D* h1=new TH1D("h1","htemp1",200,0.,1.2);
    
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
            std::vector<FairLink> links = (cluster->GetTrackEntering()).GetSortedMCTracks();
            nlnks = links.size();
            for (int j = 0; j < nlnks; j++){
            int id = links[j].GetIndex();
            PndMCTrack* track = (PndMCTrack*)fMCtrackArray->At(id);
                TLorentzVector s4momentum = track->Get4Momentum();
                Double_t EE = s4momentum.E();
                h1->Fill(EE);
                
            }
            
            
            
            
            
            
        }
    }
            
            TCanvas* c1=new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Energy(GeV)");
    h1->GetYaxis()->SetTitle("Entries");
            h1->SetLineColor(kRed);
            h1->Draw();
            
            
            
            
            
            
            
            
            
            
            
            
            
            /*
            std::map<int, std::vector<int> > table;
            for (int j = 0; j < nlnks-1; j++){
                for (int k = j+1; k < nlnks; k++){
                    int id1 = links[j].GetIndex();
                    int id2 = links[k].GetIndex();
                    TrackFinder* TF = new TrackFinder(fMCtrackArray);
                    TF->AddTrackID(id1,id2);
                    //TF->Print();
                    //TF->PrintPath();
                    std::vector <int> P1 = TF->GetPath1();
                    std::vector <int> P2 = TF->GetPath2();
                    for (int n = (P1.size()) - 1; n >= 0; n-- ) {table[in].push_back(P1[n]);}
                    in++;
                    for (int n = (P2.size()) - 1; n >= 0; n-- ) {table[in].push_back(P2[n]);}
                    in++;
                }
            }
            
            for (int i = 0;i < in-1;i++) {
                int N1 = table[i].size();
                if ( (table[i].at(N1-1))==-1 ) continue;
                for (int j =i+1;j < in;j++) {
                    int N2 = table[j].size();
                    if ( (table[j].at(N2-1))==-1 ) continue;
                    int NI = table[i].size();
                    int NJ = table[j].size();
                    for (int l = 0;l < NJ ; l++){
                        if ( table[i].at(NI-1) == table[j].at(l) && table[i].at(0) == table[j].at(l) )  {table[i].push_back(-1);break;}
                    }
                    for (int k = 0;k < NI ;k++ ){
                        if ( table[i].at(k) == table[j].at(NJ-1) && table[i].at(k) == table[j].at(0) )  {table[i].push_back(-1);break;}
                    }
                }
            }
            /
             for (int i = 0;i < in-1;i++) {
             int NN(0);
             for (int j =i+1;j < in;j++) {
             int NI = table[i].size();
             int NJ = table[j].size();
             for (int k = 0;k < NI ;k++ ){
             for (int l = 0;l < NJ ; l++){
             if ( table[i].at(k) == table[j].at(l) ){
             k = NI;
             NN++;
             }
             }
             }
             }
             }
             //
            for (int i = 0;i < in;i++) {
                int N = table[i].size();
                if ( table[i].at(N-1)==-1 ) continue;
                for (int j =0;j < N-1;j++) {
                    cout << table[i].at(j) << " -> ";
                }
                cout << table[i].at(N-1) << endl;
                
            }
            TrackFinder* TF1 = new TrackFinder(fMCtrackArray);
            std::vector <int> pr;
            pr.clear();
            for (int i = 0;i < in;i++){
                int N = table[i].size();
                if ( (table[i].at(N-1))==-1 ) continue;
                vector<int>::iterator result = find( pr.begin(), pr.end(), table[i].at(0) );
                if ( result == pr.end( ) ) { pr.push_back(table[i].at(0));cout<< "id:" << table[i].at(0) << "\t";TF1->Print(table[i].at(0)); }
            }
            std::vector <int> vect{0,1,2};
            TF1->AddTrackID(1,2);
            TF1->PrintPath();
            TF1->Print(0);
            TF1->PrintPDG(vect);
            cout << endl;
        }
    }*/
    return 0;
}

