bool cmp(const pair<Int_t, Int_t>& a, const pair<Int_t, Int_t>& b) {
    return a.second > b.second;
}
int Energy_Compare_new()
{
    int bin1(500),bin2(500);
    float tx(800),ty(600);
    double xmin(0),xmax(4),ymin(-1.2),ymax(1.2);
    
    //*****************************************************
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/Compare_n/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../data/Compare_n/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    
    TClonesArray* fHitArray = new TClonesArray("PndEmcHit");
    t->SetBranchAddress("EmcHit",&fHitArray);
    if (!fHitArray) return -1;
    
    PndEmcMapper::Init(1);
    TFile *parfile = new TFile("../data/Compare_n/evtcomplete_par.root");
    parfile->Get("FairGeoParSet");
    PndEmcStructure *fEmcStr = PndEmcStructure::Instance();
    PndEmcMapper *fMapper = PndEmcMapper::Instance();
    
    TClonesArray* fBumpArray = (TClonesArray*) ioman->GetObject("EmcBump");
    if (!fBumpArray) return -1;
    TClonesArray* fSharedDigiArray = (TClonesArray*) ioman->GetObject("EmcSharedDigi");
    if (!fSharedDigiArray) return -1;
    
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    TCanvas* c1=new TCanvas("PANDA","MCTruth",tx,ty);
    TCanvas* c2=new TCanvas("PANDA1","MCTruth1",tx,ty);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy->GetXaxis()->SetTitle("E_{turth}");
    histxy->GetYaxis()->SetTitle("w - w_{turth}");
    histxy->GetXaxis()->CenterTitle();
    histxy->GetYaxis()->CenterTitle();
    
    //TH1D* hist=new TH1D("hvx0","vx",bin2,-1.2,-0.9);
    TH1D* hist=new TH1D("hvx0","vx",100,ymin,ymax);
    hist->GetXaxis()->SetTitle("w - w_{turth}");
    hist->GetYaxis()->SetTitle("Entries");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    
    
    std::map<Int_t, Int_t>::iterator it;
    std::map<Int_t, Double_t>::iterator p;
    Int_t NN = 0;
    
    
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        
        int nbump = fBumpArray->GetEntriesFast();
        
        for (int n = 0; n < nbump ; n++){
            Int_t TruthIndex;
            std::map<Int_t, Int_t> count;
            std::vector<pair<Int_t, Int_t>> mylist;
            PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(n);
            //cout << "Bump " << n << " ---> " << " Cluster " << bump->GetClusterIndex() << endl;
            std::vector<Int_t> list = bump->DigiList();
            for (int i=0; i < list.size(); i++){
                //cout << i <<" "<<list[i]<<endl;
                PndEmcSharedDigi* sharedigi = (PndEmcSharedDigi*)fSharedDigiArray->At(list[i]);
                double w = sharedigi->weight();
                //double DetID = sharedigi->GetDetectorId();
                PndEmcTwoCoordIndex* BumpTCI = sharedigi->GetTCI();
                PndEmcHit* hit = (PndEmcHit*)fHitArray->At(sharedigi->GetHitIndex());
                std::map<Int_t, Double_t> Mclist = hit->GetDepositedEnergyMap();
                if ( Mclist.size() == 1 ) {
                    count[(Mclist.begin())->first]++;
                }
            }
            for (it = count.begin(); it !=count.end() ; it++){
                mylist.push_back( std::pair<Int_t, Int_t>(it->first,it->second) );
            }
            sort(mylist.begin(), mylist.end(), cmp);
            if ( mylist.size() == 1 ) {
                TruthIndex = mylist[0].first;
            }else if ( mylist.size() >1 && mylist[0].second / mylist[1].second > 2 ){
                    TruthIndex = mylist[0].first;
            }else {
                NN++;
                continue;
            }
            
            for (int i=0; i < list.size(); i++){
                PndEmcSharedDigi* sharedigi = (PndEmcSharedDigi*)fSharedDigiArray->At(list[i]);
                PndEmcHit* hit = (PndEmcHit*)fHitArray->At(sharedigi->GetHitIndex());
                double w = sharedigi->weight();
                std::map<Int_t, Double_t> Mclist = hit->GetDepositedEnergyMap();
                double E = hit->GetEnergy();
                double Etotal = hit->GetEnergy();
                double wt = 0;
                if ( Mclist.find(TruthIndex) != Mclist.end() ) wt = ((Mclist.find(TruthIndex))->second)/Etotal;
                //if ( E < 1 ) continue;
                //if ( w <= 0 ) cout << "########" << endl;
                histxy->Fill(E, w - wt);
                hist->Fill( w - wt);
            }
        }
    }
    cout << NN << endl;
    c1->cd();
    histxy->Draw("SCAT");
    c2->cd();
    hist->Draw();
    return 0;
}

