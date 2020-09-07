int MC_backup5()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/test/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/test/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    TClonesArray* fHitArray = new TClonesArray("PndEmcHit");
    t->SetBranchAddress("EmcHit",&fHitArray);
    if (!fHitArray) return -1;
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    int bin1(40),bin2(40);
    double xmin(0.9),xmax(1.8),ymin(-0.6),ymax(0.6);
    //TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",bin1,xmin,xmax,bin2,ymin,ymax);
    TH1D* h1 = new TH1D("h1", "hit1", 200, 0.0, 1.2);
    TH1D* h2 = new TH1D("h2", "hit2", 200, 0.0, 1.2);
    
    std::map<Int_t, Double_t>::iterator it;
    
    int N=0;
    
    //for (Int_t ievt = 0; ievt < 1; ievt++) {
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        
        int nhits = fHitArray->GetEntriesFast();
        for ( Int_t i = 0; i < nhits; i++ ){
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t> ds = hit->GetDepositedEnergyMap();
            cout << "hit: " << i << endl;
            double hitE=0;
            double E = hit->GetEnergy();
            for( it=ds.begin(); it!=ds.end(); ++it){
                //cout << it->first << endl;
                //cout << it->second << endl << endl;
                if ( it->second > 1 ) N++;
                hitE += it->second;
            }
            cout << "Hit Energy(MC): " << hitE << endl;
            cout << "Hit Energy(MN): " << E << endl;
            
            h1->Fill(hitE);
            h2->Fill(E);
        }
    }
    
    cout << N << endl;
    //TCanvas* c1=new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    /*
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->SetFillColor(46);
    histxy->Draw("BOX");
    */
    histxy1->SetFillColor(13);
    //histxy1->Draw("SAMEBOX");
    
    h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Energy(GeV)");
    h1->GetYaxis()->SetTitle("Entries");
    h1->SetLineColor(kRed);
    
    h2->SetLineWidth(2);
    h2->GetXaxis()->SetTitle("Energy(GeV)");
    h2->GetYaxis()->SetTitle("Entries");
    h2->SetLineColor(kBlue);
    
    h1->Draw();
    h2->Draw("SAME");
    
    return 0;
}

