int Energy_MCTruth(int aa=1)
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/new1/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/new1/evtcomplete_sim.root");
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
    
    int bin1(100),bin2(100);
    float tx(800),ty(600);
    float xsct(74),ysct(0);
    float Rad(15);
    double xmin(xsct-Rad*tx/ty),xmax(xsct+Rad*tx/ty),ymin(ysct-Rad),ymax(ysct+Rad);
    
    TCanvas* c1=new TCanvas("PANDA","MCTruth",tx,ty);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->Draw();
    
    std::map<Int_t, Double_t>::iterator it;
    
    for (Int_t ievt = aa; ievt < aa+1; ievt++) {
        //for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        
        int nhits = fHitArray->GetEntriesFast();
        for ( Int_t i = 0; i < nhits; i++ ){
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t> ds = hit->GetMcSourceEnergy();
            cout << "hit: " << i << endl;
            double E = hit->GetEnergy();
            double theta = hit->GetTheta();
            double phi = hit->GetPhi();
            
            for( it=ds.begin(); it!=ds.end(); ++it){
                cout << it->first << endl;
                cout << it->second << endl;
                double S = 0.8 * pow(it->second, 1 / M_E);
                TWbox *twb = new TWbox(theta-S,phi-S,theta+S,phi+S,90-20*((it->first)),-9,0);
                twb->Draw("SAME");
            }
            
            
        }
    }
    return 0;
}

