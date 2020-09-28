bool cmp(const pair<Int_t, Double_t>& a, const pair<Int_t, Double_t>& b) {
        return a.second > b.second;
}
int Energy_MCTruth_v2(int aa=1)
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/Npi0_3GeV/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/Npi0_3GeV/evtcomplete_sim.root");
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
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->GetXaxis()->CenterTitle();
    histxy->GetYaxis()->CenterTitle();
    histxy->Draw();
    
    for (Int_t ievt = aa; ievt < aa+1; ievt++) {
        //for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        
        int nhits = fHitArray->GetEntriesFast();
        for ( Int_t i = 0; i < nhits; i++ ){
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t> ds = hit->GetDepositedEnergyMap();
            vector<pair<Int_t, Double_t>> SourceEnergy(ds.begin(), ds.end());
            sort(SourceEnergy.begin(), SourceEnergy.end(), cmp);
            cout << "hit: " << i << endl;
            double E = hit->GetEnergy();
            double theta = hit->GetTheta();
            double phi = hit->GetPhi();
            
            for( Int_t j = 0; j < SourceEnergy.size(); j++){
                cout << SourceEnergy[j].first << endl;
                cout << SourceEnergy[j].second << endl;
                double S = 0.8 * pow(SourceEnergy[j].second, 1 / M_E);
                TWbox *twb = new TWbox(theta-S,phi-S,theta+S,phi+S,90-20*((SourceEnergy[j].first)),-9,0);
                twb->Draw("SAME");
            }
        }
    }
    return 0;
}

