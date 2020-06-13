bool cmp(const pair<Int_t, Double_t>& a, const pair<Int_t, Double_t>& b) {
        return a.second > b.second;
}
int Energy_Bump(int aa=1)
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
    
    PndEmcMapper::Init(1);
    TFile *parfile = new TFile("../../data/new1/evtcomplete_par.root");
    parfile->Get("FairGeoParSet");
    PndEmcStructure *fEmcStr = PndEmcStructure::Instance();
    PndEmcMapper *fMapper = PndEmcMapper::Instance();
    typedef std::map<Int_t, Float_t> mapper;
    mapper emcX = fEmcStr->GetEmcX();
    mapper emcY = fEmcStr->GetEmcY();
    mapper emcZ = fEmcStr->GetEmcZ();
    
    TClonesArray* fHitArray = new TClonesArray("PndEmcHit");
    t->SetBranchAddress("EmcHit",&fHitArray);
    if (!fHitArray) return -1;
    
    TClonesArray* fBumpArray = (TClonesArray*) ioman->GetObject("EmcBump");
    if (!fBumpArray) return -1;
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    TClonesArray* fSharedDigiArray = (TClonesArray*) ioman->GetObject("EmcSharedDigi");
    if (!fSharedDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    int bin1(100),bin2(100);
    float tx(800),ty(600);
    float xsct(74),ysct(0);
    float Rad(15);
    double xmin(xsct-Rad*tx/ty),xmax(xsct+Rad*tx/ty),ymin(ysct-Rad),ymax(ysct+Rad);
    
    TCanvas* c1=new TCanvas("PANDA","Bump",tx,ty);
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
    
    std::map<Int_t, std::vector<pair<Int_t, Double_t>> >::iterator it;
    
    for (Int_t ievt = aa; ievt < aa+1; ievt++) {
        //for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        std::map<Int_t, std::vector<pair<Int_t, Double_t>> > SourceEnergy;
        
        int nbump = fBumpArray->GetEntriesFast();
        for (int n = 0; n < nbump ; n++){
            PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(n);
            cout << "Bump " << n << " ---> " << " Cluster " << bump->GetClusterIndex() << endl;
            std::vector<Int_t> list = bump->DigiList();
            for (int i=0; i < list.size(); i++){
                cout << i <<" "<<list[i]<<endl;
                PndEmcDigi* digi = (PndEmcDigi*)fSharedDigiArray->At(list[i]);
                double E = digi->GetEnergy();
                double DetID = digi->GetDetectorId();
                SourceEnergy[DetID].push_back(make_pair(n,E));
            }
        }
        for ( it = SourceEnergy.begin(); it != SourceEnergy.end(); it++){
            sort((it->second).begin(), (it->second).end(), cmp);
            TVector3 position(emcX[it->first], emcY[it->first], emcZ[it->first]);
            double theta = position.Theta();
            double phi = position.Phi();
            theta = theta * ( 180.0 / 3.1416);
            phi = phi * ( 180.0 / 3.1416);
            for (int j = 0; j < (it->second).size(); j++ ){
                double S = 0.8 * pow((it->second)[j].second, 1 / M_E);
                TWbox *twb = new TWbox(theta-S,phi-S,theta+S,phi+S,90-20*((it->second)[j].first),-9,0);
                twb->Draw("SAME");
            }
        }
    }
    return 0;
}
