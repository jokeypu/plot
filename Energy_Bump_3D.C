bool cmp(const pair<Int_t, Double_t>& a, const pair<Int_t, Double_t>& b) {
    return a.first < b.first;
}
int Energy_Bump_3D(int evNo = 6, bool UsedGrid = false)
{
    Double_t px1(72.8607),px2(74.9418);
    Double_t py1(0.496197),py2(-1.70964);
    Double_t xstart = (px1+px2)/2, ystart = (py1+py2)/2;
    Double_t xstep = abs(px1-px2), ystep = abs(py1-py2);
    
    float tx(800),ty(600);
    float xsct(74),ysct(0);
    float Rad(12);
    double xmin(xsct-Rad*tx/ty),xmax(xsct+Rad*tx/ty),ymin(ysct-Rad),ymax(ysct+Rad);
    int bin1 = (xmax-xmin)/xstep;
    int bin2 = (ymax-ymin)/ystep;
    
    //*****************************************************
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/new1/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../data/new1/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    PndEmcMapper::Init(1);
    TFile *parfile = new TFile("../data/new1/evtcomplete_par.root");
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
    
    TCanvas* c1=new TCanvas("PANDA","Bump",tx,ty);
    c1->SetLogz();
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
    gStyle->SetTitleOffset(1.3,"xyz");
    
    TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy1->GetXaxis()->SetTitle("#theta");
    histxy1->GetYaxis()->SetTitle("#phi");
    histxy1->GetZaxis()->SetTitle("Energy (MeV)");
    histxy1->GetXaxis()->CenterTitle();
    histxy1->GetYaxis()->CenterTitle();
    histxy1->GetZaxis()->SetTitleOffset(1.0);
    
    TH2D* histxy2=new TH2D("hvx0vy02","vx vs vy2",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy2->GetXaxis()->SetTitle("#theta");
    histxy2->GetYaxis()->SetTitle("#phi");
    histxy2->GetZaxis()->SetTitle("Energy (MeV)");
    histxy2->GetXaxis()->CenterTitle();
    histxy2->GetYaxis()->CenterTitle();
    histxy2->GetZaxis()->SetTitleOffset(1.0);
    /*
    TH2D* histxy3=new TH2D("hvx0vy03","vx vs vy3",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy3->GetXaxis()->SetTitle("#theta");
    histxy3->GetYaxis()->SetTitle("#phi");
    histxy3->GetZaxis()->SetTitle("Energy (MeV)");
    histxy3->GetXaxis()->CenterTitle();
    histxy3->GetYaxis()->CenterTitle();
    histxy3->GetZaxis()->SetTitleOffset(1.0);*/
    
    std::map<Int_t, std::vector<pair<Int_t, Double_t>> >::iterator it;
    
    for (Int_t ievt = evNo; ievt < evNo+1; ievt++) {
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
                double w = (it->second)[j].second;
                //cout << "*******" << (it->second)[j].first << endl;
                if ((it->second)[j].first == 0 ) {
                    histxy1->Fill(theta,phi,1000 * w);
                }else{
                    histxy2->Fill(theta,phi,1000 * w);
                }
            }
        }
    }
    histxy1->SetFillColorAlpha(kAzure+3,0.5);
    histxy1->SetLineColor(kAzure+3);
    //histxy1->GetZaxis()->SetRangeUser(1,10000);
    //histxy1->SetFillStyle(3013);
    histxy1->Draw("LEGO1");
    
    histxy2->SetFillColorAlpha(kYellow-2,0.5);
    histxy2->SetLineColor(kYellow-2);
    //histxy2->GetZaxis()->SetRangeUser(1,10000);
    //histxy2->SetFillStyle(3013);
    histxy2->Draw("LEGO1,same");
    
    //histxy3->SetFillColorAlpha(kWhite,1);
    //histxy3->Draw("LEGO1same");
    
    return 0;
}
