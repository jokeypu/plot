bool cmp(const pair<Int_t, Double_t>& a, const pair<Int_t, Double_t>& b) {
    return a.first < b.first;
}
int Energy_Bump(int evNo = 1, bool UsedGrid = true)
{
    int bin1(100),bin2(100);
    float tx(800),ty(600);
    float xsct(74),ysct(0);
    float Rad(12);
    double xmin(xsct-Rad*tx/ty),xmax(xsct+Rad*tx/ty),ymin(ysct-Rad),ymax(ysct+Rad);
    
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
    //c1->cd()->SetGrid(1,1);
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
    
    if ( UsedGrid ) {
        Double_t px1(72.8607),px2(74.9418);
        Double_t py1(0.496197),py2(-1.70964);
        Double_t xstart = (px1+px2)/2, ystart = (py1+py2)/2;
        Double_t xstep = abs(px1-px2), ystep = abs(py1-py2);
        Int_t width(1),style(8);
        Color_t color(kGray);
        for (int n = 0; xstart + n * xstep < xmax; n++) {
            TLine* line = new TLine();
            line->SetLineColor(color);
            line->SetLineWidth(width);
            line->SetLineStyle(style);
            line->SetX1(xstart + n * xstep);
            line->SetY1(ymin);
            line->SetX2(xstart + n * xstep);
            line->SetY2(ymax);
            line->Draw("SAME");
        }
        for (int n = 1; xstart - n * xstep > xmin; n++) {
            TLine* line = new TLine();
            line->SetLineColor(color);
            line->SetLineWidth(width);
            line->SetLineStyle(style);
            line->SetX1(xstart - n * xstep);
            line->SetY1(ymin);
            line->SetX2(xstart - n * xstep);
            line->SetY2(ymax);
            line->Draw("SAME");
        }
        for (int n = 0; ystart + n * ystep < ymax; n++) {
            TLine* line = new TLine();
            line->SetLineColor(color);
            line->SetLineWidth(width);
            line->SetLineStyle(style);
            line->SetX1(xmin);
            line->SetY1(ystart + n * ystep);
            line->SetX2(xmax);
            line->SetY2(ystart + n * ystep);
            line->Draw("SAME");
        }
        for (int n = 1; ystart - n * ystep > ymin; n++) {
            TLine* line = new TLine();
            line->SetLineColor(color);
            line->SetLineWidth(width);
            line->SetLineStyle(style);
            line->SetX1(xmin);
            line->SetY1(ystart - n * ystep);
            line->SetX2(xmax);
            line->SetY2(ystart - n * ystep);
            line->Draw("SAME");
        }
    }
    
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
                double S = 0.8 * pow((it->second)[j].second, 1 / M_E);
                float alpha = 0.2 + 2 * sqrt(S) / 2;
                TWbox *twb = new TWbox(theta-S,phi-S,theta+S,phi+S,0,-9,0);
                twb->SetFillColorAlpha(90-20*((it->second)[j].first), alpha);
                twb->Draw("SAME");
            }
        }
    }
    return 0;
}
