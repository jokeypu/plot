bool cmp(const pair<Int_t, Double_t>& a, const pair<Int_t, Double_t>& b) {
        return a.second > b.second;
}
int crystal_pos(int aa=1)
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/Gamma_1G_all_all/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/Gamma_1G_all_all/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    PndEmcMapper::Init(1);
    TFile *parfile = new TFile("../../data/Gamma_1G_all_all/evtcomplete_par.root");
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
    
    TCanvas* c1=new TCanvas("PANDA","Bump",1000,800);
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
    
    TGraph *g = new TGraph();
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kBlue+1, 1);
    g->GetXaxis()->SetTitle("Theta");
    g->GetYaxis()->SetTitle("pt");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    int N = 0;
    vector<double> v;
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        
        int nhits = fHitArray->GetEntriesFast();

        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            Int_t detID = hit->GetDetectorID();
            TVector3 position(emcX[detID], emcY[detID], emcZ[detID]);
            Double_t theta = position.Theta();
            theta *= TMath::RadToDeg();
            if (theta >140 || theta <20) continue;
            if (position.Pt()>70) continue;
            g->SetPoint(N,theta,position.Pt());
            v.push_back(theta);
            N++;
        }
    }
    sort(v.begin(),v.end());
    
    map<int, double> mm;
    N=-1;
    int count=0;
    double value_std = -100;
    for (int i=0;i<v.size();i++){
        double value = v[i];
        if (fabs(value - value_std)<0.5) {mm[N] += value;count++;}
        else {
            if (N>0) mm[N] = mm[N]/count;
            N++;
            value_std = value;
            mm[N] = value;
            count = 1;
        }
    }
    
    N=0;
    for (std::map<int, double>::iterator it = mm.begin(); it != mm.end(); it++) {
        if (it->second >142) continue;
        //g->SetPoint(N,it->second,60);
        cout << N+1 << " " << it->second << endl;
        N++;
    }
    g->GetYaxis()->SetRangeUser(0,70);
    g->Draw("AP.");
    return 0;
}
