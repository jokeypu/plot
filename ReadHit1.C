int ReadHit1()
{
    int bin1(300),bin2(300);
    float tx(800),ty(600);
    double xmin(0),xmax(3),ymin(0),ymax(1);
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/Gamma1/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../data/Gamma1/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    /*
     PndEmcMapper::Init(1);
     TFile *parfile = new TFile("../data/new1/evtcomplete_par.root");
     parfile->Get("FairGeoParSet");
     PndEmcStructure *fEmcStr = PndEmcStructure::Instance();
     PndEmcMapper *fMapper = PndEmcMapper::Instance();
     typedef std::map<Int_t, Float_t> mapper;
     mapper emcX = fEmcStr->GetEmcX();
     mapper emcY = fEmcStr->GetEmcY();
     mapper emcZ = fEmcStr->GetEmcZ();
     */
    
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
    //histxy->Draw();
    int excnum(0);
    
    for (Int_t ievt = 1; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        
        int nbump = fBumpArray->GetEntriesFast();
        if ( nbump > 1 ) continue;
        
        int ndigi = fDigiArray->GetEntriesFast();
        PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(n);
        for ( int idigi = 0; idigi < ndigi ; idigi++) {
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(idigi);
            double E = digi->GetEnergy();
            double d = bump->DistanceToCentre(digi);
            histxy->Fill(d,E);
        }
        excnum++;
    }
    histxy->Draw("SCAT");
    cout << "Exc Num:" << excnum << endl;
    return 0;
}
