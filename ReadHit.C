int ReadHit()
{
    int bin1(300),bin2(300);
    float tx(800),ty(600);
    double xmin(0),xmax(0.6),ymin(0),ymax(1);
    TVector3 V;
    V.SetMagThetaPhi(68,1.29,0);
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/Gamma_1_1/evtcomplete_sim.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    PndEmcMapper::Init(1);
    TFile *parfile = new TFile("../data/Gamma_1_1/evtcomplete_par.root");
    parfile->Get("FairGeoParSet");
    PndEmcStructure *fEmcStr = PndEmcStructure::Instance();
    PndEmcMapper *fMapper = PndEmcMapper::Instance();
    typedef std::map<Int_t, Float_t> mapper;
    mapper emcX = fEmcStr->GetEmcX();
    mapper emcY = fEmcStr->GetEmcY();
    mapper emcZ = fEmcStr->GetEmcZ();
    
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
    
    TH1D* hist=new TH1D("h","v",bin1,xmin,xmax);
    hist->GetXaxis()->SetTitle("#theta");
    hist->GetYaxis()->SetTitle("Entries");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    int maxEvtNo = ioman->CheckMaxEventNo();
    
    for (int ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        int nhits = fHitArray->GetEntriesFast();
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            Int_t DetID = hit->GetDetectorID();
            double E = hit->GetEnergy();
            TVector3 position(emcX[DetID], emcY[DetID], emcZ[DetID]);
            double distance =  position.Angle(V);
            //cout << position.Mag() << " " << position.Theta() << " " << position.Phi() << endl;
            //histxy->Fill(position.Theta(),position.Phi());
            //histxy->Fill(position.Theta(),E);
            histxy->Fill(distance,E);
            //hist->Fill(position.Theta(),E/maxEvtNo);
        }
    }
    //hist->Draw();
    histxy->Draw("SCAT");
    return 0;
}


