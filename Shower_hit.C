int Shower_hit(){
    int bin1(200),bin2(200),bin3(200);
    float tx(800),ty(600);
    double xmin(0),xmax(5),ymin(0),ymax(45),zmin(0),zmax(1.1);
    TString dir_name("Gamma_one_1G");
    TVector3 vz(0, 0, 1);
    
    //******************************************//
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/"+dir_name+"/evtcomplete_sim.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TClonesArray* fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
    if (!fPointArray) return -1;
    TClonesArray* fMCTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMCTrackArray) return -1;
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    TFile* f = new TFile("../data/"+dir_name+"/evtcomplete_digi.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    
    PndEmcMapper::Init(1);
    TFile *parfile = new TFile("../data/"+dir_name+"/evtcomplete_par.root");
    parfile->Get("FairGeoParSet");
    PndEmcStructure *fEmcStr = PndEmcStructure::Instance();
    PndEmcMapper *fMapper = PndEmcMapper::Instance();
    typedef std::map<Int_t, Float_t> mapper;
    mapper emcX = fEmcStr->GetEmcX();
    mapper emcY = fEmcStr->GetEmcY();
    mapper emcZ = fEmcStr->GetEmcZ();
    
    //TCanvas* c1=new TCanvas("PANDA1","Hit1",tx,ty);
    TCanvas* c2=new TCanvas("PANDA2","Hit2",tx,ty);
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
    gStyle->SetOptFit(1111);
    
    TH1D* h1D = new TH1D("h1D","h1",200,0,15);
    h1D->SetLineColor(kBlue);
    h1D->SetLineWidth(2);
    h1D->GetXaxis()->SetTitle("distance");
    h1D->GetYaxis()->SetTitle("Energy");
    h1D->GetXaxis()->CenterTitle();
    h1D->GetYaxis()->CenterTitle();
    
    TH3D* h2D = new TH3D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax,bin3,zmin,zmax);
    //TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",50,0,5,50,0,1);
    h2D->GetYaxis()->SetTitle("angle");
    h2D->GetXaxis()->SetTitle("d (cm)");
    h2D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 0.5);
    
    /*
    TF1 *f=new TF1("f","[1]*[0]*x",2,10);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.25,1);
    f->SetParLimits(0, 0.01, 10);
    f->SetParLimits(1, 1, 1);
     */
    /*
    TF1 *f=new TF1("f","[0]*[1]*sqrt(x-1.2)",1.2,15);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.25,1);
    f->SetParLimits(0, 0.01, 10);
    f->SetParLimits(1, 1.0, 1.0);
    */
    
    int N(0);
    int num(5);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int npoints = fPointArray->GetEntriesFast();
        int ntrack = fMCTrackArray->GetEntriesFast();
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        
        if (ntrack < 2) continue;
        PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(0);
        TVector3 mom(mcTrack->GetMomentum());
        
        //Exclude events generated electron-positron
        bool Exist = false;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
            for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++)
            if (linkIter->GetIndex() == 0) Exist = true;
        }
        if (!Exist) continue;
        
        if (nclusters != 1) continue;
            
        PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(0);
        TVector3 Cent = cluster->where();
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            Int_t DetID = hit->GetDetectorID();
            Double_t E = hit->GetEnergy();
            TVector3 DetPos(emcX[DetID], emcY[DetID], emcZ[DetID]);
            Double_t angle = abs(57.29578*TMath::ATan((DetPos.Cross(vz).Unit().Dot(Cent-DetPos)) / (Cent-DetPos).Dot(DetPos.Unit())));
            Double_t distance = (DetPos - Cent).Mag();
            h2D->Fill(distance,angle,E);
        }
    N++;
    }
    //cout << test/cunt << endl;
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    //c1->cd();
    //h1D->Draw("HIST");
    c2->cd();
    //h2D->Fit(f,"R");
    h2D->Draw("HIST");
    //h2D->Draw("LEGO");
    //f->Draw("SAME");
    return 0;
}
