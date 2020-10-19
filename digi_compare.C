int Exec(TString dir_name, TH2D* h2D, Int_t NGamma=2, bool IsSplit=1);
int digi_compare()
{
    float tx(800),ty(600);
    TCanvas* c1=new TCanvas("PANDA1","c1",tx,ty);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(1);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",200,0,1,200,0,1.5);
    h2D->GetYaxis()->SetTitle("#daltaw");
    h2D->GetXaxis()->SetTitle("#daltaE");
    h2D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 1);
    
    TH1D* h1D = new TH1D("Hist1_1","h1_1", 200, 0, 1);
    h1D->SetLineColor(kBlue);
    h1D->SetLineWidth(2);
    h1D->GetXaxis()->SetTitle("Energy");
    h1D->GetYaxis()->SetTitle("Entries");
    h1D->GetXaxis()->CenterTitle();
    h1D->GetYaxis()->CenterTitle();
    
    if( Exec( "Gamma_tow_non_1G", h2D, 2, true) ) return 1;

    c1->cd();
    h2D->Draw();
   
    //TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    //leg->AddEntry(h1D1,"Bump Energy old" , "L");
    //leg->AddEntry(h1D2,"Bump Energy new", "L");
    //leg->AddEntry(h1D3,"Cluster Energy 1Gamma", "L");
    //leg->Draw();
    return 0;
}

//*****************************************************************************************//
int Exec(TString dir_name, TH2D* h2D, Int_t NGamma, bool IsSplit){
    //IsSplit: Whether shower separation is required
    //NGamma: Number of photons produced
    
    TString file_path_sim = "../data/"+dir_name+"/evtcomplete_sim.root";
    TString file_path_digi = "../data/"+dir_name+"/evtcomplete_digi.root";
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile(file_path_sim);
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TClonesArray* fMCTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMCTrackArray) return -1;
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    
    TFile* f = new TFile(file_path_digi);
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray);
    if (!fBumpArray) return -1;
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    TClonesArray* fSharedDigiArray = new TClonesArray("PndEmcSharedDigi");
    t->SetBranchAddress("EmcSharedDigi",&fSharedDigiArray);
    if (!fSharedDigiArray) return -1;
    
    int N(0);
    //Int_t maxEvtNo = ioman->CheckMaxEventNo();
    Int_t maxEvtNo = t->GetEntries();
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
    	if (ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        
        //Get the momentum of each photon
        std::vector<TVector3> Gamma_mom;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(iGamma);
            Gamma_mom.push_back(mcTrack->GetMomentum());
        }
        
        //Exclude events generated electron-positron
        std::map<Int_t, bool> Exist;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
            for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) {
                for (int iGamma = 0; iGamma < NGamma; iGamma++)
                    if (linkIter->GetIndex() == iGamma) Exist[iGamma] = true;
            }
        }
        if (Exist.size() != NGamma) continue;
        if (nclusters!=1) continue;

        //Match bump for each photon
        std::vector<Int_t> match;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            Double_t min_d(99999);
            Int_t index(-1);
            for (int i = 0; i < nbumps; i++) {
                PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(i);
                TVector3 pos = Bump->position();
                Double_t d = pos.Mag()*sin(Gamma_mom[iGamma].Angle(pos));
                if (d < min_d) { min_d = d; index = i; }
            }
            if ( index == -1 ) return 1;
            match.push_back(index);
        }
        
        if (nbumps != 2) continue;
        
        std::map<Int_t, Double_t> ww0;
        std::map<Int_t, Double_t> ww1;
        std::map<Int_t, Double_t> EE0;
        std::map<Int_t, Double_t> EE1;
        std::map<Int_t, bool> detid;
        for (int n = 0; n < nbumps ; n++){
            PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(n);
            std::vector<Int_t> list = bump->DigiList();
            for (int i=0; i < list.size(); i++){
                PndEmcSharedDigi* digi = (PndEmcSharedDigi*)fSharedDigiArray->At(list[i]);
                if (n == 0) {ww0[digi->GetDetectorId()] = digi->weight(); EE0[digi->GetDetectorId()] = digi->GetEnergy();}
                if (n == 1) {ww1[digi->GetDetectorId()] = digi->weight(); EE1[digi->GetDetectorId()] = digi->GetEnergy();}
                detid[digi->GetDetectorId()] = true;
            }
        }
        
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            int DID = hit->GetDetectorID();
            if (detid.find(DID) != detid.end()){
            std::map<Int_t, Double_t> ds = hit->GetDepositedEnergyMap();
            h2D->Fill(abs(ww0[DID] - ww1[DID]), ds[match[0]] - EE0[DID]);
            h2D->Fill(abs(ww0[DID] - ww1[DID]), ds[match[1]] - EE1[DID]);
            }
        }
        
        N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
