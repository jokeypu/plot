const Int_t Mode = 2;
int Exec(string dir_name, TH1D *h, vector<Int_t> NBin, Int_t NGamma=2, bool IsSplit=1);
int NEW_Algorithm_test1(string dir_name)
{
    int bin1(100);
    float tx(800),ty(600);
    //double xmin(0.7),xmax(1.3);
    double xmin(0),xmax(8);
    
    TH1D* h1D1 = new TH1D("Hist1_1","h1_1", bin1, xmin, xmax);
    h1D1->SetLineColor(kRed);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("distance");
    h1D1->GetYaxis()->SetTitle("ratio (%)");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    
    vector<Int_t> NBin;
    NBin.push_back(-1);
    if (Mode != 1) {
        ifstream file;
        file.open("doc/NBin.txt", ios::in);
        string str;
        while (getline(file,str)) {
            Int_t value= atof(str.c_str());
            NBin.push_back(value);
        }
        file.close();
    }
    
    if( Exec( dir_name, h1D1, NBin, 2, true) ) return 1;
    
    if (Mode == 1){
        ofstream file;
        file.open("doc/NBin.txt", ios::out);
        for (int i=1; i<= bin1; i++) file << h1D1->GetBinContent(i) << endl;
        file.close();
    }else{
        TCanvas* c1=new TCanvas("PANDA1","c1",tx,ty);
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
        h1D1->Draw("HIST");
        TLegend * leg1 = new TLegend(0.61,0.72,0.88,0.85);
        leg1->AddEntry(h1D1, "The ratio of successful splitting", "L");
        leg1->Draw("SAME");
    }
    
    return 0;
}

//*****************************************************************************************//
int Exec(string dir_name, TH1D *h, vector<Int_t> NBin, Int_t NGamma, bool IsSplit){
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
        //if (nclusters!=1) continue;
        
        //Calculate the average distance between photons
        Double_t distance(0);
        Int_t Ncunt(0);
        for (int iGamma = 0; iGamma < NGamma-1; iGamma++) {
            for (int jGamma = iGamma+1; jGamma < NGamma; jGamma++) {
                Double_t TheDistance = ((65.0/Gamma_mom[iGamma].Pt())*Gamma_mom[iGamma]-(65.0/Gamma_mom[jGamma].Pt())*Gamma_mom[jGamma]).Mag();
                //Double_t TheDistance = 2 * 65.0 * sin(Gamma_mom[iGamma].Angle(Gamma_mom[jGamma])/2.0);
                distance += TheDistance;
                Ncunt++;
            }
        }
        distance /= Ncunt;
        if (Mode == 1) h->Fill(distance);
        else if (nbumps == 2) h->Fill(distance,100*1.0/NBin[h->FindBin(distance)]);
        N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
