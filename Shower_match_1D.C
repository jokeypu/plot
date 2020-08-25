int Exec(TString dir_name, string out_name, Int_t NGamma=2, bool IsSplit=1);
int Shower_match_1D()
{
    int bin1(100),bin2(200);
    float tx(1200),ty(900);
    double xmin(0),xmax(20),ymin(0),ymax(0.6);
    string out1_name("out1.txt"), out2_name("out2.txt"), out3_name("out3.txt");
    ifstream out1, out2, out3;
    out1.open(out1_name, ios::in);
    out2.open(out2_name, ios::in);
    out3.open(out3_name, ios::in);
    
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
    
    TH1D* h1D1 = new TH1D("Hist1_1","h1_1", bin1, 0.7, 1.3);
    h1D1->SetLineColor(kBlue);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("Energy");
    h1D1->GetYaxis()->SetTitle("Entries");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    TH1D* h1D2 = new TH1D("Hist1_2","h1_2", bin1, 0.7, 1.3);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);
    h1D2->GetXaxis()->SetTitle("Energy");
    h1D2->GetYaxis()->SetTitle("Entries");
    h1D2->GetXaxis()->CenterTitle();
    h1D2->GetYaxis()->CenterTitle();
    
    TH1D* h1D3 = new TH1D("Hist1_3","h1_3", bin1, 0.7, 1.3);
    h1D3->SetLineColor(kGreen);
    h1D3->SetLineWidth(2);
    h1D3->GetXaxis()->SetTitle("Energy");
    h1D3->GetYaxis()->SetTitle("Entries");
    h1D3->GetXaxis()->CenterTitle();
    h1D3->GetYaxis()->CenterTitle();
    
    string str;
    int mode(5); //1,2,3,4,5
    int min(15000);
    if (mode == 1){
        if( Exec( "Gamma_tow_1G_o", out1_name, 2, true) ) return 1;
    }else if (mode == 2) {
        if( Exec( "Gamma_tow_1G_n", out2_name, 2, true) ) return 1;
    }else if (mode == 3) {
        if( Exec( "Gamma_one_1G", out3_name, 1, true) ) return 1;
    }else if (mode == 4) {
        int cunt;
        cunt = 0;
        while (!out1.eof()) {
            getline(out1,str);
            double value= atof(str.c_str());
            h1D1->Fill(value);
            cunt++;
        }
        cout << "N1:" << cunt << endl;
        cunt = 0;
        while (!out2.eof()) {
            getline(out2,str);
            double value= atof(str.c_str());
            h1D2->Fill(value);
            cunt++;
        }
        cout << "N2:" << cunt << endl;
        cunt = 0;
        while (!out3.eof()) {
            getline(out3,str);
            double value= atof(str.c_str());
            h1D3->Fill(value);
            cunt++;
        }
        cout << "N3:" << cunt << endl;
        cunt = 0;
    }else if (mode == 5) {
        int cunt;
        cunt = 0;
        for (int i= 0;i < min; i++) {
            getline(out1,str);
            double value= atof(str.c_str());
            h1D1->Fill(value);
        }
        for (int i= 0;i < min; i++) {
            getline(out2,str);
            double value= atof(str.c_str());
            h1D2->Fill(value);
        }
        for (int i= 0;i < min; i++) {
            getline(out3,str);
            double value= atof(str.c_str());
            h1D3->Fill(value);
        }
    }
    out1.close();
    out2.close();
    out3.close();
    c1->cd();
    h1D3->Draw();
    h1D1->Draw("SAME");
    h1D2->Draw("SAME");
    TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    leg->AddEntry(h1D1,"Bump Energy old" , "L");
    leg->AddEntry(h1D2,"Bump Energy new", "L");
    leg->AddEntry(h1D3,"Cluster Energy 1Gamma", "L");
    leg->Draw();
    return 0;
}

//*****************************************************************************************//
int Exec(TString dir_name, string out_name, Int_t NGamma, bool IsSplit){
    //IsSplit: Whether shower separation is required
    //NGamma: Number of photons produced
    ofstream out;
    out.open(out_name,ios::out);
    
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
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
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
        
        //Count the number of times that Bump is shared by the MCtrack
        std::map<Int_t, Int_t> Nshare;
        for (int i = 0; i < match.size(); i++) Nshare[match[i]]++;
        
        //Whether to skip events where the shower is not separated
        bool result(false);
        std::map<Int_t, Int_t>::iterator it;
        for ( it = Nshare.begin(); it != Nshare.end(); it++) if (it->second != 1) result = true;
        if (IsSplit && result) continue;
        
        //Calculate the error of energy and position
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(match[iGamma]);
            Double_t bump_E = Bump->energy();
            out << bump_E << endl;
        }
        N++;
    }
    out.close();
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
