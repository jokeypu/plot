int Exec(TH1D* hist, Int_t NGamma, bool IsSplit);
int Shower_instance()
{
    int bin1(100),bin2(200);
    float tx(1200),ty(900);
    double xmin(0.7),xmax(1.3);
    //double xmin(0.3),xmax(0.7);
    //string file_name("doc/standard_5.txt");
    //string file_name("doc/standard.txt");
    string file_name("doc/old.txt");
    
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
    
    TH1D* h1D1 = new TH1D("Hist1_1","h1_1", bin1, xmin, xmax);
    h1D1->SetLineColor(kGray+3);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("Energy");
    h1D1->GetYaxis()->SetTitle("Entries");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    TH1D* h1D2 = new TH1D("Hist1_2","h1_2", bin1, xmin, xmax);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);
    h1D2->GetXaxis()->SetTitle("Energy");
    h1D2->GetYaxis()->SetTitle("Entries");
    h1D2->GetXaxis()->CenterTitle();
    h1D2->GetYaxis()->CenterTitle();

    int MaxNo = Exec( h1D2, 2, true);
    if (  MaxNo == -1 || MaxNo > 38570 ) {
        std::cout << "Error OR MaxNo > 38570 !!!" << std::endl;
        return -1;
    }
    
    string str;
    ifstream file;
    file.open(file_name, ios::in);
    for (int i = 0; i < MaxNo; i++) {
        getline(file,str);
        double value= atof(str.c_str());
        h1D1->Fill(value);
    }
    
    c1->cd();
    h1D2->Draw();
    h1D1->Draw("SAME");
   
    TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    leg->AddEntry(h1D1,"Bump Energy old" , "L");
    leg->AddEntry(h1D2,"Bump Energy new", "L");
    //leg->AddEntry(h1D3,"Cluster Energy 1Gamma", "L");
    leg->Draw();
    return 0;
}

//*****************************************************************************************//
int Exec(TH1D* hist, Int_t NGamma, bool IsSplit){
    //IsSplit: Whether shower separation is required
    //NGamma: Number of photons produced
    
    TString file_path_sim = "../data/evtcomplete_sim.root";
    TString file_path_digi = "../data/evtcomplete_digi.root";
    
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
    for (Int_t ievt = 0; ievt < 0.9*maxEvtNo; ievt++) {
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
            if ( index == -1 ) return -1;
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
        
        if (distance > 5) continue;
        if (nclusters != 1 || nbumps != 2) continue;
        //Calculate the error of energy and position
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(match[iGamma]);
            Double_t bump_E = Bump->energy();
            hist->Fill(bump_E);
            N++;
        }
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return N;
}
