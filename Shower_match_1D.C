int Exec(TString dir_name, TH1D* h1D1, Int_t NGamma=2, bool IsSplit=1);
int Shower_match_1D( TString dir_name="Gamma_tow_1G_o" )
{
    int bin1(100),bin2(200);
    float tx(1200),ty(900);
    double xmin(0),xmax(20),ymin(0),ymax(0.6);
    
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
    
    TH1D* h1D1 = new TH1D("Hist1_1","h1_1", 100, 0, 1.5);
    h1D1->SetLineColor(kBlue);
    h1D1->SetLineWidth(2);

    TH1D* h1D2 = new TH1D("Hist1_2","h1_2", 100, 0, 1.5);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);

    if( Exec(dir_name, h1D1, 2, true) ) return 1;
    
    dir_name = "Gamma_tow_1G_n";
    if( Exec(dir_name, h1D2, 2, true) ) return 1;
    
    c1->cd();
    h1D1->Draw();
    h1D2->Draw("SAME");
    return 0;
}

//*****************************************************************************************//
int Exec(TString dir_name, TH1D* h1D1, Int_t NGamma, bool IsSplit){
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
        
        //Calculate the average distance between photons
        Double_t distance(0);
        Int_t Ncunt(0);
        for (int iGamma = 0; iGamma < NGamma-1; iGamma++) {
            for (int jGamma = iGamma+1; jGamma < NGamma; jGamma++) {
                Double_t TheDistance = 2 * 65.0 * sin(Gamma_mom[iGamma].Angle(Gamma_mom[jGamma])/2.0);
                distance += TheDistance;
                Ncunt++;
            }
        }
        distance /= Ncunt;
        
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
        
        //h2D1->Fill(distance,nclusters);
        //h2D2->Fill(distance, nbumps);
        
        //Get the true energy of each shower
        std::map<Int_t, Double_t> truth_E;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) truth_E[iGamma] = 0.0;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t>  dep = hit->GetDepositedEnergyMap();
            std::map<Int_t, Double_t>::iterator ptr;
            for ( ptr = dep.begin(); ptr != dep.end(); ptr++) {
                for (int iGamma = 0; iGamma < NGamma; iGamma++)
                if (ptr->first == iGamma) truth_E[iGamma] += ptr->second;
            }
        }
        
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
        Double_t delta_E(0), delta_pos(0);
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(match[iGamma]);
            Double_t bump_E = Bump->energy();
            TVector3 bump_pos = Bump->position();
            delta_E += (truth_E[iGamma] - bump_E/Nshare[match[iGamma]]) * (truth_E[iGamma] - bump_E/Nshare[match[iGamma]]);
            delta_pos += sin(Gamma_mom[iGamma].Angle(bump_pos)/2.0) * sin(Gamma_mom[iGamma].Angle(bump_pos)/2.0);
	        h1D1->Fill(bump_E);
        }
        delta_E = sqrt(delta_E/NGamma);
        delta_pos = 2.0 * 65.0 * sqrt(delta_pos/NGamma);
        
        //h2D3->Fill(distance, delta_E);
        //h2D4->Fill(distance, delta_pos);
        N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
