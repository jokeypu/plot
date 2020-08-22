int Shower_match_v1( TString dir_name="Gamma_tow_1G" )
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    double xmin(0),xmax(6),ymin(0),ymax(180);
    
    TCanvas* c1=new TCanvas("PANDA1","c1",tx,ty);
    TCanvas* c2=new TCanvas("PANDA2","c2",tx,ty);
    TCanvas* c3=new TCanvas("PANDA3","c3",tx,ty);
    TCanvas* c4=new TCanvas("PANDA4","c4",tx,ty);
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
    
    TH2D* h2D1 = new TH2D("Hist1","h1",100,0,20, 100,0,3.5);
    h2D1->SetMarkerStyle(7);
    h2D1->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D1->GetXaxis()->SetTitle("distance");
    h2D1->GetYaxis()->SetTitle("N_{cluster}");
    h2D1->GetXaxis()->CenterTitle();
    h2D1->GetYaxis()->CenterTitle();
    
    TH2D* h2D2 = new TH2D("Hist2","h2",100,0,20, 100,0,3.5);
    h2D2->SetMarkerStyle(7);
    h2D2->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D2->GetXaxis()->SetTitle("distance");
    h2D2->GetYaxis()->SetTitle("N_{bump}");
    h2D2->GetXaxis()->CenterTitle();
    h2D2->GetYaxis()->CenterTitle();
    
    TH2D* h2D3 = new TH2D("Hist3","h3",100,0,20, 100,0,5);
    h2D3->SetMarkerStyle(7);
    h2D3->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D3->GetXaxis()->SetTitle("distance");
    h2D3->GetYaxis()->SetTitle("#DeltaE");
    h2D3->GetXaxis()->CenterTitle();
    h2D3->GetYaxis()->CenterTitle();
    
    TH2D* h2D4 = new TH2D("Hist4","h4",100,0,20, 100,0,3);
    h2D4->SetMarkerStyle(7);
    h2D4->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D4->GetXaxis()->SetTitle("distance");
    h2D4->GetYaxis()->SetTitle("#Deltad");
    h2D4->GetXaxis()->CenterTitle();
    h2D4->GetYaxis()->CenterTitle();
    
    //********************************************************//
    
    TString file_path_sim = "../data/"+dir_name+"/evtcomplete_sim.root";
    TString file_path_digi = "../data/"+dir_name+"/evtcomplete_digi.root";
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile(file_path_sim);
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
    
    TFile* f = new TFile(file_path_digi);
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray);
    if (!fBumpArray) return -1;
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    
    int N(0), NGamma(2); // NGamma: Number of photons produced
    bool IsSplit(1); // Whether shower separation is required
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
        
        h2D1->Fill(distance,nclusters);
        h2D2->Fill(distance, nbumps);
        
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
        }
        delta_E = sqrt(delta_E/NGamma);
        delta_pos = 2.0 * 65.0 * sqrt(delta_pos/NGamma);
        
        h2D3->Fill(distance, delta_E);
        h2D4->Fill(distance, delta_pos);
        N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    
    //********************************************************//
    
    c1->cd();
    h2D1->Draw("SCAT");
    c2->cd();
    h2D2->Draw("SCAT");
    c3->cd();
    h2D3->Draw("SCAT");
    c4->cd();
    h2D4->Draw("SCAT");
    return 0;
}
