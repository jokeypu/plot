int Exec(string dir_name, TH1D *h_pass, TH1D *h_Get, Int_t NGamma=2, bool IsSplit=1);
int NEW_DigiVsHit()
{
    string dir_name = "WorkData_2Gamma_A7_E1.0_OR";
    
    TString file_path_sim = "../data/"+dir_name+"/evtcomplete_sim.root";
    TString file_path_digi = "../data/"+dir_name+"/evtcomplete_digi.root";
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile(file_path_digi);
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    
    TFile* f = new TFile(file_path_sim);
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fHitArray = new TClonesArray("PndEmcHit");
    t->SetBranchAddress("EmcHit",&fHitArray);
    if (!fHitArray) return -1;
    
    TCanvas* c1=new TCanvas("PANDA","Point",1200,900);
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
    
    TH1D *h1 = new TH1D("h1","hist1",100,1.7,2.1);
    h1->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Energy");
    h1->GetYaxis()->SetTitle("Entries");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();

    TH1D *h2 = new TH1D("h2","hist2",100,1.7,2.1);
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(2);
    h2->GetXaxis()->SetTitle("Energy");
    h2->GetYaxis()->SetTitle("Entries");
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();
    
    TH1D *h3 = new TH1D("h3","hist3",100,1.7,2.1);
    h3->SetLineColor(kBlack);
    h3->SetLineWidth(2);
    h3->GetXaxis()->SetTitle("Energy");
    h3->GetYaxis()->SetTitle("Entries");
    h3->GetXaxis()->CenterTitle();
    h3->GetYaxis()->CenterTitle();
    
    TH1D *h4 = new TH1D("h4","hist4",100,1.7,2.1);
    h4->SetLineColor(kGreen);
    h4->SetLineWidth(2);
    h4->GetXaxis()->SetTitle("Energy");
    h4->GetYaxis()->SetTitle("Entries");
    h4->GetXaxis()->CenterTitle();
    h4->GetYaxis()->CenterTitle();
    
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        if (ievt%(maxEvtNo/100)==0) cout << "Step1: " << 100 * (int)ievt/maxEvtNo << "%" << endl;
        int nhits = fHitArray->GetEntriesFast();
        int ndigis = fDigiArray->GetEntriesFast();
        
        Double_t Energy_hit = 0.0;
        Double_t Energy_digi = 0.0;
        
        for (int i = 0; i < nhits; i++) {
            PndEmcHit *hit = (PndEmcHit*) fHitArray->At(i);
            Energy_hit += hit->GetEnergy();
        }
        for (int i = 0; i < ndigis; i++) {
            PndEmcDigi *digi = (PndEmcDigi*) fDigiArray->At(i);
            Energy_digi += digi->GetEnergy();
        }
        h1->Fill(Energy_hit);
        h2->Fill(Energy_digi);
        
    }
    
    h1->Draw();
    h2->Draw("SAME");
    
    if( Exec( dir_name, h3, h4, 2, true) ) return 1;
    
    h3->Add(h4);
    
    h3->Draw("SAME");
    //h4->Draw("SAME");
    
    TLegend * leg1 = new TLegend(0.61,0.72,0.88,0.85);
    leg1->AddEntry(h1, "Hit sum", "L");
    leg1->AddEntry(h2, "Digi sum", "L");
    leg1->AddEntry(h3, "Bump pass", "L");
    leg1->AddEntry(h4, "Bump Get", "L");
    leg1->Draw("SAME");
    
    return 0;
}

//*****************************************************************************************//
int Exec(string dir_name, TH1D *h_pass, TH1D *h_Get, Int_t NGamma, bool IsSplit){
    //IsSplit: Whether shower separation is required
    //NGamma: Number of photons produced
    
    TString file_path_sim = "../data/"+dir_name+"/evtcomplete_sim.root";
    TString file_path_digi = "../data/"+dir_name+"/evtcomplete_digi.root";
    
    //FairRunAna *fRun = new FairRunAna();
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
        if (ievt%(maxEvtNo/100)==0) cout <<  "Step2: " << 100 * (int)ievt/maxEvtNo << "%" << endl;
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
        if (Exist.size() != NGamma) {
            Double_t E_pass = 0.0;
            for (int i = 0; i < nclusters; i++) {
                PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(i);
                E_pass += cluster->energy();
            }
            h_pass->Fill(E_pass);
            continue;
        }
        if (nclusters!=1) {
            Double_t E_pass = 0.0;
            for (int i = 0; i < nclusters; i++) {
                PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(i);
                E_pass += cluster->energy();
            }
            h_pass->Fill(E_pass);
            continue;
        }
        
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
        if (IsSplit && result) {
            Double_t E_pass = 0.0;
            for (int i = 0; i < nclusters; i++) {
                PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(i);
                E_pass += cluster->energy();
            }
            h_pass->Fill(E_pass);
            continue;
        }
        
        //if (distance > 5) continue;
        if (nclusters != 1 || nbumps != 2) {
            Double_t E_pass = 0.0;
            for (int i = 0; i < nclusters; i++) {
                PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(i);
                E_pass += cluster->energy();
            }
            h_pass->Fill(E_pass);
            continue;
        }
        
        //Calculate the error of energy and position
        Double_t bump_E = 0.0;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(match[iGamma]);
            bump_E += Bump->energy();
            N++;
        }
        h_Get->Fill(bump_E);
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}

