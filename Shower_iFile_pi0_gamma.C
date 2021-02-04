int Exec(string dir_name, string out_name_min, string out_name_max, Int_t NGamma=2, bool IsSplit=1);
int Shower_iFile_pi0_gamma(string dir_name)
{
    string out_name1 = "doc/" + dir_name + "_pi0_gamma_min.txt";
    string out_name2 = "doc/" + dir_name + "_pi0_gamma_max.txt";
    if( Exec( dir_name, out_name1, out_name2, 3, true) ) return 1;
    return 0;
}

//*****************************************************************************************//
int Exec(string dir_name, string out_name_min, string out_name_max, Int_t NGamma, bool IsSplit){
    //IsSplit: Whether shower separation is required
    //NGamma: Number of photons produced
    ofstream out_min;
    out_min.open(out_name_min,ios::out);
    
    ofstream out_max;
    out_max.open(out_name_max,ios::out);
    
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
    	if (maxEvtNo>=100 && ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        //if (nclusters != 1 || nbumps != 2) continue;
        
        //Get the momentum of each photon
        std::vector<TVector3> Gamma_mom;
        for (int iGamma = 1; iGamma < NGamma; iGamma++) {
            PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(iGamma);
            Gamma_mom.push_back(mcTrack->GetMomentum());
        }
        
        //Exclude events generated electron-positron
        std::map<Int_t, bool> Exist;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
            for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) {
                for (int iGamma = 1; iGamma < NGamma; iGamma++)
                    if (linkIter->GetIndex() == iGamma) Exist[iGamma] = true;
            }
        }
        //if (Exist.size() != 2) continue;
        
        //Calculate the average distance between photons
        Double_t distance(0);
        Int_t Ncunt(0);
        for (int iGamma = 1; iGamma < NGamma-1; iGamma++) {
            for (int jGamma = iGamma+1; jGamma < NGamma; jGamma++) {
                Double_t TheDistance = ((65.0/Gamma_mom[iGamma].Pt())*Gamma_mom[iGamma]-(65.0/Gamma_mom[jGamma].Pt())*Gamma_mom[jGamma]).Mag();
                //Double_t TheDistance = 2 * 65.0 * sin(Gamma_mom[iGamma].Angle(Gamma_mom[jGamma])/2.0);
                distance += TheDistance;
                Ncunt++;
            }
        }
        distance /= Ncunt;

        //Get the true energy of each shower
        std::map<Int_t, Double_t> truth_E;
        for (int iGamma = 1; iGamma < NGamma; iGamma++) truth_E[iGamma] = 0.0;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t>  dep = hit->GetDepositedEnergyMap();
            std::map<Int_t, Double_t>::iterator ptr;
            for ( ptr = dep.begin(); ptr != dep.end(); ptr++) {
	    //cout << ptr->first << endl;
                for (int iGamma = 1; iGamma < NGamma; iGamma++)
                    if (ptr->first == iGamma) truth_E[iGamma] += ptr->second;
            }
        }
        
        //Match bump for each photon
        std::vector<Int_t> match;
        for (int iGamma = 1; iGamma < NGamma; iGamma++) {
            Double_t min_d(99999);
            Int_t index(-1);
            for (int i = 0; i < nbumps; i++) {
                PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(i);
                TVector3 pos = Bump->where();
                Double_t d = pos.Mag()*sin(Gamma_mom[iGamma-1].Angle(pos));
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
        //if (IsSplit && result) continue;
        
        //if (distance > 5) continue;
        //Calculate the error of energy and position
        Double_t bump_E_min(-1.0), bump_E_max(-1.0);
        Double_t bump_E[2];
        for (int iGamma = 1; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(match[iGamma-1]);
            bump_E[iGamma-1] = Bump->energy();
        }
        //cout << bump_E[0] << " " << truth_E[0] << endl;
        //cout << bump_E[0] << " " << truth_E[1] << endl;
        out_min << bump_E[0] << " " << truth_E[1] << endl;
        out_max << bump_E[0] << " " << truth_E[1] << endl;
        //if (bump_E[0] < bump_E[1]) {
          //  out_min << bump_E[0] << " " << truth_E[1] << endl;
          //  out_max << bump_E[1] << " " << truth_E[2] << endl;
        //}else{
          //  out_min << bump_E[1] << " " << truth_E[2] << endl;
           // out_max << bump_E[0] << " " << truth_E[1] << endl;
        //}
        N++;
    }
    out_min.close();
    out_max.close();
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
