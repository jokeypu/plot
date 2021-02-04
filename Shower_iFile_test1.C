Double_t DD(const TVector3 *DetPos, const TVector3 *Cent, const Double_t par){
    Double_t distance(0), angle(0);
    if (*DetPos != *Cent) {
        TVector3 vz(0, 0, 1);
        TVector3 DetPos_o = (*DetPos) - 3.7*vz;
        TVector3 DetPos_n;
        DetPos_n.SetMagThetaPhi(DetPos_o.Mag(), DetPos_o.Theta(), DetPos_o.Phi()-0.06981317);
        TVector3 ey = DetPos_n.Cross(vz).Unit();
        TVector3 ex = DetPos_n.Cross(ey).Unit();
        Double_t dx = abs((*Cent-*DetPos).Dot(ex));
        Double_t dy = abs((*Cent-*DetPos).Dot(ey));
        distance = sqrt(dx*dx+dy*dy);
        angle = 57.2957*TMath::ATan(dy/dx);
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    }
    return distance;
}
Double_t AA(const TVector3 *DetPos, const TVector3 *Cent, const Double_t par){
    Double_t distance(0), angle(0);
    if (*DetPos != *Cent) {
        TVector3 vz(0, 0, 1);
        TVector3 DetPos_o = (*DetPos) - 3.7*vz;
        TVector3 DetPos_n;
        DetPos_n.SetMagThetaPhi(DetPos_o.Mag(), DetPos_o.Theta(), DetPos_o.Phi()-0.06981317);
        TVector3 ey = DetPos_n.Cross(vz).Unit();
        TVector3 ex = DetPos_n.Cross(ey).Unit();
        Double_t dx = abs((*Cent-*DetPos).Dot(ex));
        Double_t dy = abs((*Cent-*DetPos).Dot(ey));
        distance = sqrt(dx*dx+dy*dy);
        angle = 57.2957*TMath::ATan(dy/dx);
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    }
    return angle;
}
int Exec(string dir_name, string out_name_min, string out_name_max, Int_t NGamma=2, bool IsSplit=1);
int Shower_iFile_test1(string dir_name)
{
    string out_name_min = "doc/" + dir_name + "_test_cent.txt";
    string out_name_max = "doc/" + dir_name + "_test_band.txt";
    if( Exec( dir_name, out_name_min, out_name_max, 2, true) ) return 1;
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
    TClonesArray* fDigiArray = new TClonesArray("PndEmcDigi");
    t->SetBranchAddress("EmcDigi",&fDigiArray);
    if (!fDigiArray) return -1;
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
        int ndigis = fDigiArray->GetEntriesFast();
        
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
        
        //if (distance > 5) continue;
        if (nclusters != 1 || nbumps != 2) continue;
        //Calculate the error of energy and position
        PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(0);
        std::map<Int_t, Int_t> localmax = cluster->LocalMaxMap();
        std::map<Int_t, Int_t>::iterator thelocal;
        
        std::vector<Int_t> v_macth;
        for (thelocal = localmax.begin(); thelocal != localmax.end(); ++thelocal) {
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(thelocal->second);
            TVector3 Seed_pos = digi->where();
            PndEmcBump* Bump0 = (PndEmcBump*)fBumpArray->At(0);
            PndEmcBump* Bump1 = (PndEmcBump*)fBumpArray->At(1);
            TVector3 Cent_pos0 = Bump0->where();
            TVector3 Cent_pos1 = Bump1->where();
            if ((Cent_pos0-Seed_pos).Mag() < (Cent_pos1-Seed_pos).Mag()) v_macth.push_back(0);
            else v_macth.push_back(1);
        }
        if (v_macth[0] == v_macth[1]) {cout << "XXX" <<endl;continue;}
        
        Int_t NN = 0;
        for (thelocal = localmax.begin(); thelocal != localmax.end(); ++thelocal) {
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(thelocal->second);
            double Seed_Energy = digi->GetEnergy();
            TVector3 Seed_pos = digi->where();
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(v_macth[NN]);
            TVector3 Cent_pos = Bump->where();
            Double_t bump_E = Bump->energy();
            Double_t d = DD(&Seed_pos, &Cent_pos, 1.25);
            Double_t A = AA(&Seed_pos, &Cent_pos, 1.25);
            //if (d < 1.1 && fabs(A - 45)<10 ) out_min << bump_E << endl;
            //if (Seed_Energy>0.7) out_min << bump_E << endl;
            //if ((Seed_pos-Cent_pos).Mag()<0.7) out_min << bump_E << endl;
            if (distance < 5 ) out_min << bump_E << endl;
            else out_max << bump_E << endl;
            NN++;
        }
    }
    out_min.close();
    out_max.close();
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
