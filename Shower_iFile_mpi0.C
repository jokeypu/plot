int Exec(string dir_name, string out_name, Int_t NGamma=2, bool IsSplit=1);
int Shower_iFile_mpi0(string dir_name)
{
    string out_name = "doc/" + dir_name + "_mpi0.txt";
    if( Exec( dir_name, out_name, 2, true) ) return 1;
    return 0;
}

//*****************************************************************************************//
int Exec(string dir_name, string out_name, Int_t NGamma, bool IsSplit){
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
    //Int_t maxEvtNo = ioman->CheckMaxEventNo();
    Int_t maxEvtNo = t->GetEntries();
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
    	if (ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        
        //if (distance > 5) continue;
        if (nclusters != 1 || nbumps != 2) continue;
        //Calculate the error of energy and position
<<<<<<< HEAD
        Double_t PX = PY = PZ = PE = 0.0;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(match[iGamma]);
=======
        Double_t PX, PY, PZ, PE;
        PX = PY = PZ = PE = 0.0;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(iGamma);
>>>>>>> 560032398e62720358ab05de2eaee2c7b602b026
            Double_t bump_E = Bump->energy();
            TVector3 bump_mom = Bump->where();
            bump_mom.SetMag(bump_E);
            PX += bump_mom.X();
            PY += bump_mom.Y();
            PZ += bump_mom.Z();
            PE += bump_E;
            //out << bump_mom.X() << " " << bump_mom.Y() << " " << bump_mom.Z() << " " << bump_E << endl;
        }
        out << sqrt(PE*PE - PX*PX - PY*PY - PZ*PZ) << endl;
        N++;
    }
    out.close();
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
