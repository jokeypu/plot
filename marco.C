int marco(){
    TString file_path_digi = "../data/evtcomplete_digi.root";
    TString file_path_sim = "../data/evtcomplete_sim.root";
    
    TFile* f = new TFile(file_path_digi);
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray);
    if (!fBumpArray) return -1;
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile(file_path_sim);
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    
    //Int_t maxEvtNo = ioman->CheckMaxEventNo();
    Int_t maxEvtNo = t->GetEntries();
    for (Int_t ievt = 9; ievt < 10; ievt++) {
    	ioman->ReadEvent(ievt);
        t->GetEntry(ievt);
        int nclusters = fClusterArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        //cout << "ievt:" << ievt << " nbump:" << nbumps << ", ncluster:" << nclusters << endl;
        
        std::map<Int_t, Double_t>::iterator it;
        int nhits = fHitArray->GetEntriesFast();
        for ( Int_t i = 0; i < nhits; i++ ){
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t> ds = hit->GetDepositedEnergyMap();
            if (hit->GetDetectorID() == 108130002) {
                for ( it = ds.begin(); it != ds.end(); it++ )
                cout << it->first << ", "<< (it->second)/(hit->GetEnergy()) << endl;
            }
        }
    }
    //cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
