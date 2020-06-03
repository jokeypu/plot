int ReadMCTrack_cluster()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    int maxEvtNo = ioman->CheckMaxEventNo();
    int nlnks;
    cout << "maxEvtNo " << maxEvtNo << endl;
    for (int ievt = 0; ievt < 10; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nclusters = fClusterArray->GetEntriesFast();
        cout << "\n" << ievt+1 << "\tclusters is "<< nclusters << endl;
        std::map<int, std::vector<int> > Tracks;
        
        for (int i = 0; i < nclusters; i++) {
            PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(i);
            std::vector<FairLink> links = (cluster->GetTrackEntering()).GetSortedMCTracks();
            nlnks = links.size();
            std::vector<int> vect;
            vect.clear();
            for (int k = 0; k < nlnks; k++){
                int idx = links[k].GetIndex();
                vect.push_back(idx);
            }
            Tracks[i] = vect;
        }
        for (int i = 0; i < nclusters; i++) {
            for (int k = 0; k < nlnks; k++){
                int idx = Tracks[i].at(k);
                PndMCTrack* track = (PndMCTrack*)fMCtrackArray->At(idx);
                TVector3 svertex,smomentum;
                int mid = idx;
                while ( mid >=  0) {
                    PndMCTrack* mtrack = (PndMCTrack*)fMCtrackArray->At(mid);
                    if (((mid = mtrack->GetMotherID()) < 0) && (((mtrack->GetStartVertex()).Mag()) <= 0.1)) {
                        cout << "\t*signal* " ;
                        svertex = mtrack->GetStartVertex();
                        smomentum = mtrack->GetMomentum();
                        cout << "StartVertex:" << "("<< svertex.x() << "," << svertex.y() << "," << svertex.z() << ")\t" << "mag:" << svertex.Mag();
                        cout << "\ttheta,phi:" << "(" << smomentum.Theta() << "," << smomentum.Phi() << ")" << endl;
                        k = nlnks;
                    }
                }
            }
        }
    }
    return 0;
}

