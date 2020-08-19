int ReadMCTrack_hit()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/evtcomplete_sim.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    TClonesArray* fMCtrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMCtrackArray) return -1;
    int maxEvtNo = ioman->CheckMaxEventNo();
    int nlnks;
    cout << "maxEvtNo " << maxEvtNo << endl;
    
    for (int ievt = 0; ievt < 3; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        int nhits = fHitArray->GetEntriesFast();
        cout << "\n" << ievt+1 << "\thits is " << nhits << endl;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::vector<FairLink> links = (hit->GetTrackEntering()).GetSortedMCTracks();
            nlnks = links.size();
	    cout << "-------- hit: " << i << " --------" << endl;
            for (int k = 0; k < nlnks; k++){
                int idx = links[k].GetIndex();
                PndMCTrack* track = (PndMCTrack*)fMCtrackArray->At(idx);
                TVector3 svertex, smomentum;
                int mid = idx;
                while ( mid >=  0) {
                    PndMCTrack* mtrack = (PndMCTrack*)fMCtrackArray->At(mid);
                    if (((mid = mtrack->GetMotherID())<0) && (((mtrack->GetStartVertex()).Mag()) <=  0.1)) {
                        cout << "\t*signal* " ;
                        svertex = mtrack->GetStartVertex();
                        smomentum = mtrack->GetMomentum();
                        cout << "StartVertex:" << "(" << svertex.x() << "," << svertex.y() << "," << svertex.z() << ")\t" << "mag:" << svertex.Mag();
                        cout << "\ttheta,phi:" << "(" << smomentum.Theta() << "," << smomentum.Phi() << ")\t";
                        k = nlnks;
                    }
                }
		cout << "id:" << idx << endl;
            }
	    cout << endl;
        }
    }
    return 0;
}

