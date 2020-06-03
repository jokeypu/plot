int MC_backup4()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/Gamma/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/Gamma/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    PndEmcMapper::Init(1);
    TFile *parfile = new TFile("../../data/Gamma/evtcomplete_par.root");
    parfile->Get("FairGeoParSet");
    PndEmcStructure *fEmcStr = PndEmcStructure::Instance();
    PndEmcMapper *fMapper = PndEmcMapper::Instance();
    
    TClonesArray* fPointArray = new TClonesArray("PndEmcPoint");
    t->SetBranchAddress("EmcPoint",&fPointArray);
    if (!fPointArray) return -1;
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    int bin1(40),bin2(40);
    double xmin(0.9),xmax(1.8),ymin(-0.6),ymax(0.6);
    //TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",bin1,xmin,xmax,bin2,ymin,ymax);
    TH1D* h1 = new TH1D("h1", "hit", 500, 0.0, 1.2);
    
    double fDistanceMotherCutF(50.0);
    double fDistanceMotherCutB(80.0);
    for (Int_t ievt = 0; ievt < 1; ievt++) {
    //for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        int DetId;
        
        std::vector<Int_t> isAdded;
        std::vector<Int_t> smid;
        std::map<int, std::vector<int> > shower;
        std::map<int, std::vector<int> > mod;
        std::vector<std::map<int, double> > ShowerHitEnergy;
        std::vector<std::map<int, std::vector<int> > > fPointMatch;
        std::vector<std::vector<int> > fMatch;
        std::vector<std::vector<double> > fEnergy;

        isAdded.clear();
        smid.clear();
        shower.clear();
        mod.clear();
        ShowerHitEnergy.clear();
        fPointMatch.clear();
        fMatch.clear();
        fEnergy.clear();
        
        std::map<int, std::vector<int> >::iterator it;
        std::map<int, double>::iterator pp;
        
        int npoints = fPointArray->GetEntriesFast();
        for (Int_t i=0;i<npoints;i++){
            isAdded.push_back(i);
        }
        for (Int_t j=0;j<npoints;j++){
            if (isAdded[j] == -1) continue;
            FairMCPoint* point = (FairMCPoint* ) fPointArray->At(j);
            int mid;
            int id = point->GetTrackID();
            double rho;
            do {
                PndMCTrack* track = (PndMCTrack*)fMCtrackArray->At(id);
                TVector3 StartVertex = track->GetStartVertex();
                rho = sqrt(StartVertex.x()*StartVertex.x()+StartVertex.y()*StartVertex.y());
                mid = id;
                id = track->GetMotherID();
            } while ( rho > fDistanceMotherCutF && rho < fDistanceMotherCutB );
                shower[mid].push_back(j);
        }
        
        int k=0;
        for( it=shower.begin(); it!=shower.end(); ++it){
            smid.push_back(it->first);
            mod[k]=it->second;
            k++;
        }
        
        for (int i=0;i<mod.size();i++){
            std::map<int, double> he;
            std::map<int, std::vector<int> > mth;
            for (int j=0;j<mod[i].size();j++){
                FairMCPoint* point = (FairMCPoint* ) fPointArray->At(mod[i].at(j));
                DetId = point->GetDetectorID();
                if ( he.find(DetId) != he.end() ){
                    he[DetId] += point->GetEnergyLoss();
                }else{
                    he[DetId] = point->GetEnergyLoss();
                }
                mth[DetId].push_back(mod[i].at(j));
            }
            ShowerHitEnergy.push_back(he);
            fPointMatch.push_back(mth);
        }
       
        for (int i=0;i<ShowerHitEnergy.size();i++){
            std::vector<int> idm;
            std::vector<double> ide;
            for( pp=ShowerHitEnergy[i].begin(); pp!=ShowerHitEnergy[i].end(); ++pp){
                idm.push_back(pp->first);
                ide.push_back(pp->second);
            }
            fMatch.push_back(idm);
            fEnergy.push_back(ide);
        }
        
        //******************************************************************************//
        double E=0;
        for (int i=0;i<fEnergy.size();i++){
            for (int j=0;j<fEnergy[i].size();j++){
                E+=fEnergy[i].at(j);
            }
        }
        
        typedef std::map<Int_t, Float_t> mapper;
        mapper emcX = fEmcStr->GetEmcX();
        mapper emcY = fEmcStr->GetEmcY();
        mapper emcZ = fEmcStr->GetEmcZ();
        
        int cc=0;
        for (int i=0;i<fMatch.size();i++){
            TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
            PndMCTrack* track = (PndMCTrack*)fMCtrackArray->At(smid[i]);
            bool islow=false;
            if ( track->Get4Momentum().E() < 0.00 ) { islow = true; cc--; }
            cc++;
            for (int j=0;j<fMatch[i].size();j++){
                int idx = fMatch[i].at(j);
                double e = fEnergy[i].at(j);
                double S = 7000 * pow(e, 0.4);
                TVector3 position(emcX[idx], emcY[idx], emcZ[idx]);
                double theta = position.Theta();
                double phi = position.Phi();
                for (int k = 0; k < S; k++) {
                    if ( e > 0.00 ){
                        if ( !islow ) histxy->Fill(theta,phi);
                        else histxy1->Fill(theta,phi);
                    }
                }
                
            }
            if ( !islow ){
            histxy->GetXaxis()->SetTitle("#theta");
            histxy->GetYaxis()->SetTitle("#phi");
            histxy->SetFillColor(30+4*(i%30));
            histxy->Draw("SAMEBOX");
            }
            
        }
        cout << "Number of Clusters: " << fMatch.size() << endl;
        //h1->Fill(E);
    }
    
    //TCanvas* c1=new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    /*
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->SetFillColor(46);
    histxy->Draw("BOX");
    */
    histxy1->SetFillColor(13);
    //histxy1->Draw("SAMEBOX");
    
    h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Energy(GeV)");
    h1->GetYaxis()->SetTitle("Entries");
    h1->SetLineColor(kRed);
    //h1->Draw();
    
    return 0;
}

