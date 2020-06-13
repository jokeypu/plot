int Energy_Compare()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/new1/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/new1/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    TClonesArray* fHitArray = new TClonesArray("PndEmcHit");
    t->SetBranchAddress("EmcHit",&fHitArray);
    if (!fHitArray) return -1;
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fBumpArray = (TClonesArray*) ioman->GetObject("EmcBump");
    if (!fBumpArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    TClonesArray* fSharedDigiArray = (TClonesArray*) ioman->GetObject("EmcSharedDigi");
    if (!fSharedDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    int bin1(100),bin2(100);
    float tx(800),ty(600);
    double xmin(0),xmax(3),ymin(0),ymax(1.2);
    
    TCanvas* c1=new TCanvas("PANDA","MCTruth",tx,ty);
    //TCanvas* c2=new TCanvas("PANDA","MCTruth",tx,ty);
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
    
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",bin1,xmin,xmax,bin2,ymin,ymax);
    
    std::map<Int_t, Double_t>::iterator it;
    
    int aa=2;
    int Nmiss = 0;
    int Nh1(0),Nh2(0),N11(0);
    //for (Int_t ievt = aa; ievt < aa+1; ievt++) {
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
    //for (Int_t ievt = 0; ievt < 100; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        int noverlap_M = 0;
        int nhit = fHitArray->GetEntriesFast();
        for (int n = 0; n < nhit; n++){
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(n);
            std::map<Int_t, Double_t> ds = hit->GetMcSourceEnergy();
            if ( ds.size() > noverlap_M ) noverlap_M = ds.size() - 1;
        }
        int nbump = fBumpArray->GetEntriesFast();
        int ncluster = fClusterArray->GetEntriesFast();
        int noverlap_B = -1;
        if ( nbump == 2 &&  ncluster == 1 ) noverlap_B = 1;
        else if ( nbump == ncluster ) noverlap_B = 0;
        if ( noverlap_M == 1 && noverlap_B == 1 ) {
            std::map<Int_t, bool> overlap_list;
            std::map<Int_t, std::vector<Double_t>> B_overlap_energy;
            std::map<Int_t, std::vector<Double_t>> M_overlap_energy;
            std::map<Int_t, Double_t> B_others_energy;
            std::map<Int_t, Double_t> M_others_energy;
            int nShared = fSharedDigiArray->GetEntriesFast();
            for (int i=0; i < nShared-1; i++){
                for (int j=i+1; j < nShared; j++){
                    PndEmcDigi* digi1 = (PndEmcDigi*)fSharedDigiArray->At(i);
                    PndEmcDigi* digi2 = (PndEmcDigi*)fSharedDigiArray->At(j);
                    if ( digi1->GetDetectorId() == digi2->GetDetectorId() ) {
                        overlap_list[i] = 1;
                        overlap_list[j] = 1;
                        B_overlap_energy[digi1->GetDetectorId()].push_back(digi1->GetEnergy());
                        B_overlap_energy[digi1->GetDetectorId()].push_back(digi2->GetEnergy());
                    }
                }
            }
            for (int i=0; i < nShared; i++){
                if ( overlap_list.find(i) == overlap_list.end() ) {
                    PndEmcDigi* digi = (PndEmcDigi*)fSharedDigiArray->At(i);
                    B_others_energy[digi->GetDetectorId()] = digi->GetEnergy();
                }
            }
            
            for (int n =0; n < nhit; n++){
                PndEmcHit* hit = (PndEmcHit*)fHitArray->At(n);
                std::map<Int_t, Double_t> ds = hit->GetMcSourceEnergy();
                if ( ds.size() == 2 ) {
                    for( it=ds.begin(); it!=ds.end(); ++it){
                        M_overlap_energy[hit->GetDetectorID()].push_back(it->second);
                    }
                }else if ( ds.size() == 1 ) M_others_energy[hit->GetDetectorID()] = (ds.begin())->second;
            }
            std::map<Int_t, std::vector<Double_t>>::iterator p;
            for ( p = M_overlap_energy.begin(); p != M_overlap_energy.end(); p++){
                std::map<Int_t, std::vector<Double_t>>::iterator finder = B_overlap_energy.find(p->first);
                Double_t M_energy1 = (p->second)[0];
                Double_t M_energy2 = (p->second)[1];
                Double_t B_energy1(0);
                Double_t B_energy2(0);
                if ( finder != B_overlap_energy.end() ){
                    B_energy1 = (finder->second)[0];
                    B_energy2 = (finder->second)[2];
                }else continue;
                if(B_energy1/M_energy1 > B_energy2/M_energy1){
                    //histxy->Fill(M_energy1,B_energy1/M_energy1);
                    //histxy->Fill(M_energy2,B_energy2/M_energy2);
                    histxy->Fill(M_energy1,B_energy1);
                    histxy->Fill(M_energy2,B_energy2);
                    Nh1+=2;
                }else{
                    //histxy->Fill(M_energy1,B_energy2/M_energy1);
                    //histxy->Fill(M_energy2,B_energy1/M_energy2);
                    histxy->Fill(M_energy1,B_energy2);
                    histxy->Fill(M_energy2,B_energy1);
                    Nh1+=2;
                }
            }
            std::map<Int_t, Double_t>::iterator ptr;
            for ( ptr = M_others_energy.begin(); ptr != M_others_energy.end(); ptr++){
                std::map<Int_t, Double_t>::iterator fr = B_others_energy.find(ptr->first);
                Double_t M_energy = (ptr->second);
                Double_t B_energy(0);
                if ( fr != B_others_energy.end() ){
                    B_energy = (fr->second);
                }else continue;
                //histxy1->Fill(M_energy,B_energy/M_energy);
                histxy1->Fill(M_energy,B_energy);
                Nh2++;
            }
            
        }else if ( noverlap_M == 0 && noverlap_B == 0 ) {
            N11++;
        }else if ( noverlap_M == 1 && noverlap_B == 0 ) Nmiss++;
        else continue;
    }
    
    histxy->GetXaxis()->SetTitle("Energy");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->GetXaxis()->CenterTitle();
    histxy->GetYaxis()->CenterTitle();
    
    histxy->Draw("SCAT");
    
    histxy1->GetXaxis()->SetTitle("Energy");
    histxy1->GetYaxis()->SetTitle("#phi");
    histxy1->GetXaxis()->CenterTitle();
    histxy1->GetYaxis()->CenterTitle();
    
   // histxy1->Draw("SCAT");
    
    cout << "miss: " << Nmiss << endl;
    cout << "Nhist1: " << Nh1 << endl;
    cout << "Nhist2: " << Nh2 << endl;
    cout << "N11: " << N11 << endl;
    return 0;
}

