int MC_backup7()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/Npi0_3GeV/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../../data/Npi0_3GeV/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    TClonesArray* fHitArray = new TClonesArray("PndEmcHit");
    t->SetBranchAddress("EmcHit",&fHitArray);
    if (!fHitArray) return -1;
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    int bin1(100),bin2(100);
    double xmin(60),xmax(90),ymin(-15),ymax(15);
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    
    TCanvas* c1=new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->Draw();
    
    std::map<Int_t, Double_t>::iterator it;
    
    int aa=7;
    for (Int_t ievt = aa; ievt < aa+1; ievt++) {
    //for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        
        int ncluster = fClusterArray->GetEntriesFast();
        
        cout << ncluster << "*******" << endl;
	for (int n = 0; n < ncluster ; n++){
        PndEmcCluster* cluster = (PndEmcCluster*)fClusterArray->At(n);
        std::vector<Int_t> list = cluster->DigiList();
        for (int i=0; i < list.size(); i++){
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(list[i]);
            double E = digi->GetEnergy();
            double theta = digi->GetTheta();
            double phi = digi->GetPhi();
            double S = 1.2 * pow(E, 0.3);
            theta = theta * ( 180.0 / 3.1416);
            phi = phi * ( 180.0 / 3.1416);
            TWbox *twb = new TWbox(theta-S,phi-S,theta+S,phi+S,66+5*n,-2,-3);
            twb->Draw("SAME");
            
        }
	}
        /*for ( Int_t i = 0; i < ncluster; i++ ){
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t> ds = hit->GetMcSourceEnergy();
            cout << "hit: " << i << endl;
            double E = hit->GetEnergy();
            double theta = hit->GetTheta();
            double phi = hit->GetPhi();

            for( it=ds.begin(); it!=ds.end(); ++it){
                cout << it->first << endl;
                cout << it->second << endl;
                double S = 1.2 * pow(it->second, 0.3);
                TWbox *twb = new TWbox(theta-S,phi-S,theta+S,phi+S,46+4*((it->first)%6),-2,-3);
                twb->Draw("SAME");
            }
            
            
        }*/
    }
    return 0;
}

