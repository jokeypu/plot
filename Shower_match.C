int Shower_match( TString dir_name="Gamma_tow_1G" )
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    double xmin(0),xmax(6),ymin(0),ymax(180);
    //******************************************//
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/"+dir_name+"/evtcomplete_sim.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TClonesArray* fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
    if (!fPointArray) return -1;
    TClonesArray* fMCTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMCTrackArray) return -1;
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    TFile* f = new TFile("../../data/"+dir_name+"/evtcomplete_digi.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray);
    if (!fBumpArray) return -1;
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    
    TCanvas* c1=new TCanvas("PANDA1","c1",tx,ty);
    TCanvas* c2=new TCanvas("PANDA2","c2",tx,ty);
    TCanvas* c3=new TCanvas("PANDA3","c3",tx,ty);
    TCanvas* c4=new TCanvas("PANDA4","c4",tx,ty);
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
    
    TH2D* h2D1 = new TH2D("Hist1","h1",100,0,20, 100,0,3.5);
    h2D1->SetMarkerStyle(7);
    h2D1->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D1->GetXaxis()->SetTitle("distance");
    h2D1->GetYaxis()->SetTitle("N_{cluster}");
    h2D1->GetXaxis()->CenterTitle();
    h2D1->GetYaxis()->CenterTitle();
    
    TH2D* h2D2 = new TH2D("Hist2","h2",100,0,20, 100,0,3.5);
    h2D2->SetMarkerStyle(7);
    h2D2->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D2->GetXaxis()->SetTitle("distance");
    h2D2->GetYaxis()->SetTitle("N_{bump}");
    h2D2->GetXaxis()->CenterTitle();
    h2D2->GetYaxis()->CenterTitle();
    
    TH2D* h2D3 = new TH2D("Hist3","h3",100,0,20, 100,0,0.5);
    h2D3->SetMarkerStyle(7);
    h2D3->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D3->GetXaxis()->SetTitle("distance");
    h2D3->GetYaxis()->SetTitle("#delta");
    h2D3->GetXaxis()->CenterTitle();
    h2D3->GetYaxis()->CenterTitle();
    
    TH2D* h2D4 = new TH2D("Hist4","h4",100,0,20, 100,0,3);
    h2D4->SetMarkerStyle(7);
    h2D4->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D4->GetXaxis()->SetTitle("distance");
    h2D4->GetYaxis()->SetTitle("#Deltad");
    h2D4->GetXaxis()->CenterTitle();
    h2D4->GetYaxis()->CenterTitle();
    
    int N(0);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int npoints = fPointArray->GetEntriesFast();
        int ntrack = fMCTrackArray->GetEntriesFast();
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        PndMCTrack *mcTrack0 = (PndMCTrack *)fMCTrackArray->At(0);
        TVector3 mom0(mcTrack0->GetMomentum());
        PndMCTrack *mcTrack1 = (PndMCTrack *)fMCTrackArray->At(1);
        TVector3 mom1(mcTrack1->GetMomentum());
        Double_t angle = mom0.Angle(mom1);
        
        std::vector<Int_t> seed0;
        std::vector<Int_t> seed1;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
            for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) {
                if (linkIter->GetIndex() == 0) seed0.push_back( hit->GetDetectorID() );
                if (linkIter->GetIndex() == 1) seed1.push_back( hit->GetDetectorID() );
            }
        }
        
        if ((seed0.size() == 0) || (seed1.size() == 0)) continue;
        Double_t truth_E0 = 0;
        Double_t truth_E1 = 0;
        
        for (int i = 0; i < nhits; i++) {
            //********** truth energy **********//
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t>  dep = hit->GetDepositedEnergyMap();
            std::map<Int_t, Double_t>::iterator ptr;
            for ( ptr = dep.begin(); ptr != dep.end(); ptr++){
                if (ptr->first == 0) truth_E0 += ptr->second;
                if (ptr->first == 1) truth_E1 += ptr->second;
            }
        }
        
        Double_t distance = 2 * 65.0 * sin(angle/2.0);
        h2D1->Fill(distance,nclusters);
        h2D2->Fill(distance, nbumps);
        
        if (nbumps != 2) continue;
        std::vector<Int_t> bumpmatch0;
        std::vector<Int_t> bumpmatch1;
        for (int i = 0; i < nbumps; i++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(i);
            TVector3 pb = Bump->position();
            /*if (pb.Mag()*sin(mom0.Angle(pb)) < pb.Mag()*sin(mom1.Angle(pb))) bumpmatch0.push_back(i);
            else bumpmatch1.push_back(i);*///method1
            if (((pb.Mag()*sin(mom0.Angle(pb))) / (pb.Mag()*sin(mom1.Angle(pb)))) < 0.33 ) bumpmatch0.push_back(i);
            else if (((pb.Mag()*sin(mom1.Angle(pb))) / (pb.Mag()*sin(mom0.Angle(pb)))) < 0.33 ) bumpmatch1.push_back(i);
            else continue;
        }
        if ((bumpmatch0.size() != 1) || (bumpmatch1.size() != 1)) continue;
        PndEmcBump* Bump0 = (PndEmcBump*)fBumpArray->At(bumpmatch0[0]);
        PndEmcBump* Bump1 = (PndEmcBump*)fBumpArray->At(bumpmatch1[0]);
        Double_t bump_E0 = Bump0->energy();
        Double_t bump_E1 = Bump1->energy();
        
        Double_t delta =sqrt( (truth_E0 - bump_E0)*(truth_E0 - bump_E0) + (truth_E1 - bump_E1)*(truth_E1 - bump_E1) );
        h2D3->Fill(distance, sqrt(delta));
        
        TVector3 bump_p0 = Bump0->position();
        TVector3 bump_p1 = Bump1->position();
        Double_t distance_e = 2 * 65.0 * sqrt( sin(mom0.Angle(bump_p0)/2.0)*sin(mom0.Angle(bump_p0)/2.0) + sin(mom1.Angle(bump_p1)/2.0)*sin(mom1.Angle(bump_p1)/2.0) );
        h2D4->Fill(distance,distance_e);
        N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    c1->cd();
    h2D1->Draw("SCAT");
    c2->cd();
    h2D2->Draw("SCAT");
    c3->cd();
    h2D3->Draw("SCAT");
    c4->cd();
    h2D4->Draw("SCAT");
    return 0;
}
