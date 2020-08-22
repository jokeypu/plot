int Shower_match( TString dir_name="Gamma_tow_1G_n" )
{
    int bin1(100),bin2(200);
    float tx(1200),ty(900);
    double xmin(0),xmax(20),ymin(0),ymax(0.6);
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
    
    TH2D* h2D1 = new TH2D("Hist1","h1",bin1,xmin,xmax, 10,0.5,2.5);
    h2D1->SetMarkerStyle(7);
    h2D1->SetMarkerColorAlpha(kRed+3, 0.5);
    h2D1->GetXaxis()->SetTitle("distance");
    h2D1->GetYaxis()->SetTitle("N_{cluster}");
    h2D1->GetXaxis()->CenterTitle();
    h2D1->GetYaxis()->CenterTitle();
    h2D1->GetXaxis()->SetTitleSize(0);
    h2D1->GetYaxis()->SetTitleSize(0.2);
    h2D1->GetXaxis()->SetTitleOffset(0.5);
    h2D1->GetYaxis()->SetTitleOffset(0.15);
    h2D1->GetXaxis()->SetLabelSize(0);
    h2D1->GetYaxis()->SetLabelSize(0.2);
    h2D1->GetYaxis()->SetNdivisions(503);

    
    TH2D* h2D2 = new TH2D("Hist2","h2",bin1,xmin,xmax, 10,0.5,2.5);
    h2D2->SetMarkerStyle(7);
    h2D2->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D2->GetXaxis()->SetTitle("distance");
    h2D2->GetYaxis()->SetTitle("N_{bump}");
    h2D2->GetXaxis()->CenterTitle();
    h2D2->GetYaxis()->CenterTitle();
    h2D2->GetXaxis()->SetTitleSize(0.15);
    h2D2->GetYaxis()->SetTitleSize(0.15);
    h2D2->GetXaxis()->SetTitleOffset(0.8);
    h2D2->GetYaxis()->SetTitleOffset(0.2);
    h2D2->GetXaxis()->SetLabelSize(0.15);
    h2D2->GetYaxis()->SetLabelSize(0.15);
    h2D2->GetYaxis()->SetNdivisions(503);
    
    TH2D* h2D3 = new TH2D("Hist3","h3",bin1,xmin,xmax, bin2,ymin,ymax);
    h2D3->SetMarkerStyle(7);
    h2D3->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D3->GetXaxis()->SetTitle("distance");
    h2D3->GetYaxis()->SetTitle("#delta");
    h2D3->GetXaxis()->CenterTitle();
    h2D3->GetYaxis()->CenterTitle();
    h2D3->GetXaxis()->SetTitleSize(0);
    h2D3->GetYaxis()->SetTitleSize(0.1);
    h2D3->GetYaxis()->SetTitleOffset(0.45);
    h2D3->GetXaxis()->SetLabelSize(0);
    h2D3->GetYaxis()->SetLabelSize(0.07);

    int N(0);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
    //for (Int_t ievt = 0; ievt < 1000; ievt++) {
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
            if (pb.Mag()*sin(mom0.Angle(pb)) < pb.Mag()*sin(mom1.Angle(pb))) bumpmatch0.push_back(i);
            else bumpmatch1.push_back(i);
        }
        if ((bumpmatch0.size() != 1) || (bumpmatch1.size() != 1)) continue;
        PndEmcBump* Bump0 = (PndEmcBump*)fBumpArray->At(bumpmatch0[0]);
        PndEmcBump* Bump1 = (PndEmcBump*)fBumpArray->At(bumpmatch1[0]);
        Double_t bump_E0 = Bump0->energy();
        Double_t bump_E1 = Bump1->energy();
        
        Double_t delta_2 = (truth_E0 - bump_E0)*(truth_E0 - bump_E0) + (truth_E1 - bump_E1)*(truth_E1 - bump_E1);
        h2D3->Fill(distance, sqrt(delta_2));
        N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    c1->Divide(1, 4);
    c1->GetPad(1)->SetPad(0,1,1,0.45);
    c1->GetPad(2)->SetPad(0,0.475,1,0.275);
    c1->GetPad(3)->SetPad(0,0.3,1,0.025);
    c1->GetPad(4)->SetPad(0,0.025,1,0);
    c1->GetPad(1)->SetGridx();
    c1->GetPad(2)->SetGridx();
    c1->GetPad(3)->SetGridx();
    c1->cd(1);
    h2D3->Draw("SCAT");
    c1->cd(2);
    h2D1->Draw("CONT");
    c1->cd(3);
    h2D2->Draw("CONT");
    return 0;
}
