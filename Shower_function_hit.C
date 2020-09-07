int Shower_function_hit()
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    double xmin(-0.1),xmax(5),ymin(0),ymax(1.1);
    TString dir_name("Gamma_1G_all");
    
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
    
    //TCanvas* c1=new TCanvas("PANDA1","Hit1",tx,ty);
    TCanvas* c2=new TCanvas("PANDA2","Hit2",tx,ty);
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
    gStyle->SetOptFit(1111);
    
    TH1D* h1D = new TH1D("h1D","h1",200,0,15);
    h1D->SetLineColor(kBlue);
    h1D->SetLineWidth(2);
    h1D->GetXaxis()->SetTitle("distance");
    h1D->GetYaxis()->SetTitle("Energy");
    h1D->GetXaxis()->CenterTitle();
    h1D->GetYaxis()->CenterTitle();
    
    TH3D* h2D = new TH3D("hvx0vy0","vx vs vy",bin1,xmin,xmax,600,0,2.5,bin2,ymin,ymax);
    h2D->GetXaxis()->SetTitle("d (cm)");
    h2D->GetYaxis()->SetTitle("d (cm)");
    h2D->GetZaxis()->SetTitle("E/E_{0}");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 0.5);
    
    /*
    TF1 *f=new TF1("f","[1]*[0]*x",2,10);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.25,1);
    f->SetParLimits(0, 0.01, 10);
    f->SetParLimits(1, 1, 1);
     */
    /*
    TF1 *f=new TF1("f","[0]*[1]*sqrt(x-1.2)",1.2,15);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.25,1);
    f->SetParLimits(0, 0.01, 10);
    f->SetParLimits(1, 1.0, 1.0);
    */
    
    
    Double_t fitmin(2);
    TF1 *f=new TF1("f","exp(-1*[0]*sqrt(x-[1]))",fitmin,14);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.25,1.4);
    f->SetParLimits(0, 0, 200);
    f->SetParLimits(1, 0, fitmin);
    
    //double test(0.0);
    //int cunt(0);
    
    int N(0);
    int num(5);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        int npoints = fPointArray->GetEntriesFast();
        int ntrack = fMCTrackArray->GetEntriesFast();
        int nhits = fHitArray->GetEntriesFast();
        
        if (ntrack < 2) continue;
        PndMCTrack *mcTrack_1 = (PndMCTrack *)fMCTrackArray->At(1);
        TVector3 mcStartPos(mcTrack_1->GetStartVertex());
        if (sqrt(mcStartPos.X()*mcStartPos.X()+mcStartPos.Y()*mcStartPos.Y()) < 56) continue;
        PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(0);
        TVector3 mom(mcTrack->GetMomentum());
        
        PndEmcPoint* point_0 = (PndEmcPoint*)fPointArray->At(0);
        Int_t seedID = point_0->GetDetectorID();
        Int_t seedHit = -1;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            Int_t DetID = hit->GetDetectorID();
            if (DetID == seedID) {
                seedHit = i;
                break;
            }
        }
        
        Double_t maxE = -1.0;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            if ((hit->GetEnergy()) > maxE) {maxE=(hit->GetEnergy());seedHit=i;seedID = hit->GetDetectorID();}
        }
        
        
        if ( seedHit == -1 ) continue;
        if ( npoints == 0 ) continue;
        PndEmcHit* hit_0 = (PndEmcHit*)fHitArray->At(seedHit);
        Double_t E_0 = hit_0->GetEnergy();
        TVector3 pos_0(hit_0->GetX(),hit_0->GetY(),hit_0->GetZ());
        //if ( E_0 < 0.5 ) continue;
        for (int i = 0; i < nhits; i++) {
            // computing distance from each hit to track
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            Int_t DetID = hit->GetDetectorID();
            TVector3 pos(hit->GetX(),hit->GetY(),hit->GetZ());
            Double_t distance = pos.Mag() * sin(pos.Angle(mom));
            Double_t distance_0 = pos_0.Mag() * sin(pos_0.Angle(mom));
            Double_t E = hit->GetEnergy();
            h1D->Fill(distance,E/E_0);
            //h2D->Fill(distance,-1*log(E/E_0));
            if (DetID == seedID) continue;
            h2D->Fill(distance,distance_0,E/E_0);
	    //if ( distance < 1.7 ) {test+=E/E_0;cunt++;}
        }
    N++;
    }
    //cout << test/cunt << endl;
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    //c1->cd();
    //h1D->Draw("HIST");
    c2->cd();
    //h2D->Fit(f,"R");
    h2D->Draw("HIST");
    //f->Draw("SAME");
    return 0;
}
