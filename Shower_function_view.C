int Shower_function_view()
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    double xmin(0),xmax(6),ymin(0),ymax(180);
    TString dir_name("Gamma_1G_90d");
    
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
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    TCanvas* c1=new TCanvas("PANDA","Point",tx,ty);
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
    
    //TH3D* h3D = new TH3D("h3D","3D",200,5,35,200,-15,15,200,55,80);
    TH3D* h3D = new TH3D("h3D","3D",200,-100,100,200,-15,15,200,55,80);
    h3D->SetMarkerStyle(1);
    //h3D->SetMarkerStyle(6);
    h3D->SetMarkerColorAlpha(kAzure+3, 0.7);
    //h3D->SetMarkerColorAlpha(kRed, 0.7);
    h3D->GetXaxis()->SetTitle("x");
    h3D->GetYaxis()->SetTitle("y");
    h3D->GetZaxis()->SetTitle("z");
    h3D->GetXaxis()->CenterTitle();
    h3D->GetYaxis()->CenterTitle();
    
    int N(0);
    int num(2);
    for (Int_t ievt = num; ievt < num+100; ievt++) {
    //for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        int npoints = fPointArray->GetEntriesFast();
        int ntrack = fMCTrackArray->GetEntriesFast();
        
        if (ntrack < 2) continue;
        PndMCTrack *mcTrack_1 = (PndMCTrack *)fMCTrackArray->At(1);
        TVector3 mcStartPos(mcTrack_1->GetStartVertex());
        if (sqrt(mcStartPos.X()*mcStartPos.X()+mcStartPos.Y()*mcStartPos.Y()) < 56) continue;
        PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(0);
        TVector3 mom(mcTrack->GetMomentum());
        
        if (npoints == 0 ) continue;
        PndEmcPoint* point_0 = (PndEmcPoint*)fPointArray->At(0);
        Double_t E_0 = point_0->GetEnergyLoss();
        for (int i = 0; i < npoints; i++) {
            // computing distance from each point to track
            PndEmcPoint* point = (PndEmcPoint*)fPointArray->At(i);
            Double_t x = point->GetX();
            Double_t y = point->GetY();
            Double_t z = point->GetZ();
            TVector3 pos(x, y, z);
            Double_t distance = pos.Mag()*sin(mom.Angle(pos));
            Double_t E = point->GetEnergyLoss();
            h3D->Fill(z, y, x);
            N++;
        }
    }
    c1->cd();
    h3D->Draw("SCAT");
    return 0;
}
