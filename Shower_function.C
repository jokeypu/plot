int Shower_function()
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    double xmin(0),xmax(6),ymin(0),ymax(180);
    TString dir_name("Gamma_1G");
    
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
    gStyle->SetOptStat(1);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    gStyle->SetOptFit(1111);
    
    TH1D* h1D = new TH1D("h1D","h1",50,0,5);
    h1D->SetLineColor(kBlue);
    h1D->SetLineWidth(2);
    h1D->GetXaxis()->SetTitle("distance");
    h1D->GetYaxis()->SetTitle("Energy");
    h1D->GetXaxis()->CenterTitle();
    h1D->GetYaxis()->CenterTitle();
    
    /*
    TF1 *f=new TF1("f","[0]*exp(-1*[1]*pow(x-[3],[2]))",0,5);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(300,1.25,0.8,0);
    f->SetParLimits(0, 0, 500);
    f->SetParLimits(1, 0.5, 8);
    f->SetParLimits(2, 0.1, 1.5);
    f->SetParLimits(3, -1.5, 1.5);
    */
    
    TF1 *f=new TF1("f","[0]/(pow(x-[1],2)-[2])",0,5);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1,0,0);
    f->SetParLimits(0, 0, 100);
    f->SetParLimits(1, -3, 3);
    f->SetParLimits(2, -3, 3);
    
    int N(0);
    int num(5);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
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
            h1D->Fill(distance,E);
        }
    N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    c1->cd();
    h1D->Fit(f,"R");
    h1D->Draw("HIST");
    f->Draw("SAME");
    return 0;
}
