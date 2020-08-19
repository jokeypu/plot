int PlotHit()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../../data/evtcomplete_sim.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TCanvas* c1=new TCanvas("PANDA","Bump",800,600);
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
    
    TH1D* hist=new TH1D("h","v",1000,0,0.001);
    hist->GetXaxis()->SetTitle("#theta");
    hist->GetYaxis()->SetTitle("Entries");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    int maxEvtNo = ioman->CheckMaxEventNo();
    
    for (int ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        int nhits = fHitArray->GetEntriesFast();
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            double E = hit->GetEnergy();
            hist->Fill(E);
	    if ( E < 0.00001 ) cout << E << endl;
        }
    }
    hist->Draw();
    return 0;
}


