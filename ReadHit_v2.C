int ReadHit_v2()
{
    int bin1(300),bin2(300);
    float tx(800),ty(600);
    double xmin(0),xmax(6),ymin(0),ymax(1.05);
    TString dir_name("Gamma_1G");
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/"+dir_name+"/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TClonesArray* fBumpArray = (TClonesArray*) ioman->GetObject("EmcBump");
    if (!fBumpArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    TCanvas* c1=new TCanvas("PANDA","Bump",tx,ty);
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
    
    TH2D* histxy = new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy->GetXaxis()->SetTitle("distance");
    histxy->GetYaxis()->SetTitle("E/E_{0}");
    histxy->GetXaxis()->CenterTitle();
    histxy->GetYaxis()->CenterTitle();
    
    TF1 *f=new TF1("f","exp(-2.5*x/2.0)",xmin,xmax);
    
    
    int excnum(0);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt);
        int nbump = fBumpArray->GetEntriesFast();
        if ( nbump > 1 ) continue;
        int ndigi = fDigiArray->GetEntriesFast();
        PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(0);
        PndEmcDigi* mdigi = bump->Maxima(fDigiArray);
        double E_0 = mdigi->GetEnergy();
        
        for ( int idigi = 0; idigi < ndigi ; idigi++) {
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(idigi);
            double E = digi->GetEnergy();
            double d = bump->DistanceToCentre(digi);
            //double d = ((digi->where()) - (mdigi->where())).Mag();
            histxy->Fill(d,E/E_0);
        }
        
        excnum++;
    }
    histxy->SetMarkerStyle(7);
    histxy->SetMarkerColorAlpha(kAzure+3, 0.5);
    histxy->Draw("SCAT");
    f->Draw("SAME");
    cout << "Exc Num:" << excnum << endl;
    return 0;
}
