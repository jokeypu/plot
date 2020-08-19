int ReadHit_v3()
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    double xmin(0),xmax(360),ymin(0),ymax(5);
    //TString dir_name("Gamma_6G");
    //TString dir_name("Gamma_0.1to6G");
    TString dir_name("Gamma_0.1to6G_all");
    
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
    histxy->GetXaxis()->SetTitle("#phi");
    histxy->GetYaxis()->SetTitle("#xi");
    histxy->GetXaxis()->CenterTitle();
    histxy->GetYaxis()->CenterTitle();
    
    int excnum(0);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt);
        int nbump = fBumpArray->GetEntriesFast();
        if ( nbump != 1 ) continue;
        int ndigi = fDigiArray->GetEntriesFast();
        PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(0);
        PndEmcDigi* mdigi = bump->Maxima(fDigiArray);
        double E_0 = mdigi->GetEnergy();
        
        //**********************************
        double theta = 180*(bump->theta())/3.14;
        double phi = 180 + 180*(bump->phi())/3.14;
        //int aa = 5
        //if ( E_0 < aa ||  E_0 > aa+1 ) continue;
        //if ( theta < 0 || theta > 180 ) continue;
        //if ( phi < 0 || phi > 360 ) continue;
        //**********************************
        
        for ( int idigi = 0; idigi < ndigi ; idigi++) {
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(idigi);
            double E = digi->GetEnergy();
            if ( E/E_0 > 0.99 ) continue;
            double d = bump->DistanceToCentre(digi);
            double dE = -2*log(E/E_0)/d;
            //histxy->Fill(E_0,dE);
            histxy->Fill(phi,dE);
        }
        excnum++;
    }
    c1->cd();
    histxy->SetMarkerStyle(7);
    histxy->SetMarkerColorAlpha(kAzure+3, 0.5);
    histxy->Draw("SCAT");
    cout << "Exc Num:" << excnum << endl;
    return 0;
}
