int ReadHit_v4()
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    //double xmin(0),xmax(10),ymin(0),ymax(1.1);
    double xmin(0),xmax(6),ymin(0),ymax(180);
    TString dir_name("Gamma_1G");
    //TString dir_name("Gamma_0.1to6G");
    //TString dir_name("Gamma_0.1to6G_all");
    
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
    TCanvas* c2=new TCanvas("PANDA1","Bump1",tx,ty);
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
    
    //TH1D* hist = new TH1D("h","hist",bin1,-.1,.1);
    TH1D* hist = new TH1D("h","hist",50,0,0.1);
    hist->GetXaxis()->SetTitle("E_{1}");
    hist->GetYaxis()->SetTitle("Entries");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    //TH1D* hist1 = new TH1D("h1","hist1",bin1,0,3);
    
    TF1 *f=new TF1("f","exp(-2.5*x/2.0)",xmin,xmax);
    //TF1 *f=new TF1("f","exp(-[0]*x)+[1]",xmin,xmax);
    //f->SetParameters(1.25,0.0025);
    //f->SetParameters(1.38,0.0025);
    
    int excnum(0);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt);
        int nbump = fBumpArray->GetEntriesFast();
        if ( nbump != 1 ) continue;
        int ndigi = fDigiArray->GetEntriesFast();
        PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(0);
        PndEmcDigi* mdigi = bump->Maxima(fDigiArray);
        double E_0 = mdigi->GetEnergy();
        //cout << mdigi->GetDetectorId() << endl;
        
        //**********************************
        double theta = 180*(bump->theta())/3.14;
        double phi = 180 + 180*(bump->phi())/3.14;
        //int aa = 5
        if ( E_0 < 0.8 ||  E_0 > 0.9 ) continue;
        //if ( theta < 0 || theta > 180 ) continue;
        //if ( phi < 0 || phi > 360 ) continue;
        //**********************************
        
        for ( int idigi = 0; idigi < ndigi ; idigi++) {
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(idigi);
            double E = digi->GetEnergy();
            if ( E/E_0 > 0.99 ) continue;
            //double d = bump->DistanceToCentre(digi);
            double d = (digi->where()-mdigi->where()).Mag();
            //double dE = E/E_0 - exp(-2.5*d/2.0);
            double dE = -1*log(E/E_0)/d;
            Int_t DetID = digi->GetDetectorId();
            //if ( E < 0.15 ) continue;
            //cout << DetID << endl;
            if ( DetID == 106130005 ) {
                hist->Fill(E);
                excnum++;
                cout << d << endl;
                continue;
            }
            //hist1->Fill(-1*log(E/E_0-0.0025)/d);
            //histxy->Fill(d,E/E_0);
        }
        
    }
    c1->cd();
    histxy->SetMarkerStyle(7);
    histxy->SetMarkerColorAlpha(kAzure+3, 0.5);
    histxy->Draw("SCAT");
    f->Draw("SAME");
    c2->cd();
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue);
    hist->Draw();
    //hist1->SetLineColor(kRed);
    //hist1->Draw("SAME");
    cout << "Exc Num:" << excnum << endl;
    return 0;
}
