int ReadHit_v8()
{
    int bin1(300),bin2(300);
    float tx(800),ty(600);
    double xmin(0),xmax(20),ymin(0),ymax(20);
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
    gStyle->SetOptFit(1111);
    
    TH2D* histxy = new TH2D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax);
    histxy->GetXaxis()->SetTitle("d (cm)");
    histxy->GetYaxis()->SetTitle("-ln(E/E_{0})");
    histxy->GetXaxis()->CenterTitle();
    histxy->GetYaxis()->CenterTitle();
    
    //TH1D* hist = new TH1D("h","hist",bin1,-.1,.1);
    TH1D* hist = new TH1D("h","hist",100,0,5);
    hist->GetXaxis()->SetTitle("#xi");
    hist->GetYaxis()->SetTitle("Entries");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    //TH1D* hist1 = new TH1D("h1","hist1",bin1,0,3);
    
    TGraph* g = new TGraph();
    //TF1 *f=new TF1("f","2.5*x/2.0",xmin,xmax);
    //TF1 *f=new TF1("f","2.5*pow(x-1.4,0.5)",1.4,xmax);
    /*TF1 *f=new TF1("f","pow([0]*log(x-[1]),[2])",1.34,9);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(6.0,0.2,0.8);
    f->SetParLimits(0, 5.8, 6.8);
    f->SetParLimits(1, 0.0, 1.0);
    f->SetParLimits(2, 0.5, 1.34);*/
    //f->SetParLimits(3, 0.1, 5.0);
    //f->SetParameters(1.25,0.0025);
    //f->SetParameters(1.38,0.0025);
    
    
    /*TF1 *f=new TF1("f","[0]*pow((x-[1]),[2])",1.3,9);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1,1.2,1);
    f->SetParLimits(0, 0, 3);
    f->SetParLimits(1, 0, 1.3);
    f->SetParLimits(2, 0, 1.0);
    */

    TF1 *f=new TF1("f","[0]-[1]*exp(-[2]*x)",1.2,16);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(7,20,1);
    f->SetParLimits(0, 6, 8);
    f->SetParLimits(1, 10, 100);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    f->SetParLimits(2, 0.1, 5);
    
    int excnum(0),N(0);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt);
        int nbump = fBumpArray->GetEntriesFast();
        if ( nbump != 1 ) continue;
        int ndigi = fDigiArray->GetEntriesFast();
        PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(0);
        PndEmcDigi* mdigi = bump->Maxima(fDigiArray);
        double E_0 = mdigi->GetEnergy();
        
        for ( int idigi = 0; idigi < ndigi ; idigi++) {
            PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(idigi);
            double E = digi->GetEnergy();
            if ( E/E_0 > 0.99 ) continue;
            double d = bump->DistanceToCentre(digi);
            //double dE = E/E_0 - exp(-2.5*d/2.0);
            double dE = -2*log(E/E_0)/d;
            hist->Fill(dE);
            //hist1->Fill(-1*log(E/E_0-0.0025)/d);
            //histxy->Fill(d,-1*log(E/E_0));
            histxy->Fill(d,-(1/0.48454)*log((6+log(E/E_0))/10));
            //histxy->Fill(-0.8*log(E/E_0),0.1*pow(5.15*log(d-0.2),0.8)*d);
            if ( d > 16 ) continue;
            g->SetPoint(N,d,-1*log(E/E_0));
            //g->SetPoint(N,2.865*pow(d-1.3,0.44),-1*log(E/E_0));
            //g->SetPoint(N,pow(6.4*log(d-0.33),0.7),-1*log(E/E_0));
            N++;
        }
        excnum++;
    }
    c2->cd();
    histxy->SetMarkerStyle(7);
    histxy->SetMarkerColorAlpha(kAzure+3, 0.5);
    //histxy->Draw("SCAT");
    //histxy->Fit(f);
    //f->Draw("SAME");
    
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
    g->GetXaxis()->SetTitle("d (cm)");
    g->GetYaxis()->SetTitle("-ln(E/E_{0})");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->Fit(f,"R");
    g->Draw("AP");
    
    c1->cd();
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue);
    hist->Draw();
    //hist1->SetLineColor(kRed);
    //hist1->Draw("SAME");
    cout << "Exc Num:" << excnum << endl;
    return 0;
}
