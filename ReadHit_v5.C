int ReadHit_v5()
{
    int bin1(600),bin2(600);
    float tx(800),ty(600);
    //double xmin(0),xmax(10),ymin(0),ymax(1.1);
    double xmin(0),xmax(6),ymin(0),ymax(180);
    TString dir_name("aaa");
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
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fSharedDigiArray = (TClonesArray*) ioman->GetObject("EmcSharedDigi");
    if (!fSharedDigiArray) return -1;
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
    
    //TH1D* hist = new TH1D("h","hist",bin1,-.1,.1);
    TH1D* hist = new TH1D("h","hist",50,0,6);
    hist->GetXaxis()->SetTitle("E_{seed}");
    hist->GetYaxis()->SetTitle("Entries");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    //TH1D* hist1 = new TH1D("h1","hist1",bin1,0,3);
    
    int excnum(0);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt);
        int nbump = fBumpArray->GetEntriesFast();
        int ncluster = fClusterArray->GetEntriesFast();
        if ( ncluster != 1 ) continue;
        if ( nbump != 2 ) continue;
        PndEmcBump* bump1 = (PndEmcBump*)fBumpArray->At(0);
        //cout << bump1->energy() << endl;
        //PndEmcBump* bump2 = (PndEmcBump*)fBumpArray->At(1);
        PndEmcDigi* mdigi = bump1->Maxima(fSharedDigiArray);
        double E_0 = mdigi->GetEnergy();
        
        hist->Fill(E_0);
        //hist->Fill(bump2->energy());
        
        excnum++;
    }
    
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue);
    hist->Draw();
    //hist1->SetLineColor(kRed);
    //hist1->Draw("SAME");
    cout << "Exc Num:" << excnum << endl;
    return 0;
}
