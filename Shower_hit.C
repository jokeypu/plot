Double_t func_h(Double_t d,Double_t p0,Double_t p1,Double_t p2,Double_t p3,Double_t p4,Double_t p5,Double_t p6,Double_t p7,Double_t p8,Double_t p9 ){
    if (d < 1.4 ) return p0*pow(d,p1)+p2;
    else if ( d < 3.5) return p3*TMath::Landau((d-p4),p5,p6);
    else return p7*exp(p8*d+p9);
}
Double_t newfunc1(Double_t distance){
    if (distance < 1.4 ) return -0.19945*pow(distance,2.90303)+0.855836;
    else if ( distance < 3.5) return 2.14778*TMath::Landau(distance-0.625135,0.690762,0.121077);
    else return 0.0590966*exp(-0.428928*distance - 0.473589);
}
int Shower_hit(){
    int bin1(50),bin2(50),bin3(50);
    float tx(800),ty(600);
    double xmin(0),xmax(3),ymin(0),ymax(46),zmin(0),zmax(1.01);
    TString dir_name("Gamma_one_1G");
    TVector3 vz(0, 0, 1);
    
    //******************************************//
    ofstream out;
    //out.open("doc/Shower_hit_90.txt",ios::out);
    out.open("doc/Shower_hit_1.txt",ios::out);


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
    
    TFile* ff = new TFile("../data/"+dir_name+"/evtcomplete_digi.root");
    TTree* t = (TTree*)ff->Get("pndsim");
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray);
    if (!fBumpArray) return -1;
    
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
    
    //TH3D* h2D = new TH3D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax,bin3,zmin,zmax);
    TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",200,0,5,200,0,1);
    h2D->GetYaxis()->SetTitle("angle");
    h2D->GetXaxis()->SetTitle("d (cm)");
    h2D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 1);
    
    
    TF1 *f=new TF1("f","func_h(x,[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])",0,5);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(-0.0849149, 2.58396, 0.419646, 10, 0.0124098, 0.0780432, 0.182612, 0.144939, -0.435278, 0.0642399);
    /*f->SetParLimits(0, -2, 0);
     f->SetParLimits(1, 0, 10);
     f->SetParLimits(2, 0, 1);
     f->SetParLimits(3, 0, 50);
     f->SetParLimits(4, 0, 1);
     f->SetParLimits(5, 0, 1);
     f->SetParLimits(6, 0, 1);
     f->SetParLimits(7, 0, 1);
     f->SetParLimits(8, -3, 0);
     f->SetParLimits(9, 0, 1);*/
    
    TF1 *f1=new TF1("f","newfunc1(x)",0,5);
    f1->SetLineWidth(2);
    f1->SetLineColor(kRed);
    /*
     TF1 *f=new TF1("f","[0]*[1]*sqrt(x-1.2)",1.2,15);
     f->SetLineWidth(2);
     f->SetLineColor(kRed);
     f->SetParameters(1.25,1);
     f->SetParLimits(0, 0.01, 10);
     f->SetParLimits(1, 1.0, 1.0);
     */
    
    TF1 *f2=new TF1("f2","[0]*exp(-1*[1]*x)+[2]*exp(-1*[3]*x)",0,15);
    f2->SetLineWidth(2);
    f2->SetLineColor(kRed);
    
    int N(0);
    int num(5);
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int npoints = fPointArray->GetEntriesFast();
        int ntrack = fMCTrackArray->GetEntriesFast();
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        
        if (ntrack < 2) continue;
        PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(0);
        TVector3 mom(mcTrack->GetMomentum());
        
        //Exclude events generated electron-positron
        bool Exist = false;
        int seedid(-1);
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
            for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++)
                if (linkIter->GetIndex() == 0) { Exist = true; seedid = i; }
        }
        if (!Exist) continue;
        
        if (nclusters != 1) continue;
        
        PndEmcHit* seedhit = (PndEmcHit*)fHitArray->At(seedid);
        Double_t E0 = seedhit->GetEnergy();
        if (E0 < 0.7) continue;
        
        PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(0);
        TVector3 Cent = bump->where();
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            Double_t E = hit->GetEnergy();
            TVector3 DetPos(hit->GetX(), hit->GetY(), (hit->GetZ()));
            TVector3 DetPos_o = (DetPos) - 3.7*vz;
            //TVector3 DetPos_n = DetPos;
            TVector3 DetPos_n;
            DetPos_n.SetMagThetaPhi(DetPos_o.Mag(), DetPos_o.Theta(), DetPos_o.Phi()-0.06981317);
            TVector3 ey = DetPos_n.Cross(vz).Unit();
            TVector3 ex = DetPos_n.Cross(ey).Unit();
            Double_t dx = abs((Cent-DetPos).Dot(ex));
            Double_t dy = abs((Cent-DetPos).Dot(ey));
            Double_t angle = 57.29578*TMath::ATan(dy/dx);
            angle = abs(fmod(angle,45.0) - 45*(((int)(angle/45.0))%2));
            //angle = abs(fmod(angle,90.0) - 90*(((int)(angle/90.0))%2));
            Double_t distance = sqrt(dx*dx+dy*dy);
            //if (angle>1 && angle <2)
            //h2D->Fill(distance,angle,E);
            h2D->Fill(distance,E/E0);
            out << distance << endl;
	        out << angle << endl;
            out << E/E0 << endl;
        }
        N++;
    }
    //cout << test/cunt << endl;
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    //c1->cd();
    //h1D->Draw("HIST");
    c2->cd();
    h2D->Fit(f2,"R");
    h2D->Draw("HIST");
    //h2D->Draw("LEGO");
    f2->Draw("SAME");
    out.close();
    return 0;
}
