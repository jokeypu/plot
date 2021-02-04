Double_t DD(const TVector3 *DetPos, const TVector3 *Cent, const Double_t par){
    Double_t distance(0), angle(0);
    if (*DetPos != *Cent) {
        TVector3 vz(0, 0, 1);
        TVector3 DetPos_o = (*DetPos) - 3.7*vz;
        TVector3 DetPos_n;
        DetPos_n.SetMagThetaPhi(DetPos_o.Mag(), DetPos_o.Theta(), DetPos_o.Phi()-0.06981317);
        TVector3 ey = DetPos_n.Cross(vz).Unit();
        TVector3 ex = DetPos_n.Cross(ey).Unit();
        Double_t dx = abs((*Cent-*DetPos).Dot(ex));
        Double_t dy = abs((*Cent-*DetPos).Dot(ey));
        distance = sqrt(dx*dx+dy*dy);
        angle = 57.2957*TMath::ATan(dy/dx);
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    }
    return distance;
}

Double_t AA(const TVector3 *DetPos, const TVector3 *Cent, const Double_t par){
    Double_t distance(0), angle(0);
    if (*DetPos != *Cent) {
        TVector3 vz(0, 0, 1);
        TVector3 DetPos_o = (*DetPos) - 3.7*vz;
        TVector3 DetPos_n;
        DetPos_n.SetMagThetaPhi(DetPos_o.Mag(), DetPos_o.Theta(), DetPos_o.Phi()-0.06981317);
        TVector3 ey = DetPos_n.Cross(vz).Unit();
        TVector3 ex = DetPos_n.Cross(ey).Unit();
        Double_t dx = abs((*Cent-*DetPos).Dot(ex));
        Double_t dy = abs((*Cent-*DetPos).Dot(ey));
        distance = sqrt(dx*dx+dy*dy);
        angle = 57.2957*TMath::ATan(dy/dx);
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        //if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    }
    return angle;
}

TVector2 Pos2(const TVector3 *DetPos, const TVector3 *Cent){
    Double_t distance(0), angle(0);
    Double_t dx(0), dy(0);
    if (*DetPos != *Cent) {
        TVector3 vz(0, 0, 1);
        TVector3 DetPos_o = (*DetPos) - 3.7*vz;
        TVector3 DetPos_n;
        DetPos_n.SetMagThetaPhi(DetPos_o.Mag(), DetPos_o.Theta(), DetPos_o.Phi()-0.06981317);
        TVector3 ey = DetPos_n.Cross(vz).Unit();
        TVector3 ex = DetPos_n.Cross(ey).Unit();
        dx = (*DetPos-*Cent).Dot(ex);
        dy = (*DetPos-*Cent).Dot(ey);
        distance = sqrt(dx*dx+dy*dy);
        angle = 57.2957*TMath::ATan(dy/dx);
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        //if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    }
    TVector2 Pos2D(dx,dy);
    return Pos2D;
}
Double_t myfunc(Double_t x, Double_t p0, Double_t p1, Double_t p2, Double_t x0);
int Shower_function_fix()
{
    int bin1(600),bin2(600);
    float tx(800),ty(800);
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
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    TFile* f1 = new TFile("../data/"+dir_name+"/evtcomplete_digi.root");
    TTree* t = (TTree*)f1->Get("pndsim");
    TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray);
    if (!fBumpArray) return -1;
    
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
    
    //TH2D* h2D = new TH2D("h1D","h1",100,-3.5,3.5,100,-3.5,3.5);
    TH2D* h2D = new TH2D("h1D","h1",80,0.5,2.5,80,0,0.8);
    h2D->SetLineColor(kBlue);
    h2D->SetLineWidth(2);
    h2D->GetXaxis()->SetTitle("x");
    h2D->GetYaxis()->SetTitle("y");
    h2D->GetZaxis()->SetTitle("Energy");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    
    TH1D* h1D = new TH1D("h1D","h1",100,0,20);
    h1D->SetLineColor(kBlue);
    h1D->SetLineWidth(2);
    h1D->GetXaxis()->SetTitle("distance");
    h1D->GetYaxis()->SetTitle("Energy");
    h1D->GetXaxis()->CenterTitle();
    h1D->GetYaxis()->CenterTitle();
    
    /*TF1 *f=new TF1("f","myfunc(x,[0],[1],[2],[3])",0,9);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(210,-0.2,-0.1,4);
    f->SetParLimits(0, 100, 300);
    f->SetParLimits(1, -0.5, 0.5);
    f->SetParLimits(2, -0.5, 0.5);
    f->SetParLimits(3, 0, 8);*/
    
    TF1 *f=new TF1("f","[0]*exp(-1*[1]*x)+[2]*exp(-1*[3]*x)",0,4);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(28.5983,0.619019,225.041,3.13661);
    /*f->SetParameters(29,0.62,225,3.14);
    f->SetParLimits(0, 0, 1000);
    f->SetParLimits(1, 0.1, 0.7);
    f->SetParLimits(2, 0, 2000);
    f->SetParLimits(3, 0, 5);*/

    int N(0);
    int num(5);
    //maxEvtNo /= 10;
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        if (maxEvtNo>=100 && ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int npoints = fPointArray->GetEntriesFast();
        int ntrack = fMCTrackArray->GetEntriesFast();
        int nhits = fHitArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        
        //Exclude events generated electron-positron
        bool Exist(false);
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
            for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++)
            if (linkIter->GetIndex() == 0) Exist = true;
        }
        if (!Exist) continue;
        if (nbumps!=1) continue;
        
        if (ntrack < 2) continue;
        PndMCTrack *mcTrack_1 = (PndMCTrack *)fMCTrackArray->At(1);
        TVector3 mcStartPos(mcTrack_1->GetStartVertex());
        if (sqrt(mcStartPos.X()*mcStartPos.X()+mcStartPos.Y()*mcStartPos.Y()) < 50) continue;
        PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(0);
        TVector3 mom(mcTrack->GetMomentum());
        
        if (npoints == 0 ) continue;
        PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(0);
        TVector3 cent = Bump->where();
            
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
            TVector2 Pos2D = Pos2(&pos,&cent);
            //h2D->Fill(Pos2D.X(),Pos2D.Y(),E);
            Double_t seed1 = 1.0/TMath::Tan(mom.Theta());
            Double_t seed2 = mom.Phi();
            Double_t cood1 = 1.0/TMath::Tan(point->GetTheta());
            Double_t cood2 = point->GetPhi();
            h1D->Fill(sqrt((cood1-seed1)*(cood1-seed1)+(cood2-seed2)*(cood2-seed2)),E);
            //h2D->Fill(1/TMath::Tan(point->GetTheta()),point->GetPhi(),E);
            //cout << point->GetTheta() << ", " << point->GetPhi() << endl;
        }
    N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    c1->cd();
    h1D->Draw("HIST");
    //h2D->Draw("LEGO");
    //h2D->Draw("cont");
    //h1D->Fit(f,"R");
    //f->Draw("SAME");
    return 0;
}
