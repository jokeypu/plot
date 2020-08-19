int Energy_Compare_t()
{
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile("../data/Compare_n/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("../data/Compare_n/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    TClonesArray* fHitArray = new TClonesArray("PndEmcHit");
    t->SetBranchAddress("EmcHit",&fHitArray);
    if (!fHitArray) return -1;
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fBumpArray = (TClonesArray*) ioman->GetObject("EmcBump");
    if (!fBumpArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    TClonesArray* fSharedDigiArray = (TClonesArray*) ioman->GetObject("EmcSharedDigi");
    if (!fSharedDigiArray) return -1;
    Int_t maxEvtNo = ioman->CheckMaxEventNo();
    
    int bin1(1000),bin2(1000);
    float tx(800),ty(580);
    double xmin(0),xmax(4.3),ymin(-1),ymax(1); //double xmin(0),xmax(3),ymin(0),ymax(2);
    
    TCanvas* c1=new TCanvas("PANDA1","MCTruth1",tx,ty);
    TCanvas* c2=new TCanvas("PANDA2","MCTruth2",tx,ty);
    //gStyle->SetOptTitle(0);
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
    
    TH2D* hxy1 = new TH2D("hist2D1", "hxy1" , bin1, xmin, xmax, bin2, ymin, ymax);
    TH2D* hxy2 = new TH2D("hist2D2", "hxy2" , bin1, xmin, xmax, bin2, ymin, ymax);
    TH1D* h1 = new TH1D("hist","",300,-0.5,0.5);
    TH1D* h2 = new TH1D("hist2","",300,-0.5,0.5);
    
    std::map<Int_t, Double_t>::iterator it;
    
    int aa=2;
    int Nmiss = 0;
    int Nh1(0),Nh2(0),N11(0);
    //for (Int_t ievt = aa; ievt < aa+1; ievt++) {
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
    //    for (Int_t ievt = 0; ievt < 350; ievt++) {
        t->GetEntry(ievt);
        ioman->ReadEvent(ievt);
        std::map<Int_t, Double_t> TEnergy_M;
        TEnergy_M.clear();
        std::map<Int_t, Double_t> TEnergy_B;
        TEnergy_B.clear();
        int noverlap_M = 0;
        int nhit = fHitArray->GetEntriesFast();
        for (int n = 0; n < nhit; n++){
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(n);
            std::map<Int_t, Double_t> ds = hit->GetDepositedEnergyMap();
            for ( it = ds.begin(); it != ds.end(); it++ ) {
                TEnergy_M[it->first] += it->second;
            }
        }
        if ( TEnergy_M.size() != 2 ) continue;
        int nbump = fBumpArray->GetEntriesFast();
        if ( nbump == 2 ) {
            PndEmcBump* bump1 = (PndEmcBump*)fBumpArray->At(0);
            std::vector<Int_t> list1 = bump1->DigiList();
            PndEmcBump* bump2 = (PndEmcBump*)fBumpArray->At(1);
            std::vector<Int_t> list2 = bump2->DigiList();
            for (int i=0; i < list1.size(); i++){
                PndEmcDigi* digi1 = (PndEmcDigi*)fSharedDigiArray->At(list1[i]);
                TEnergy_B[1] += digi1->GetEnergy();
            }
            for (int i=0; i < list2.size(); i++){
                PndEmcDigi* digi2 = (PndEmcDigi*)fSharedDigiArray->At(list2[i]);
                TEnergy_B[2] += digi2->GetEnergy();
            }
            double B_Energy1 = TEnergy_B[1];
            double B_Energy2 = TEnergy_B[2];
            it = TEnergy_M.begin();
            double M_Energy1 = it->second;
            it++;
            double M_Energy2 = it->second;
            if((abs(B_Energy1-M_Energy1)+abs(B_Energy2-M_Energy2)) < (abs(B_Energy2-M_Energy1)+abs(B_Energy1-M_Energy2))){
                hxy1->Fill(M_Energy1,B_Energy1-M_Energy1);  //hxy1->Fill(M_Energy1,B_Energy1/M_Energy1);
                hxy1->Fill(M_Energy2,B_Energy2-M_Energy2);  //hxy1->Fill(M_Energy2,B_Energy2/M_Energy2);
                h1->Fill(B_Energy1-M_Energy1);
                h1->Fill(B_Energy2-M_Energy2);
                Nh1+=2;
            }else{
                hxy1->Fill(M_Energy1,B_Energy2-M_Energy1);  //hxy1->Fill(M_Energy1,B_Energy2/M_Energy1);
                hxy1->Fill(M_Energy2,B_Energy1-M_Energy2);  //hxy1->Fill(M_Energy2,B_Energy1/M_Energy2);
                h1->Fill(B_Energy2-M_Energy1);
                h1->Fill(B_Energy1-M_Energy2);
                Nh1+=2;
            }
        
        }else if ( nbump == 1 ) {
            PndEmcBump* bump = (PndEmcBump*)fBumpArray->At(0);
            std::vector<Int_t> list = bump->DigiList();
            for (int i=0; i < list.size(); i++){
                PndEmcDigi* digi = (PndEmcDigi*)fSharedDigiArray->At(list[i]);
                TEnergy_B[1] += digi->GetEnergy();
            }
            double B_Energy = TEnergy_B[1];
            it = TEnergy_M.begin();
            double M_Energy = it->second;
            it++;
            M_Energy += it->second;
            hxy2->Fill(M_Energy,B_Energy-M_Energy); //hxy2->Fill(M_Energy,B_Energy/M_Energy);
            h2->Fill(B_Energy-M_Energy);
            Nh2++;
            
        }else continue;
    }
    
    TLegend * leg1 = new TLegend(0.61,0.72,0.88,0.85);
    TLegend * leg2 = new TLegend(0.7,0.75,0.88,0.85);
    
    hxy1->SetMarkerStyle(8);
    hxy1->SetMarkerColorAlpha(kRed-7, 0.4);
    //hxy1->SetMarkerColorAlpha(kAzure+3, 0.4);
    hxy1->SetMarkerSize(0.8);
    hxy1->GetXaxis()->SetTitle("E_{truth} (GeV)");
    hxy1->GetYaxis()->SetTitle("E_{bump} - E_{truth} (GeV)");
    hxy1->GetXaxis()->CenterTitle();
    hxy1->GetYaxis()->CenterTitle();
    hxy1->GetXaxis()->SetRangeUser(xmin,xmax);
    hxy1->GetYaxis()->SetRangeUser(ymin,ymax);
    hxy1->GetXaxis()->SetLabelSize(0.038);
    hxy1->GetYaxis()->SetLabelSize(0.038);
    hxy1->GetXaxis()->SetTitleSize(0.04);
    hxy1->GetYaxis()->SetTitleSize(0.04);
    hxy1->GetXaxis()->SetTitleOffset(1.2);
    hxy1->GetYaxis()->SetTitleOffset(1.2);
    hxy1->SetTitle("Bump and Truth Energy");
    
    hxy2->SetMarkerStyle(22);   //45 22
    hxy2->SetMarkerColorAlpha(kAzure+3, 0.5);
    //hxy2->SetMarkerColorAlpha(kRed-7, 0.5);
    hxy2->GetXaxis()->SetTitle("E_{truth} (GeV)");
    hxy2->GetYaxis()->SetTitle("E_{bump} - E_{truth} (GeV)");
    hxy2->GetXaxis()->CenterTitle();
    hxy2->GetYaxis()->CenterTitle();
    hxy2->GetXaxis()->SetRangeUser(xmin,xmax);
    hxy2->GetYaxis()->SetRangeUser(ymin,ymax);
    hxy2->GetXaxis()->SetLabelSize(0.038);
    hxy2->GetYaxis()->SetLabelSize(0.038);
    hxy2->GetXaxis()->SetTitleSize(0.04);
    hxy2->GetYaxis()->SetTitleSize(0.04);
    hxy2->GetXaxis()->SetTitleOffset(1.2);
    hxy2->GetYaxis()->SetTitleOffset(1.2);
    hxy2->SetTitle("Bump and Truth Energy");
    
    c1->cd();
    //c1->SetLogx();
    hxy1->Draw("SCAT");
    hxy2->Draw("SCATsame");
    
    leg1->AddEntry(hxy1, "Overlapping Entries", "P");
    leg1->AddEntry(hxy2, "Non-overlapping Entries", "P");
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.03);
    leg1->Draw("SAME");
    
    h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("E_{bump} - E_{truth} (GeV)");
    h1->GetYaxis()->SetTitle("Entries");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->SetLineColor(kRed);
    h1->GetXaxis()->SetLabelSize(0.038);
    h1->GetYaxis()->SetLabelSize(0.038);
    h1->GetXaxis()->SetTitleSize(0.04);
    h1->GetYaxis()->SetTitleSize(0.04);
    h1->GetXaxis()->SetTitleOffset(1.2);
    h1->GetYaxis()->SetTitleOffset(1.2);
    h1->SetTitle("Bump and Truth Energy");
    
    h2->SetLineWidth(2);
    h2->GetXaxis()->SetTitle("E_{bump} - E_{truth} (GeV)");
    h2->GetYaxis()->SetTitle("Entries");
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();
    h2->SetLineColor(kBlue);
    h2->GetXaxis()->SetLabelSize(0.038);
    h2->GetYaxis()->SetLabelSize(0.038);
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetYaxis()->SetTitleSize(0.04);
    h2->GetXaxis()->SetTitleOffset(1.2);
    h2->GetYaxis()->SetTitleOffset(1.2);
    h2->SetTitle("Bump and Truth Energy");
    
    c2->cd();
    h1->Draw();
    h2->Draw("SAME");
    leg2->AddEntry(h1, "Overlapping Entries", "L");
    leg2->AddEntry(h2, "Non-overlapping Entries", "L");
    leg2->SetTextFont(42);
    leg2->Draw("SAME");
    
    cout << "miss: " << Nmiss << endl;
    cout << "Nhist1: " << Nh1 << endl;
    cout << "Nhist2: " << Nh2 << endl;
    cout << "N11: " << N11 << endl;
    return 0;
}
/*TEST/
 std::map<Int_t, Double_t>::iterator pg;
 for ( pg = M_others_energy.begin(); pg != M_others_energy.end(); ++pg){
 cout << "DetID: " << pg->first << " Energy: " << pg->second << endl;
 }
 std::map<Int_t, std::vector<Double_t>>::iterator ps;
 for ( ps = M_overlap_energy.begin(); ps != M_overlap_energy.end(); ++ps){
 cout << "DetID: " << ps->first << " Energy: " << (ps->second)[0] << ", " << (ps->second)[1] << "    N:" << (ps->second).size() << endl;
 }
 cout << "************************************" << endl;
 for ( ps = B_overlap_energy.begin(); ps != B_overlap_energy.end(); ++ps){
 cout << "DetID: " << ps->first << " Energy: " << (ps->second)[0] << ", " << (ps->second)[1] << "    N:" << (ps->second).size() << endl;
 }
 /TEST*/

