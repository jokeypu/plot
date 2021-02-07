int Exec(string dir_name, string out_name, Int_t NGamma=2, bool IsSplit=1);
Double_t C_min = -1, C_max = -1;
double Z_MyTool_GetDeltaTheta_mpi0(double Theta_cent = 181, double Phi_Range = 9){
    double Phi_cent = 0;
    Phi_Range = fabs(Phi_Range);
    if (Theta_cent < 0 || fabs(Theta_cent)>180.001){
        cout << "Input theta center:(DEG)" << endl;
        cin >> Theta_cent;
        cout << "Input Phi range:(DEG)   default: 9 DEG" << endl;
        cin >> Phi_Range;
        if (Theta_cent <= 0) {cout << "-E  ERROR!!" << endl;return "XXXXXXXXX";}
        if (Phi_Range <= 0) Phi_Range = 9;
    }
    cout << "-INFO  Theta Center:  " << Theta_cent << " DEG" << endl;
    cout << "-INFO  Phi Center:  " << Phi_cent << " DEG" << endl;
    cout << "-INFO  Phi Range:  " << Phi_Range << " DEG" << endl;
    double Theta_cent_save = Theta_cent;
    double Phi_Range_save = Phi_Range;
    if (Theta_cent == 90) Theta_cent = 89.99999;
    Theta_cent *= TMath::DegToRad();
    Phi_Range *= TMath::DegToRad();
    TF1 *f = new TF1("func","tan([0]+x)-tan([0]-x)",0,Phi_Range);
    f->SetParameter(0, fabs(90*TMath::DegToRad() - Theta_cent));
    double DeltaTheta = f->GetX(Phi_Range);
    DeltaTheta *= TMath::RadToDeg();
    return DeltaTheta;
}
void Set_Angle_Cut(double min, double max){C_min = min; C_max = max;}
int Shower_iFile_mpi0(string dir_name)
{
    int NO_Angle = 15;
    
    if (NO_Angle == 1) {t_min = 23.8514; t_max = 24.6978;}
    if (NO_Angle == 2) {t_min = 26.4557; t_max = 27.3781;}
    if (NO_Angle == 3) {t_min = 29.4579; t_max = 30.4916;}
    if (NO_Angle == 4) {t_min = 32.6536; t_max = 33.7759;}
    if (NO_Angle == 5) {t_min = 36.1172; t_max = 37.3507;}
    if (NO_Angle == 6) {t_min = 39.9051; t_max = 41.2390;}
    if (NO_Angle == 7) {t_min = 44.2385; t_max = 45.7355;}
    if (NO_Angle == 8) {t_min = 48.8451; t_max = 50.4459;}
    if (NO_Angle == 9) {t_min = 53.7548; t_max = 55.4790;}
    if (NO_Angle == 10) {t_min = 59.0059; t_max = 60.8229;}
    if (NO_Angle == 11) {t_min = 64.7855; t_max = 66.7591;}
    if (NO_Angle == 12) {t_min = 70.8088; t_max = 72.8652;}
    if (NO_Angle == 13) {t_min = 77.0506; t_max = 79.1942;}
    if (NO_Angle == 14) {t_min = 83.4997; t_max = 85.6749;}
    if (NO_Angle == 15) {t_min = 90.2068; t_max = 92.4062;}
    if (NO_Angle == 16) {t_min = 96.8200; t_max = 99.0099;}
    if (NO_Angle == 17) {t_min = 103.361; t_max = 105.534;}
    if (NO_Angle == 18) {t_min = 109.793; t_max = 111.893;}
    if (NO_Angle == 19) {t_min = 116.067; t_max = 118.019;}
    if (NO_Angle == 20) {t_min = 121.838; t_max = 123.686;}
    if (NO_Angle == 21) {t_min = 127.273; t_max = 129.033;}
    if (NO_Angle == 22) {t_min = 132.400; t_max = 134.031;}
    if (NO_Angle == 23) {t_min = 137.230; t_max = 138.679;}
    double delta_cut = Z_MyTool_GetDeltaTheta_mpi0(t_min,6.75);
    Set_Angle_Cut(t_min-delta_cut, t_min+delta_cut);
    string out_name = "doc/" + dir_name + "_mpi0.txt";
    if( Exec( dir_name, out_name, 2, true) ) return 1;
    return 0;
}

//*****************************************************************************************//
int Exec(string dir_name, string out_name, Int_t NGamma, bool IsSplit){
    //IsSplit: Whether shower separation is required
    //NGamma: Number of photons produced
    ofstream out;
    out.open(out_name,ios::out);
    
    TString file_path_sim = "../data/"+dir_name+"/evtcomplete_sim.root";
    TString file_path_digi = "../data/"+dir_name+"/evtcomplete_digi.root";
    
    FairRunAna *fRun = new FairRunAna();
    TFile* file = new TFile(file_path_sim);
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TClonesArray* fMCTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMCTrackArray) return -1;
    TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
    if (!fHitArray) return -1;
    
    TFile* f = new TFile(file_path_digi);
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
    t->SetBranchAddress("EmcBump",&fBumpArray);
    if (!fBumpArray) return -1;
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    
    int N(0);
    //Int_t maxEvtNo = ioman->CheckMaxEventNo();
    Int_t maxEvtNo = t->GetEntries();
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
    	if (maxEvtNo>=100 && ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        
        //if (distance > 5) continue;
        if (nclusters != 1 || nbumps != 2) continue;
        //Calculate the error of energy and position
        
        PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(0);
        TVector3 pi0_mom = mcTrack->GetMomentum();
        double pi0_theta = pi0_mom.Theta();
        if ((C_min >0 && C_max >0)&&(pi0_theta < C_min || pi0_theta >C_max)) continue;

        Double_t PX, PY, PZ, PE;
        PX = PY = PZ = PE = 0.0;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(iGamma);
            Double_t bump_E = Bump->energy();
            TVector3 bump_mom = Bump->where();
            bump_mom.SetMag(bump_E);
            PX += bump_mom.X();
            PY += bump_mom.Y();
            PZ += bump_mom.Z();
            PE += bump_E;
            //out << bump_mom.X() << " " << bump_mom.Y() << " " << bump_mom.Z() << " " << bump_E << endl;
        }
        out << sqrt(PE*PE - PX*PX - PY*PY - PZ*PZ) << endl;
        N++;
    }
    out.close();
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
