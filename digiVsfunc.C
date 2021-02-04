struct INTEGRAL {
    Double_t y, ci, L0;
    void SetPar(Int_t index, Double_t par){
        if (index == 1) L0 = par;
        else if (index == 2) ci = par;
        else if (index == 3) y = par;
        else std::cout << "Parameter error!!" << std::endl;
    }
    Double_t func_Int(Double_t a, Double_t b){
        float value(0);
        if (y < 0.1){
            Int_t N = 3+100*(b-a);
            float step = (b-a)/N, k = a+step, m = b-step/2;
            for (float i = k; i < m; i+=step) value += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            return step*((y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y))/2+value);
        }else{
            Int_t n = 4;
            float step = (b-a)/(2*n), Twostep = 2*step, StepOver2 = step/2;
            float sum1(0), sum2(0), k1 = a + step, k2 = k1 + step, m1 = b - StepOver2, m2 = m1 - step;
            for (float i = k1; i < m1; i += Twostep) sum1 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            for (float i = k2; i < m2; i += Twostep) sum2 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            return step/3*(y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y)+4*sum1+2*sum2);
        }
    }
    Double_t line_Int(Double_t x_n, Double_t y_n){
        SetPar(3,y_n);
        if ( x_n-L0 > 0) return func_Int(x_n-L0,x_n+L0);
        else return 2*func_Int(0,L0-x_n)+func_Int(L0-x_n,x_n+L0);
    }
    Double_t block_Int(Double_t distance, Double_t angle, Double_t ci){
        Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
        Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
        SetPar(2,ci);
        Int_t w1(1),w3(1);
        if (xm>0) { w1 = -1; if (ym>0) w3 = -1; }
        return w1*line_Int(y0,fabs(xm))+line_Int(y0,xp)+w3*line_Int(x0,fabs(ym))+line_Int(x0,yp);
    }
    Double_t shower_Digi(Double_t distance,Double_t angle){
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        else if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
        angle *= TMath::DegToRad();
        //Double_t p[9] = {L_1, L_2, L_3, p0, p1, p2, p3, p4, p5};
        Double_t p[9] = {1.57344, 0.943614, 1.20918, 0.139749, 4.94142, 1.18306, 0.494592, 3.37466, 1.89294};
        SetPar(1,p[0]);
        Double_t T1 = block_Int(distance,angle,p[4])/p[4];
        SetPar(1,p[1]);
        Double_t T2 = p[5]*block_Int(distance,angle,p[6])/p[6];
        SetPar(1,p[2]);
        Double_t T3 = p[7]*block_Int(distance,angle,p[8])/p[8];
        return p[3]*(T1+T2+T3)/3;
    }
}Shower_Function;




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
        Double_t dx = (*DetPos-*Cent).Dot(ex);
        Double_t dy = (*DetPos-*Cent).Dot(ey);
        TVector2 vv(dx,dy);
        distance = sqrt(dx*dx+dy*dy);
        angle = fabs(TMath::RadToDeg()*vv.Phi_mpi_pi(vv.Phi()));
        //if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        //if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    }
    return angle;
}

Double_t Getfunc(Double_t angle);

int digiVsfunc(std::string dir_name){
    //std::string dir_name="Gamma_one_1G";
    //File_out << "Index " << "Distance " << "Angle " << "Energy" << endl;
    
    Int_t NGamma(1);
    
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
    TClonesArray* fDigiArray = new TClonesArray("PndEmcDigi");
    t->SetBranchAddress("EmcDigi",&fDigiArray);
    if (!fBumpArray) return -1;
    TClonesArray* fClusterArray = new TClonesArray("PndEmcCluster");
    t->SetBranchAddress("EmcCluster",&fClusterArray);
    if (!fClusterArray) return -1;
    
    TH1D *h1 = new TH1D("h1","hist1",100,0,0.7);
    h1->SetLineColor(kRed);
    
    TH1D *h2 = new TH1D("h2","hist2",100,0,0.7);
    h2->SetLineColor(kBlue);
    
    
    int N(0);
    Int_t maxEvtNo = t->GetEntries();
    //maxEvtNo /= 10;
    for (Int_t ievt = 0; ievt < maxEvtNo; ievt++) {
        ioman->ReadEvent(ievt); // read event by event
        t->GetEntry(ievt);
        if (maxEvtNo>=100 && ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        int nhits = fHitArray->GetEntriesFast();
        int nclusters = fClusterArray->GetEntriesFast();
        int nbumps = fBumpArray->GetEntriesFast();
        int ndigis = fDigiArray->GetEntriesFast();
        
        //Get the momentum of each photon
        std::vector<TVector3> Gamma_mom;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            PndMCTrack *mcTrack = (PndMCTrack *)fMCTrackArray->At(iGamma);
            Gamma_mom.push_back(mcTrack->GetMomentum());
        }
        
        //Calculate the average distance between photons
        Double_t distance(0);
        Int_t Ncunt(0);
        for (int iGamma = 0; iGamma < NGamma-1; iGamma++) {
            for (int jGamma = iGamma+1; jGamma < NGamma; jGamma++) {
                Double_t TheDistance = ((65.0/Gamma_mom[iGamma].Pt())*Gamma_mom[iGamma]-(65.0/Gamma_mom[jGamma].Pt())*Gamma_mom[jGamma]).Mag();
                distance += TheDistance;
                Ncunt++;
            }
        }
        distance /= Ncunt;
        
        //Exclude events generated electron-positron
        std::map<Int_t, bool> Exist;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
            for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) {
                for (int iGamma = 0; iGamma < NGamma; iGamma++)
                    if (linkIter->GetIndex() == iGamma) Exist[iGamma] = true;
            }
        }
        if (Exist.size() != NGamma) continue;
        
        //Get the true energy of each shower
        std::map<Int_t, Double_t> truth_E;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) truth_E[iGamma] = 0.0;
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            std::map<Int_t, Double_t>  dep = hit->GetDepositedEnergyMap();
            std::map<Int_t, Double_t>::iterator ptr;
            for ( ptr = dep.begin(); ptr != dep.end(); ptr++) {
                for (int iGamma = 0; iGamma < NGamma; iGamma++)
                    if (ptr->first == iGamma) truth_E[iGamma] += ptr->second;
            }
        }
        
        //Match bump for each photon
        std::vector<Int_t> match;
        for (int iGamma = 0; iGamma < NGamma; iGamma++) {
            Double_t min_d(99999);
            Int_t index(-1);
            for (int i = 0; i < nbumps; i++) {
                PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(i);
                TVector3 pos = Bump->position();
                Double_t d = pos.Mag()*sin(Gamma_mom[iGamma].Angle(pos));
                if (d < min_d) { min_d = d; index = i; }
            }
            if ( index == -1 ) return 1;
            match.push_back(index);
        }
        
        PndEmcBump* Bump = (PndEmcBump*)fBumpArray->At(0);
        PndEmcCluster *theCluster = (PndEmcCluster *)fClusterArray->At(0);
        std::map<Int_t, Int_t> theMaximaDigis = theCluster->LocalMaxMap();
        
        std::map<Int_t, Int_t>::iterator p;
        Int_t digi_seed(-1),digi_seed_id(-1);
        int c(0);
        for (p = theMaximaDigis.begin(); p != theMaximaDigis.end(); p++){
            digi_seed = p->second;
            digi_seed_id = p->first;
            c++;
        }
        if ( c!=1 ) continue;
        
        PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(digi_seed);
        double Seed_Energy = digi->GetEnergy();
        TVector3 Seed_pos = digi->where();
        TVector3 Cent_pos = Bump->where();
        
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            //if (hit->GetDetectorID() == digi_seed_id) continue;
            TVector3 Det_Pos;
            double Truth_Energy = hit->GetEnergy();
            double Digi_Energy = -1;
            for (int j = 0 ; j < ndigis ; j++){
                PndEmcDigi* idigi = (PndEmcDigi*)fDigiArray->At(j);
                if (idigi->GetDetectorId() == hit->GetDetectorID()) {
                    Digi_Energy = idigi->GetEnergy();
                    Det_Pos = idigi->where();
                }
            }
            if (Digi_Energy == -1) continue;
            double Distance = DD(&Det_Pos, &Cent_pos, 1.25);
            double angle = AA(&Det_Pos, &Cent_pos, 1.25);
            if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
            if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
            if (Distance > 14) continue;
            N++;
            if (Distance >1.45 && Distance<1.55) {
                h1->Fill(Digi_Energy);
                //h2->Fill(Shower_Function.shower_Digi(Distance,angle));
                h2->Fill(Getfunc(angle));
            }
        }
    }
    h2->Draw();
    h1->Draw("same");
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}

Double_t Getfunc(Double_t angle){
    Double_t a1 =  1.125, b1 = 0.183206,
    a2 = 3.375, b2 =  0.201581,
    a3 =5.625, b3 =  0.189013,
    a4 =7.875, b4 =  0.218617,
    a5 =10.125, b5 =  0.204281,
    a6 =12.375, b6 =  0.215612,
    a7 =14.625, b7 =  0.232127,
    a8 =16.875, b8 =  0.228904,
    a9 =19.125, b9 =  0.249902,
    a10 =21.375, b10 = 0.263043,
    a11 =23.625, b11 =  0.240984,
    a12 =25.875, b12 =  0.267332,
    a13 =28.125, b13 =  0.311402,
    a14 =30.375, b14 =  0.314986,
    a15 =32.625, b15 =  0.313625,
    a16 =34.875, b16 =  0.32179,
    a17 =37.125, b17 =  0.324177,
    a18 =39.375, b18 =  0.38823,
    a19 =41.625, b19 =  0.37198,
    a20 =43.875, b20 =  0.351777;
    
    if (angle >=0 && angle < a1) return b1;
    else if ( angle >=a1 && angle < a2 ) return b1+(b2-b1)*((angle-a1)/(a2-a1));
    else if ( angle >=a2 && angle < a3 ) return b2+(b3-b2)*((angle-a2)/(a3-a2));
    else if ( angle >=a3 && angle < a4 ) return b3+(b4-b3)*((angle-a3)/(a4-a3));
    else if ( angle >=a4 && angle < a5 ) return b4+(b5-b4)*((angle-a4)/(a5-a4));
    else if ( angle >=a5 && angle < a6 ) return b5+(b6-b5)*((angle-a5)/(a6-a5));
    else if ( angle >=a6 && angle < a7 ) return b6+(b7-b6)*((angle-a6)/(a7-a6));
    else if ( angle >=a7 && angle < a8 ) return b7+(b8-b7)*((angle-a7)/(a8-a7));
    else if ( angle >=a8 && angle < a9 ) return b8+(b9-b8)*((angle-a8)/(a9-a8));
    else if ( angle >=a9 && angle < a10 ) return b9+(b10-b9)*((angle-a9)/(a10-a9));
    else if ( angle >=a10 && angle < a11 ) return b10+(b11-b10)*((angle-a10)/(a11-a10));
    else if ( angle >=a11 && angle < a12 ) return b11+(b12-b11)*((angle-a11)/(a12-a11));
    else if ( angle >=a12 && angle < a13 ) return b12+(b13-b12)*((angle-a12)/(a13-a12));
    else if ( angle >=a13 && angle < a14 ) return b13+(b14-b13)*((angle-a13)/(a14-a13));
    else if ( angle >=a14 && angle < a15 ) return b14+(b15-b14)*((angle-a14)/(a15-a14));
    else if ( angle >=a15 && angle < a16 ) return b15+(b16-b15)*((angle-a15)/(a16-a15));
    else if ( angle >=a16 && angle < a17 ) return b16+(b17-b16)*((angle-a16)/(a17-a16));
    else if ( angle >=a17 && angle < a18 ) return b17+(b18-b17)*((angle-a17)/(a18-a17));
    else if ( angle >=a18 && angle < a19 ) return b18+(b19-b18)*((angle-a18)/(a19-a18));
    else if ( angle >=a19 && angle < a20 ) return b19+(b20-b19)*((angle-a19)/(a20-a19));
    else return b20;
    
}

