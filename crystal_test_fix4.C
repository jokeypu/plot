struct myfunc {
    std::vector<double> par1 = {8.62369, 5.32678, 12.0, 0.696003, 22.8543};
    std::vector<double> par2 = {83.4757, -0.171518, 1.88588, 0.199225, 23.0441};
    std::vector<double> par3 = {0.0227033, 4.97094, 0.772748, 4.0, 0.0018315};
    std::vector<double> par4 = {0.00790974, 1.56139, 1.39, 0.428328, 0.00115601};
    std::vector<double> par5 = {-0.0849149, 2.58396, 0.419646};
    std::vector<double> par6 = {10, 0.0124098, 0.0780432, 0.182612};
    std::vector<double> par7 = {0.144939, -0.435278, 0.0642399};
    Double_t func_x0(Double_t d){
        if (d < 1.7) return par1[0]*TMath::Vavilov(d - par1[1], par1[2], par1[3]) + par1[4];
        else return par2[0]*TMath::Landau(d - par2[1], par2[2], par2[3]) + par2[4];
    }
    Double_t func_a(Double_t d){
        if (d < 1.39) return par3[0]*TMath::Poisson(par3[1]*(d-par3[2]), par3[3])+par3[4];
        else return par4[0]*TMath::Poisson(par4[1]*(d-par4[2]), par4[3])+par4[4];
    }
    Double_t func_h(Double_t d){
        if (d < 1.4 ) return par5[0]*pow(d,par5[1])+par5[2];
        else if ( d < 3.5) return par6[0]*TMath::Landau(d-par6[1],par6[2],par6[3]);
        else return par7[0]*exp(par7[1]*d+par7[2]);
    }
    Double_t m(Double_t d, Double_t angle){
        Double_t h = func_h(d);
        if ((d < 0.8) || (d > 2.8)) return h;
        else{
            Double_t a = func_a(d);
            Double_t x0 = func_x0(d);
            a *= a;
            angle = abs(fmod(angle,45.0) - 45*(((int)(angle/45.0))%2));
            if (angle<x0) return a*angle*angle+h;
            else {cout<< "Angle:" << angle << ", " << a*x0*x0+(a/(45/x0-1))*(x0-45)*(x0-45)-(a/(45/x0-1))*(angle-45)*(angle-45)<<
                "+" << h << endl;return a*x0*x0+(a/(45/x0-1))*(x0-45)*(x0-45)-(a/(45/x0-1))*(angle-45)*(angle-45)+h;}
        }
    }
    Double_t m(Double_t distance, Double_t angle, Double_t par){
        if ( angle > 90 && angle <= 180 ) angle = 180 -angle;
        if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
        if ( distance < 4 ) {
            Double_t value = (1.585 - 0.0237*angle)*pow(distance,5) - (1.371/(angle + 24.672))*(pow(distance,4)+22.845*pow(distance,3) - 3.794*distance*distance) + 0.485*distance + 2.787;
            return 2.456/value;
        }else return exp(-1* par * distance);
    }
}func;

std::vector<double> par1 = {1.47994, 1.40911, 1.5216, 1.58767, 1.41823, 1.20107, 1.40021, 1.61076, 1.55799, 1.69814, 1.44313, 1.36737, 1.44847, 1.50126, 1.30191, 1.50295, 1.68052, 1.53869, 1.32379, 1.44949, 1.50036, 1.31779, 1.28991, 1.55822, 1.37487, 1.11861, 1.59666, 1.35665, 1.42687, 1.45577, 1.45591, 1.21475, 1.00299, 1.30838, 1.05776, 1.46028, 1.20945, 1.12975, 1.29005, 1.12541, 1.30164, 0.89004, 1.03526, 1.04499, 1.42351, 1.2598, 1.15151, 1.03511, 1.35549, 0.858096, 1.10612, 0.965125, 1.10845, 0.89562, 1.03077, 1.13563, 1.2431, 1.12137, 0.936256, 0.68764, 1.19257, 0.78799, 0.954587, 0.830432, 1.04553, 1.13091, 0.89148, 0.896716, 1.04182, 0.999725, 0.827203, 0.66813, 0.768195, 0.789239, 0.880294, 0.748822, 0.749049, 0.898385, 0.917733, 0.657576, 0.76885, 1.00577, 1.06092, 0.509508, 0.783106, 0.973449, 0.847542, 0.758128, 0.721299, 0.631344, 0.874684, 0.83243, 0.855176, 0.782924, 0.834866, 0.834289, 0.772338, 0.848204, 0.804116, 0.765132, 0.716634, 0.720747, 0.905901, 0.669254, 0.781137, 0.995052, 0.908544, 0.619489, 0.725752, 0.841246, 0.773001, 0.782085, 0.520848, 0.682983, 0.678196, 0.635669, 0.632214, 0.783848, 0.729467, 0.820897, 0.728379, 0.695053, 0.722519, 0.580005, 0.661338, 0.653723, 0.594421, 0.468908, 0.438184, 0.535993, 0.572852, 0.742575, 0.579848, 0.495763, 0.601748, 0.539577, 0.475001, 0.616259, 0.478793, 0.499877, 0.722398, 0.495558, 0.3571, 0.441352, 0.628777, 0.609481, 0.513991, 0.588428, 0.601055, 0.446958, 0.622156, 0.49616, 0.643373, 0.458008, 0.656871, 0.625916, 0.608822, 0.562916, 0.671652, 0.582634, 0.539649, 0.787993, 0.579912, 0.627778, 0.583095, 0.526202, 0.710461, 0.635143, 0.585342, 0.635929, 0.639902, 0.690184, 0.724497, 0.55786, 0.688477, 0.630092, 0.549355, 0.643222, 0.721614, 0.640634, 0.600949, 0.68907, 0.674864, 0.759479, 0.655113, 0.828516, 0.673711, 0.944664, 0.767714, 0.720903, 0.701229, 0.850507, 0.677267, 0.694624, 0.8888, 0.762769, 0.658335, 0.791415, 0.78603, 0.627016, 0.651832, 0.873996, 0.882778, 0.722487, 0.756483, 0.607489, 0.840972, 0.679022, 1.00984, 0.766365, 0.762409, 0.851544, 1.05027, 0.871671, 0.909959, 1.05971, 0.873106, 0.78498, 0.85526, 0.854952, 0.893633, 0.759977, 1.11885, 0.726261, 1.10145, 1.00594, 0.936423, 0.898264, 0.566903, 1.01108, 0.978162, 0.97062, 0.972432, 0.790211, 0.949527, 0.677132, 0.849114, 1.20114, 1.02572, 0.829375, 1.01202, 0.856844, 1.00655, 0.908292, 1.10332, 1.0529, 1.01916, 1.13544, 0.833919, 1.05235, 0.950042, 1.24953, 0.837041, 1.34852, 1.24002, 0.927093, 1.34141, 1.26009, 1.38737, 1.37204, 1.05675, 1.21017, 1.07819, 1.33668, 1.20768, 1.26512, 1.4366, 1.21145, 1.31817, 1.25951, 1.31265, 1.18841, 1.41653, 1.12476, 1.17037, 1.31719, 1.2444, 1.24971, 1.26595, 1.42063, 1.32398, 1.37531, 1.41557, 1.48662, 1.42628, 1.4907, 1.47549, 1.47494, 1.3504, 1.34299, 1.51453, 1.41799, 1.37948, 1.3809, 1.30753, 1.51892, 1.56084, 1.29576, 1.47021, 1.52589, 1.52589};
std::vector<double> par2 = {2.66896, 2.49278, 2.80345, 2.89587, 2.53116, 2.03032, 2.43975, 2.95661, 2.77252, 3.15394, 2.55249, 2.42622, 2.65767, 2.72616, 2.2347, 2.73141, 3.05106, 2.86981, 2.20703, 2.63139, 2.70011, 2.28799, 2.30774, 2.83522, 2.47116, 1.8769, 2.91516, 2.45296, 2.65155, 2.65813, 2.78106, 2.14554, 1.55177, 2.2633, 1.7169, 2.78115, 2.04959, 2.00416, 2.30041, 1.95049, 2.31355, 1.22492, 1.56356, 1.79032, 2.63897, 2.28102, 2.0282, 1.67633, 2.45868, 1.35526, 1.82565, 1.61954, 2.09909, 1.51749, 1.82267, 2.06525, 2.39333, 2.03185, 1.64258, 0.95157, 2.18628, 1.2758, 1.62803, 1.32005, 1.86366, 2.11769, 1.40808, 1.4002, 1.76142, 1.76, 1.36814, 1.00798, 1.15659, 1.35867, 1.468, 1.31071, 1.13268, 1.52369, 1.51906, 0.939691, 1.28152, 1.83448, 2.05368, 0.637378, 1.35827, 1.90374, 1.4623, 1.23614, 1.15256, 0.95972, 1.66724, 1.40318, 1.54213, 1.36875, 1.56704, 1.51898, 1.40027, 1.56842, 1.53894, 1.21049, 1.33156, 1.24108, 1.62466, 1.06124, 1.47042, 1.91208, 1.81024, 1.03995, 1.36047, 1.58137, 1.47129, 1.45371, 0.827003, 1.12369, 1.27271, 1.11212, 1.17264, 1.55889, 1.45651, 1.63882, 1.3695, 1.34007, 1.40464, 1.01638, 1.27249, 1.20306, 1.13775, 0.840592, 0.743743, 0.94896, 0.978785, 1.41735, 1.10671, 0.927937, 1.15159, 1.0082, 0.870135, 1.16996, 0.869439, 0.924373, 1.46248, 0.926641, 0.465128, 0.73435, 1.325, 1.185, 0.953477, 1.13338, 1.1283, 0.785783, 1.20098, 0.853101, 1.34885, 0.754249, 1.29577, 1.27105, 1.22525, 1.1137, 1.37496, 1.04011, 0.967997, 1.6301, 1.05949, 1.25039, 1.13691, 0.992427, 1.48312, 1.1932, 1.08105, 1.2027, 1.2287, 1.34198, 1.4273, 1.0085, 1.35113, 1.2724, 0.949598, 1.26141, 1.43119, 1.13857, 1.1127, 1.31641, 1.26746, 1.46623, 1.2242, 1.63115, 1.18876, 1.90111, 1.48521, 1.39906, 1.3169, 1.63189, 1.18844, 1.23203, 1.71187, 1.34889, 1.15075, 1.46681, 1.47676, 0.959881, 1.04886, 1.60121, 1.62965, 1.31811, 1.37114, 0.941688, 1.53951, 1.12537, 1.97132, 1.31065, 1.28464, 1.58935, 2.00845, 1.5215, 1.67351, 1.92489, 1.62972, 1.39437, 1.48894, 1.52074, 1.57092, 1.19209, 2.12267, 1.11734, 2.16532, 1.68846, 1.66906, 1.4942, 0.745874, 1.77681, 1.7075, 1.79041, 1.71171, 1.20515, 1.47032, 0.974393, 1.38426, 2.1419, 1.75069, 1.4106, 1.67161, 1.32163, 1.70253, 1.57238, 1.98548, 1.82058, 1.75653, 2.05942, 1.38215, 1.89048, 1.54624, 2.25449, 1.15548, 2.40144, 2.29742, 1.47754, 2.45524, 2.20288, 2.51807, 2.51233, 1.72627, 2.18322, 1.80719, 2.45504, 2.11215, 2.24465, 2.6528, 2.10092, 2.28634, 2.17699, 2.26368, 2.04988, 2.53467, 1.91532, 2.07669, 2.30477, 2.16611, 2.21703, 2.11469, 2.61702, 2.29666, 2.38737, 2.49955, 2.71136, 2.58961, 2.748, 2.67885, 2.72404, 2.40714, 2.41498, 2.78601, 2.5474, 2.45347, 2.44354, 2.23669, 2.67769, 2.81414, 2.29152, 2.67774, 2.76884, 2.76884};
Double_t newfunc3(const TVector3 *DetPos, const TVector3 *Cent, const Double_t par){
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
    //distance = (distance-1.66);
    //if (distance<0) distance = 0;
    
    if ( distance < 2 ) return func.func_h(distance);
    else return exp(-1* par * distance);  //0.2*
}

Double_t rat(const TVector3 *DetPos_i, const TVector3 *DetPos_0, const TVector3 *Cent, const Double_t par) {
    Double_t distance(0);
    if (*DetPos_i != *Cent) {
        TVector3 vz(0, 0, 1);
        TVector3 DetPos_o = (*DetPos_i) - 3.7*vz;
        TVector3 DetPos_n;
        DetPos_n.SetMagThetaPhi(DetPos_o.Mag(), DetPos_o.Theta(), DetPos_o.Phi()-0.06981317);
        TVector3 ey = DetPos_n.Cross(vz).Unit();
        TVector3 ex = DetPos_n.Cross(ey).Unit();
        Double_t dx = abs((*Cent-*DetPos_i).Dot(ex));
        Double_t dy = abs((*Cent-*DetPos_i).Dot(ey));
        distance = sqrt(dx*dx+dy*dy);
    }
    Double_t value;
    if (distance < 2)
    value = newfunc3(DetPos_i, Cent, par)/newfunc3(DetPos_0, Cent, par);
    else value = newfunc3(DetPos_i, Cent, par);
    //std::cout << "value:" << value << std::endl;
    //if ( value >= 1.0 ) return 0.99;
    return value;
}
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


int Exec(TString dir_name, TH2D *h, Int_t NGamma=1);
int crystal_test_fix4( TString dir_name="Gamma_one_1G" )
{
    int bin1(400),bin2(400),bin3(150);
    float tx(800),ty(600);
    //double xmin(0),xmax(3.5),ymin(-0.6),ymax(0.6),zmin(0),zmax(0.1);
    double xmin(0),xmax(1),ymin(0),ymax(1),zmin(0),zmax(0.1);
    
    TCanvas* c1=new TCanvas("PANDA1","c1",tx,ty);
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
    
    //TH3D* h2D1 = new TH3D("Hist1","h1",bin1,xmin,xmax, bin2,ymin,ymax, bin3,zmin,zmax);
    TH2D* h2D = new TH2D("Hist1","h1",bin1,xmin,xmax,bin2,ymin,ymax);
    h2D->SetMarkerStyle(22);
    h2D->SetMarkerColorAlpha(kAzure+3, 0.5);
    h2D->GetYaxis()->SetTitle("E_{ci}");h2D->GetXaxis()->SetTitle("E_{truth}");
    //h2D->GetYaxis()->SetTitle("E_{ci}-E_{truth}");h2D->GetXaxis()->SetTitle("distance");
    //h1D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    
    if( Exec(dir_name, h2D, 1) ) return 1;
    
    c1->cd();
    c1->SetGridy();
    h2D->Draw("SCAT");
    
    TF1* f1=new TF1("f1","x",0,1);
    f1->Draw("SAME");
    
    TLegend * leg1 = new TLegend(0.61,0.72,0.88,0.85);
    leg1->AddEntry(h2D, "Crystal calculated", "P");
    leg1->Draw("SAME");
    
    return 0;
}

//*******************************************************************************************************//
int Exec(TString dir_name, TH2D *h, Int_t NGamma){
    //NGamma: Number of photons produced
    
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
    
    int N(0);
    Int_t maxEvtNo = t->GetEntries();
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
        //cout <<  "c:" << c << endl;
        
        PndEmcDigi* digi = (PndEmcDigi*)fDigiArray->At(digi_seed);
        double Seed_Energy = digi->GetEnergy();
        TVector3 Seed_pos = digi->where();
        TVector3 Cent_pos = Bump->where();
        //TVector3 Cent_pos = (65.0/Gamma_mom[0].Pt())*Gamma_mom[0];
        
        for (int i = 0; i < nhits; i++) {
            PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
            //if (hit->GetDetectorID() == digi_seed_id) continue;
            TVector3 Det_Pos(hit->GetX(), hit->GetY(), (hit->GetZ()));
            double Truth_Energy = hit->GetEnergy();
            double Digi_Energy = -1;
            for (int j = 0 ; j < ndigis ; j++){
                PndEmcDigi* idigi = (PndEmcDigi*)fDigiArray->At(j);
                if (idigi->GetDetectorId() == hit->GetDetectorID()) Digi_Energy = idigi->GetEnergy();
            }
            double Eci = Seed_Energy * rat(&Det_Pos, &Seed_pos, &Cent_pos, 1.25);
            double Distance = DD(&Det_Pos, &Cent_pos, 1.25);
            //if ((Eci - Truth_Energy) < 0.02) continue;
            h->Fill(Distance,Eci - Truth_Energy);
            //h->Fill(Truth_Energy,Eci);
            //h->Fill(Truth_Energy,Eci - Truth_Energy);
            if (abs(Eci - Truth_Energy) < 0.03) N++;
            //if (Digi_Energy >= 0) h->Fill(Eci - Digi_Energy);
        }
        cout <<  "c:" << c << endl;
        //N++;
    }
    cout << "Max Event Nomber:" << maxEvtNo << ", " << "Passed:" << N << endl;
    return 0;
}
