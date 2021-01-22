/*Double_t FABC(Double_t x,Double_t A, Double_t B, Double_t C, Double_t p1, Double_t p2){
    //p2 *= (1.0-(1.0-A)*exp(-pow(x/B,C)));
    p2 *= (1.0-A*exp(-pow(x/B,C)));
    //c2 *= 4*(1-exp(-A*pow(x,3)));
    return p1*exp(-p2*x);
}*/
const Double_t X0 = 0.89;
const Double_t RM = 2.00;
Double_t FABC(Double_t t,Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4){
    Double_t xi = t - p2*t*exp(-pow(X0*t/p3/RM,p4));
    return p0*exp(-p1*xi*X0/RM);
}

int Fit_DigiEnergy_test6( Int_t NO_Angle = 7, Double_t Energy = 1.0){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    
    TCanvas* c1=new TCanvas("PANDA1","test1",800,600);
    //TCanvas* c2=new TCanvas("PANDA2","test2",800,600);
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
    
    //std::string file_name(title);
    std::string file_name = "doc/WorkData_1Gamma_A7_E1.0_OR_R.txt";
    std::ifstream in_file;
    in_file.open(file_name,std::ios::in);
    
    TGraph *g = new TGraph();
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
    g->GetYaxis()->SetTitle("E_{digi}");
    g->GetXaxis()->SetTitle("distance");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    Double_t distance_cut = 3;
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        //if (angle>10 || angle<0) continue;
        if (distance > distance_cut) continue;
        g->SetPoint(N,distance,(energy));
        //g->SetPoint(N,distance,energy);
        N++;
    }
    
    //TF1* f=new TF1("f1","TMath::Log(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    TF1* f=new TF1("f1","(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    f->SetParameters(Energy, 2.5, 0.9, 0.7, 3);
    //f->SetParLimits(0, 0.5, 1);
    f->SetParLimits(2, 0, 1);
    //f->SetParLimits(2, 1.01, 25);

    //f->SetParameters(1, 0.6, 5, 3.37, 1.45);
    c1->cd();
    g->Draw("AP.");
    g->Fit(f,"R");
    
    in_file.close();
    return 0;
}
