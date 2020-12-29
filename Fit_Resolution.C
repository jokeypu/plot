int Fit_Resolution(Int_t NO_Angle){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1;
    out1 << NO_Angle;
    string str_NO_Angle = out1.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_resolution_par.txt";
    std::ifstream par_file;
    par_file.open(out_name,std::ios::in);
    
    TString ts = "Angle "+str_NO_Angle+"0 deg";
    TCanvas* c1=new TCanvas("PANDA1",ts,1000,700);
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
    //gStyle->SetOptFit(1111);
    
    TGraph *g1 = new TGraph();
    g1->SetMarkerStyle(20);
    g1->SetMarkerColorAlpha(kRed-3, 1);
    g1->GetXaxis()->SetTitle("Energy");
    g1->GetYaxis()->SetTitle("Resolution");
    g1->GetXaxis()->CenterTitle();
    g1->GetYaxis()->CenterTitle();
    
    TGraph *g2 = new TGraph();
    g2->SetMarkerStyle(22);
    g2->SetMarkerColorAlpha(kGreen+1, 1);
    g2->GetXaxis()->SetTitle("Energy");
    g2->GetYaxis()->SetTitle("Resolution");
    g2->GetXaxis()->CenterTitle();
    g2->GetYaxis()->CenterTitle();
    
    double R1_x[12], R1_y[12], R1_Ex[12], R1_Ey[12];
    double R2_x[12], R2_y[12], R2_Ex[12], R2_Ey[12];
    string str;
    Int_t N = 0;
    while (std::getline(par_file, str)) {
        std::stringstream strStream(str);
        float energy, R1, ER1, R2, ER2;
        strStream >> energy >> R1 >> ER1 >> R2 >> ER2;
        R1_x[N] = energy;
        R1_y[N] = R1;
        R1_Ex[N] = 0.0;
        R1_Ey[N] = ER1;
        R2_x[N] = energy;
        R2_y[N] = R2;
        R2_Ex[N] = 0.0;
        R2_Ey[N] = ER2;
        N++;
    }
    
    TGraphErrors *gr1 = new TGraphErrors(12,R1_x,R1_y,R1_Ex,R1_Ey);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerColorAlpha(kGreen+1, 1);
    gr1->GetXaxis()->SetTitle("Energy");
    gr1->GetYaxis()->SetTitle("Resolution");
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
    
    TGraphErrors *gr2 = new TGraphErrors(12,R2_x,R2_y,R2_Ex,R2_Ey);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColorAlpha(kRed-3, 1);
    gr2->GetXaxis()->SetTitle("Energy");
    gr2->GetYaxis()->SetTitle("Resolution");
    gr2->GetXaxis()->CenterTitle();
    gr2->GetYaxis()->CenterTitle();
    
    gr2->Draw("LAP");
    gr1->Draw("LPsame");
    
    TLegend * leg1 = new TLegend(0.61,0.72,0.88,0.85);
    leg1->AddEntry(gr1, "Old algorithm", "P");
    leg1->AddEntry(gr2, "New algorithm", "P");
    leg1->Draw("SAME");
    
    return 0;
}
