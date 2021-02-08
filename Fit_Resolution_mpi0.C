int Fit_Resolution_mpi0(Int_t NO_Angle){
    ostringstream out1;
    out1 << NO_Angle;
    string str_NO_Angle = out1.str();
    std::string out_name = "doc/mpi0_A"+str_NO_Angle+"_resolution_par.txt";
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
    
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(-0.15);
    c1->SetTopMargin(-0.15);
    c1->SetBottomMargin(0.15);
    
    TGraph *g1 = new TGraph();
    g1->SetMarkerStyle(20);
    g1->SetMarkerColorAlpha(kRed-3, 1);
    g1->GetXaxis()->SetTitle("Energy");
    g1->GetYaxis()->SetTitle("Resolution");
    g1->GetXaxis()->CenterTitle();
    g1->GetYaxis()->CenterTitle();
    
    TGraph *g2 = new TGraph();
    g2->SetMarkerStyle(22);
    g2->SetMarkerColorAlpha(kBlack+1, 1);
    g2->GetXaxis()->SetTitle("Energy");
    g2->GetYaxis()->SetTitle("Resolution");
    g2->GetXaxis()->CenterTitle();
    g2->GetYaxis()->CenterTitle();
    
    int const N_point = 36;
    double R1_x[N_point], R1_y[N_point], R1_Ex[N_point], R1_Ey[N_point], R1_Ey_m[N_point], R1_Ey_p[N_point];
    double R2_x[N_point], R2_y[N_point], R2_Ex[N_point], R2_Ey[N_point], R2_Ey_m[N_point], R2_Ey_p[N_point];
    string str;
    Int_t N = 0;
    while (std::getline(par_file, str)) {
        std::stringstream strStream(str);
        float x_1, E_x1, R1, ER1_m, ER1_p, x_2, E_x2, R2, ER2_m, ER2_p;
        strStream >> x_1 >> R1 >> ER1_m >> x_2 >> R2 >> ER2_m;
        ER1_p = ER1_m;
        ER2_p = ER2_m;
        R1_x[N] = x_1;
        R1_y[N] = 100*R1;
        R1_Ex[N] = 0.0;
        //R1_Ex[N] = E_x1;
        R1_Ey[N] = 100*(ER1_p > ER1_m ? ER1_p : ER1_m);
        R2_x[N] = x_2;
        R2_y[N] = 100*R2;
        R2_Ex[N] = 0.0;
        //R2_Ex[N] = E_x2;
        R2_Ey[N] = 100*(ER2_p > ER2_m ? ER2_p : ER2_m);
        N++;
    }
    
    TGraphErrors *gr1 = new TGraphErrors(N_point,R1_x,R1_y,R1_Ex,R1_Ey);
    gr1->SetMarkerStyle(22);
    gr1->SetMarkerColorAlpha(kBlack, 1);
    gr1->GetXaxis()->SetTitle("distance   [cm]");
    gr1->GetYaxis()->SetTitle("#sigma(E_{#gamma})/E_{#gamma}   [%]");
    gr1->GetYaxis()->SetTitleOffset(1.0);
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
    
    TGraphErrors *gr2 = new TGraphErrors(N_point,R2_x,R2_y,R2_Ex,R2_Ey);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColorAlpha(kRed, 1);
    gr2->GetXaxis()->SetTitle("distance   [cm]");
    gr2->GetYaxis()->SetTitle("#sigma(E_{#gamma})/E_{#gamma}   [%]");
    gr2->GetYaxis()->SetTitleOffset(1.0);
    gr2->GetXaxis()->CenterTitle();
    gr2->GetYaxis()->CenterTitle();
    
    //gr1->GetYaxis()->SetRangeUser(1.4,5);
    //gr2->GetYaxis()->SetRangeUser(1.4,5);
    
    gr1->Draw("AP");
    gr2->Draw("Psame");
    
    TLegend * leg1 = new TLegend(0.554108,0.76,0.824649,0.89037);
    leg1->AddEntry(gr1, "Raw algorithm", "EP");
    leg1->AddEntry(gr2, "New algorithm", "EP");
    leg1->SetBorderSize(1);
    leg1->Draw("SAME");
    
    return 0;
}
