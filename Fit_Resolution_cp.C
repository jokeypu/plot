int Fit_Resolution_cp(Int_t NO_Angle = 7){
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
    g2->SetMarkerColorAlpha(kBlack+1, 1);
    g2->GetXaxis()->SetTitle("Energy");
    g2->GetYaxis()->SetTitle("Resolution");
    g2->GetXaxis()->CenterTitle();
    g2->GetYaxis()->CenterTitle();
    
    double R1_x[12], R1_y[12], R1_Ex[12], R1_Ey[12], R1_Ey_m[12], R1_Ey_p[12];
    double R2_x[12], R2_y[12], R2_Ex[12], R2_Ey[12], R2_Ey_m[12], R2_Ey_p[12];
    string str;
    Int_t N = 0;
    while (std::getline(par_file, str)) {
        std::stringstream strStream(str);
        float energy_1, R1, ER1_m, ER1_p, energy_2, R2, ER2_m, ER2_p;
        strStream >> energy_1 >> R1 >> ER1_m >> ER1_p >> energy_2 >> R2 >> ER2_m >> ER2_p;
        R1_x[N] = energy_1;
        R1_y[N] = R1;
        R1_Ex[N] = 0.0;
        R1_Ey[N] = (ER1_p > ER1_m ? ER1_p : ER1_m);
        R2_x[N] = energy_2;
        R2_y[N] = R2;
        R2_Ex[N] = 0.0;
        R2_Ey[N] = (ER2_p > ER2_m ? ER2_p : ER2_m);
        N++;
    }
    
    TGraphErrors *gr1 = new TGraphErrors(12,R1_x,R1_y,R1_Ex,R1_Ey);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColorAlpha(kBlack, 1);
    gr1->GetXaxis()->SetTitle("E_{#gamma}  (GeV/c^{2})");
    gr1->GetYaxis()->SetTitle("#sigma(E_{#gamma})/E_{#gamma}");
    gr1->GetYaxis()->SetTitleOffset(1.4);
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
    
    TGraphErrors *gr2 = new TGraphErrors(12,R2_x,R2_y,R2_Ex,R2_Ey);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColorAlpha(kRed, 1);
    gr2->GetXaxis()->SetTitle("E_{#gamma}  (GeV/c^{2})");
    gr2->GetYaxis()->SetTitle("#sigma(E_{#gamma})/E_{#gamma}");
    gr2->GetYaxis()->SetTitleOffset(1.4);
    gr2->GetXaxis()->CenterTitle();
    gr2->GetYaxis()->CenterTitle();
    
    gr1->GetYaxis()->SetRangeUser(0.014,0.05);
    gr2->GetYaxis()->SetRangeUser(0.014,0.05);
    
    gr1->Draw("AP");
    gr2->Draw("Psame");
    
    TF1 *f1 = new TF1("f1","[0]/sqrt(x) + [1]",0,6.5);
    f1->SetLineColor(kBlack);
    f1->SetLineWidth(2);
    f1->SetLineStyle(2);
    f1->SetParName(0, "#sigma_{0}");
    f1->SetParName(1, "#sigma_{1}");
    
    TF1 *f2 = new TF1("f2","[0]/sqrt(x) + [1]",0,6.5);
    f2->SetLineColor(kRed);
    f2->SetLineWidth(2);
    f2->SetLineStyle(2);
    f2->SetParName(0, "#sigma_{0}");
    f2->SetParName(1, "#sigma_{1}");
    
    gr1->Fit(f1,"R");
    gr2->Fit(f2,"R");
    
    ostringstream os1, os2, os3, os4;
    os1 << f1->GetParameter(0);
    os2 << f1->GetParameter(1);
    os3 << f2->GetParameter(0);
    os4 << f2->GetParameter(1);
    string par1_str = os1.str(), par2_str = os2.str(), par3_str = os3.str(), par4_str = os4.str();
    string sss1 = par1_str+"/#sqrt{E} + "+par2_str;
    string sss2 = par3_str+"/#sqrt{E} + "+par4_str;
    cout << sss1 << endl;
    cout << sss2 << endl;
    TPaveText *pt = new TPaveText(0.609218, 0.437037, 0.878758, 0.691852,"NDC");
    //TText *tt1 = pt->AddText("Raw");
    TText *t1 = pt->AddText("Raw:  #frac{#sigma(E_{#gamma})}{E_{#gamma}}  =  #frac{2.50%}{#sqrt{E_{#gamma}}} + 1.39%");
    pt->AddLine(.0,.5,1.,.5);
    //TText *tt2 = pt->AddText("New");
    TText *t2 = pt->AddText("New: #frac{#sigma(E_{#gamma})}{E_{#gamma}}  =  #frac{2.16%}{#sqrt{E_{#gamma}}} + 0.79%");
    
    //tt1->SetTextColor(kBlack);
    t1->SetTextColor(kBlack);
    //tt2->SetTextColor(kRed);
    t2->SetTextColor(kRed);
    //tt1->SetTextSize(0.04);
    //tt2->SetTextSize(0.04);
    t1->SetTextSize(0.028);
    t2->SetTextSize(0.028);
    pt->Draw();

    TLegend * leg1 = new TLegend(0.61,0.72,0.88,0.85);
    leg1->AddEntry(gr1, "Raw algorithm", "EP");
    leg1->AddEntry(gr2, "New algorithm", "EP");
    leg1->Draw("SAME");
    
    return 0;
}
