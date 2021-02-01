int Fit_Resolution_range(Int_t NO_Angle = 12, double Energy = 1.0, const int NBins = 100){
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/range_A"+str_NO_Angle+ "_E" + str_Energy + "_resolution_par.txt";
    std::ifstream par_file;
    par_file.open(out_name,std::ios::in);
    
    TString ts = "Angle "+str_NO_Angle;
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
    
    double cut_min[NBins], cut_max[NBins];
    double mean_OR[NBins], Resolution_OR[NBins], D_Resolution_OR[NBins];
    double mean_fix[NBins], Resolution_fix[NBins], D_Resolution_fix[NBins];
    double cut[NBins], zero[NBins];
    string str;
    Int_t N = 0;
    while (std::getline(par_file, str)) {
        std::stringstream strStream(str);
        double tmp_cut_min, tmp_cut_max, tmp_mean_OR, tmp_Resolution_OR, tmp_D_Resolution_OR, tmp_mean_fix, tmp_Resolution_fix, tmp_D_Resolution_fix;
        strStream >> tmp_cut_min >> tmp_cut_max >> tmp_mean_OR >> tmp_Resolution_OR >> tmp_D_Resolution_OR >> tmp_mean_fix >> tmp_Resolution_fix >> tmp_D_Resolution_fix;
        cut_min[N] = tmp_cut_min;
        cut_max[N] = tmp_cut_max;
        mean_OR[N] = tmp_mean_OR;
        Resolution_OR[N] = 100.0*tmp_Resolution_OR;
        D_Resolution_OR[N] = 100.0*tmp_D_Resolution_OR;
        mean_fix[N] = tmp_mean_fix;
        Resolution_fix[N] = 100.0*tmp_Resolution_fix;
        D_Resolution_fix[N] = 100.0*tmp_D_Resolution_fix;
        cut[N] = (tmp_cut_min+tmp_cut_max)/2.0;
        zero[N] = 0.0;
        N++;
    }
    
    TGraphErrors *gr1 = new TGraphErrors(NBins,cut,Resolution_OR,zero,D_Resolution_OR);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColorAlpha(kBlack, 1);
    gr1->GetXaxis()->SetTitle("d   [cm]");
    gr1->GetYaxis()->SetTitle("#sigma(E_{#gamma})/E_{#gamma}   [%]");
    gr1->GetYaxis()->SetTitleOffset(1.0);
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
    
    TGraphErrors *gr2 = new TGraphErrors(NBins,cut,Resolution_fix,zero,D_Resolution_fix);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColorAlpha(kRed, 1);
    gr2->GetXaxis()->SetTitle("d   [cm]");
    gr2->GetYaxis()->SetTitle("#sigma(E_{#gamma})/E_{#gamma}   [%]");
    gr2->GetYaxis()->SetTitleOffset(1.0);
    gr2->GetXaxis()->CenterTitle();
    gr2->GetYaxis()->CenterTitle();
    
    //gr1->GetYaxis()->SetRangeUser(1.4,5);
    //gr2->GetYaxis()->SetRangeUser(1.4,5);
    
    gr1->Draw("AP");
    gr2->Draw("Psame");
    
    /*
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
    */
    TLegend * leg1 = new TLegend(0.61,0.72,0.88,0.85);
    leg1->AddEntry(gr1, "Raw algorithm", "EP");
    leg1->AddEntry(gr2, "New algorithm", "EP");
    leg1->Draw("SAME");
    
    return 0;
}
