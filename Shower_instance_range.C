Double_t Novosibirsk(Double_t x,Double_t peak=0.,Double_t width=0.,Double_t tail=0.)
{
  if (TMath::Abs(tail) < 1.e-7) {
    return TMath::Exp( -0.5 * TMath::Power( ( (x - peak) / width ), 2 ));
  }

  Double_t arg = 1.0 - ( x - peak ) * tail / width;

  if (arg < 1.e-7) {
    //Argument of logarithm negative. Real continuation -> function equals zero
    return 0.0;
  }

  Double_t log = TMath::Log(arg);
  static const Double_t xi = 2.3548200450309494; // 2 Sqrt( Ln(4) )

  Double_t width_zero = ( 2.0 / xi ) * TMath::ASinH( tail * xi * 0.5 );
  Double_t width_zero2 = width_zero * width_zero;
  Double_t exponent = ( -0.5 / (width_zero2) * log * log ) - ( width_zero2 * 0.5 );

  return TMath::Exp(exponent);
}

Double_t par[5];

void SetPar(Double_t par0, Double_t par1, Double_t par2, Double_t par3, Double_t par4){
    par[0] = par0;
    par[1] = par1;
    par[2] = par2;
    par[3] = par3;
    par[4] = par4;
}
Double_t fit_func(Double_t x){return par[0]*TMath::Gaus(x,par[1],par[2])+par[3]*TMath::Gaus(x,par[1],par[4]);}

Double_t func_Int(Double_t mean, Double_t sigma){
    ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
    ROOT::Math::Functor1D wf(&fit_func);
    //ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kLEGENDRE,10,10,2);
    ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kGAUSS);
    ig.SetFunction(wf);
    Double_t value = ig.Integral(mean-sigma,mean+sigma);
    return value;
}

Double_t finding_sigma(Double_t mean, Double_t Init_sigma){
    double p = 0.0001;
    double cut = 0.6826 * func_Int(mean, 0.5);
    double sigma_shift = 9999;
    double sigma_max = 2*Init_sigma;
    double sigma_min = 0.0001;
    while (sigma_shift > 0.000001) {
        double value_min = func_Int(mean, sigma_min);
        double value_max = func_Int(mean, sigma_max);
        double sigma_mid = (sigma_max+sigma_min)/2.0;
        double value_mid = func_Int(mean, sigma_mid);
        if (cut < value_mid) sigma_max = sigma_mid;
        else sigma_min = sigma_mid;
        sigma_shift = fabs(sigma_max-sigma_min);
    }
    return (sigma_max+sigma_min)/2.0;
}

int Shower_instance_range(const char old_file[30], const char new_file[30], double Energy , int NO_Angle, double cut_min, double cut_max)
{
    const Double_t X0 = 0.89;
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/range_A"+str_NO_Angle+ "_E" + str_Energy + "_resolution_par.txt";
    std::ofstream par_file;
    par_file.open(out_name,std::ios::app);
    
    double cut_mean = fabs(cut_max + cut_min)/2.0;
    double dex = 22.709*exp(-1.11254*cut_mean)+0.139594;
    if (dex > 1.0) dex = 1.0;

    int bin1 = 120, bin2 = 200;
    float tx(1200),ty(900);
    double xmin = 0.98*Energy - dex*Energy, xmax = 0.98*Energy + dex*Energy;
    //double xmin(0),xmax(2.0);
    cout << "-INFO  Old File : " << old_file << endl;
    cout << "-INFO  New File : " << new_file << endl;
    
    string file_str1(old_file), file_str2(new_file);
    
    string file_name1 = "doc/" + file_str1 + "_range.txt";
    string file_name2 = "doc/" + file_str2 + "_range.txt";
    
    TCanvas* c1=new TCanvas("PANDA1","c1",tx,ty);
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
    
    TH1D* h1D1 = new TH1D("Energy_Raw","h1_1", bin1, xmin, xmax);
    h1D1->SetLineColor(kGray+3);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("E_{#gamma}   [GeV]");
    h1D1->GetYaxis()->SetTitle("Entries");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    TH1D* h1D2 = new TH1D("Energy_New","h1_2", bin1, xmin, xmax);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);
    h1D2->GetXaxis()->SetTitle("E_{#gamma}   [GeV]");
    h1D2->GetYaxis()->SetTitle("Entries");
    h1D2->GetXaxis()->CenterTitle();
    h1D2->GetYaxis()->CenterTitle();

    string str1, str2;

    ifstream file1;
    file1.open(file_name1, ios::in);
    Int_t num1(0);
    while (getline(file1,str1)) num1++;
    cout << num1 << endl;
    file1.clear();
    file1.seekg(0, ios::beg);
    
    ifstream file2;
    file2.open(file_name2, ios::in);
    Int_t num2(0);
    while (getline(file2,str2)) num2++;
    cout << num2 << endl;
    file2.clear();
    file2.seekg(0, ios::beg);
    
    Int_t MaxNo = num1 < num2 ? num1 : num2;
    Int_t N = 0;
    for (int i = 0; i < MaxNo; i++) {
        if (!getline(file1,str1)) continue;
        if (!getline(file2,str2)) continue;
        std::stringstream strStream1(str1);
        std::stringstream strStream2(str2);
        double distance1, value1_1, value1_2;
        double distance2, value2_1, value2_2;
        strStream1 >> distance1 >> value1_1 >> value1_2;
        strStream2 >> distance2 >> value2_1 >> value2_2;
        if (distance1 > cut_min && distance1 < cut_max) {
            h1D1->Fill(value1_1);
            h1D1->Fill(value1_2);
        }
        if (distance2 > cut_min && distance2 < cut_max) {
            h1D2->Fill(value2_1);
            h1D2->Fill(value2_2);
        }
        N++;
    }
    file1.clear();
    file2.clear();
    cout << "Entries : " << N << endl;
    
    double NewRange_min = h1D2->GetMean()-(4.4 - 0.4*(Energy))*(h1D2->GetStdDev());
    double NewRange_max = h1D2->GetMean()+(4.4 - 0.4*(Energy))*(h1D2->GetStdDev());
                                            
    //h1D2->SetAxisRange(NewRange_min, NewRange_max);
    //h1D1->SetAxisRange(NewRange_min, NewRange_max);
    
    int N1_max(-1), N2_max(-1);
    for (int i = 1; i <= bin1; i++){
        int N1 = h1D1->GetBinContent(i);
        int N2 = h1D2->GetBinContent(i);
    if (N1 > N1_max) N1_max = N1;
    if (N2 > N2_max) N2_max = N2;
    }
    
    Double_t set_p00 = N1_max, set_p01 = h1D1->GetMean(), set_p02 = dex*Energy/5.0;
    TF1 *f1=new TF1("f1","[0]*TMath::Gaus(x,[1],[2])",NewRange_min+0.01*Energy, NewRange_max-0.01*Energy);
    f1->SetLineColor(kBlack);
    f1->SetLineStyle(2);
    f1->SetParameters(set_p00,set_p01,set_p02);
    
    Double_t set_p10 = N2_max, set_p11 = h1D2->GetMean(), set_p12 = dex*Energy/5;
    TF1 *f2=new TF1("f2","[0]*TMath::Gaus(x,[1],[2])",NewRange_min+0.01*Energy, NewRange_max-0.01*Energy);
    f2->SetLineColor(kRed);
    f2->SetLineStyle(2);
    f2->SetParameters(set_p10,set_p11,set_p12);
    
    h1D1->Fit(f1,"R");
    h1D2->Fit(f2,"R");

    c1->cd();
    if (N1_max > N2_max){
	h1D1->Draw();
	h1D2->Draw("SAME");
    }else{
    	h1D2->Draw();
    	h1D1->Draw("SAME");
   }

    TLegend * leg = new TLegend(0.68, 0.76, 0.88, 0.86);
    //leg->AddEntry(h1D1, old_file, "L");
    //leg->AddEntry(h1D2, new_file, "L");
    leg->AddEntry(h1D1,"Bump Energy Raw" , "L");
    leg->AddEntry(h1D2,"Bump Energy New", "L");
    leg->Draw();
    
    TPaveText *pt = new TPaveText(0.68, 0.65, 0.88, 0.73,"NDC");
    pt->AddText(Form("CUT : %g ~ %g (X_{0})",cut_min/X0,cut_max/X0));
    pt->AddLine(.0,.5,1.,.5);
    pt->AddText(Form("d : %g (X_{0})",cut_mean/X0));
    pt->Draw();
    
    Double_t mean_OR = f1->GetParameter(1);
    Double_t mean_fix = f2->GetParameter(1);
    Double_t D_mean_OR = f1->GetParError(1);
    Double_t D_mean_fix = f2->GetParError(1);
    
    Double_t sigma_OR = f1->GetParameter(2);
    Double_t sigma_fix = f2->GetParameter(2);
    Double_t D_sigma_OR = f1->GetParError(2);
    Double_t D_sigma_fix = f2->GetParError(2);

    Double_t Resolution_OR = sigma_OR/mean_OR;
    Double_t Resolution_fix = sigma_fix/mean_fix;
    Double_t D_Resolution_OR = sqrt(  pow(D_sigma_OR,2) + pow(sigma_OR,2)*pow(D_mean_OR,2)/pow(mean_OR,2)  )/mean_OR;
    Double_t D_Resolution_fix = sqrt(  pow(D_sigma_fix,2) + pow(sigma_fix,2)*pow(D_mean_fix,2)/pow(mean_fix,2)  )/mean_fix;

    par_file << cut_min << " " << cut_max << " " << mean_OR << " " << Resolution_OR << " " << D_Resolution_OR << " " << mean_fix << " " << Resolution_fix << " " << D_Resolution_fix << endl;
    
    cout << cut_min << " " << cut_max << " " << mean_OR << " " << Resolution_OR << " " << D_Resolution_OR << " " << mean_fix << " " << Resolution_fix << " " << D_Resolution_fix << endl;
    
    ostringstream oo;
    //oo << (cut_min+cut_max)/2.0;
    oo << fixed << setprecision(2) << (cut_min+cut_max)/2.0;
    string str_cut = oo.str();
    TString picture_name= "doc/range_A"+str_NO_Angle+ "_E" + str_Energy +"_resolution_Picture/A"+str_NO_Angle+"_E"+str_Energy+"_CUT"+str_cut+"_resolution_Picture.png";
    c1->Print(picture_name);
    
    par_file.close();
    
    return 0;
}
