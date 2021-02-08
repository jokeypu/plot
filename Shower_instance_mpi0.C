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

int Shower_instance_mpi0(const char old_file[30], const char new_file[30], double Energy = 1.0 , int NO_Angle = 7)
{
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/mpi0_A"+str_NO_Angle+"_resolution_par.txt";
    std::ofstream par_file;
    par_file.open(out_name,std::ios::app);
    
    int bin1 = 100,bin2(200);
    float tx(1200),ty(900);
    double xmin(60),xmax(200);
    //double xmin(0),xmax(2.0);
    cout << "-INFO  Old File : " << old_file << endl;
    cout << "-INFO  New File : " << new_file << endl;
    
    string file_str1(old_file), file_str2(new_file);
    
    string file_name1 = "doc/" + file_str1 + "_mpi0.txt";
    string file_name2 = "doc/" + file_str2 + "_mpi0.txt";
    
    TCanvas* c2=new TCanvas("PANDA2","c2",tx,ty);
    TCanvas* c1=new TCanvas("PANDA1","c1",tx,ty);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(1);
    //gStyle->SetOptStat(111110);
    //gStyle->SetOptStat(1001);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(-0.13);
    c1->SetTopMargin(-0.13);
    c1->SetBottomMargin(0.13);
    
    c2->SetLeftMargin(0.13);
    c2->SetRightMargin(-0.13);
    c2->SetTopMargin(-0.13);
    c2->SetBottomMargin(0.13);
    
    TH1D* h1D1 = new TH1D("m_{#pi^{0}} Raw","h1_1", bin1, xmin, xmax);
    h1D1->SetLineColor(kGray+3);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("m_{#pi0}  [MeV]");
    h1D1->GetYaxis()->SetTitle("Entries");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    TH1D* h1D2 = new TH1D("m_{#pi^{0}} New","h1_2", bin1, xmin, xmax);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);
    h1D2->GetXaxis()->SetTitle("m_{#pi0}  [MeV]");
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
        double value1= atof(str1.c_str());
        double value2= atof(str2.c_str());
        h1D1->Fill(1000*value1);
        h1D2->Fill(1000*value2);
        N++;
    }
    file1.clear();
    file2.clear();
    cout << "Entries : " << N << endl;
    
    c1->cd();
    //h1D2->Draw();
    //h1D1->Draw("SAME");

    h1D1->Draw();
    gPad->Update();
    h1D2->Draw();
    gPad->Update();
    
    //THStack* hs= new THStack("hs","m_{#pi^{0}} histograms  (Raw vs. New)");
    //hs->Add(h1D1);
    //hs->Add(h1D2);
    
    TPaveStats *ps1 = (TPaveStats*)h1D1->GetListOfFunctions()->FindObject("stats");
    TPaveStats *ps2 = (TPaveStats*)h1D2->GetListOfFunctions()->FindObject("stats");
    
    ps1->SetX1NDC(0.132832);
    ps1->SetY1NDC(0.773913);
    ps1->SetX2NDC(0.333333);
    ps1->SetY2NDC(0.874783);
    
    ps2->SetX1NDC(0.132832);
    ps2->SetY1NDC(0.624348);
    ps2->SetX2NDC(0.333333);
    ps2->SetY2NDC(0.725217);
    
    //c2->cd();
    //hs->Draw();
    ps1->Draw();
    ps2->Draw();
    //hs->Draw();
    //h1D1->Draw("SAME");
    
    gPad->Update();
    h1D1->Draw("SAME");
    
    TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    leg->AddEntry(h1D1,"m_{#pi^{0}} Raw" , "L");
    leg->AddEntry(h1D2,"m_{#pi^{0}} New", "L");
    leg->Draw();
    
    double NewRange_min = 130;
    double NewRange_max = 140;
                                            
    h1D2->SetAxisRange(NewRange_min, NewRange_max);
    h1D1->SetAxisRange(NewRange_min, NewRange_max);
    
    int N1_max(-1), N2_max(-1);
     for (int i = 1; i <= bin1; i++){
         int N1 = h1D1->GetBinContent(i);
         int N2 = h1D2->GetBinContent(i);
     if (N1 > N1_max) N1_max = N1;
     if (N2 > N2_max) N2_max = N2;
     }
     
     Double_t set_p00 = N1_max, set_p01 = h1D1->GetMean(), set_p02 = 1.0;
     TF1 *f1=new TF1("f1","[0]*Novosibirsk(x,[1],[2],[3])",NewRange_min, NewRange_max);
     f1->SetLineColor(kBlack);
     f1->SetLineStyle(2);
     f1->SetLineWidth(2);
    f1->SetParameters(set_p00,set_p01,set_p02,0.0);
     
     Double_t set_p10 = N2_max, set_p11 = h1D2->GetMean(), set_p12 = 1.0;
     TF1 *f2=new TF1("f2","[0]*Novosibirsk(x,[1],[2],[3])",NewRange_min, NewRange_max);
     f2->SetLineColor(kRed);
     f2->SetLineStyle(2);
     f2->SetLineWidth(2);
    f2->SetParameters(set_p10,set_p11,set_p12,0.0);
     
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

   
    TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    leg->AddEntry(h1D1, old_file, "L");
    leg->AddEntry(h1D2, new_file, "L");
    //leg->AddEntry(h1D1,"Bump Energy old" , "L");
    //leg->AddEntry(h1D2,"Bump Energy new", "L");
    leg->Draw();
    
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

    par_file << str_Energy << " " << Resolution_OR << " " << D_Resolution_OR << " " << str_Energy << " " << Resolution_fix << " " << D_Resolution_fix << endl;
    
    cout << str_Energy << " " << Resolution_OR << " " << D_Resolution_OR << " " << str_Energy << " " << Resolution_fix << " " << D_Resolution_fix << endl;
    
    TString picture_name= "doc/mpi0_A"+str_NO_Angle+"_resolution_Picture/mpi0_A"+str_NO_Angle+"_E"+str_Energy+"_resolution_Picture.png";
    c1->Print(picture_name);
    
    par_file.close();
    
    return 0;
}
