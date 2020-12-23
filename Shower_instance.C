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

int Shower_instance(const char old_file[20], const char new_file[20], double Energy = 1.0 , int NO_Angle = 7)
{
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_resolution_par.txt";
    std::ofstream par_file;
    par_file.open(out_name,std::ios::app);
    
    int bin1 = 100*(15.0/(0.6*Energy)),bin2(200);
    float tx(1200),ty(900);
    double xmin(0),xmax(15);
    //double xmin(0),xmax(2.0);
    cout << "-INFO  Old File : " << old_file << endl;
    cout << "-INFO  New File : " << new_file << endl;
    
    string file_str1(old_file), file_str2(new_file);
    
    string file_name1 = "doc/" + file_str1 + ".txt";
    string file_name2 = "doc/" + file_str2 + ".txt";
    
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
    
    TH1D* h1D1 = new TH1D("Hist1_1","h1_1", bin1, xmin, xmax);
    h1D1->SetLineColor(kGray+3);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("Energy");
    h1D1->GetYaxis()->SetTitle("Entries");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    TH1D* h1D2 = new TH1D("Hist1_2","h1_2", bin1, xmin, xmax);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);
    h1D2->GetXaxis()->SetTitle("Energy");
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
        h1D1->Fill(value1);
        h1D2->Fill(value2);
        N++;
    }
    file1.clear();
    file2.clear();
    cout << "Entries : " << N << endl;
    
    double NewRange_min = h1D2->GetMean()-(4.4 - 0.4*(Energy))*(h1D2->GetStdDev());
    double NewRange_max = h1D2->GetMean()+(4.4 - 0.4*(Energy))*(h1D2->GetStdDev());
                                            
    h1D2->SetAxisRange(NewRange_min, NewRange_max);
    h1D1->SetAxisRange(NewRange_min, NewRange_max);
    
    
    TF1 *f1=new TF1("f1","[0]*TMath::Gaus(x,[1],[2])",NewRange_min, NewRange_max);
    f1->SetLineColor(kBlack);
    f1->SetParameters(200,1,0.2);
    
    TF1 *f2=new TF1("f2","[0]*TMath::Gaus(x,[1],[2])",NewRange_min, NewRange_max);
    f2->SetParameters(200,1,0.2);
    f2->SetLineColor(kRed);
    
    //h1D1->Fit(f1,"R");
    //h1D2->Fit(f2,"R");
    
    c1->cd();
    h1D2->Draw();
    h1D1->Draw("SAME");
   
    TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    leg->AddEntry(h1D1, old_file, "L");
    leg->AddEntry(h1D2, new_file, "L");
    //leg->AddEntry(h1D1,"Bump Energy old" , "L");
    //leg->AddEntry(h1D2,"Bump Energy new", "L");
    leg->Draw();
    
    par_file << str_Energy << " " << f1->GetParameter(2)/(f1->GetParameter(1)) << " " << f2->GetParameter(2)/(f2->GetParameter(1)) << endl;
    
    TString picture_name= "doc/A"+str_NO_Angle+"_resolution_Picture/A"+str_NO_Angle+"_E"+str_Energy+"_resolution_Picture.png";
    c1->Print(picture_name);
    
    par_file.close();
    
    return 0;
}
