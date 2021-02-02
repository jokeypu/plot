int Shower_instance_pi0_gamma(const char old_file[30], const char new_file[30], double Energy = 1.0 , int NO_Angle = 7)
{
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/pi0_gamma_A"+str_NO_Angle+"_resolution_par.txt";
    std::ofstream par_file;
    par_file.open(out_name,std::ios::app);
    
    int bin1(100),bin2(200);
    float tx(1200),ty(900);
    double xmin(0),xmax(Energy);
    //double xmin(0),xmax(2.0);
    cout << "-INFO  Old File : " << old_file << endl;
    cout << "-INFO  New File : " << new_file << endl;
    
    string file_str1(old_file), file_str2(new_file);
    
    string file_name1_min = "doc/" + file_str1 + "_pi0_gamma_min.txt";
    string file_name1_max = "doc/" + file_str1 + "_pi0_gamma_max.txt";
    string file_name2_min = "doc/" + file_str2 + "_pi0_gamma_min.txt";
    string file_name2_max = "doc/" + file_str2 + "_pi0_gamma_max.txt";
    
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
    
    TH1D* h1D1 = new TH1D("Hist_Raw","h_raw", bin1, xmin, xmax);
    h1D1->SetLineColor(kGray+3);
    h1D1->SetLineWidth(2);
    h1D1->GetXaxis()->SetTitle("E_{#gamma}   [GeV]");
    h1D1->GetYaxis()->SetTitle("Entries");
    h1D1->GetXaxis()->CenterTitle();
    h1D1->GetYaxis()->CenterTitle();
    
    TH1D* h1D1_truth = new TH1D("Hist_Truth","h_truth", bin1, xmin, xmax);
    //TH1D* h1D1_truth = new TH1D("Hist_Truth_Raw","h_truth_raw", bin1, xmin, xmax);
    h1D1_truth->SetLineColor(kGreen+1);
    h1D1_truth->SetLineWidth(2);
    //h1D1_truth->SetLineStyle(2);
    h1D1_truth->GetXaxis()->SetTitle("E_{#gamma}   [GeV]");
    h1D1_truth->GetYaxis()->SetTitle("Entries");
    h1D1_truth->GetXaxis()->CenterTitle();
    h1D1_truth->GetYaxis()->CenterTitle();
    
    TH1D* h1D2 = new TH1D("Hist_New","h_new", bin1, xmin, xmax);
    h1D2->SetLineColor(kRed);
    h1D2->SetLineWidth(2);
    h1D2->GetXaxis()->SetTitle("E_{#gamma}   [GeV]");
    h1D2->GetYaxis()->SetTitle("Entries");
    h1D2->GetXaxis()->CenterTitle();
    h1D2->GetYaxis()->CenterTitle();
    
    TH1D* h1D2_truth = new TH1D("Hist_Truth_New","h_truth_new", bin1, xmin, xmax);
    h1D2_truth->SetLineColor(kGreen+1);
    h1D2_truth->SetLineWidth(2);
    //h1D2_truth->SetLineStyle(2);
    h1D2_truth->GetXaxis()->SetTitle("E_{#gamma}   [GeV]");
    h1D2_truth->GetYaxis()->SetTitle("Entries");
    h1D2_truth->GetXaxis()->CenterTitle();
    h1D2_truth->GetYaxis()->CenterTitle();

    string str1_min, str1_max, str2_min, str2_max;

    ifstream file1_min;
    file1_min.open(file_name1_min, ios::in);
    Int_t num1(0);
    while (getline(file1_min,str1_min)) num1++;
    cout << num1 << endl;
    file1_min.clear();
    file1_min.seekg(0, ios::beg);
    
    ifstream file2_min;
    file2_min.open(file_name2_min, ios::in);
    Int_t num2(0);
    while (getline(file2_min,str2_min)) num2++;
    cout << num2 << endl;
    file2_min.clear();
    file2_min.seekg(0, ios::beg);
    
    ifstream file1_max;
    file1_max.open(file_name1_max, ios::in);
    
    ifstream file2_max;
    file2_max.open(file_name2_max, ios::in);
    
    
    Int_t MaxNo = num1 < num2 ? num1 : num2;
    Int_t N = 0;
    for (int i = 0; i < MaxNo; i++) {
        if (!getline(file1_min,str1_min)) continue;
        if (!getline(file1_max,str1_max)) continue;
        if (!getline(file2_min,str2_min)) continue;
        if (!getline(file2_max,str2_max)) continue;
        std::stringstream strStream1_min(str1_min);
        std::stringstream strStream1_max(str1_max);
        std::stringstream strStream2_min(str2_min);
        std::stringstream strStream2_max(str2_max);
        double value1_min, value1_min_truth;
        double value1_max, value1_max_truth;
        double value2_min, value2_min_truth;
        double value2_max, value2_max_truth;
        strStream1_min >> value1_min >> value1_min_truth;
        strStream1_max >> value1_max >> value1_max_truth;
        strStream2_min >> value2_min >> value2_min_truth;
        strStream2_max >> value2_max >> value2_max_truth;
        
        h1D1->Fill(value1_min);
        h1D1->Fill(value1_max);
        h1D2->Fill(value2_min);
        h1D2->Fill(value2_max);
        
        h1D1_truth->Fill(value1_min_truth);
        h1D1_truth->Fill(value1_max_truth);
        //h1D2_truth->Fill(value2_min_truth);
        //h1D2_truth->Fill(value2_max_truth);
        
        N++;
    }
    file1_min.clear();
    file1_max.clear();
    file2_min.clear();
    file2_max.clear();
    cout << "Entries : " << N << endl;
    
    c1->cd();
    h1D1_truth->Draw();
    h1D2->Draw("SAME");
    h1D1->Draw("SAME");
   
    TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    leg->AddEntry(h1D1_truth, "Bump Energy Truth", "L");
    leg->AddEntry(h1D1,"Bump Energy Raw" , "L");
    leg->AddEntry(h1D2,"Bump Energy New", "L");
    leg->Draw();
    
    TString picture_name= "doc/pi0_gamma_A"+str_NO_Angle+"_resolution_Picture/pi0_gamma_A"+str_NO_Angle+"_E"+str_Energy+"_resolution_Picture.png";
    c1->Print(picture_name);
    
    par_file.close();
    
    return 0;
}
