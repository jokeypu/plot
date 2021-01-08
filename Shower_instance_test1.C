int Shower_instance_test1(const char new_file[30], double Energy = 1.0 , int NO_Angle = 7)
{
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_resolution_par_test.txt";
    
    int bin1 = 100*(15.0/(0.6*Energy)),bin2(200);
    float tx(1200),ty(900);
    double xmin(0),xmax(15);
    //double xmin(0),xmax(2.0);
    cout << "-INFO  New File : " << new_file << endl;
    
    string file_str(new_file);
    
    string file_name = "doc/" + file_str + ".txt";
    string file_name_min = "doc/" + file_str + "_test_cent.txt";
    string file_name_max = "doc/" + file_str + "_test_band.txt";
    
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
    
    TH1D* h1D3 = new TH1D("Hist1_3","h1_3", bin1, xmin, xmax);
    h1D3->SetLineColor(kGreen);
    h1D3->SetLineWidth(2);
    h1D3->GetXaxis()->SetTitle("Energy");
    h1D3->GetYaxis()->SetTitle("Entries");
    h1D3->GetXaxis()->CenterTitle();
    h1D3->GetYaxis()->CenterTitle();

    string str;

    ifstream file1;
    file1.open(file_name, ios::in);
    while (getline(file1,str)) {
        double value= atof(str.c_str());
        h1D1->Fill(value);
    }
    
    ifstream file2;
    file2.open(file_name_min, ios::in);
    while (getline(file2,str)) {
        double value= atof(str.c_str());
        h1D2->Fill(value);
    }
    
    ifstream file3;
    file3.open(file_name_max, ios::in);
    while (getline(file3,str)) {
        double value= atof(str.c_str());
        h1D3->Fill(value);
    }
    
    file1.close();
    file2.close();
    file3.close();
    
    double NewRange_min = h1D1->GetMean()-(4.4 - 0.4*(Energy))*(h1D1->GetStdDev());
    double NewRange_max = h1D1->GetMean()+(4.4 - 0.4*(Energy))*(h1D1->GetStdDev());
                                            
    h1D2->SetAxisRange(NewRange_min, NewRange_max);
    h1D1->SetAxisRange(NewRange_min, NewRange_max);
    h1D3->SetAxisRange(NewRange_min, NewRange_max);
    
    TF1 *f = new TF1("f","[0]*TMath::Gaus(x,[1],[2])",0.6,1.3);
    f->SetParameters(200,0.98,0.023);

    c1->cd();
    h1D1->Draw();
    h1D2->Draw("SAME");
    h1D3->Draw("SAME");
    h1D3->Fit(f,"R");
   
    TLegend * leg = new TLegend(0.7,0.7 , 0.9, 0.8);
    leg->AddEntry(h1D1, "Bump Energy", "L");
    leg->AddEntry(h1D2, "Bump Energy cent", "L");
    leg->AddEntry(h1D3, "Bump Energy boundary", "L");
    //leg->AddEntry(h1D1,"Bump Energy old" , "L");
    //leg->AddEntry(h1D2,"Bump Energy new", "L");
    leg->Draw("SAME");
    
    return 0;
}
