int Fit_Read_Par_one(Int_t NO_Angle){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1;
    out1 << NO_Angle;
    string str_NO_Angle = out1.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_par_cp.txt";
    std::ifstream par_file;
    par_file.open(out_name,std::ios::in);
    std::ofstream AllPar_file;
    AllPar_file.open("doc/AllPar_cp.txt",std::ios::app);
    
    TString ts = "Angle "+str_NO_Angle;
    TCanvas* c2=new TCanvas("PANDA2",ts,800,600);
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
    
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(-0.13);
    c1->SetTopMargin(-0.13);
    c1->SetBottomMargin(0.13);
    
    TGraph *g0 = new TGraph();
    g0->SetMarkerStyle(20);
    g0->SetMarkerColorAlpha(kGray-1, 1);
    g0->GetXaxis()->SetTitle("Energy");
    g0->GetYaxis()->SetTitle("p0");
    g0->GetXaxis()->CenterTitle();
    g0->GetYaxis()->CenterTitle();
    
    TGraph *g1 = new TGraph();
    g1->SetMarkerStyle(21);
    g1->SetMarkerColorAlpha(kRed-3, 1);
    g1->GetXaxis()->SetTitle("Energy");
    g1->GetYaxis()->SetTitle("p1");
    g1->GetXaxis()->CenterTitle();
    g1->GetYaxis()->CenterTitle();
    
    TGraph *g2 = new TGraph();
    g2->SetMarkerStyle(22);
    g2->SetMarkerColorAlpha(kGreen+1, 1);
    g2->GetXaxis()->SetTitle("Energy");
    g2->GetYaxis()->SetTitle("p2");
    g2->GetXaxis()->CenterTitle();
    g2->GetYaxis()->CenterTitle();
    
    TGraph *g3 = new TGraph();
    g3->SetMarkerStyle(33);
    g3->SetMarkerColorAlpha(kAzure+3, 1);
    g3->GetXaxis()->SetTitle("Energy");
    g3->GetYaxis()->SetTitle("p3");
    g3->GetXaxis()->CenterTitle();
    g3->GetYaxis()->CenterTitle();
    
    TGraph *g4 = new TGraph();
    g4->SetMarkerStyle(34);
    g4->SetMarkerColorAlpha(kBlue+1, 1);
    g4->GetXaxis()->SetTitle("Energy");
    g4->GetYaxis()->SetTitle("p4");
    g4->GetXaxis()->CenterTitle();
    g4->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    Double_t Max_Energy(0);
    while (std::getline(par_file, str)) {
        std::stringstream strStream(str);
        float energy, p0, p1, p2, p3, p4;
        strStream >> energy >> p0 >> p1 >> p2 >> p3 >> p4;
        if (energy < 0.2) continue;
        g0->SetPoint(N,energy,p0);
        g1->SetPoint(N,energy,p1);
        g2->SetPoint(N,energy,p2);
        g3->SetPoint(N,energy,p3);
        g4->SetPoint(N,energy,p4);
        if (energy>Max_Energy) Max_Energy = energy;
        N++;
    }
    
    TF1* f0=new TF1("f0","[0]*x+[1]",0,Max_Energy);
    f0->SetLineColor(kRed-7);
    f0->SetParameters(-0.004,1.5,1.38);
    
    TF1* f1=new TF1("f1","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f1->SetLineColor(kRed-7);
    f1->SetParameters(-0.18, 1.4, 2.7);
    f1->SetParLimits(0,-0.01,1.0);
    f1->SetParLimits(1, 0.13.0);
    f1->SetParLimits(2,2.0,3.5);

    TF1* f2=new TF1("f2","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f2->SetLineColor(kRed-7);
    f2->SetParameters(-0.23,2.8,0.92);
    f2->SetParLimits(0,-0.01,1.0);
    f2->SetParLimits(1,1.0,5.0);
    f2->SetParLimits(2,0.5,1.5);
    
    TF1* f3=new TF1("f3","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f3->SetLineColor(kRed-7);
    f3->SetParameters(0.62,8.6,0.77);
    f3->SetParLimits(0,-0.01,1.0);
    f3->SetParLimits(1,1.0,15.0);
    f3->SetParLimits(2,0.3,1.5);
    
    TF1* f4=new TF1("f4","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f4->SetLineColor(kRed-7);
    f4->SetParameters(-2.8,1.8,5.27);
    f4->SetParLimits(0,-0.1,5.0);
    f4->SetParLimits(1,1.0,5.0);
    f4->SetParLimits(2,2.1,7.0);
    
    
    g0->Fit(f0,"R");
    g1->Fit(f1,"R");
    g2->Fit(f2,"R");
    g3->Fit(f3,"R");
    g4->Fit(f4,"R");
    
    //g0->GetYaxis()->SetRangeUser(0,6.5);
    //g1->GetYaxis()->SetRangeUser(0,3.5);
    //g2->GetYaxis()->SetRangeUser(0.87,0.98);
    //g3->GetYaxis()->SetRangeUser(0.6,0.8);
    //g4->GetYaxis()->SetRangeUser(3.5,5.5);
    
    g0->GetYaxis()->SetRangeUser(0,6.5);
    //g1->GetYaxis()->SetRangeUser(1,4);
    g1->GetYaxis()->SetRangeUser(2.5,3.1);
    //g2->GetYaxis()->SetRangeUser(0.5,1.5);
    g2->GetYaxis()->SetRangeUser(0.5,1.1);
    //g3->GetYaxis()->SetRangeUser(0.2,1.2);
    g3->GetYaxis()->SetRangeUser(0.5,1.1);
    //g4->GetYaxis()->SetRangeUser(2,7);
    g4->GetYaxis()->SetRangeUser(2.5,7.5);
    
    c1->Divide(2, 2);
    c1->cd(1);
    //c1->cd(1)->SetGridx();
    g1->Draw("AP.");
    c1->cd(2);
    //c1->cd(2)->SetGridx();
    g2->Draw("AP.");
    c1->cd(3);
    //c1->cd(3)->SetGridx();
    g3->Draw("AP.");
    //c1->cd(1);
    //c1->cd(4)->SetGridx();
    //g3->Draw("AP.");
    c1->cd(4);
    g4->Draw("AP.");

    c2->cd();
    g0->Draw("AP.");
    
    AllPar_file << str_NO_Angle << " "
    << f1->GetParameter(0) << " " << f1->GetParameter(1) << " " << f1->GetParameter(2) << " "
    << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << " "
    << f3->GetParameter(0) << " " << f3->GetParameter(1) << " " << f3->GetParameter(2) << " "
    << f4->GetParameter(0) << " " << f4->GetParameter(1) << " " << f4->GetParameter(2) << " "
    << endl;
    AllPar_file.close();
    par_file.close();
    
    TString picture_name= "doc/AllPar_FitPicture_cp/A"+str_NO_Angle+"_FitPar_cp.png";
    c1->Print(picture_name);
    return 0;
}
