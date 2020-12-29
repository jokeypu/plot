int Fit_Read_Par_cp(Int_t NO_Angle){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1;
    out1 << NO_Angle;
    string str_NO_Angle = out1.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_par_cp.txt";
    std::ifstream par_file;
    par_file.open(out_name,std::ios::in);
    std::ofstream AllPar_file;
    AllPar_file.open("doc/AllPar_cp.txt",std::ios::app);
    
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
        g0->SetPoint(N,energy,p0);
        g1->SetPoint(N,energy,p1);
        g2->SetPoint(N,energy,p2);
        g3->SetPoint(N,energy,p3);
        g4->SetPoint(N,energy,p4);
        Max_Energy = energy;
        N++;
    }

    TF1* f0=new TF1("f0","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f0->SetLineColor(kRed-7);
    f0->SetParameters(-0.16,0.24,0.28);
    
    TF1* f1=new TF1("f1","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f1->SetLineColor(kRed-7);
    //f1->SetParameters(-0.16,0.24,0.28);
    
    TF1* f2=new TF1("f2","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f2->SetLineColor(kRed-7);
    f2->SetParameters(0.3,0.6,1.4);
    
    TF1* f3=new TF1("f3","[0]*exp(-[1]*x)+[2]",0,Max_Energy);
    f3->SetLineColor(kRed-7);
    f3->SetParameters(0.04,1.5,0.004);
    
    TF1* f4=new TF1("f4","[0]*x+[1]",0,Max_Energy);
    f4->SetLineColor(kRed-7);
    f4->SetParameters(0.015,0.051);
    
    g0->Fit(f0,"R");
    g1->Fit(f1,"R");
    g2->Fit(f2,"R");
    g3->Fit(f3,"R");
    g4->Fit(f4,"R");
    
    c1->Divide(2, 3);
    c1->cd(1);
    //c1->cd(1)->SetGridx();
    g1->Draw("AP.");
    c1->cd(2);
    //c1->cd(2)->SetGridx();
    g2->Draw("AP.");
    c1->cd(3);
    //c1->cd(3)->SetGridx();
    g3->Draw("AP.");
    c1->cd(4);
    //c1->cd(4)->SetGridx();
    g4->Draw("AP.");
    
    AllPar_file << str_NO_Angle << " "
    << f0->GetParameter(0) << " " << f0->GetParameter(1) << " " << f0->GetParameter(2) << " "
    << f1->GetParameter(0) << " " << f1->GetParameter(1) << " " << f1->GetParameter(2) << " " 
    << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << " " 
    << f3->GetParameter(0) << " " << f3->GetParameter(1) << " " << f3->GetParameter(2) << " " 
    << f4->GetParameter(0) << " " << f4->GetParameter(1) << endl;
    AllPar_file.close();
    par_file.close();
    
    TString picture_name= "doc/AllPar_FitPicture_cp/A"+str_NO_Angle+"_FitPar_cp.png";
    c1->Print(picture_name);
    return 0;
}
