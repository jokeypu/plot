int Fit_Read_Par(Int_t NO_Angle){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1;
    out1 << NO_Angle;
    string str_NO_Angle = out1.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_par.txt";
    std::ifstream par_file;
    par_file.open(out_name,std::ios::in);
    
    TCanvas* c1=new TCanvas("PANDA1","test1",800,600);
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
    
    TGraph *g1 = new TGraph();
    g1->SetMarkerStyle(20);
    g1->SetMarkerColorAlpha(kRed-3, 0.5);
    g1->GetXaxis()->SetTitle("Energy");
    g1->GetYaxis()->SetTitle("p1");
    g1->GetXaxis()->CenterTitle();
    g1->GetYaxis()->CenterTitle();
    
    TGraph *g2 = new TGraph();
    g2->SetMarkerStyle(22);
    g2->SetMarkerColorAlpha(kGreen+1, 0.5);
    g2->GetXaxis()->SetTitle("Energy");
    g2->GetYaxis()->SetTitle("p2");
    g2->GetXaxis()->CenterTitle();
    g2->GetYaxis()->CenterTitle();
    
    TGraph *g3 = new TGraph();
    g3->SetMarkerStyle(33);
    g3->SetMarkerColorAlpha(kAzure+3, 0.5);
    g3->GetXaxis()->SetTitle("Energy");
    g3->GetYaxis()->SetTitle("p3");
    g3->GetXaxis()->CenterTitle();
    g3->GetYaxis()->CenterTitle();
    
    TGraph *g4 = new TGraph();
    g4->SetMarkerStyle(34);
    g4->SetMarkerColorAlpha(kBlue+1, 0.5);
    g4->GetXaxis()->SetTitle("Energy");
    g4->GetYaxis()->SetTitle("p1");
    g4->GetXaxis()->CenterTitle();
    g4->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    while (std::getline(par_file, str)) {
        std::stringstream strStream(str);
        float energy, p1, p2, p3, p4;
        strStream >> energy >> p1 >> p2 >> p3 >> p4;
        g1->SetPoint(N,energy,p1);
        g2->SetPoint(N,energy,p2);
        g3->SetPoint(N,energy,p3);
        g4->SetPoint(N,energy,p4);
        N++;
    }

    c1->Divide(2, 2);
    c1->cd(1);
    g1->Draw("AP.");
    c1->cd(2);
    g2->Draw("AP.");
    c1->cd(3);
    g3->Draw("AP.");
    c1->cd(4);
    g4->Draw("AP.");
    
    par_file.close();
    return 0;
}
