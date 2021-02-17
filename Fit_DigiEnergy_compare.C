const Double_t RM = 2.00;
Double_t FABC(Double_t r,Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4){
    Double_t xi = r - p2*r*exp(-pow(r/p3/RM,p4));
    return p0*exp(-p1*xi/RM);
}

int Fit_DigiEnergy_compare(std::string dir_name, Int_t NO_Angle, Double_t Energy){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    //root Fit_DigiEnergy_compare.C'("WorkData_1Gamma_A12_E1.0_OR",12,1.0)';
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    //std::string out_name = "doc/compare_A"+str_NO_Angle+"_par_cp.txt";
    //std::ofstream par_file;
    //par_file.open(out_name,std::ios::app);
    
    TCanvas* c1=new TCanvas("PANDA1","test1",800,600);
    TCanvas* c2=new TCanvas("PANDA2","test2",800,600);
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
    gStyle->SetTitleOffset(1.2,"xyz");
    //gStyle->SetPalette(1);
    
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.13);
    c1->SetTopMargin(-0.13);
    c1->SetBottomMargin(0.13);
    
    c2->SetLeftMargin(0.13);
    c2->SetRightMargin(-0.13);
    c2->SetTopMargin(-0.13);
    c2->SetBottomMargin(0.13);
    
    std::string file_name = "doc/"+ dir_name +".txt";
    std::ifstream in_file;
    in_file.open(file_name,std::ios::in);
    
    TGraph *g = new TGraph();
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kYellow-4);
    g->GetYaxis()->SetTitle("E_{truth}   [GeV]");
    g->GetXaxis()->SetTitle("t   [X_{0}]");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    Double_t distance_cut = 3.5;
    Int_t binx = 200, biny = 200;
    
    TH2D* h = new TH2D("Hist", "h", binx, 0, distance_cut, biny, 0, Energy);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 1);
    h->GetYaxis()->SetTitle("E_{truth}   [GeV]");
    h->GetXaxis()->SetTitle("r   [cm]");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        if (distance > distance_cut) continue;
        h->Fill(distance,energy);
        N++;
    }
    
    for (int i = 1; i < 201; i++){
        for (int j = 1; j < 201; j++){
            if (N < 30000) {if (h->GetBinContent(i,j)<2) h->SetBinContent(i,j,0);}
            else if (h->GetBinContent(i,j)<2) h->SetBinContent(i,j,0);
        }
    }
    
    int cunt = 0;
    int step = 9;
    for (int i = 5; i < binx+1; i+=step){
        Double_t mean = h->ProfileY("px",i,i+step-1)->GetMean();
        Double_t wx = distance_cut/binx;
        Double_t nx =( (i+(step-1)/2)*wx - wx/2  );
        g->SetPoint(cunt, nx, mean);
        cunt++;
    }
    
    TF1* f=new TF1("f1","(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    f->SetParameters(0.8*Energy, 2.5, 0.9, 1.4, 3);
    f->SetLineWidth(3);
    f->SetLineColor(kRed);
    
    TF1* f_cp=new TF1("f_cp","exp(-1.25*x)",0,distance_cut);
    f_cp->SetLineWidth(3);
    f_cp->SetLineStyle(2);
    f_cp->SetLineColor(kCyan);
    
    TF1* f_test=new TF1("f_test","exp(-[0]*pow(x-[1],[2]))",0,distance_cut);
    f_test->SetLineWidth(3);
    //f_test->SetLineStyle(2);
    f_test->SetLineColor(kYellow);
    
    c1->cd();
    h->Draw("PCOLZ");
    g->Draw("Psame");
    g->Fit(f,"R");
    f->Draw("SAME");
    //f_cp->Draw("SAME");
    //par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << endl;
    
    //TString picture_name= "doc/compare_A"+str_NO_Angle+"_FitPicture_cp/A"+str_NO_Angle+"_E"+str_Energy+"_FitPar_cp.png";
    //c1->Print(picture_name);
    
    TGraph *g_Error = new TGraph();
    g_Error->SetMarkerStyle(7);
    g_Error->SetMarkerColorAlpha(kRed, 0.5);
    g_Error->GetYaxis()->SetTitle("E_{digi}/E_{seed}");
    g_Error->GetXaxis()->SetTitle("r   [cm]");
    g_Error->GetXaxis()->CenterTitle();
    g_Error->GetYaxis()->CenterTitle();
    
    TGraph *g_Error_cp = new TGraph();
    g_Error_cp->SetMarkerStyle(7);
    g_Error_cp->SetMarkerColorAlpha(kAzure+3, 0.5);
    g_Error_cp->GetYaxis()->SetTitle("E_{digi}/E_{seed}");
    g_Error_cp->GetXaxis()->SetTitle("r   [cm]");
    g_Error_cp->GetXaxis()->CenterTitle();
    g_Error_cp->GetYaxis()->CenterTitle();
    
    TH2D* h_Error = new TH2D("Hist_Error","h_Error",200,0,distance_cut,200,0,1.1);
    h_Error->SetMarkerStyle(7);
    h_Error->SetMarkerColorAlpha(kAzure+3, 0.5);
    h_Error->GetYaxis()->SetTitle("E_{digi}/E_{seed}");
    h_Error->GetXaxis()->SetTitle("r   [cm]");
    h_Error->GetXaxis()->CenterTitle();
    h_Error->GetYaxis()->CenterTitle();
    h_Error->GetZaxis()->CenterTitle();
    
    std::string file_name_n = "doc/compare_"+ dir_name +".txt";
    std::ifstream in_file_n;
    in_file_n.open(file_name_n,std::ios::in);
    
    //in_file_n.clear();
    //in_file_n.seekg(0, ios::beg);
    N = 0;
    while (std::getline(in_file_n, str)) {
        std::stringstream strStream(str);
        float dseed, distance, Eseed, energy;
        strStream >> dseed >> distance >> Eseed >> energy;
        if (distance > distance_cut) continue;
        h_Error->Fill(distance,energy/Eseed);
        g_Error_cp->SetPoint(N,distance,energy/Eseed);
        g_Error->SetPoint(N,distance,FABC(distance,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4))/FABC(dseed,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4)));
        N++;
    }
    
    c2->cd();
    g_Error->GetYaxis()->SetRangeUser(0,1);
    g_Error_cp->GetYaxis()->SetRangeUser(0,1);
    g_Error->GetXaxis()->SetRangeUser(0,distance_cut);
    g_Error_cp->GetXaxis()->SetRangeUser(0,distance_cut);
    //h_Error->Draw("PCOLZ");
    g_Error_cp->Draw("AP.");
    g_Error->Draw("Psame");
    g_Error->Fit(f_test,"R");
    //f->Draw("SAME");
    f_cp->Draw("SAME");
    
    TLegend * leg = new TLegend(0.625, 0.6, 0.88, 0.86);
    leg->AddEntry(f_cp,"Raw: exp(-#epsilonr/R_{M})" , "L");
    leg->AddEntry(g_Error,"New: f(r)/f(r_{seed})" , "p");
    leg->AddEntry(g_Error_cp,"data: E_{digi}/E_{seed}", "p");
    leg->Draw();
    
    //TString picture_name_error= "doc/compare_A"+str_NO_Angle+"_FitPicture_cp/Error_A"+str_NO_Angle+"_E"+str_Energy+"_FitPar_cp.png";
    //c2->Print(picture_name_error);
    
    in_file.close();
    in_file_n.close();
    //par_file.close();
    return 0;
}
